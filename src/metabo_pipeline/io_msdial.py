from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

# Optional pandas import for faster merging
try:
    import pandas as pd  # type: ignore
    import numpy as np  # type: ignore

    HAS_PANDAS = True
except Exception:
    HAS_PANDAS = False

# Shared thresholds/constants
from .constants import (
    MSMS_MIN_IONS,
    SNR_MIN,
    ISOMER_RT_WINDOW_MIN,
    BLANK_FOLD_MIN,
    PRESENT_PERCENT_MIN,
    CV_PERCENT_MAX,
)


@dataclass
class MSDialMetadata:
    classes: Dict[str, str]
    file_types: Dict[str, str]
    injection_order: Dict[str, int]
    batch_id: Dict[str, str]


@dataclass
class MSDialSummary:
    path: Path
    n_features: int
    samples: List[str]
    metadata: MSDialMetadata


def _read_first_rows(path: Path, n: int = 6, encoding: str = "utf-8-sig") -> List[List[str]]:
    rows: List[List[str]] = []
    with path.open("r", encoding=encoding, newline="") as f:
        r = csv.reader(f)
        for i, row in enumerate(r):
            rows.append(row)
            if i + 1 >= n:
                break
    return rows


def _parse_embedded_metadata(rows: List[List[str]]) -> Tuple[List[str], MSDialMetadata]:
    # Expect first 4 lines to be: Class, File type, Injection order, Batch ID
    # Values start after some leading empty columns; sample names appear in the main header.
    class_row = rows[0] if len(rows) > 0 else []
    filetype_row = rows[1] if len(rows) > 1 else []
    injection_row = rows[2] if len(rows) > 2 else []
    batch_row = rows[3] if len(rows) > 3 else []
    header_row = rows[4] if len(rows) > 4 else []

    # Sample columns start after the fixed MS-DIAL feature metadata columns.
    # Find the index of the first sample column by locating a known feature header then moving right.
    # We assume header contains "Alignment ID" and later specific sample names; sample names in header align with metadata rows.
    try:
        align_idx = header_row.index("Alignment ID")
    except ValueError:
        # Fallback: assume first 30 columns are feature metadata; this is conservative.
        align_idx = 0

    # Heuristic: sample headers begin after the fixed columns; find the last known fixed column index, then take the rest as samples.
    KNOWN_LAST_FIXED = "MS/MS spectrum"
    try:
        fixed_last_idx = header_row.index(KNOWN_LAST_FIXED)
        sample_start = fixed_last_idx + 1
    except ValueError:
        # Fallback if not present
        sample_start = max(align_idx + 1, 30)

    sample_names = header_row[sample_start:]

    def _to_map(name: str, row: List[str]) -> Dict[str, str]:
        return {s: v for s, v in zip(sample_names, row[sample_start:])}

    def _to_int_map(row: List[str]) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for s, v in zip(sample_names, row[sample_start:]):
            try:
                out[s] = int(v)
            except Exception:
                # leave missing or non-int as 0
                out[s] = 0
        return out

    metadata = MSDialMetadata(
        classes=_to_map("Class", class_row),
        file_types=_to_map("File type", filetype_row),
        injection_order=_to_int_map(injection_row),
        batch_id=_to_map("Batch ID", batch_row),
    )

    return sample_names, metadata


def summarize_alignment_table(path: Path) -> MSDialSummary:
    rows = _read_first_rows(path)
    if len(rows) < 5:
        raise ValueError(f"Not enough header rows for MS-DIAL file: {path}")
    sample_names, metadata = _parse_embedded_metadata(rows)

    # Count features by streaming the file once; skip top 5 header rows
    n_features = 0
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        for i, _ in enumerate(r):
            if i >= 5:
                n_features += 1

    return MSDialSummary(path=path, n_features=n_features, samples=sample_names, metadata=metadata)


# --- Merging utilities ---

def _infer_chrom_from_name(name: str) -> str:
    n = name.lower()
    if "hilic" in n:
        return "HILIC"
    if "c18" in n:
        return "C18"
    if "lipid" in n:
        return "Lipidomics"
    return "unknown"


def _normalize_sample_id_core(name: str) -> str:
    v = str(name or "").strip().lower()
    v = re.sub(r"[\s\-]+", "_", v)
    token_pat = re.compile(r"(?i)(^|[\W_])(lipidomics|lipids|lipid|hilic|c18|pos|neg|ms1)(?=($|[\W_]))")
    v = token_pat.sub(lambda m: "_" if m.group(1) else "", v)
    v = re.sub(r"_+", "_", v).strip("_")
    return v


def iter_alignment_long_rows(path: Path) -> Iterable[Dict[str, str]]:
    """Yield long-format rows from an MS-DIAL alignment table.

    Columns yielded per row:
      - source_file, chrom, polarity
      - alignment_id, rt_min, mz, adduct, metabolite_name
      - sample_id, intensity
    """
    first_rows = _read_first_rows(path)
    sample_names, _md = _parse_embedded_metadata(first_rows)

    header_row = first_rows[4]
    # Find indices for feature metadata columns we care about
    def col_idx(col: str) -> int:
        try:
            return header_row.index(col)
        except ValueError:
            return -1

    idx_alignment = col_idx("Alignment ID")
    idx_rt = col_idx("Average Rt(min)")
    idx_mz = col_idx("Average Mz")
    idx_name = col_idx("Metabolite name")
    idx_adduct = col_idx("Adduct type")
    idx_msms = col_idx("MS/MS spectrum")
    idx_msms_assigned = col_idx("MS/MS assigned")
    idx_weighted_dot = col_idx("Weighted dot product")
    idx_reverse_dot = col_idx("Reverse dot product")
    idx_matched_peaks = col_idx("Matched peaks count")

    # Determine where sample intensities start
    try:
        fixed_last_idx = header_row.index("MS/MS spectrum")
        sample_start = fixed_last_idx + 1
    except ValueError:
        sample_start = max(30, idx_alignment + 1)

    chrom = _infer_chrom_from_name(path.name)

    def _count_msms_ions(msms: str) -> int:
        if not msms:
            return 0
        s = msms.strip().lower()
        if s in {"", "null", "na", "none"}:
            return 0
        # Expect space-separated mz:intensity pairs
        count = 0
        for tok in s.split():
            if ":" not in tok:
                continue
            parts = tok.split(":", 1)
            if len(parts) != 2:
                continue
            try:
                inten = float(parts[1])
            except Exception:
                continue
            if inten > 0:
                count += 1
        return count

    def _to_float(val: str) -> float | None:
        try:
            if val is None:
                return None
            v = str(val).strip()
            if v == "" or v.lower() in {"null", "na", "none"}:
                return None
            return float(v)
        except Exception:
            return None

    def _to_int(val: str) -> int | None:
        f = _to_float(val)
        if f is None:
            return None
        try:
            return int(round(f))
        except Exception:
            return None

    def _assign_annotation_level(row: List[str], metabolite_name: str) -> str:
        # Determine if annotated: metabolite name not unknown
        annotated = bool(metabolite_name and metabolite_name.strip().lower() != "unknown")
        if not annotated:
            # As a fallback, consider MS/MS assigned flag
            if 0 <= idx_msms_assigned < len(row):
                annotated = str(row[idx_msms_assigned]).strip().upper() == "TRUE"

        if annotated:
            wdot = _to_float(row[idx_weighted_dot]) if 0 <= idx_weighted_dot < len(row) else None
            rdot = _to_float(row[idx_reverse_dot]) if 0 <= idx_reverse_dot < len(row) else None
            mcount = _to_int(row[idx_matched_peaks]) if 0 <= idx_matched_peaks < len(row) else None
            if (
                wdot is not None
                and rdot is not None
                and mcount is not None
                and wdot > 750
                and abs(wdot - rdot) < 200
                and mcount >= 3
            ):
                return "1"
            return "2"
        else:
            return "3"

    # Prepare sample ID normalizer: remove assay/polarity tokens and tidy
    token_pat = re.compile(r"(?i)(^|[\W_])(lipidomics|lipids|lipid|hilic|c18|pos|neg|ms1)(?=($|[\W_]))")

    def _normalize_sample_id(s: str) -> str:
        v = str(s or "").strip().lower()
        v = re.sub(r"[\s\-]+", "_", v)
        v = token_pat.sub(lambda m: "_" if m.group(1) else "", v)
        v = re.sub(r"_+", "_", v).strip("_")
        return v

    # Stream remaining lines
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        for i, row in enumerate(r):
            if i < 5:
                continue
            # Guard against short rows
            if len(row) <= sample_start:
                continue
            # Filter by MS/MS presence and minimum ions (>=3)
            msms_val = row[idx_msms] if idx_msms >= 0 and idx_msms < len(row) else ""
            if _count_msms_ions(msms_val) < MSMS_MIN_IONS:
                continue
            alignment_id = row[idx_alignment] if idx_alignment >= 0 else ""
            rt_min = row[idx_rt] if idx_rt >= 0 else ""
            mz = row[idx_mz] if idx_mz >= 0 else ""
            name = row[idx_name] if idx_name >= 0 else ""
            adduct = row[idx_adduct] if idx_adduct >= 0 else ""
            ann_level = _assign_annotation_level(row, name)

            for sample, val in zip(sample_names, row[sample_start:]):
                norm_sample = _normalize_sample_id(sample)
                yield {
                    "chrom": chrom,
                    "annotation_level": ann_level,
                    "alignment_id": alignment_id,
                    "rt_min": rt_min,
                    "mz": mz,
                    "metabolite_name": name,
                    "adduct": adduct,
                    "sample_id": norm_sample,
                    "intensity": val,
                }


def list_alignment_files(input_dir: Path, recursive: bool = False) -> List[Path]:
    patterns = ["*.csv", "*.txt"]
    files: List[Path] = []
    if recursive:
        for pat in patterns:
            files.extend(input_dir.rglob(pat))
    else:
        for pat in patterns:
            files.extend(input_dir.glob(pat))
    # Basic filter to avoid obviously non-MS-DIAL CSVs if desired could be added
    return sorted({p.resolve() for p in files})


def merge_folder_to_long_csv(input_dir: Path, output_csv: Path, recursive: bool = False, engine: str = "pandas") -> Dict[str, int]:
    """Merge all MS-DIAL tables in folder to one long-format CSV.

    Returns a summary dict with counts.
    """
    files = list_alignment_files(input_dir, recursive=recursive)
    # Use pandas path if available and selected
    if engine == "pandas" and HAS_PANDAS:
        return _merge_folder_to_long_csv_pandas(files, output_csv)
    counts = {"files": len(files), "rows": 0}
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", encoding="utf-8", newline="") as f_out:
        w = csv.DictWriter(
            f_out,
            fieldnames=[
                "chrom",
                "annotation_level",
                "alignment_id",
                "rt_min",
                "mz",
                "metabolite_name",
                "adduct",
                "sample_id",
                "intensity",
            ],
        )
        w.writeheader()
        for p in files:
            for row in iter_alignment_long_rows(p):
                w.writerow(row)
                counts["rows"] += 1

    return counts


def _merge_folder_to_long_csv_pandas(files: List[Path], output_csv: Path) -> Dict[str, int]:
    rows_total = 0
    frames: List["pd.DataFrame"] = []

    # Helper functions reused from csv path
    def count_msms_ions(msms: str) -> int:
        if not isinstance(msms, str):
            return 0
        s = msms.strip().lower()
        if s in {"", "null", "na", "none"}:
            return 0
        cnt = 0
        for tok in s.split():
            if ":" not in tok:
                continue
            try:
                inten = float(tok.split(":", 1)[1])
            except Exception:
                continue
            if inten > 0:
                cnt += 1
        return cnt

    def assign_level(row: "pd.Series") -> str:
        # Annotated if name not Unknown or MS/MS assigned TRUE
        name = str(row.get("Metabolite name", "")).strip()
        annotated = bool(name and name.lower() != "unknown")
        if not annotated:
            ms2 = str(row.get("MS/MS assigned", "")).strip().upper()
            annotated = (ms2 == "TRUE")
        if annotated:
            try:
                wdot = float(row.get("Weighted dot product", float("nan")))
            except Exception:
                wdot = float("nan")
            try:
                rdot = float(row.get("Reverse dot product", float("nan")))
            except Exception:
                rdot = float("nan")
            try:
                mcount = int(round(float(row.get("Matched peaks count", float("nan")))))
            except Exception:
                mcount = -1
            if (pd.notna(wdot) and pd.notna(rdot) and mcount >= 3 and (wdot > 750) and (abs(wdot - rdot) < 200)):
                return "1"
            return "2"
        return "3"

    for p in files:
        # Read with header row at index 4 (0-based), so header=4
        try:
            df = pd.read_csv(p, header=4, encoding="utf-8-sig")
        except Exception:
            # Fallback to python engine for odd quoting
            df = pd.read_csv(p, header=4, encoding="utf-8-sig", engine="python")

        # Determine sample columns: after 'MS/MS spectrum'
        cols = list(df.columns)
        if "MS/MS spectrum" in cols:
            sample_start = cols.index("MS/MS spectrum") + 1
        else:
            sample_start = max(30, cols.index("Alignment ID") + 1 if "Alignment ID" in cols else 30)
        sample_cols = cols[sample_start:]

        # Filter features based on MS/MS ions count >= constant
        df["_msms_ion_count"] = df.get("MS/MS spectrum", "").map(count_msms_ions)
        df = df[df["_msms_ion_count"] >= MSMS_MIN_IONS].copy()
        # Additional pre-merge filter: S/N average threshold
        if "S/N average" in df.columns:
            sn = pd.to_numeric(df["S/N average"], errors="coerce")
            df = df[sn >= SNR_MIN].copy()

        # Assign annotation level per feature
        df["annotation_level"] = df.apply(assign_level, axis=1)

        # Melt to long
        id_cols = [
            "Alignment ID",
            "Average Rt(min)",
            "Average Mz",
            "Metabolite name",
            "Adduct type",
            "annotation_level",
        ]
        present_id_cols = [c for c in id_cols if c in df.columns]
        long_df = df.melt(id_vars=present_id_cols, value_vars=sample_cols, var_name="sample_id", value_name="intensity")

        # Normalize sample names to common core across assays
        token_pat = re.compile(r"(?i)(^|[\W_])(lipidomics|lipids|lipid|hilic|c18|pos|neg|ms1)(?=($|[\W_]))")
        def _norm_sample(s: str) -> str:
            v = str(s or "").strip().lower()
            v = re.sub(r"[\s\-]+", "_", v)
            v = token_pat.sub(lambda m: "_" if m.group(1) else "", v)
            v = re.sub(r"_+", "_", v).strip("_")
            return v
        long_df["sample_id"] = long_df["sample_id"].map(_norm_sample)

        # Add chrom column at front
        chrom = _infer_chrom_from_name(p.name)
        long_df.insert(0, "chrom", chrom)

        # Rename columns to normalized names
        rename_map = {
            "Alignment ID": "alignment_id",
            "Average Rt(min)": "rt_min",
            "Average Mz": "mz",
            "Metabolite name": "metabolite_name",
            "Adduct type": "adduct",
        }
        long_df = long_df.rename(columns=rename_map)

        # Keep exact output column order
        out_cols = [
            "chrom",
            "annotation_level",
            "alignment_id",
            "rt_min",
            "mz",
            "metabolite_name",
            "adduct",
            "sample_id",
            "intensity",
        ]
        # Some columns might be missing in rare files; ensure existence
        for c in out_cols:
            if c not in long_df.columns:
                long_df[c] = ""
        long_df = long_df[out_cols]

        rows_total += len(long_df)
        frames.append(long_df)

    if frames:
        merged = pd.concat(frames, ignore_index=True)
    else:
        merged = pd.DataFrame(columns=[
            "chrom","annotation_level","alignment_id","rt_min","mz","metabolite_name","adduct","sample_id","intensity"
        ])

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_csv, index=False)
    return {"files": len(files), "rows": int(rows_total)}


def merge_folder_to_wide_csv(input_dir: Path, output_csv: Path, recursive: bool = False) -> Dict[str, int]:
    """Merge MS-DIAL tables into a wide table: one feature per row, normalized sample columns."""
    if not HAS_PANDAS:
        raise RuntimeError("Pandas is required for wide-format merging. Install pandas and retry.")

    files = list_alignment_files(input_dir, recursive=recursive)
    frames: List["pd.DataFrame"] = []
    # Stats
    totals = {"raw": 0, "after_msms": 0, "after_snr": 0, "after_pass_all": 0}
    per_file_stats: List[Dict[str, int]] = []

    def count_msms_ions(msms: str) -> int:
        if not isinstance(msms, str):
            return 0
        s = msms.strip().lower()
        if s in {"", "null", "na", "none"}:
            return 0
        cnt = 0
        for tok in s.split():
            if ":" not in tok:
                continue
            try:
                inten = float(tok.split(":", 1)[1])
            except Exception:
                continue
            if inten > 0:
                cnt += 1
        return cnt

    def assign_level(row: "pd.Series") -> str:
        name = str(row.get("Metabolite name", "")).strip()
        annotated = bool(name and name.lower() != "unknown")
        if not annotated:
            ms2 = str(row.get("MS/MS assigned", "")).strip().upper()
            annotated = (ms2 == "TRUE")
        if annotated:
            try:
                wdot = float(row.get("Weighted dot product", float("nan")))
            except Exception:
                wdot = float("nan")
            try:
                rdot = float(row.get("Reverse dot product", float("nan")))
            except Exception:
                rdot = float("nan")
            try:
                mcount = int(round(float(row.get("Matched peaks count", float("nan")))))
            except Exception:
                mcount = -1
            if (pd.notna(wdot) and pd.notna(rdot) and mcount >= 3 and (wdot > 750) and (abs(wdot - rdot) < 200)):
                return "1"
            return "2"
        return "3"

    # Mapping of additional MS-DIAL columns to normalized names
    meta_rename = {
        "Formula": "formula",
        "Ontology": "ontology",
        "INCHIKEY": "inchi_key",
        "SMILES": "smiles",
        "Annotation tag (VS1.0)": "annotation_tag",
        "Comment": "comment",
        "MS/MS assigned": "msms_assigned",
        "RT matched": "rt_matched",
        "m/z matched": "mz_matched",
        "MS/MS matched": "msms_matched",
        "Weighted dot product": "weighted_dot",
        "Reverse dot product": "reverse_dot",
        "Simple dot product": "simple_dot",
        "Matched peaks count": "matched_peaks_count",
        "Matched peaks percentage": "matched_peaks_pct",
        "Total score": "total_score",
        "S/N average": "snr_avg",
        "Reference RT": "reference_rt",
        "Reference m/z": "reference_mz",
        "Spectrum reference file name": "spectrum_ref",
    }

    for p in files:
        try:
            df = pd.read_csv(p, header=4, encoding="utf-8-sig")
        except Exception:
            df = pd.read_csv(p, header=4, encoding="utf-8-sig", engine="python")

        cols = list(df.columns)
        if "MS/MS spectrum" in cols:
            sample_start = cols.index("MS/MS spectrum") + 1
        else:
            sample_start = max(30, cols.index("Alignment ID") + 1 if "Alignment ID" in cols else 30)
        sample_cols = cols[sample_start:]

        # Stats: raw feature count
        raw_count = int(df.shape[0])

        df["_msms_ion_count"] = df.get("MS/MS spectrum", "").map(count_msms_ions)
        df = df[df["_msms_ion_count"] >= MSMS_MIN_IONS].copy()
        after_msms = int(df.shape[0])

        df["annotation_level"] = df.apply(assign_level, axis=1)

        # Normalize sample col names
        rename_samples = {c: _normalize_sample_id_core(c) for c in sample_cols}
        df = df.rename(columns=rename_samples)

        # Build feature frame
        id_rename = {
            "Alignment ID": "alignment_id",
            "Average Rt(min)": "rt_min",
            "Average Mz": "mz",
            "Metabolite name": "metabolite_name",
            "Adduct type": "adduct",
        }
        present_ids = [c for c in id_rename.keys() if c in df.columns]
        present_meta = [c for c in meta_rename.keys() if c in df.columns]
        base_cols = present_ids + ["annotation_level"] + present_meta
        feat_df = df[base_cols + list(rename_samples.values())].copy()
        feat_df = feat_df.rename(columns={**id_rename, **meta_rename})

        # Isomer labeling by RT clustering within metabolite_name + adduct
        RT_CLUSTER_WINDOW = ISOMER_RT_WINDOW_MIN  # minutes
        if all(c in feat_df.columns for c in ("metabolite_name", "adduct", "rt_min")):
            feat_df["_rt_num"] = pd.to_numeric(feat_df["rt_min"], errors="coerce")
            def _label_group(g):
                g = g.sort_values("_rt_num")
                rts = g["_rt_num"].to_numpy()
                labels = []
                cluster = 1
                prev = None
                for rt in rts:
                    if prev is None or pd.isna(prev) or pd.isna(rt) or abs(rt - prev) > RT_CLUSTER_WINDOW:
                        cluster = cluster if prev is None else cluster + 1
                    labels.append(f"isomer_{cluster}")
                    prev = rt
                out = g.copy()
                out["_isomer_label_tmp"] = labels
                return out
            feat_df = feat_df.groupby(["metabolite_name", "adduct"], dropna=False, group_keys=False).apply(_label_group)
            feat_df["isomer_label"] = feat_df.pop("_isomer_label_tmp")
            feat_df = feat_df.drop(columns=["_rt_num"], errors="ignore")
        else:
            feat_df["isomer_label"] = ""

        # Insert leading columns in desired order: isomer_label, chrom
        feat_df.insert(0, "chrom", _infer_chrom_from_name(p.name))
        feat_df.insert(0, "isomer_label", feat_df.pop("isomer_label"))

        # --- Per-group QC metrics ---
        # Identify blank column if present
        blank_col = "blank" if "blank" in feat_df.columns else None
        if blank_col is not None:
            blank_series = pd.to_numeric(feat_df[blank_col], errors="coerce")
        else:
            blank_series = None

        # Build mapping of group -> replicate columns for this file
        group_cols: Dict[str, List[str]] = {}
        sample_norm_cols = list(rename_samples.values())
        pat = re.compile(r"^m2_([a-z0-9]+)_(.+)$")
        for c in sample_norm_cols:
            if c in ("blank",):
                continue
            m = pat.match(c)
            if not m:
                continue
            grp = m.group(1)
            if grp in ("pool", "qc"):
                continue
            group_cols.setdefault(grp, []).append(c)

        # Thresholds imported from constants

        pass_cols_this_file: List[str] = []

        # Compute metrics per group and append as trailing columns
        for grp, cols_grp in group_cols.items():
            vals = feat_df[cols_grp].apply(pd.to_numeric, errors="coerce")
            nrep = max(1, len(cols_grp))
            # Presence: intensity > 0 treated as detected
            present = vals.gt(0)
            present_frac = present.sum(axis=1) / float(nrep)
            present_percent = (present_frac * 100.0).astype(float)

            # Blank 7x using max across replicates vs blank
            max_val = vals.max(axis=1, skipna=True)
            if blank_col is not None:
                # Use numpy.nan instead of pd.NA to avoid astype issues on older pandas
                denom = blank_series.replace(0, np.nan)
                blank_fold = (max_val / denom)
            else:
                # No blank available → leave as NaN
                blank_fold = pd.Series([float('nan')]*len(feat_df), index=feat_df.index)

            # CV <= 40% among present replicates (need at least 2 present)
            vals_present = vals.where(present)
            mean = vals_present.mean(axis=1, skipna=True)
            std = vals_present.std(axis=1, ddof=1, skipna=True)
            # Use numpy.nan instead of pd.NA to ensure float-friendly ops
            denom_mean = mean.replace(0, np.nan)
            cv_percent = (std / denom_mean) * 100.0
            # If fewer than 2 replicates present, set NaN
            cv_percent = cv_percent.where(present.sum(axis=1) >= 2)

            # Append columns
            feat_df[f"blank_fold_{grp}"] = blank_fold
            feat_df[f"present_percent_{grp}"] = present_percent
            feat_df[f"cv_percent_{grp}"] = cv_percent

            # Pass flag for this group
            grp_pass = (
                (blank_fold >= BLANK_FOLD_MIN) &
                (present_percent >= PRESENT_PERCENT_MIN) &
                (cv_percent <= CV_PERCENT_MAX)
            )
            feat_df[f"pass_{grp}"] = grp_pass.fillna(False)
            pass_cols_this_file.append(f"pass_{grp}")

        # Overall pass: require ALL groups (present in this file) to pass
        if pass_cols_this_file:
            feat_df["pass_all_groups"] = feat_df[pass_cols_this_file].all(axis=1).fillna(False)
            # Filter to passing features only
            feat_df = feat_df[feat_df["pass_all_groups"]]

        # Count after pass_all gating
        after_pass_all = int(feat_df.shape[0])

        # Update stats
        per_file_stats.append({
            "file": p.name,  # type: ignore
            "raw": raw_count,
            "after_msms": after_msms,
            "after_snr": int(df.shape[0]),
            "after_pass_all": after_pass_all,
        })
        totals["raw"] += raw_count
        totals["after_msms"] += after_msms
        totals["after_snr"] += int(df.shape[0])
        totals["after_pass_all"] += after_pass_all

        frames.append(feat_df)

    if frames:
        # unify columns
        id_preferred = [
            "isomer_label",
            "chrom",
            "annotation_level",
            "alignment_id",
            "rt_min",
            "mz",
            "metabolite_name",
            "adduct",
            # additional metadata (only those present will be kept in order)
            "formula",
            "ontology",
            "inchi_key",
            "smiles",
            "annotation_tag",
            "msms_assigned",
            "rt_matched",
            "mz_matched",
            "msms_matched",
            "weighted_dot",
            "reverse_dot",
            "simple_dot",
            "matched_peaks_count",
            "matched_peaks_pct",
            "total_score",
            "snr_avg",
            "reference_rt",
            "reference_mz",
            "spectrum_ref",
            "comment",
        ]
        all_cols = set().union(*[set(fr.columns) for fr in frames])
        id_order = [c for c in id_preferred if c in all_cols]
        # Identify metric columns (place at end)
        metric_cols = sorted(
            c for c in all_cols if str(c).startswith(("blank_fold_", "present_percent_", "cv_percent_"))
        )
        pass_cols = sorted(c for c in all_cols if str(c).startswith("pass_") or c == "pass_all_groups")
        # Identify sample columns by pattern and exclude filters/ids
        sample_pat = re.compile(r"^m2_[a-z0-9]+_.+")
        sample_cols = sorted(
            c for c in all_cols
            if (c not in set(id_order)) and (c not in set(metric_cols)) and sample_pat.match(str(c))
        )
        # Everything else (rare leftovers) goes between id and samples
        other_cols = sorted(c for c in all_cols if c not in set(id_order + sample_cols + metric_cols + pass_cols))
        out_cols = id_order + other_cols + sample_cols + metric_cols + pass_cols
        frames = [fr.reindex(columns=out_cols) for fr in frames]
        merged = pd.concat(frames, ignore_index=True)
        rows_pre_dedup = int(merged.shape[0])

        # After merging all files: append isomer label to metabolite_name and deduplicate by name
        if "isomer_label" in merged.columns and "metabolite_name" in merged.columns:
            merged["metabolite_name"] = merged.apply(
                lambda r: f"{r['metabolite_name']}_{r['isomer_label']}" if isinstance(r.get("isomer_label"), str) and r["isomer_label"] else r["metabolite_name"],
                axis=1,
            )

        # Compute cv median percent from per-group metrics
        cv_cols_all = [c for c in merged.columns if str(c).startswith("cv_percent_")]
        if cv_cols_all:
            merged["cv_median_percent"] = merged[cv_cols_all].median(axis=1, skipna=True)
        else:
            merged["cv_median_percent"] = np.nan

        # Compute average intensity across normalized sample columns
        sample_pat_all = re.compile(r"^m2_[a-z0-9]+_.+")
        sample_cols_all = [c for c in merged.columns if sample_pat_all.match(str(c))]
        if sample_cols_all:
            merged["_avg_intensity"] = merged[sample_cols_all].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True)
        else:
            merged["_avg_intensity"] = np.nan

        # Deduplicate rows with the same metabolite_name (post-isomer suffix):
        # keep the one with lowest cv_median_percent, then highest average intensity
        if "metabolite_name" in merged.columns:
            def _pick_best(g: "pd.DataFrame") -> "pd.DataFrame":
                g = g.copy()
                cv = pd.to_numeric(g["cv_median_percent"], errors="coerce")
                g["_cv_sort"] = cv.fillna(np.inf)
                g = g.sort_values(["_cv_sort", "_avg_intensity"], ascending=[True, False])
                return g.head(1)

            merged = merged.groupby("metabolite_name", dropna=False, group_keys=False).apply(_pick_best)
            merged = merged.drop(columns=["_cv_sort"], errors="ignore")
        merged = merged.drop(columns=["_avg_intensity"], errors="ignore")
        rows_post_dedup = int(merged.shape[0])
    else:
        merged = pd.DataFrame(columns=["chrom", "annotation_level", "alignment_id", "rt_min", "mz", "metabolite_name", "adduct"])  # type: ignore
        rows_pre_dedup = 0
        rows_post_dedup = 0

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_csv, index=False)
    return {
        "files": len(files),
        "rows": int(len(merged)),
        "rows_pre_dedup": rows_pre_dedup,
        "dedup_dropped": rows_pre_dedup - rows_post_dedup,
        "totals": totals,
        "per_file": per_file_stats,
    }
