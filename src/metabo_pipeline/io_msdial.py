from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

# Optional pandas import for faster merging
try:
    import pandas as pd  # type: ignore

    HAS_PANDAS = True
except Exception:
    HAS_PANDAS = False


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
            if _count_msms_ions(msms_val) < 3:
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

        # Filter features based on MS/MS ions count >= 3
        df["_msms_ion_count"] = df.get("MS/MS spectrum", "").map(count_msms_ions)
        df = df[df["_msms_ion_count"] >= 3].copy()

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

        df["_msms_ion_count"] = df.get("MS/MS spectrum", "").map(count_msms_ions)
        df = df[df["_msms_ion_count"] >= 3].copy()

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
        feat_df.insert(0, "chrom", _infer_chrom_from_name(p.name))
        frames.append(feat_df)

    if frames:
        # unify columns
        id_preferred = [
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
        sample_cols = sorted(c for c in all_cols if c not in set(id_order))
        out_cols = id_order + sample_cols
        frames = [fr.reindex(columns=out_cols) for fr in frames]
        merged = pd.concat(frames, ignore_index=True)
    else:
        merged = pd.DataFrame(columns=["chrom", "annotation_level", "alignment_id", "rt_min", "mz", "metabolite_name", "adduct"])  # type: ignore

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_csv, index=False)
    return {"files": len(files), "rows": int(len(merged))}
