from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


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

def _infer_chrom_and_polarity_from_name(name: str) -> Tuple[str, str]:
    n = name.lower()
    chrom = "unknown"
    if "hilic" in n:
        chrom = "HILIC"
    elif "c18" in n:
        chrom = "C18"
    elif "lipid" in n:
        chrom = "Lipidomics"

    polarity = "unknown"
    if "pos" in n:
        polarity = "POS"
    if "neg" in n:
        polarity = "NEG"
    return chrom, polarity


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

    # Determine where sample intensities start
    try:
        fixed_last_idx = header_row.index("MS/MS spectrum")
        sample_start = fixed_last_idx + 1
    except ValueError:
        sample_start = max(30, idx_alignment + 1)

    chrom, polarity = _infer_chrom_and_polarity_from_name(path.name)

    # Stream remaining lines
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        for i, row in enumerate(r):
            if i < 5:
                continue
            # Guard against short rows
            if len(row) <= sample_start:
                continue
            alignment_id = row[idx_alignment] if idx_alignment >= 0 else ""
            rt_min = row[idx_rt] if idx_rt >= 0 else ""
            mz = row[idx_mz] if idx_mz >= 0 else ""
            name = row[idx_name] if idx_name >= 0 else ""
            adduct = row[idx_adduct] if idx_adduct >= 0 else ""

            for sample, val in zip(sample_names, row[sample_start:]):
                yield {
                    "source_file": path.name,
                    "chrom": chrom,
                    "polarity": polarity,
                    "alignment_id": alignment_id,
                    "rt_min": rt_min,
                    "mz": mz,
                    "metabolite_name": name,
                    "adduct": adduct,
                    "sample_id": sample,
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


def merge_folder_to_long_csv(input_dir: Path, output_csv: Path, recursive: bool = False) -> Dict[str, int]:
    """Merge all MS-DIAL tables in folder to one long-format CSV.

    Returns a summary dict with counts.
    """
    files = list_alignment_files(input_dir, recursive=recursive)
    counts = {"files": len(files), "rows": 0}
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", encoding="utf-8", newline="") as f_out:
        w = csv.DictWriter(
            f_out,
            fieldnames=[
                "source_file",
                "chrom",
                "polarity",
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
