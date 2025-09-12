from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


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

