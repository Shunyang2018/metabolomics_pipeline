"""
Thin compatibility wrapper for MS-DIAL routines (merge + validate summary).

The main merge implementation lives in src/metabo_pipeline/merge.py.
This module also provides summarize_alignment_table for the CLI validate command.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Callable

from .merge import merge_folder_to_wide_csv  # re-export for CLI compatibility


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


def summarize_alignment_table(path: Path) -> MSDialSummary:
    rows = _read_first_rows(path)
    if len(rows) < 5:
        raise ValueError(f"Not enough header rows for MS-DIAL file: {path}")

    class_row = rows[0]
    filetype_row = rows[1]
    injection_row = rows[2]
    batch_row = rows[3]
    header_row = rows[4]

    # Determine sample start after the last fixed column (MS/MS spectrum)
    try:
        sample_start = header_row.index("MS/MS spectrum") + 1
    except ValueError:
        # fallback heuristic
        try:
            align_idx = header_row.index("Alignment ID")
        except ValueError:
            align_idx = 0
        sample_start = max(align_idx + 1, 30)

    sample_names = header_row[sample_start:]

    def _map_from(row: List[str]) -> Dict[str, str]:
        return {s: v for s, v in zip(sample_names, row[sample_start:])}

    def _int_map_from(row: List[str]) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for s, v in zip(sample_names, row[sample_start:]):
            try:
                out[s] = int(v)
            except Exception:
                out[s] = 0
        return out

    meta = MSDialMetadata(
        classes=_map_from(class_row),
        file_types=_map_from(filetype_row),
        injection_order=_int_map_from(injection_row),
        batch_id=_map_from(batch_row),
    )

    # Count feature rows (skip 5 header rows)
    n_features = 0
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        r = csv.reader(f)
        for i, _ in enumerate(r):
            if i >= 5:
                n_features += 1

    return MSDialSummary(path=path, n_features=n_features, samples=sample_names, metadata=meta)


__all__ = ["merge_folder_to_wide_csv", "summarize_alignment_table", "MSDialSummary", "MSDialMetadata"]
