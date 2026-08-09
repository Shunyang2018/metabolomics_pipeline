"""
Thin compatibility wrapper for MS-DIAL routines (merge + validate summary).

The main merge implementation lives in src/metabo_pipeline/merge.py.
This module also provides summarize_alignment_table for the CLI validate command.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .merge import merge_folder_to_wide_csv  # re-export for CLI compatibility

#: The last fixed MS-DIAL column; sample columns start immediately after it.
_MSMS_COLUMN = "MS/MS spectrum"


@dataclass
class MSDialMetadata:
    classes: dict[str, str]
    file_types: dict[str, str]
    injection_order: dict[str, int]
    batch_id: dict[str, str]


@dataclass
class MSDialSummary:
    path: Path
    n_features: int
    samples: list[str]
    metadata: MSDialMetadata


def _read_first_rows(path: Path, n: int = 6, encoding: str = "utf-8-sig") -> list[list[str]]:
    """Read the first few header rows from an MS-DIAL export."""
    rows: list[list[str]] = []
    with path.open("r", encoding=encoding, newline="") as f:
        r = csv.reader(f)
        for i, row in enumerate(r):
            rows.append(row)
            if i + 1 >= n:
                break
    return rows


def summarize_alignment_table(path: Path) -> MSDialSummary:
    """Summarize key metadata from an MS-DIAL alignment table."""
    rows = _read_first_rows(path)
    if len(rows) < 5:
        raise ValueError(f"Not enough header rows for MS-DIAL file: {path}")

    class_row = rows[0]
    filetype_row = rows[1]
    injection_row = rows[2]
    batch_row = rows[3]
    header_row = rows[4]

    # Determine sample start after the last fixed column (MS/MS spectrum)
    # Sample columns begin after the last fixed MS-DIAL column.
    if _MSMS_COLUMN not in header_row:
        raise ValueError(
            f"Missing the {_MSMS_COLUMN!r} column in {path}; sample columns are "
            "located relative to it. Is this an MS-DIAL alignment export?"
        )
    sample_start = header_row.index(_MSMS_COLUMN) + 1
    sample_names = header_row[sample_start:]

    # A malformed export can carry metadata rows shorter than the header row;
    # zip truncates to the shorter of the two rather than raising.
    def _map_from(row: list[str]) -> dict[str, str]:
        return dict(zip(sample_names, row[sample_start:], strict=False))

    def _int_map_from(row: list[str]) -> dict[str, int]:
        return {
            s: int(v)
            for s, v in zip(sample_names, row[sample_start:], strict=False)
            if str(v).strip().lstrip("-").isdigit()
        }

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


__all__ = [
    "merge_folder_to_wide_csv",
    "summarize_alignment_table",
    "MSDialSummary",
    "MSDialMetadata",
]
