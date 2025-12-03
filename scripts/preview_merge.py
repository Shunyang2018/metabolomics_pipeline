from __future__ import annotations

import csv
import sys
from pathlib import Path


def main(input_dir: str, output_csv: str):
    # Import from src without installing the package
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from metabo_pipeline.io_msdial import merge_folder_to_wide_csv  # type: ignore

    in_dir = Path(input_dir)
    out_path = Path(output_csv)
    summary = merge_folder_to_wide_csv(in_dir, out_path, recursive=False)
    print(f"Merged files: {summary['files']}, rows: {summary['rows']}")

    # Print first row with first few sample columns
    print("\nFirst row (truncated sample intensities):")
    with out_path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        row = next(r)
        fixed = {
            k: row[k]
            for k in ("chrom", "annotation_level", "alignment_id", "metabolite_name")
        }
        # show up to 5 sample columns
        sample_cols = [
            c
            for c in r.fieldnames
            if c
            not in (
                "chrom",
                "annotation_level",
                "alignment_id",
                "rt_min",
                "mz",
                "metabolite_name",
                "adduct",
            )
        ]
        preview_samples = {c: row[c] for c in sample_cols[:5]}
        print({**fixed, **preview_samples})

    # Print a couple rows for specific sample_ids if present
    # Show a couple of normalized sample columns that should exist
    print("\nExample normalized sample columns present:")
    print(sample_cols[:10])


if __name__ == "__main__":
    in_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    out_csv = sys.argv[2] if len(sys.argv) > 2 else "outputs/merged_preview.csv"
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    main(in_dir, out_csv)
