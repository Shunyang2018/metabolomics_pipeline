from __future__ import annotations

import csv
import sys
from pathlib import Path


def main(input_dir: str, output_csv: str):
    # Import from src without installing the package
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from metabo_pipeline.io_msdial import merge_folder_to_long_csv  # type: ignore

    in_dir = Path(input_dir)
    out_path = Path(output_csv)
    summary = merge_folder_to_long_csv(in_dir, out_path, recursive=False)
    print(f"Merged files: {summary['files']}, rows: {summary['rows']}")

    # Print first 5 rows
    print("\nFirst 5 rows:")
    with out_path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for i, row in enumerate(r):
            print({k: row[k] for k in (
                "chrom","source_file","annotation_level","alignment_id","metabolite_name","sample_id","intensity"
            )})
            if i >= 4:
                break

    # Print a couple rows for specific sample_ids if present
    targets = {"M2_Lipids_TV_5+6_pos", "M2_C18_TV_5+6_POS"}
    print("\nRows for selected sample_id(s):")
    count = 0
    with out_path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get("sample_id") in targets:
                print({k: row[k] for k in (
                    "chrom","source_file","annotation_level","metabolite_name","sample_id","intensity"
                )})
                count += 1
                if count >= 6:
                    break


if __name__ == "__main__":
    in_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    out_csv = sys.argv[2] if len(sys.argv) > 2 else "outputs/merged_preview.csv"
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    main(in_dir, out_csv)

