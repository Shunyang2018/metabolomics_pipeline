"""Reproduce the public synthetic example without network access or SIRIUS."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import pandas as pd

import metabo_pipeline
from metabo_pipeline import constants
from metabo_pipeline.bioactivity import match_bioactives
from metabo_pipeline.classify import classify_level12_with_classyfire
from metabo_pipeline.merge import merge_folder_to_wide_csv
from metabo_pipeline.reporting import write_report

HERE = Path(__file__).resolve().parent


def run_demo(output_dir: Path) -> dict:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Choose an empty output directory to keep this run independent.")
    output_dir.mkdir(parents=True, exist_ok=True)
    input_dir = HERE / "synthetic" / "alignments"
    merged = output_dir / "merged.csv"
    summary = merge_folder_to_wide_csv(input_dir, merged)
    # Copy the explicitly synthetic cache; never modify repository fixtures.
    cache = output_dir / "taxonomy_cache.json"
    cache.write_bytes((HERE / "synthetic" / "taxonomy_cache.json").read_bytes())
    final = output_dir / "annotated.csv"
    taxonomy = classify_level12_with_classyfire(
        str(merged), str(final), cache_path=str(cache), offline=True
    )
    bioactivity = match_bioactives(
        final, HERE / "synthetic" / "reference.csv", output_dir / "reference_matches.csv"
    )
    inputs = sorted(p for p in (HERE / "synthetic").rglob("*") if p.is_file())
    package_root = Path(metabo_pipeline.__file__).parent
    manifest = {
        "schema_version": 1,
        "dataset": "Synthetic MS-DIAL demonstration; not experimental observations",
        "external_services": {"classyfire": "synthetic cache only", "sirius": "export only"},
        "versions": {
            "python": platform.python_version(),
            **{name: version(name) for name in ("metabolomics-pipeline", "pandas", "numpy")},
        },
        "inputs": [
            {
                "path": p.relative_to(HERE).as_posix(),
                "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
            }
            for p in inputs
        ],
        "parameters": {
            name: getattr(constants, name)
            for name in (
                "MSMS_MIN_IONS",
                "SNR_MIN",
                "BLANK_FOLD_MIN",
                "PRESENT_PERCENT_MIN",
                "CV_PERCENT_MAX",
                "DEDUP_RT_WINDOW_MIN",
                "DEDUP_MZ_PPM",
            )
        },
        "real_sample_tokens": sorted(constants.REAL_SAMPLE_TOKENS or []),
        "source_sha256": {
            p.relative_to(package_root).as_posix(): hashlib.sha256(
                p.read_text(encoding="utf-8").encode("utf-8")
            ).hexdigest()
            for p in sorted(package_root.rglob("*.py"))
        },
        "demo_script_sha256": hashlib.sha256(
            Path(__file__).read_text(encoding="utf-8").encode("utf-8")
        ).hexdigest(),
        "merge": summary,
        "taxonomy": taxonomy,
        "bioactivity": bioactivity,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    write_report(pd.read_csv(final), manifest, output_dir / "report.html", synthetic=True)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/demo"))
    args = parser.parse_args()
    try:
        result = run_demo(args.output_dir)
    except ValueError as exc:
        parser.exit(2, f"{exc}\n")
    print(f"Retained {result['merge']['rows']} features; report: {args.output_dir / 'report.html'}")
