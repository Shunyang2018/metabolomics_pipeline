"""Bioactivity command - match MS-DIAL features to a bioactivity reference database."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from .. import constants
from ..bioactivity import match_bioactives
from ..logging import get_logger

log = get_logger()


def bioactivity(
    merged_csv: Optional[str] = typer.Argument(
        None,
        help=(
            "Merged CSV to match (default: looks in OUTPUT_DIR for "
            "merged_classyfire_final.csv, then merged_classyfire.csv, then merged.csv)"
        ),
    ),
    output_dir: str = typer.Option(
        None, "--output-dir", "-o", help="Output directory (default: INPUT_DIR/outputs)"
    ),
    db: str = typer.Option(
        None,
        "--db",
        help="Path to the bioactivity database CSV (default: constants.BIOACTIVITY_DB_PATH)",
    ),
    output: str = typer.Option(
        None, "--output", help="Output CSV path (default: OUTPUT_DIR/bioactives.csv)"
    ),
):
    """Match features against a bioactivity database via InChIKey skeleton.

    Matches on the first block of the InChIKey (molecular skeleton, ignoring
    stereochemistry/protonation state) and writes one row per feature x
    sample-it-was-detected-in x bioactivity hit.
    """
    out_dir = Path(output_dir) if output_dir else Path(constants.INPUT_DIR) / "outputs"

    if merged_csv is None:
        for name in ("merged_classyfire_final.csv", "merged_classyfire.csv", "merged.csv"):
            candidate = out_dir / name
            if candidate.exists():
                merged_csv = str(candidate)
                break
        if merged_csv is None:
            log.error(f"No merged CSV found in {out_dir}")
            raise typer.Exit(code=1)

    db_path = Path(db) if db else Path(constants.BIOACTIVITY_DB_PATH)
    if not db_path.exists():
        log.error(f"Bioactivity database not found: {db_path}")
        raise typer.Exit(code=1)

    output_csv = Path(output) if output else out_dir / "bioactives.csv"

    log.info(f"Matching {merged_csv} against {db_path.name} ...")
    summary = match_bioactives(Path(merged_csv), db_path, output_csv)
    log.ok(f"✓ Bioactivity match complete: {output_csv}")
    log.info(f"  • Features with InChIKey: {summary['features_with_inchikey']}")
    log.info(f"  • Unique matched molecular skeletons: {summary['unique_matched_skeletons']}")
    log.info(f"  • Output rows (feature x sample x hit): {summary['output_rows']}")
    log.info(f"  • Unique features with a bioactivity hit: {summary['unique_features_matched']}")
