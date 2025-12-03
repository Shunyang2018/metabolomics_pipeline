"""
Command to collect/export SIRIUS results from .sirius database files into CSV.
This allows inspection of SIRIUS annotations before merging into the final output.
"""

from pathlib import Path

import typer

from .. import constants
from ..logging import get_logger
from ..sirius_collect import collect_sirius_results
from ..sirius_export import export_tsv_summaries

log = get_logger()


def sirius_collect(
    output_dir: str = typer.Argument(
        None,
        help="Directory containing sirius_pos.sirius and sirius_neg.sirius files (default: INPUT_DIR/outputs from constants.py)",
    ),
    output_csv: str = typer.Option(
        None,
        "--output",
        "-o",
        help="Path to write SIRIUS results CSV (default: OUTPUT_DIR/sirius_identifications.csv)",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Force re-scan of SIRIUS databases (ignore cached extractions)",
    ),
):
    """
    Collect and export SIRIUS results from .sirius database files into CSV.

    This command extracts SIRIUS annotations (formulas, structures, CANOPUS classifications)
    from the binary .sirius database files and writes them to a human-readable CSV file.
    This allows you to inspect what SIRIUS found before merging into the final output.

    The output CSV will contain:
    - feature_id: Feature ID from the original MS-DIAL data
    - SIRIUS formula predictions
    - Structure predictions (if available)
    - CANOPUS classifications (superclass, class, most specific)

    Examples:
        # Collect SIRIUS results from default location
        metabo sirius-collect

        # Collect from specific directory
        metabo sirius-collect /path/to/outputs

        # Force re-scan (ignore cache)
        metabo sirius-collect --force
    """
    # Default paths
    if output_dir is None:
        output_dir = str(Path(constants.INPUT_DIR) / "outputs")

    out_dir = Path(output_dir)

    if output_csv is None:
        output_csv = str(out_dir / "sirius_identifications.csv")

    output_csv_path = Path(output_csv)

    # SIRIUS 6 creates TSV summary files in directories (not .sirius database files)
    # These directories are created by SIRIUS with the write-summaries command
    pos_tsv_dir = out_dir / "sirius_pos"
    neg_tsv_dir = out_dir / "sirius_neg"

    # Check if TSV summary directories exist
    has_pos = pos_tsv_dir.exists() and pos_tsv_dir.is_dir()
    has_neg = neg_tsv_dir.exists() and neg_tsv_dir.is_dir()

    # If TSV directories don't exist, try to export from .sirius database files
    if not has_pos or not has_neg:
        pos_db = out_dir / "sirius_pos.sirius"
        neg_db = out_dir / "sirius_neg.sirius"

        # Find SIRIUS executable
        sirius_exe = constants.SIRIUS_EXECUTABLE
        if not Path(sirius_exe).exists():
            # Try common locations
            for candidate in [
                "/Applications/sirius.app/Contents/MacOS/sirius",
                "/usr/local/bin/sirius",
                str(
                    Path.home()
                    / "Downloads"
                    / "sirius.app"
                    / "Contents"
                    / "MacOS"
                    / "sirius"
                ),
            ]:
                if Path(candidate).exists():
                    sirius_exe = candidate
                    break

        # Try to export POS if needed
        if not has_pos and pos_db.exists() and pos_db.stat().st_size > 1024:
            log.info("TSV summaries not found for POS. Attempting to export...")
            if export_tsv_summaries(pos_db, pos_tsv_dir, sirius_exe):
                has_pos = True

        # Try to export NEG if needed
        if not has_neg and neg_db.exists() and neg_db.stat().st_size > 1024:
            log.info("TSV summaries not found for NEG. Attempting to export...")
            if export_tsv_summaries(neg_db, neg_tsv_dir, sirius_exe):
                has_neg = True

    # Final check
    if not has_pos and not has_neg:
        log.error(f"No SIRIUS results found in: {out_dir}")
        log.error("  Looked for:")
        log.error(f"    - TSV directories: {pos_tsv_dir.name}/, {neg_tsv_dir.name}/")
        log.error("    - Database files: sirius_pos.sirius, sirius_neg.sirius")
        log.error("\n  Run 'metabo sirius' first to generate SIRIUS results")
        raise typer.Exit(code=1)

    log.info("Collecting SIRIUS results...")
    log.info(f"  POS summaries: {pos_tsv_dir if has_pos else 'N/A'}")
    log.info(f"  NEG summaries: {neg_tsv_dir if has_neg else 'N/A'}")
    log.info(f"  Output CSV: {output_csv_path}")

    # Collect results
    summary = collect_sirius_results(
        pos_dir=pos_tsv_dir if has_pos else None,
        neg_dir=neg_tsv_dir if has_neg else None,
        output_csv=output_csv_path,
        join_with_merged=None,  # Don't join, just export
        force_rescan=force,
    )

    # Report statistics
    log.ok(f"✓ SIRIUS results collected: {output_csv_path.name}")
    log.info(f"  Total compounds: {summary.get('rows', 0)}")

    n_pos = summary.get("n_pos", 0)
    n_neg = summary.get("n_neg", 0)
    log.info(f"  POS: {n_pos}, NEG: {n_neg}")

    # Structure and formula counts
    has_struct = summary.get("n_has_struct", 0)
    has_formula = summary.get("n_has_formula", 0)
    has_canopus = summary.get("n_has_canopus", 0)

    if summary.get("rows", 0) > 0:
        log.info(
            f"  With structures: {has_struct} ({has_struct / summary['rows'] * 100:.1f}%)"
        )
        log.info(
            f"  With formulas: {has_formula} ({has_formula / summary['rows'] * 100:.1f}%)"
        )
        log.info(
            f"  With CANOPUS: {has_canopus} ({has_canopus / summary['rows'] * 100:.1f}%)"
        )

    log.info("")
    log.info("Next steps:")
    log.info("  1. Inspect SIRIUS results:")
    log.info(f"     - Open {output_csv_path.name} in Excel or a CSV viewer")
    log.info("  2. Run final merge to combine with MS-DIAL data:")
    log.info("     - metabo final")
