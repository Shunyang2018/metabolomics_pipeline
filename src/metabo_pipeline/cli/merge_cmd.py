"""Merge command - merge MS-DIAL files."""

from __future__ import annotations

from pathlib import Path

import typer

from .. import constants
from ..io_msdial import merge_folder_to_wide_csv
from ..logging import get_logger

log = get_logger()


def merge(
    input_dir: str = typer.Argument(
        None,
        help="Folder containing MS-DIAL CSV/TXT files (default: INPUT_DIR from constants.py)",
    ),
    output_csv: str = typer.Option(
        None,
        "--output",
        "--output-csv",
        "-o",
        help="Path to write merged CSV (default: INPUT_DIR/outputs/merged.csv)",
    ),
    recursive: bool = typer.Option(False, help="Recurse into subfolders"),
):
    """Merge HILIC/C18/Lipidomics files (wide format, one feature per row)."""
    # Use defaults from constants if not provided
    if input_dir is None:
        input_dir = constants.INPUT_DIR

    if output_csv is None:
        output_csv = str(Path(input_dir) / "outputs" / "merged.csv")

    in_dir = Path(input_dir)
    if not in_dir.exists() or not in_dir.is_dir():
        log.error(f"Input directory not found: {in_dir}")
        raise typer.Exit(code=2)
    out = Path(output_csv)
    try:

        def _progress(rec: dict):
            log.info(
                f"[file] {rec.get('file')}: raw={rec.get('raw')}, MS/MS={rec.get('after_msms')}, "
                f"S/N={rec.get('after_snr')}, pass_any={rec.get('after_pass')}"
            )

        summary = merge_folder_to_wide_csv(
            in_dir, out, recursive=recursive, progress=_progress
        )
    except Exception as e:
        log.error(f"Failed to merge files from {in_dir}: {e}")
        raise typer.Exit(code=1)

    log.ok(f"Merged {summary['files']} files → {out}")
    # Print stats if available (wide mode)
    totals = summary.get("totals") if isinstance(summary, dict) else None
    if totals:
        log.info(
            f"Totals — raw: {totals.get('raw')}, after MS/MS: {totals.get('after_msms')}, "
            f"after S/N: {totals.get('after_snr')}, after pass_any: {totals.get('after_pass')}"
        )
        log.info(
            f"Annotation levels after pass_any — L1: {totals.get('ann1')}, L2: {totals.get('ann2')}, L3: {totals.get('ann3')}"
        )
        pre = summary.get("rows_pre_dedup")
        post = summary.get("rows")
        drop = summary.get("dedup_dropped")
        log.info(f"Dedup — pre: {pre}, post: {post}, dropped: {drop}")
        ann_post = summary.get("ann_post", {}) or {}
        if ann_post:
            log.info(
                f"Annotation levels after dedup — L1: {ann_post.get('1', 0)}, L2: {ann_post.get('2', 0)}, L3: {ann_post.get('3', 0)}"
            )
        iso_post = summary.get("isomer_post", {}) or {}
        if iso_post:
            log.info(
                f"Isomers after dedup — metabolites with isomers: {iso_post.get('metabolites_with_isomers', 0)}, "
                f"total isomers: {iso_post.get('total_isomers', 0)}, max per metabolite: {iso_post.get('max_isomers_per_metabolite', 0)}"
            )
        # SIRIUS export counts
        sp = summary.get("sirius_pos_count", 0)
        sn = summary.get("sirius_neg_count", 0)
        total_l3_sirius = summary.get("sirius_l3_total", 0)
        total_l3_merged = ann_post.get("3", 0)
        l3_excluded = total_l3_merged - total_l3_sirius
        log.info(
            f"SIRIUS export — L3 in merged.csv: {total_l3_merged}, exported to SIRIUS (m/z 150-800): {total_l3_sirius}, excluded: {l3_excluded}"
        )
        log.info(f"  POS entries: {sp}, NEG entries: {sn}, total: {sp + sn}")
        # Per-file breakdown
        per_file = summary.get("per_file", [])
        if per_file:
            log.info("Per-file stats:")
            for rec in per_file:
                log.info(
                    f" - {rec.get('file')}: raw={rec.get('raw')}, "
                    f"MS/MS={rec.get('after_msms')}, S/N={rec.get('after_snr')}, pass_any={rec.get('after_pass')}, "
                    f"ann: L1={rec.get('ann1', 0)}, L2={rec.get('ann2', 0)}, L3={rec.get('ann3', 0)}"
                )
    # For wide: 'rows' are features
    log.info(f"Rows written: {summary['rows']}")
