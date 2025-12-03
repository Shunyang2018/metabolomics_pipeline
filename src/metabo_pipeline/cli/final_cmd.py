"""Final command - create final merged table."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from .. import constants
from ..logging import get_logger
from ..sirius_collect import collect_sirius_results

log = get_logger()


def final(
    output_dir: str = typer.Option(
        None, "--output-dir", "-o", help="Output directory (default: INPUT_DIR/outputs)"
    ),
    rescan: bool = typer.Option(
        False, help="Force rescan of SIRIUS folders and refresh the extract cache"
    ),
):
    """Create the final merged table combining MS-DIAL, ClassyFire, and SIRIUS results.

    This is the last step in the pipeline. It merges SIRIUS identifications and CANOPUS
    classifications into the merged table (optionally with ClassyFire taxonomy if available).

    Uses automatic defaults:
    - Input: INPUT_DIR/outputs/sirius_pos, sirius_neg, merged_classyfire.csv
    - Output: INPUT_DIR/outputs/merged_classyfire_final.csv

    Just run: metabo final
    """
    # Determine output directory
    if output_dir is None:
        out_dir = Path(constants.INPUT_DIR) / "outputs"
    else:
        out_dir = Path(output_dir)

    # Derive all paths from output directory
    pos_dir = str(out_dir / "sirius_pos")
    neg_dir = str(out_dir / "sirius_neg")
    join_merged = str(out_dir / "merged.csv")
    extract_csv = str(out_dir / "sirius_identifications.csv")
    output_csv = str(out_dir / "merged_final.csv")  # Will be updated below

    j: Optional[Path] = None
    if join_merged and str(join_merged).strip():
        j = Path(join_merged)
        # Prefer classified output if present: <stem>_classyfire.<suffix>
        classified = j.with_name(j.stem + "_classyfire").with_suffix(j.suffix)
        if classified.exists():
            log.info(
                f"Detected classified table: {classified.name}; using it for merge"
            )
            j = classified
            output_csv = str(out_dir / "merged_classyfire_final.csv")
        else:
            output_csv = str(out_dir / "merged_final.csv")

    # Progress bar (rich) if available; otherwise fallback
    pd_path, nd_path = Path(pos_dir), Path(neg_dir)
    # If an extract CSV exists and no rescan requested, skip progress UI and reuse it
    reuse_extract = False
    if extract_csv and not rescan and Path(extract_csv).exists():
        reuse_extract = True

    if reuse_extract:
        summary = collect_sirius_results(
            pd_path,
            nd_path,
            Path(output_csv),
            join_with_merged=j,
            progress=None,
            extract_cache=Path(extract_csv),
            force_rescan=False,
        )
        log.info(f"Reused SIRIUS identifications from {extract_csv}")
    else:
        from rich.progress import (  # type: ignore
            BarColumn,
            Progress,
            TextColumn,
            TimeRemainingColumn,
        )

        total_dirs = 0
        if pd_path.exists():
            total_dirs += sum(1 for _ in pd_path.iterdir() if _.is_dir())
        if nd_path.exists():
            total_dirs += sum(1 for _ in nd_path.iterdir() if _.is_dir())
        prog = Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TimeRemainingColumn(),
        )
        task_id = prog.add_task("Collecting SIRIUS IDs", total=total_dirs or 1)

        def _cb(done: int, total: int, path: str):
            prog.update(task_id, completed=done)

        with prog:
            summary = collect_sirius_results(
                pd_path,
                nd_path,
                Path(output_csv),
                join_with_merged=j,
                progress=_cb,
                extract_cache=(Path(extract_csv) if extract_csv else None),
                force_rescan=bool(rescan),
            )
    log.ok(f"Wrote identifications/merge: {output_csv} (rows: {summary['rows']})")
    log.info(
        f"POS: total={summary.get('pos_compounds', 0)}, identified={summary.get('pos_identified', 0)}, "
        f"miss={summary.get('miss_pos', 0)}; "
        f"NEG: total={summary.get('neg_compounds', 0)}, identified={summary.get('neg_identified', 0)}, miss={summary.get('miss_neg', 0)}"
    )

    # Annotation level distribution after merge
    ann = summary.get("ann_after") or {}
    if ann:
        log.info("Annotation levels after SIRIUS merge:")
        log.info(f"  • L1: {ann.get('1', 0)} (MS-DIAL library hits)")
        log.info(f"  • L2: {ann.get('2', 0)} (Probable structure)")
        log.info(f"  • L3: {ann.get('3', 0)} (Putative structure - SIRIUS)")
        log.info(f"  • L4: {ann.get('4', 0)} (Molecular formula/class)")
        log.info(f"  • L5: {ann.get('5', 0)} (Exact mass only)")
