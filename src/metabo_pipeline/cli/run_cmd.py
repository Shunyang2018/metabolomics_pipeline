"""Run command - full pipeline orchestration."""

from __future__ import annotations

import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import typer

from .. import constants
from ..bioactivity import match_bioactives
from ..classify import classify_level12_with_classyfire
from ..logging import get_logger
from ..merge import merge_folder_to_wide_csv
from ..sirius_collect import collect_sirius_results

log = get_logger()


def _get_metabo_executable() -> str:
    """Resolve the current metabo executable path for subprocess calls."""
    candidates = []
    try:
        candidates.append(Path(sys.argv[0]).resolve())
    except Exception:
        pass
    which = shutil.which("metabo")
    if which:
        candidates.append(Path(which))
    for c in candidates:
        if c and c.exists():
            return str(c)
    return "metabo"


def _run_classify_task(csv_path: str, output_dir: str):
    """Run ClassyFire classification task."""
    log.info("  → Starting ClassyFire classification...")

    output_path = str(Path(csv_path).parent / f"{Path(csv_path).stem}_classyfire.csv")
    summary = classify_level12_with_classyfire(
        input_csv=csv_path,
        output_csv=output_path,
        sleep_sec=1.0,
        results_only=False,
    )

    # Enhanced logging
    unique_keys = summary.get("unique_keys", 0)
    hits = summary.get("hits", 0)
    miss = summary.get("miss", 0)
    cached = hits - summary.get("api_calls", 0)
    rows = summary.get("rows", 0)
    added_cols = summary.get("added_cols", 0)

    log.info(f"  ClassyFire: {unique_keys} unique InChIKeys processed")
    log.info(f"    • Cached: {cached}, API calls: {summary.get('api_calls', 0)}")
    log.info(f"    • Hits: {hits}, Misses: {miss}")
    log.info(f"    • Output: {rows} rows, +{added_cols} classification columns")
    return ("classify", True, summary)


def _run_sirius_task(output_dir: str):
    """Run SIRIUS analysis task."""
    log.info("  → Starting SIRIUS analysis...")
    metabo_exe = _get_metabo_executable()
    result = subprocess.run(
        [metabo_exe, "sirius", "--output-dir", output_dir],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return ("sirius", False, result.stderr or result.stdout)
    return ("sirius", True, None)


def _dispatch_task(task_info):
    """Dispatch a task based on task type."""
    task_name, csv_path, out_dir = task_info
    if task_name == "classify":
        return _run_classify_task(csv_path, out_dir)
    elif task_name == "sirius":
        return _run_sirius_task(out_dir)


def _run_parallel_tasks(tasks, task_names):
    """Run tasks in parallel and return results."""
    log.info(f"Running: {', '.join(task_names)}")

    results = {}
    with ThreadPoolExecutor(max_workers=2) as executor:
        future_to_task = {
            executor.submit(_dispatch_task, task): task[0] for task in tasks
        }

        for future in as_completed(future_to_task):
            result_name, success, error = future.result()
            results[result_name] = (success, error)
            if success:
                log.ok(f"  ✓ {result_name.title()} completed")
            else:
                log.error(f"  ✗ {result_name.title()} failed: {error}")

    return results


def _create_final_merge(merged_csv: str, output_dir: str):
    """Create final merged table with SIRIUS results."""
    log.info("\n[3/4] CREATING final merged table...")

    # Determine output CSV path
    j = Path(merged_csv)
    # Prefer classified output if present
    classified = j.with_name(j.stem + "_classyfire").with_suffix(j.suffix)
    if classified.exists():
        log.info(f"Using classified table: {classified.name}")
        j = classified

    output_csv_path = j.with_name(j.stem + "_final").with_suffix(j.suffix)

    summary = collect_sirius_results(
        pos_dir=Path(f"{output_dir}/sirius_pos"),
        neg_dir=Path(f"{output_dir}/sirius_neg"),
        output_csv=output_csv_path,
        join_with_merged=j,
        progress=None,
        extract_cache=Path(f"{output_dir}/sirius_identifications.csv"),
        force_rescan=False,
    )

    log.ok(f"✓ Final merge complete: {output_csv_path.name}")
    log.info(f"  • Total rows: {summary['rows']}")

    # Show SIRIUS statistics
    log.info(
        f"  • SIRIUS compounds: {summary.get('n_pos', 0)} POS, {summary.get('n_neg', 0)} NEG"
    )

    # Show annotation level distribution after SIRIUS
    ann_counts = summary.get("ann_after", {}) or summary.get(
        "annotation_level_counts", {}
    )
    if ann_counts:
        log.info("  • Final annotation levels:")
        log.info(f"      L1: {ann_counts.get('1', 0)} (MS-DIAL library hits)")
        log.info(f"      L2: {ann_counts.get('2', 0)} (Probable structure)")
        log.info(f"      L3: {ann_counts.get('3', 0)} (Putative structure - SIRIUS)")
        log.info(f"      L4: {ann_counts.get('4', 0)} (Molecular formula/class)")
        log.info(f"      L5: {ann_counts.get('5', 0)} (Exact mass only)")


def _run_bioactivity(output_dir: str, bioactivity_db: str | None) -> None:
    """Match the final merged table against the bioactivity database."""
    log.info("\n[4/4] MATCHING bioactivity database...")

    out_dir = Path(output_dir)
    merged_csv = None
    for name in ("merged_classyfire_final.csv", "merged_classyfire.csv", "merged.csv"):
        candidate = out_dir / name
        if candidate.exists():
            merged_csv = candidate
            break
    if merged_csv is None:
        log.warn(f"No merged CSV found in {out_dir}; skipping bioactivity match")
        return

    db_path = Path(bioactivity_db) if bioactivity_db else Path(constants.BIOACTIVITY_DB_PATH)
    if not db_path.exists():
        log.warn(f"Bioactivity database not found: {db_path}; skipping bioactivity match")
        return

    output_csv = out_dir / "bioactives.csv"
    summary = match_bioactives(merged_csv, db_path, output_csv)
    log.ok(f"✓ Bioactivity match complete: {output_csv.name}")
    log.info(f"  • Unique matched molecular skeletons: {summary['unique_matched_skeletons']}")
    log.info(f"  • Output rows (feature x sample x hit): {summary['output_rows']}")


def run(
    input_dir: str = typer.Option(
        None,
        help="Input directory with MS-DIAL files (default: INPUT_DIR from constants.py)",
    ),
    output_dir: str = typer.Option(
        None, help="Output directory (default: INPUT_DIR/outputs)"
    ),
    skip_merge: bool = typer.Option(
        False, help="Skip merge step (use existing merged.csv)"
    ),
    skip_classify: bool = typer.Option(False, help="Skip ClassyFire classification"),
    skip_sirius: bool = typer.Option(False, help="Skip SIRIUS analysis"),
    skip_bioactivity: bool = typer.Option(
        False, help="Skip bioactivity database matching"
    ),
    bioactivity_db: str = typer.Option(
        None, help="Path to bioactivity database CSV (default: constants.BIOACTIVITY_DB_PATH)"
    ),
):
    """Run the complete pipeline: merge → (classify + sirius in parallel) → final → bioactivity.

    This command orchestrates the full workflow:
    1. Merge MS-DIAL files (unless --skip-merge)
    2. Run ClassyFire and SIRIUS in parallel (independent operations)
    3. Create final merged table with all results
    4. Match features against a bioactivity database (unless --skip-bioactivity)

    ClassyFire and SIRIUS run concurrently to save time!
    """
    # Determine directories
    if input_dir is None:
        input_dir = constants.INPUT_DIR
    if output_dir is None:
        output_dir = str(Path(input_dir) / "outputs")

    merged_csv = str(Path(output_dir) / "merged.csv")

    log.info("=" * 60)
    log.info("METABOLOMICS PIPELINE - FULL RUN")
    log.info("=" * 60)

    # Step 1: Merge (unless skipped)
    if not skip_merge:
        log.info("\n[1/4] MERGING MS-DIAL files...")
        # Call merge directly
        summary = merge_folder_to_wide_csv(
            input_dir=Path(input_dir),
            output_csv=Path(merged_csv),
            recursive=False,
        )
        log.ok(f"✓ Merge complete: {merged_csv}")
        log.info(
            f"  • {summary['rows']} features across {summary.get('n_samples', 'N/A')} samples"
        )
        log.info(f"  • Files merged: {summary.get('n_files_merged', 'N/A')}")

        # Show annotation level distribution
        ann_post = summary.get("ann_post", {})
        if ann_post:
            log.info(
                f"  • Annotation levels: L1={ann_post.get('1', 0)}, L2={ann_post.get('2', 0)}, L3={ann_post.get('3', 0)}"
            )

        # Show SIRIUS export stats
        if summary.get("sirius_exported"):
            log.info(
                f"  • SIRIUS inputs: {summary.get('sirius_pos_compounds', 0)} POS, {summary.get('sirius_neg_compounds', 0)} NEG compounds"
            )
    else:
        log.info(f"\n[1/4] SKIPPING merge (using existing: {merged_csv})")
        if not Path(merged_csv).exists():
            log.error(f"Merged file not found: {merged_csv}")
            raise typer.Exit(code=1)

    # Step 2: Run ClassyFire and SIRIUS in parallel
    log.info("\n[2/4] RUNNING ClassyFire and SIRIUS in parallel...")

    tasks = []
    task_names = []

    if not skip_classify:
        tasks.append(("classify", merged_csv, output_dir))
        task_names.append("ClassyFire")

    if not skip_sirius:
        tasks.append(("sirius", None, output_dir))
        task_names.append("SIRIUS")

    if not tasks:
        log.warn("Both ClassyFire and SIRIUS skipped. Moving to final merge...")
    else:
        results = _run_parallel_tasks(tasks, task_names)

        # Check if any failed
        failed = [name for name, (success, _) in results.items() if not success]
        if failed:
            log.error(f"\nFailed tasks: {', '.join(failed)}")
            log.warn("Proceeding to final merge with available results...")

    # Step 3: Create final merged table
    _create_final_merge(merged_csv, output_dir)

    # Step 4: Match against bioactivity database
    if not skip_bioactivity:
        _run_bioactivity(output_dir, bioactivity_db)
    else:
        log.info("\n[4/4] SKIPPING bioactivity matching")

    log.info("\n" + "=" * 60)
    log.ok("✓ PIPELINE COMPLETE!")
    log.info("=" * 60)

    # Determine final output filename
    final_name = (
        "merged_classyfire_final.csv" if not skip_classify else "merged_final.csv"
    )
    log.info(f"📁 Final output: {output_dir}/{final_name}")
    if not skip_bioactivity:
        log.info(f"📁 Bioactivity matches: {output_dir}/bioactives.csv")
    log.info("")
    log.info("Next steps:")
    log.info("  • Review annotation level distribution above")
    log.info("  • Check ClassyFire taxonomies (cf_* columns)")
    log.info("  • Verify SIRIUS identifications (SIRIUS_* columns)")
    log.info("=" * 60)
