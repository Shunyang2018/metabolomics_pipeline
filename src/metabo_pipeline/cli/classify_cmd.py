"""Classify commands - ClassyFire classification and API check."""

from __future__ import annotations

from pathlib import Path

import typer

from .. import constants
from ..classify import classify_level12_with_classyfire, probe_classyfire
from ..logging import get_logger

log = get_logger()


def _progress_callback(i: int, n: int, key: str, cached: bool, status: str):
    """Log progress for ClassyFire classification."""
    src = "cache" if cached else "net"
    log.info(f"[{i}/{n}] {key} -> {src}/{status}")


def _adjust_output_for_join_merged(
    output_csv: str, input_csv: str, join_merged: bool
) -> tuple[bool, str]:
    """Adjust results_only and output_csv if join_merged is enabled.

    Returns (results_only, output_csv).
    """
    if not join_merged:
        return True, output_csv

    results_only = False
    if (
        str(output_csv).replace("/", "\\") == "outputs\\merged_classified.csv"
        and input_csv
    ):
        output_csv = input_csv

    return results_only, output_csv


def classify(
    input_csv: str = typer.Argument(
        None,
        help="Merged CSV to classify (default: INPUT_DIR/outputs/merged.csv from constants.py)",
    ),
    output_csv: str = typer.Option(
        None,
        help="Where to write classified CSV (default derives from input)",
    ),
    provider: str = typer.Option(
        "classyfire", help="Classification provider (only 'classyfire' supported)"
    ),
    base_url: str = typer.Option(
        "http://classyfire.wishartlab.com", help="ClassyFire base URL"
    ),
    sleep_sec: float = typer.Option(10.0, help="Sleep between API calls (seconds)"),
    results_only: bool = typer.Option(
        True, help="Write only identifier + classification columns"
    ),
    id_column: str = typer.Option(
        "feature_id", help="Identifier column to include (default: feature_id)"
    ),
    cache_path: str = typer.Option(
        None,
        help="Cache file for InChIKey->taxonomy (default: package cache with 1500+ compounds)",
    ),
    offline: bool = typer.Option(False, help="Use cache only; skip network calls"),
    timeout: float = typer.Option(15.0, help="Per-request timeout (seconds)"),
    join_merged: bool = typer.Option(
        False,
        help="Append cf_* columns into the merged table (writes back to input CSV if output path is default).",
    ),
    force_api: bool = typer.Option(
        False, help="Force API calls even if a classified CSV already exists"
    ),
):
    """Classify Level 1/2 rows via Wishart lab's ClassyFire API using unique InChIKeys and merge back.

    Prints how many unique keys were queried and how many columns were added; you can then decide which
    classification columns to keep in your table.
    """
    if provider.lower() != "classyfire":
        log.error("Only 'classyfire' provider is supported currently.")
        raise typer.Exit(code=2)

    # Use default paths from constants if not provided
    if input_csv is None:
        input_csv = str(Path(constants.INPUT_DIR) / "outputs" / "merged.csv")

    if output_csv is None:
        output_csv = str(
            Path(input_csv).parent / f"{Path(input_csv).stem}_classified.csv"
        )

    # Adjust output settings if joining with merged table
    results_only, output_csv = _adjust_output_for_join_merged(
        output_csv, input_csv, join_merged
    )

    # Run classification (cache mechanism is built-in via cache_path)
    summary = classify_level12_with_classyfire(
        input_csv,
        output_csv,
        sleep_sec=sleep_sec,
        base_url=base_url,
        results_only=results_only,
        id_column=id_column,
        cache_path=cache_path,
        progress=_progress_callback,
        offline=offline,
        timeout=timeout,
    )

    log.ok(
        f"Classified: {summary['unique_keys']} unique InChIKeys (hits: {summary['hits']}, miss: {summary['miss']}, skipped_offline: {summary.get('skipped_offline', 0)})"
    )
    log.info(
        f"Output: {output_csv} (rows: {summary['rows']}, classification columns added: {summary['added_cols']})"
    )


def classify_check(
    base_url: str = typer.Option(
        "http://classyfire.wishartlab.com", help="ClassyFire base URL to probe"
    ),
    timeout: float = typer.Option(10.0, help="Timeout per request (seconds)"),
):
    """Probe the ClassyFire API endpoint for availability and latency."""
    res = probe_classyfire(base_url, timeout=timeout)

    ok = bool(res.get("ok"))
    log.info(f"Base URL: {res.get('base_url')}")
    checks = res.get("checks", []) or []
    for c in checks:
        log.info(
            f" - {c.get('url')}: ok={c.get('ok')} status={c.get('status')} ms={c.get('ms')} error={c.get('error')}"
        )
    if not ok:
        log.error("ClassyFire API not reachable or returned non-200 status.")
        raise typer.Exit(code=1)
    log.ok("ClassyFire API reachable.")
