from __future__ import annotations

import json
import os
import platform
from pathlib import Path
from typing import Optional

import typer

from .io_msdial import summarize_alignment_table
from .logging import get_logger
from .io_msdial import merge_folder_to_long_csv, merge_folder_to_wide_csv

app = typer.Typer(add_completion=False, help="Metabolomics pipeline for MS-DIAL outputs.")
log = get_logger()


@app.command()
def init(
    output: str = typer.Option(".", help="Where to scaffold config/templates"),
):
    """Scaffold a minimal config and folders."""
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    (out / "configs").mkdir(exist_ok=True)
    (out / "templates").mkdir(exist_ok=True)

    cfg = out / "configs" / "config.yml"
    if not cfg.exists():
        cfg.write_text(
            """
# Minimal config example
inputs:
  files: []  # paths to MS-DIAL alignment CSV/TXT
processing:
  missingness_threshold: 0.2
  rsd_qc_threshold: 0.3
  normalization: median  # options: median, tss, pqn, is
  imputation: half_min
report:
  enabled: true
  output_dir: report
            """.strip()
        )
        log.ok(f"Created {cfg}")
    else:
        log.warn(f"Exists: {cfg}")


@app.command()
def validate(
    input: Optional[str] = typer.Option(None, "--input", "-i", help="MS-DIAL alignment CSV/TXT"),
):
    """Validate MS-DIAL file(s), summarize samples/features, and embedded metadata rows."""
    if not input:
        log.error("Please provide --input path to an MS-DIAL alignment table.")
        raise typer.Exit(code=2)

    path = Path(input)
    if not path.exists():
        log.error(f"File not found: {path}")
        raise typer.Exit(code=2)

    try:
        summary = summarize_alignment_table(path)
    except Exception as e:
        log.error(f"Failed to parse {path}: {e}")
        raise typer.Exit(code=1)

    log.ok(f"Parsed: {summary.path}")
    log.info(f"Features: {summary.n_features}")
    log.info(f"Samples: {len(summary.samples)}")

    # Show a quick JSON block with first few samples and their metadata
    preview = {
        "samples_preview": summary.samples[:5],
        "classes_preview": {k: summary.metadata.classes.get(k) for k in summary.samples[:5]},
        "file_types_preview": {k: summary.metadata.file_types.get(k) for k in summary.samples[:5]},
        "injection_order_preview": {k: summary.metadata.injection_order.get(k) for k in summary.samples[:5]},
        "batch_id_preview": {k: summary.metadata.batch_id.get(k) for k in summary.samples[:5]},
    }
    print(json.dumps(preview, indent=2))


@app.command()
def run(
    config: Optional[str] = typer.Option(None, "--config", "-c", help="Path to config.yml"),
):
    """Run the pipeline (placeholder)."""
    log.warn("Pipeline execution is not implemented yet. See PLAN.md milestones M3–M7.")
    if config:
        log.info(f"Config: {config}")


@app.command()
def diag():
    """Print environment diagnostics."""
    info = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "cwd": str(Path.cwd()),
        "env": {k: v for k, v in os.environ.items() if k in ("CONDA_DEFAULT_ENV", "VIRTUAL_ENV")},
    }
    print(json.dumps(info, indent=2))


@app.command()
def merge(
    input_dir: str = typer.Argument(..., help="Folder containing MS-DIAL CSV/TXT files"),
    output_csv: str = typer.Option("outputs/merged.csv", help="Path to write merged CSV"),
    recursive: bool = typer.Option(False, help="Recurse into subfolders"),
    format: str = typer.Option("wide", help="Output format: wide or long", case_sensitive=False),
    engine: str = typer.Option("pandas", help="Engine for parsing (long mode only): pandas or csv", case_sensitive=False),
):
    """Merge HILIC/C18/Lipidomics files. Default is wide format (one feature per row)."""
    in_dir = Path(input_dir)
    if not in_dir.exists() or not in_dir.is_dir():
        log.error(f"Input directory not found: {in_dir}")
        raise typer.Exit(code=2)
    out = Path(output_csv)
    try:
        if format.lower() == "wide":
            summary = merge_folder_to_wide_csv(in_dir, out, recursive=recursive)
        else:
            summary = merge_folder_to_long_csv(in_dir, out, recursive=recursive, engine=engine.lower())
    except Exception as e:
        log.error(f"Failed to merge files from {in_dir}: {e}")
        raise typer.Exit(code=1)

    log.ok(f"Merged {summary['files']} files → {out}")
    # For wide: 'rows' are features; for long: 'rows' are long rows
    log.info(f"Rows written: {summary['rows']}")


if __name__ == "__main__":
    app()
