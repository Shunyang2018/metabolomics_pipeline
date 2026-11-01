"""Inspect the MS-DIAL header contract before processing a dataset."""

import json
from dataclasses import asdict
from pathlib import Path
from typing import Annotated

import typer

from ..io_msdial import summarize_alignment_table


def validate(path: Annotated[Path, typer.Argument(exists=True, dir_okay=False)]):
    """Print sample metadata and feature counts for a CSV/TXT alignment export."""
    try:
        summary = summarize_alignment_table(path)
    except (OSError, ValueError) as exc:
        typer.echo(f"Invalid MS-DIAL export: {exc}", err=True)
        raise typer.Exit(code=2) from exc
    typer.echo(json.dumps(asdict(summary), default=str, indent=2))
