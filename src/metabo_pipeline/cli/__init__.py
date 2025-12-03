"""Command-line interface for metabolomics pipeline."""

from __future__ import annotations

import typer

from .classify_cmd import classify, classify_check
from .final_cmd import final
from .merge_cmd import merge
from .run_cmd import run
from .sirius_cmd import sirius
from .sirius_collect_cmd import sirius_collect

# Create the main CLI app
app = typer.Typer(
    add_completion=False, help="Metabolomics pipeline for MS-DIAL outputs."
)

# Register commands with the app
app.command()(run)
app.command()(merge)
app.command()(classify)
app.command("classify-check")(classify_check)
app.command()(sirius)
app.command("sirius-collect")(sirius_collect)
app.command()(final)

__all__ = ["app"]
