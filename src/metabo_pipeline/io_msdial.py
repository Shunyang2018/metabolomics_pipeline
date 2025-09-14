"""
Thin compatibility wrapper for MS-DIAL merge routines.

The implementation has been moved to src/metabo_pipeline/merge.py
to keep the codebase modular and review-friendly.
CLI and external imports should continue using:
    from metabo_pipeline.io_msdial import merge_folder_to_wide_csv
"""

from pathlib import Path
from typing import Dict, Optional, Callable

from .merge import merge_folder_to_wide_csv  # re-export for CLI compatibility

__all__ = ["merge_folder_to_wide_csv"]
