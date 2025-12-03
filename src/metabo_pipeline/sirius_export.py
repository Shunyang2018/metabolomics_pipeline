from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import List, Tuple

import pandas as pd

from .logging import get_logger
from .utils import parse_spectrum

log = get_logger()


def build_ms_entries(l3_df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """Build SIRIUS .ms entries in classic CLI format using feature_id.

    Format per entry:
      >compound <feature_id>
      >parentmass <Average Mz>
      >retentiontime <Average Rt(min)>
      >ionization <Adduct type>
      >feature_id <feature_id>

      >ms1
      <mz> <intensity>
      ...
      >collision
      <mz> <intensity>
      ...
    """
    pos_entries: List[str] = []
    neg_entries: List[str] = []
    for _, r in l3_df.iterrows():
        fid = r.get("feature_id")
        precursor = r.get("Average Mz", "")
        rt = r.get("Average Rt(min)", "")
        ion = r.get("Adduct type", "")
        ms1 = parse_spectrum(r.get("MS1 isotopic spectrum", ""))
        ms2 = parse_spectrum(r.get("MS/MS spectrum", ""))

        # Compose a readable compound name; mandatory in .ms
        comp_name = str(r.get("Metabolite name") or "").strip() or f"Unknown_{fid}"

        block: List[str] = []
        block.append(f">compound\t{comp_name}")
        block.append(f">parentmass\t{precursor}")
        block.append(f">retentiontime\t{rt}")
        block.append(f">ionization\t{ion}")
        # Add a comment with feature_id so downstream tools can recover it without triggering SIRIUS warnings
        block.append(f"#feature_id {fid}")
        block.append("")
        block.append(">ms1")
        for mz, inten in ms1:
            block.append(f"{mz} {inten}")
        block.append("")
        block.append(">collision")
        for mz, inten in ms2:
            block.append(f"{mz} {inten}")
        block.append("")

        entry = "\n".join(block)
        pol = (
            "POS"
            if "+" in str(ion or "")
            else ("NEG" if "-" in str(ion or "") else "UNK")
        )
        if pol == "POS":
            pos_entries.append(entry)
        elif pol == "NEG":
            neg_entries.append(entry)
    return pos_entries, neg_entries


def write_ms_files(
    pos_entries: List[str], neg_entries: List[str], out_dir: Path
) -> Tuple[int, int]:
    """Write POS/NEG .ms files and return the counts of entries written."""
    out_dir.mkdir(parents=True, exist_ok=True)
    pos_path = out_dir / "sirius_unknown_pos.ms"
    neg_path = out_dir / "sirius_unknown_neg.ms"
    if pos_entries:
        pos_path.write_text("\n".join(pos_entries), encoding="utf-8")
    if neg_entries:
        neg_path.write_text("\n".join(neg_entries), encoding="utf-8")
    return len(pos_entries), len(neg_entries)


def export_tsv_summaries(
    sirius_db: Path, output_dir: Path, sirius_exe: str, timeout: int = 300
) -> bool:
    """
    Export TSV summaries from a .sirius database file.

    Args:
        sirius_db: Path to .sirius database file
        output_dir: Directory to write TSV summaries
        sirius_exe: Path to SIRIUS executable
        timeout: Timeout in seconds (default: 300 = 5 minutes)

    Returns:
        True if successful, False otherwise
    """
    if not sirius_db.exists():
        log.warn(f"Database file not found: {sirius_db}")
        return False

    if sirius_db.stat().st_size < 1024:
        log.warn(f"Database file too small (< 1KB): {sirius_db}")
        return False

    log.info(f"Exporting TSV summaries from: {sirius_db.name}")
    log.info(f"  Output directory: {output_dir}")

    cmd = [
        sirius_exe,
        "--input",
        str(sirius_db),
        "write-summaries",
        "--output",
        str(output_dir),
    ]

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
            env=os.environ.copy(),
        )

        # Check for login errors in output (SIRIUS returns 0 even on login failure!)
        if "Login ERROR" in result.stdout or "not logged in" in result.stdout.lower():
            log.error("✗ Export failed: SIRIUS requires login")
            log.error("  Please login via SIRIUS GUI:")
            log.error("    1. Open SIRIUS application")
            log.error("    2. Go to Account → Login")
            log.error("    3. Enter credentials and wait for 'Logged in' confirmation")
            log.error("    4. Close SIRIUS and retry this command")
            return False

        # Check if output directory was created
        if output_dir.exists() and output_dir.is_dir():
            tsv_files = list(output_dir.glob("*.tsv"))
            if tsv_files:
                log.ok(
                    f"✓ TSV summaries exported to: {output_dir.name}/ ({len(tsv_files)} files)"
                )
                return True

        # If we get here, something went wrong
        if result.returncode != 0:
            log.error(f"✗ Export failed (exit code {result.returncode})")
            if result.stdout:
                error_lines = [
                    line
                    for line in result.stdout.split("\n")
                    if "ERROR" in line or "SEVERE" in line
                ]
                if error_lines:
                    log.error("  Error messages:")
                    for line in error_lines[:5]:
                        log.error(f"    {line}")
        else:
            log.error("✗ Export failed: No TSV files created")
        return False

    except subprocess.TimeoutExpired:
        log.error(f"✗ Export timed out (>{timeout}s)")
        return False
    except Exception as e:
        log.error(f"✗ Export failed: {e}")
        return False
