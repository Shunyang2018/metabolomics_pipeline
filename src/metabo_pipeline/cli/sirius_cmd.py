"""SIRIUS command - run SIRIUS analysis."""

from __future__ import annotations

import os
import signal
import subprocess
import sys
from pathlib import Path

import typer
import yaml

from .. import constants
from ..logging import get_logger
from ..sirius_export import export_tsv_summaries
from ..sirius_utils import guess_sirius_executable

log = get_logger()

# Global variable to track running SIRIUS processes
_sirius_processes = []


def _cleanup_sirius_processes(signum=None, frame=None):
    """Terminate all running SIRIUS processes."""
    global _sirius_processes
    for proc in _sirius_processes:
        if proc and proc.poll() is None:  # Process still running
            log.warn(f"Terminating SIRIUS process (PID: {proc.pid})...")
            try:
                proc.terminate()
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                log.warn(f"Force killing SIRIUS process (PID: {proc.pid})...")
                proc.kill()
            except Exception as e:
                log.warn(f"Error terminating SIRIUS process: {e}")
    _sirius_processes.clear()
    if signum is not None:
        sys.exit(1)


# Register signal handlers for graceful shutdown
signal.signal(signal.SIGINT, _cleanup_sirius_processes)
signal.signal(signal.SIGTERM, _cleanup_sirius_processes)


def sirius(
    output_dir: str = typer.Option(
        None,
        "--output-dir",
        "-o",
        help="Output directory containing SIRIUS .ms files (default: INPUT_DIR/outputs)",
    ),
    sirius_exe: str = typer.Option(
        None,
        help="Path to SIRIUS executable (auto-detected if not specified; override in constants.py: SIRIUS_EXECUTABLE)",
    ),
    cores: int = typer.Option(
        None, help="CPU cores for SIRIUS (default from constants.py: SIRIUS_CORES)"
    ),
    mzdev_ppm: float = typer.Option(
        None,
        help="MS2 mass deviation in ppm (default from constants.py: SIRIUS_MZDEV_PPM)",
    ),
    config: str = typer.Option(
        "configs/sirius.yml",
        help="Optional YAML with SIRIUS settings (adducts, DBs, tasks). CLI flags override YAML.",
    ),
    dry_run: bool = typer.Option(
        False,
        "--dry-run",
        "--test",
        help="Test mode: create fake SIRIUS output directories without running SIRIUS",
    ),
    test_real_case: bool = typer.Option(
        False,
        "--test-real-case",
        help="Run SIRIUS on small test .ms files (4 compounds, ~2-3 min) for quick end-to-end validation",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Force re-run even if SIRIUS output files already exist",
    ),
):
    """Run SIRIUS on the exported POS/NEG .ms files with identification and classification tasks.

    Runs: formulas, fingerprints, structures, canopus. Uses Orbitrap profile and your adduct lists.

    Quick test options:
    - --dry-run: Create fake outputs without running SIRIUS (instant)
    - --test-real-case: Run SIRIUS on 4 small test compounds (~2-3 min, real SIRIUS execution)
    """

    # Handle test-real-case mode
    if test_real_case:
        log.info("=" * 60)
        log.info("TEST MODE: Running SIRIUS on small test dataset")
        log.info("=" * 60)
        log.info("Test dataset: 4 compounds (2 POS, 2 NEG)")
        log.info("Expected runtime: ~2-3 minutes")
        log.info("")

        # Use test_data directory
        from pathlib import Path as P

        # Navigate from cli/sirius_cmd.py up to project root: cli -> metabo_pipeline -> src -> root
        test_dir = P(__file__).parent.parent.parent.parent / "test_data"
        if not test_dir.exists():
            log.error(f"Test data directory not found: {test_dir}")
            log.error("Run this command from the project root directory")
            raise typer.Exit(code=1)

        # Override paths to use test data
        if output_dir is None:
            output_dir = "test_output"

        # Create output directory (SIRIUS requires it to exist)
        P(output_dir).mkdir(parents=True, exist_ok=True)

        pos_ms = str(test_dir / "test_pos.ms")
        neg_ms = str(test_dir / "test_neg.ms")
        # SIRIUS 6 expects .sirius file extension (it's a database file, not a directory)
        out_pos = f"{output_dir}/test_sirius_pos.sirius"
        out_neg = f"{output_dir}/test_sirius_neg.sirius"

        # Verify test files exist
        if not P(pos_ms).exists() or not P(neg_ms).exists():
            log.error(f"Test .ms files not found: {test_dir}")
            log.error("Expected files: test_pos.ms, test_neg.ms")
            raise typer.Exit(code=1)

        log.info(f"Input: {test_dir}/test_*.ms")
        log.info(f"Output: {output_dir}/")
        log.info("")
    else:
        # Determine output directory for normal mode
        pass  # Will be handled below

    # Determine output directory and file paths
    if not test_real_case:
        if output_dir is None:
            out_dir = Path(constants.INPUT_DIR) / "outputs"
        else:
            out_dir = Path(output_dir)

        # Fixed file/folder names within output directory
        pos_ms = str(out_dir / "sirius_unknown_pos.ms")
        neg_ms = str(out_dir / "sirius_unknown_neg.ms")
        # SIRIUS 6 expects .sirius file extension (it's a database file, not a directory)
        out_pos = str(out_dir / "sirius_pos.sirius")
        out_neg = str(out_dir / "sirius_neg.sirius")
    # else: test_real_case paths already set above

    # Load YAML config if present
    cfg_yaml = {}
    cfg_path = Path(config)
    if cfg_path.exists():
        try:
            cfg_yaml = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
            log.ok(f"Loaded SIRIUS config: {cfg_path}")
        except Exception as e:
            log.warn(f"Failed to parse {cfg_path}: {e}")

    # Priority: CLI > YAML > constants.py
    if sirius_exe is None:
        sirius_exe = (
            cfg_yaml.get("sirius", {}).get("executable") or constants.SIRIUS_EXECUTABLE
        )
    if cores is None:
        cores = cfg_yaml.get("sirius", {}).get("cores") or constants.SIRIUS_CORES
    if mzdev_ppm is None:
        mzdev_ppm = (
            cfg_yaml.get("sirius", {}).get("mzdev_ppm") or constants.SIRIUS_MZDEV_PPM
        )

    sirius_exe, exe_note = guess_sirius_executable(sirius_exe)
    note = " (source: %s)" % exe_note if exe_note else ""
    log.info(f"Using SIRIUS executable: {sirius_exe}{note}")
    if not Path(sirius_exe).exists():
        log.warn(
            "SIRIUS executable not found on disk; update --sirius-exe or configs/sirius.yml if this fails at runtime."
        )

    # Check if SIRIUS is logged in (required for SIRIUS 6+)
    log.info("Checking SIRIUS login status...")
    try:
        login_check = subprocess.run(
            [sirius_exe, "login", "--show"], capture_output=True, text=True, timeout=30
        )
        if (
            "not logged in" in login_check.stdout.lower()
            or "please login" in login_check.stderr.lower()
        ):
            log.error("SIRIUS requires authentication.")
            log.error("Please login via SIRIUS GUI:")
            log.error("  1. Open SIRIUS application")
            log.error("  2. Login with your account")
            log.error("  3. Wait for 'Logged in' confirmation")
            log.error("  4. Close SIRIUS and retry this command")
            log.error("")
            log.error(
                "Get a free account at: https://bio.informatik.uni-jena.de/sirius-register/"
            )
            raise typer.Exit(code=2)
        log.ok("SIRIUS login verified")
    except subprocess.TimeoutExpired:
        log.warn("Login check timed out, proceeding anyway...")

    # SIRIUS 6: Most config options are not available as CLI flags
    # Use SIRIUS defaults which are already optimized for metabolomics
    # SIRIUS 6: Use formulas only by default (others require project-space workflow)
    tasks_cfg = cfg_yaml.get("tasks", ["formulas"]) or ["formulas"]

    def run_sirius(
        ms_path: Path, out_dir: Path, mode: str, force_rerun: bool = False
    ) -> int:
        if not ms_path.exists():
            log.warn(f"Missing .ms file: {ms_path}")
            return 0

        # Resolve absolute paths to avoid SIRIUS writing under Program Files
        ms_path = ms_path.resolve()
        out_dir = out_dir.resolve()

        # Check if output already exists (skip if present, unless --force)
        if not force_rerun and out_dir.exists():
            # Check if it's a valid SIRIUS output (has some content)
            if out_dir.stat().st_size > 1024:  # At least 1KB
                log.info(f"SIRIUS output already exists ({mode}): {out_dir}")
                log.info(f"  Size: {out_dir.stat().st_size / (1024**2):.1f} MB")
                log.info("  Skipping SIRIUS run. Use --force to re-run.")
                return 1  # Return 1 to indicate success (output exists)

        # SIRIUS 6: Command-line options are very limited
        # Use SIRIUS's built-in intelligent defaults (they are already optimized!)
        # Adduct inference is automatic based on ionization mode in the .ms file

        # Simple command structure: sirius [global-options] [tasks...]
        cmd = [
            sirius_exe,
            "-i",
            str(ms_path),
            "-o",
            str(out_dir),
            "--cores",
            str(cores),
            *tasks_cfg,  # Tasks to run (e.g., "formulas")
        ]
        log.info(f"Running SIRIUS ({mode}) → {out_dir}")
        log.info(f"Command: {' '.join(cmd[:5])}... (+ {len(cmd) - 5} more args)")

        try:
            # Use Popen to track the process for proper cleanup
            global _sirius_processes

            # Ensure SIRIUS credentials are passed to subprocess

            env = os.environ.copy()

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # Merge stderr into stdout for unified output
                text=True,
                bufsize=1,  # Line-buffered
                env=env,  # Pass environment variables including SIRIUS_USER/PASSWORD
            )
            _sirius_processes.append(process)

            # Stream output in real-time
            stdout_lines = []
            log.info(f"SIRIUS {mode} started (PID: {process.pid}). Streaming output...")
            for line in process.stdout:
                line = line.rstrip()
                if line:
                    # Show relevant progress lines (filter out too much noise)
                    if any(
                        keyword in line
                        for keyword in [
                            "Processing",
                            "Analyzing",
                            "Computing",
                            "Writing",
                            "Finished",
                            "Done",
                            "compounds",
                            "%",
                            "progress",
                            "WARNING",
                            "ERROR",
                            "SEVERE",
                        ]
                    ):
                        log.info(f"  SIRIUS [{mode}]: {line}")
                    stdout_lines.append(line)

            # Wait for process to finish
            process.wait()
            stdout = "\n".join(stdout_lines)

            # Remove from tracking list after completion
            _sirius_processes.remove(process)

            # Check exit code
            if process.returncode != 0:
                log.error(f"SIRIUS ({mode}) failed with exit code {process.returncode}")
                if stdout:
                    log.error("Last 20 lines of output:")
                    for line in stdout.split("\n")[-20:]:
                        if line.strip():
                            log.error(f"  {line}")
                return 0

            # SIRIUS 6: Export TSV summaries from the .sirius database
            # This creates a directory with TSV files for easy inspection
            tsv_dir = out_dir.parent / out_dir.stem.replace(".sirius", "")
            export_success = export_tsv_summaries(
                out_dir, tsv_dir, sirius_exe, timeout=300
            )
            if not export_success:
                log.warn(
                    "  TSV summary export failed - you can still inspect the .sirius database directly"
                )
                log.warn(
                    "  Or run 'metabo sirius-collect' after logging in via SIRIUS GUI"
                )

            # Store result for error checking
            result = type(
                "Result",
                (),
                {"stdout": stdout, "stderr": "", "returncode": process.returncode},
            )()

            # Check for Java errors in output even if exit code is 0
            # Be more specific to avoid false positives from normal log messages
            stdout_text = result.stdout

            # Real error patterns (case-sensitive to avoid false positives)
            has_java_error = False
            error_lines = []
            all_lines = stdout_text.split("\n")

            for i, line in enumerate(all_lines):
                line_lower = line.lower()
                # Look for actual error indicators, not just the word "error"
                if any(
                    [
                        line.startswith("Exception in thread"),
                        line.startswith("Error:"),
                        "java.lang.Exception" in line,
                        "java.lang.Error" in line,
                        "Fatal error" in line,
                        "Could not create the Java Virtual Machine" in line,
                        line_lower.startswith("severe:"),
                        "failed with code" in line_lower and "error" in line_lower,
                    ]
                ):
                    has_java_error = True
                    # Include context: 2 lines before and 5 lines after
                    start = max(0, i - 2)
                    end = min(len(all_lines), i + 6)
                    context = all_lines[start:end]
                    error_lines.extend(context)
                    error_lines.append("---")  # Separator between errors

            if has_java_error:
                log.error(f"SIRIUS ({mode}) completed but Java errors detected:")
                for line in error_lines[:50]:  # Show first 50 lines with context
                    if line.strip():
                        log.error(f"  {line}")
                return 0

            # Validate output: SIRIUS 6 creates a .sirius project file
            sirius_project = out_dir.with_suffix(".sirius")
            if sirius_project.exists():
                size_mb = sirius_project.stat().st_size / (1024 * 1024)
                log.ok(
                    f"SIRIUS ({mode}) completed successfully: {sirius_project.name} ({size_mb:.1f} MB)"
                )
            else:
                log.error(
                    f"SIRIUS ({mode}) finished but no .sirius project file created: {sirius_project}"
                )
                # Check if directory was created (might be old output)
                if out_dir.exists() and out_dir.is_dir():
                    dirs = [d for d in out_dir.iterdir() if d.is_dir()]
                    log.warn(
                        f"Found directory output (SIRIUS 5 style?): {len(dirs)} items"
                    )
                return 0

            return 1

        except KeyboardInterrupt:
            log.warn(f"SIRIUS ({mode}) interrupted by user")
            _cleanup_sirius_processes()
            raise
        except Exception as e:
            log.error(f"SIRIUS ({mode}) failed: {e}")
            # Make sure to clean up the process if it's still in the list
            if process in _sirius_processes:
                _sirius_processes.remove(process)
            return 0

    ran = 0
    ran += run_sirius(Path(neg_ms), Path(out_neg), "NEG", force_rerun=force)
    ran += run_sirius(Path(pos_ms), Path(out_pos), "POS", force_rerun=force)
    if ran == 0:
        raise typer.Exit(code=1)
    log.ok("SIRIUS runs finished. Check outputs for identification and classification.")
