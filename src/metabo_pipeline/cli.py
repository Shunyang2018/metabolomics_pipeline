from __future__ import annotations

import json
import os
import platform
import shutil
from pathlib import Path
from typing import List, Optional

import typer

from .io_msdial import summarize_alignment_table
from .logging import get_logger
from .io_msdial import merge_folder_to_wide_csv
import subprocess
import yaml
from .classify import classify_level12_with_classyfire
from .sirius_collect import collect_sirius_results
from .sirius_utils import guess_sirius_executable

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
def clean(
    outputs_dir: str = typer.Option("outputs", help="Root folder for generated artifacts."),
    prune_tables: bool = typer.Option(False, help="Also delete merged/classified CSV tables."),
    prune_all: bool = typer.Option(False, help="Delete the entire outputs directory."),
    dry_run: bool = typer.Option(False, help="Show what would be deleted without removing anything."),
    extra: List[str] = typer.Option([], "--extra", help="Additional files or directories to delete."),
):
    """Remove generated SIRIUS exports, caches, and optional merged tables."""
    root = Path(outputs_dir)
    collected: List[Path] = []
    seen: set[str] = set()

    def queue(path: Path) -> None:
        try:
            resolved = path.resolve()
        except Exception:
            resolved = path
        key = str(resolved)
        if key in seen or not path.exists():
            return
        seen.add(key)
        collected.append(path)

    if prune_all:
        if root.exists():
            queue(root)
    else:
        if root.exists():
            patterns = [
                "sirius_unknown_*.ms",
                "sirius_pos*",
                "sirius_neg*",
                "sirius_identifications.csv",
                "classyfire_cache.*",
            ]
            for pattern in patterns:
                for cand in root.glob(pattern):
                    queue(cand)
            if prune_tables:
                for pattern in ("merged*.csv", "*_classyfire*.csv", "*_ultimate*.csv"):
                    for cand in root.glob(pattern):
                        queue(cand)
        else:
            log.info(f"Outputs directory not found: {root}")

    for item in extra:
        queue(Path(item))

    if not collected:
        log.info("Nothing to clean.")
        return

    collected.sort(key=lambda p: str(p).lower())
    action = "Would remove" if dry_run else "Removing"
    for target in collected:
        log.info(f"{action}: {target}")

    if dry_run:
        return

    removed = 0
    for target in collected:
        try:
            if target.is_dir() and not target.is_symlink():
                shutil.rmtree(target)
            else:
                target.unlink(missing_ok=True)
            removed += 1
        except Exception as exc:
            log.warn(f"Failed to delete {target}: {exc}")
    log.ok(f"Removed {removed} item(s).")


@app.command()
def merge(
    input_dir: str = typer.Argument(..., help="Folder containing MS-DIAL CSV/TXT files"),
    output_csv: str = typer.Option(
        "outputs/merged.csv",
        "--output",
        "--output-csv",
        "-o",
        help="Path to write merged CSV",
    ),
    recursive: bool = typer.Option(False, help="Recurse into subfolders"),
):
    """Merge HILIC/C18/Lipidomics files (wide format, one feature per row)."""
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
        summary = merge_folder_to_wide_csv(in_dir, out, recursive=recursive, progress=_progress)
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
        total_l3 = summary.get("sirius_l3_total", 0)
        log.info(f"SIRIUS export — L3 total (merged): {total_l3}, POS entries: {sp}, NEG entries: {sn}, unique: {sp + sn}")
        # Per-file breakdown
        per_file = summary.get("per_file", [])
        if per_file:
            log.info("Per-file stats:")
            for rec in per_file:
                log.info(
                    f" - {rec.get('file')}: raw={rec.get('raw')}, "
                    f"MS/MS={rec.get('after_msms')}, S/N={rec.get('after_snr')}, pass_any={rec.get('after_pass')}, "
                    f"ann: L1={rec.get('ann1',0)}, L2={rec.get('ann2',0)}, L3={rec.get('ann3',0)}"
                )
    # For wide: 'rows' are features
    log.info(f"Rows written: {summary['rows']}")


@app.command()
def sirius(
    sirius_exe: str = typer.Option(
        "sirius",
        help="Path to SIRIUS executable (e.g., 'sirius' or 'C:\\Program Files\\sirius\\sirius.exe')",
    ),
    pos_ms: str = typer.Option("outputs/sirius_unknown_pos.ms", help="Path to POS .ms file"),
    neg_ms: str = typer.Option("outputs/sirius_unknown_neg.ms", help="Path to NEG .ms file"),
    out_pos: str = typer.Option("outputs/sirius_pos", help="Output folder for POS results"),
    out_neg: str = typer.Option("outputs/sirius_neg", help="Output folder for NEG results"),
    cores: int = typer.Option(8, help="CPU cores for SIRIUS"),
    mzdev_ppm: float = typer.Option(10.0, help="MS2 mass deviation (ppm)"),
    config: str = typer.Option(
        "configs/sirius.yml",
        help="Optional YAML with SIRIUS settings (adducts, DBs, tasks). CLI flags override YAML."
    ),
):
    """Run SIRIUS on the exported POS/NEG .ms files with identification and classification tasks.

    Runs: formula, fingerprint, structure, canopus. Uses Orbitrap profile and your adduct lists.
    """

    # Load YAML config if present
    cfg_yaml = {}
    cfg_path = Path(config)
    if cfg_path.exists():
        try:
            cfg_yaml = yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}
            log.ok(f"Loaded SIRIUS config: {cfg_path}")
        except Exception as e:
            log.warn(f"Failed to parse {cfg_path}: {e}")

    # If YAML provides defaults, apply them unless CLI overrides
    sirius_exe = cfg_yaml.get("sirius", {}).get("executable", sirius_exe)
    if cores == 8:
        cores = int(cfg_yaml.get("sirius", {}).get("cores", cores))
    if mzdev_ppm == 10.0:
        mzdev_ppm = float(cfg_yaml.get("sirius", {}).get("mzdev_ppm", mzdev_ppm))

    sirius_exe, exe_note = guess_sirius_executable(sirius_exe)
    note = ' (source: %s)' % exe_note if exe_note else ''
    log.info(f"Using SIRIUS executable: {sirius_exe}{note}")
    if not Path(sirius_exe).exists():
        log.warn('SIRIUS executable not found on disk; update --sirius-exe or configs/sirius.yml if this fails at runtime.')

    # Common config lines (exact strings)
    common_cfg = cfg_yaml.get("common_config", []) or [
        "--IsotopeSettings.filter=true",
        "--FormulaSearchDB=",
        "--Timeout.secondsPerTree=0",
        "--FormulaSettings.enforced=HCNOPS",
        "--Timeout.secondsPerInstance=0",
        "--UseHeuristic.mzToUseHeuristicOnly=650",
        "--AlgorithmProfile=orbitrap",
        "--IsotopeMs2Settings=IGNORE",
        f"--MS2MassDeviation.allowedMassDeviation={mzdev_ppm}ppm",
        "--NumberOfCandidatesPerIon=1",
        "--UseHeuristic.mzToUseHeuristic=300",
        "--FormulaSettings.detectable=,",
        "--NumberOfCandidates=10",
        "--AdductSettings.enforced=,",
        "--FormulaResultThreshold=true",
        "--InjectElGordoCompounds=true",
        "--StructureSearchDB=BIO,METACYC,CHEBI,COCONUT,ECOCYCMINE,GNPS,HMDB,HSDB,KEGG,KEGGMINE,KNAPSACK,MACONDA,MESH,NORMAN,UNDP,PLANTCYC,PUBCHEM,PUBMED,YMDB,YMDBMINE,ZINCBIO",
        "--RecomputeResults",
        "true",
    ]

    adducts_cfg = cfg_yaml.get("adducts", {})
    tasks_cfg = cfg_yaml.get("tasks", ["formula", "fingerprint", "structure", "canopus"]) or ["formula", "fingerprint", "structure", "canopus"]

    def run_sirius(ms_path: Path, out_dir: Path, mode: str) -> int:
        if not ms_path.exists():
            log.warn(f"Missing .ms file: {ms_path}")
            return 0

        # Resolve absolute paths to avoid SIRIUS writing under Program Files
        ms_path = ms_path.resolve()
        out_dir = out_dir.resolve()

        cfg = ["config", *common_cfg]

        if mode.upper() == "NEG":
            det = adducts_cfg.get("neg_detectable", "[[M - H]-, [M - H2O - H]-]")
            fb = adducts_cfg.get("neg_fallback", "[[M]-, [M + CH2O2 - H]-, [M + C2HF3O2 - H]-, [M - H + Na - H]-, [M - H + K - H]-, [M + C2H4O2 - H]-, [M + C2H3N - H]-, [M - H]-, [M - H2O - H]-]")
            cfg.extend([f"--AdductSettings.detectable={det}", f"--AdductSettings.fallback={fb}"])
        else:
            det = adducts_cfg.get("pos_detectable", "[[M + H]+, [M + Na]+]")
            fb = adducts_cfg.get("pos_fallback", "[[M]+, [M + NH4]+, [M + Na]+, [M + K]+, [M + CH3OH + H]+, [M + ACN + H]+, [M + 2Na - H]+]")
            cfg.extend([f"--AdductSettings.detectable={det}", f"--AdductSettings.fallback={fb}"])

        cmd = [
            sirius_exe,
            "-i",
            str(ms_path),
            "-o",
            str(out_dir),
            "--cores",
            str(cores),
            *cfg,
            *tasks_cfg,
        ]
        log.info(f"Running SIRIUS ({mode}) → {out_dir}")
        try:
            subprocess.run(cmd, check=True)
            return 1
        except subprocess.CalledProcessError as e:
            log.error(f"SIRIUS ({mode}) failed: {e}")
            return 0

    ran = 0
    ran += run_sirius(Path(neg_ms), Path(out_neg), "NEG")
    ran += run_sirius(Path(pos_ms), Path(out_pos), "POS")
    if ran == 0:
        raise typer.Exit(code=1)
    log.ok("SIRIUS runs finished. Check outputs for identification and classification.")


@app.command()
def sirius_init(path: str = typer.Option("configs/sirius.yml", help="Where to write the template YAML")):
    """Write a template YAML with SIRIUS settings (adducts, DBs, tasks)."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tpl = {
        "sirius": {"executable": "sirius", "cores": 8, "mzdev_ppm": 10.0},
        "adducts": {
            "neg_detectable": "[[M - H]-, [M - H2O - H]-]",
            "neg_fallback": "[[M]-, [M + CH2O2 - H]-, [M + C2HF3O2 - H]-, [M - H + Na - H]-, [M - H + K - H]-, [M + C2H4O2 - H]-, [M + C2H3N - H]-, [M - H]-, [M - H2O - H]-]",
            "pos_detectable": "[[M + H]+, [M + Na]+]",
            "pos_fallback": "[[M]+, [M + NH4]+, [M + Na]+, [M + K]+, [M + CH3OH + H]+, [M + ACN + H]+, [M + 2Na - H]+]",
        },
        "common_config": [
            "--IsotopeSettings.filter=true",
            "--FormulaSearchDB=",
            "--Timeout.secondsPerTree=0",
            "--FormulaSettings.enforced=HCNOPS",
            "--Timeout.secondsPerInstance=0",
            "--UseHeuristic.mzToUseHeuristicOnly=650",
            "--AlgorithmProfile=orbitrap",
            "--IsotopeMs2Settings=IGNORE",
            "--MS2MassDeviation.allowedMassDeviation=10.0ppm",
            "--NumberOfCandidatesPerIon=1",
            "--UseHeuristic.mzToUseHeuristic=300",
            "--FormulaSettings.detectable=,",
            "--NumberOfCandidates=10",
            "--AdductSettings.enforced=,",
            "--FormulaResultThreshold=true",
            "--InjectElGordoCompounds=true",
            "--StructureSearchDB=BIO,METACYC,CHEBI,COCONUT,ECOCYCMINE,GNPS,HMDB,HSDB,KEGG,KEGGMINE,KNAPSACK,MACONDA,MESH,NORMAN,UNDP,PLANTCYC,PUBCHEM,PUBMED,YMDB,YMDBMINE,ZINCBIO",
            "--RecomputeResults",
            "true",
        ],
        "tasks": ["formula", "fingerprint", "structure", "canopus"],
    }
    p.write_text(yaml.safe_dump(tpl, sort_keys=False), encoding="utf-8")
    log.ok(f"Wrote SIRIUS config template to {p}")


@app.command()
def sirius_export(
    pos_dir: str = typer.Option("outputs/sirius_pos", help="SIRIUS POS results folder"),
    neg_dir: str = typer.Option("outputs/sirius_neg", help="SIRIUS NEG results folder"),
    output_csv: str = typer.Option("outputs/sirius_identifications.csv", help="Where to write aggregated identifications or merged ultimate table"),
    join_merged: str = typer.Option("outputs/merged.csv", help="Merged CSV to join on feature_id; if a classified sibling (e.g., *_classyfire.csv) exists, it will be used automatically"),
    extract_csv: str = typer.Option("outputs/sirius_identifications.csv", help="Save/load extracted SIRIUS hits here; if present (and no --rescan), re-use instead of rescanning folders"),
    rescan: bool = typer.Option(False, help="Force rescan of SIRIUS folders and refresh the extract cache"),
):
    """Aggregate top CSI:FingerID candidates per SIRIUS compound and write a CSV.

    Parses each compound directory's spectrum and fingerid results to extract the best
    candidate and maps it back to the stable feature_id embedded in the .ms input.
    """
    j: Optional[Path] = None
    if join_merged and str(join_merged).strip():
        j = Path(join_merged)
        # Prefer classified output if present: <stem>_classyfire.<suffix>
        try:
            classified = j.with_name(j.stem + "_classyfire").with_suffix(j.suffix)
            if classified.exists():
                log.info(f"Detected classified table: {classified}; using it for merge")
                j = classified
        except Exception:
            pass

    # If user kept the default output path, derive an 'ultimate' path from the join target
    try:
        default_out = str(output_csv).replace("\\", "/").endswith("outputs/sirius_identifications.csv")
        if default_out and j is not None:
            outp = j.with_name(j.stem + "_ultimate").with_suffix(j.suffix)
            output_csv = str(outp)
            log.info(f"Writing ultimate results to: {output_csv}")
    except Exception:
        pass
    # Progress bar (rich) if available; otherwise fallback
    pd_path, nd_path = Path(pos_dir), Path(neg_dir)
    # If an extract CSV exists and no rescan requested, skip progress UI and reuse it
    reuse_extract = False
    if extract_csv and not rescan and Path(extract_csv).exists():
        reuse_extract = True
    if reuse_extract:
        try:
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
        except Exception as e:
            log.error(f"Failed to reuse identifications from {extract_csv}: {e}")
            raise typer.Exit(code=1)
    else:
        try:
            from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn  # type: ignore
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
                try:
                    prog.update(task_id, completed=done)
                except Exception:
                    pass
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
        except Exception:
            # Fallback without progress UI
            try:
                summary = collect_sirius_results(
                    pd_path,
                    nd_path,
                    Path(output_csv),
                    join_with_merged=j,
                    progress=None,
                    extract_cache=(Path(extract_csv) if extract_csv else None),
                    force_rescan=bool(rescan),
                )
            except Exception as e:
                log.error(f"Failed to export SIRIUS identifications: {e}")
                raise typer.Exit(code=1)
    log.ok(f"Wrote identifications/merge: {output_csv} (rows: {summary['rows']})")
    log.info(
        f"POS: total={summary.get('pos_compounds',0)}, identified={summary.get('pos_identified',0)}, "
        f"miss={summary.get('miss_pos',0)}; "
        f"NEG: total={summary.get('neg_compounds',0)}, identified={summary.get('neg_identified',0)}, miss={summary.get('miss_neg',0)}"
    )
    # Annotation level distribution after merge (if available)
    ann = summary.get('ann_after') or {}
    if ann:
        log.info(
            f"Annotation levels after merge: L1={ann.get('1',0)}, L2={ann.get('2',0)}, "
            f"L3={ann.get('3',0)}, L4={ann.get('4',0)}, L5={ann.get('5',0)}"
        )
@app.command()
def classify(
    input_csv: str = typer.Argument("outputs/merged.csv", help="Merged CSV to classify"),
    output_csv: str = typer.Option("outputs/merged_classified.csv", help="Where to write classified CSV (default derives from input)"),
    provider: str = typer.Option("classyfire", help="Classification provider (only 'classyfire' supported)"),
    base_url: str = typer.Option("http://classyfire.wishartlab.com", help="ClassyFire base URL"),
    sleep_sec: float = typer.Option(1.0, help="Sleep between API calls (seconds)"),
    results_only: bool = typer.Option(True, help="Write only identifier + classification columns"),
    id_column: str = typer.Option("feature_id", help="Identifier column to include (default: feature_id)"),
    cache_path: str = typer.Option("outputs/classyfire_cache.json", help="Cache file for InChIKey→taxonomy (json/yaml)"),
    offline: bool = typer.Option(False, help="Use cache only; skip network calls"),
    timeout: float = typer.Option(15.0, help="Per-request timeout (seconds)"),
    retries: int = typer.Option(2, help="Retries per InChIKey on transient errors"),
    backoff: float = typer.Option(1.5, help="Exponential backoff base between retries"),
    join_merged: bool = typer.Option(False, help="Append cf_* columns into the merged table (writes back to input CSV if output path is default)."),
    force_api: bool = typer.Option(False, help="Force API calls even if a classified CSV already exists"),
):
    """Classify Level 1/2 rows via Wishart lab's ClassyFire API using unique InChIKeys and merge back.

    Prints how many unique keys were queried and how many columns were added; you can then decide which
    classification columns to keep in your table.
    """
    if provider.lower() != "classyfire":
        log.error("Only 'classyfire' provider is supported currently.")
        raise typer.Exit(code=2)
    try:
        if join_merged:
            results_only = False
            if str(output_csv).replace("/", "\\") == "outputs\\merged_classified.csv" and input_csv:
                output_csv = input_csv
        # Compute derived classified path next to input for existence checks
        derived_out = None
        try:
            if str(output_csv).replace("/", "\\") == "outputs\\merged_classified.csv" and input_csv:
                ip = Path(input_csv)
                suffix = ip.suffix or ".csv"
                derived_out = ip.with_name(ip.stem + "_classyfire").with_suffix(suffix)
        except Exception:
            derived_out = None
        # If a classified CSV already exists and not forced, avoid API calls (offline)
        if not force_api:
            try:
                existing = Path(output_csv)
                if not existing.exists() and derived_out is not None and derived_out.exists():
                    existing = derived_out
                if existing.exists():
                    offline = True
                    log.info(f"Reusing existing classification: {existing} (offline)")
            except Exception:
                pass
        def _progress(i: int, n: int, key: str, cached: bool, status: str):
            src = "cache" if cached else "net"
            log.info(f"[{i}/{n}] {key}  {src}/{status}")
        summary = classify_level12_with_classyfire(
            input_csv,
            output_csv,
            sleep_sec=sleep_sec,
            base_url=base_url,
            results_only=results_only,
            id_column=id_column,
            cache_path=cache_path,
            progress=_progress,
            offline=offline,
            timeout=timeout,
            retries=retries,
            backoff=backoff,
        )
    except Exception as e:
        log.error(f"Classification failed: {e}")
        raise typer.Exit(code=1)
    log.ok(f"Classified: {summary['unique_keys']} unique InChIKeys (hits: {summary['hits']}, miss: {summary['miss']}, skipped_offline: {summary.get('skipped_offline', 0)})")
    log.info(f"Output: {output_csv} (rows: {summary['rows']}, classification columns added: {summary['added_cols']})")

if __name__ == "__main__":
    app()


@app.command()
def classify_check(
    base_url: str = typer.Option("http://classyfire.wishartlab.com", help="ClassyFire base URL to probe"),
    timeout: float = typer.Option(10.0, help="Timeout per request (seconds)"),
):
    """Probe the ClassyFire API endpoint for availability and latency."""
    try:
        res = probe_classyfire(base_url, timeout=timeout)
    except Exception as e:
        log.error(f"Probe failed: {e}")
        raise typer.Exit(code=1)
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
