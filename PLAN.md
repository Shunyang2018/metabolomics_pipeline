# MS-DIAL Post-Processing Pipeline — Project Plan

This document outlines a cross‑platform metabolomics pipeline (macOS + Windows) to process MS‑DIAL alignment outputs and generate cleaned data matrices, statistics, and an HTML report.

## Goals
- Robustly parse MS‑DIAL outputs and sample metadata.
- Provide configurable QC, filtering, normalization, imputation, and batch/drift correction.
- Export clean matrices and stats; generate an HTML report with plots.
- Deliver a simple CLI with reproducible runs and clear logs.

## Inputs
- MS‑DIAL aligned peak table: `AlignmentResult*.txt/.csv` (e.g., `*_nistannotated.csv`).
- Optional: aligned spectra `*.msp` for context.
- Sample metadata `samples.csv` with `SampleID, Group, Batch, InjectionOrder, IsQC, ...` (optional if embedded rows in MS‑DIAL are used).
- Optional: internal standards list `internal_standards.csv`.

## Outputs
- `matrix_clean.csv`: features × samples intensity matrix.
- `features.csv`: feature metadata (m/z, RT, adduct, annotation, flags).
- `samples_validated.csv`: validated sample metadata.
- `stats_univariate.csv`, `stats_multivariate.json`.
- `report/index.html`: QC, PCA, volcano, drift, missingness, summaries.
- `logs/` and `run_manifest.json` for reproducibility.

## Tech Stack
- Language: Python 3.10+ (cross‑platform).
- CLI: `typer` (+ `rich` for logging/progress).
- Data: `pandas`, `numpy`, `scipy`.
- Stats: `statsmodels`, `scikit-learn` (PCA), optional `pingouin`.
- Plots: `matplotlib`/`seaborn` (optional `plotly` for interactive sections).
- Batch correction: QC‑LOESS; ComBat (via `pycombat` or equivalent) as optional.
- Config/schema: `pydantic` for validation; YAML; `jinja2` for HTML templates.
- Packaging: `pipx`/`uv` for install; optional PyInstaller one‑file exe for Windows.

## Pipeline Stages
1) Ingest
- Parse MS‑DIAL alignment table (handles variant headers and encodings, e.g., utf‑8‑sig).
- Extract embedded metadata from top rows (Class, File type, Injection order, Batch ID) or merge with `samples.csv`.
- Build long and wide representations; map measurement columns to `SampleID`.

2) QC & Filtering
- Missingness thresholds: overall and per‑group.
- QC‑based feature RSD filter (e.g., RSD_QC ≤ 30%).
- Minimum detection within biological groups; blank filtering if present.
- Optional deduplication hints (isotopes/adducts) respecting MS‑DIAL annotations.

3) Normalization
- Options: total sum scaling, median, PQN, internal‑standard based.
- Log transform (log2/log10) with small offset; scaling (Pareto/auto) as configured.

4) Imputation
- Options: half‑min (per feature/group), KNN, left‑censored aware (QRILC‑like) as optional.

5) Batch/Drift Correction
- QC‑LOESS drift correction by injection order per feature.
- ComBat with optional covariates; allow disabling if metadata insufficient.

6) Statistics
- Univariate: t‑test/ANOVA with Benjamini–Hochberg FDR; effect sizes.
- Multivariate: PCA (scores/loadings); PLS‑DA flagged as exploratory.

7) Reporting
- QC plots: TIC proxy, missingness heatmap, RSD distributions, drift curves.
- PCA (colored by Group/Batch), volcano plots, top features table.
- Summary of thresholds, counts kept/removed, and warnings.

8) Export
- Save cleaned artifacts, report, and a JSON manifest capturing config and versions.

## Configuration
- Single YAML file: `config.yml`
  - Paths: input/output, file patterns, mode (pos/neg) if needed.
  - Thresholds: missingness, `RSD_QC`, min group coverage, blank fold‑change.
  - Methods: normalization, imputation, drift/batch correction.
  - Stats: group variable, covariates, multiple testing method.
  - Report switches: which plots/sections to include.

## CLI Design
- `metabo init` → scaffold `config.yml`, `samples.csv` template, `templates/`.
- `metabo validate` → check inputs, headers, metadata integrity; dry run.
- `metabo run` → execute full pipeline; `--steps` to run subset; `--config` path.
- `metabo report` → regenerate report from saved artifacts.
- `metabo diag` → emit environment info and versions.

## Data Model
- Stable `FeatureID` from MS‑DIAL alignment rows.
- Feature metadata: `mz, rt_min, adduct, name, msms_score, annotation_source, isotope, comment, ...`.
- Measurement columns mapped 1:1 to `SampleID`; long format available internally.

## Cross‑Platform Considerations
- Pure‑Python; avoid shell‑specific commands and R/Bioconductor hard deps.
- Robust file handling: `pathlib`, tolerant encodings, Windows path lengths, CRLF/UTF‑8‑SIG.
- Optional one‑file Windows build via PyInstaller.

## Project Structure
```
src/metabo_pipeline/
  cli.py          # Typer CLI entry points
  config.py       # Pydantic models + YAML load/validate
  io_msdial.py    # MS‑DIAL parser and harmonizer
  qc.py           # Missingness, RSD, blank, flags
  normalize.py    # Normalization options
  impute.py       # Imputation methods
  batch.py        # QC‑LOESS + ComBat wrappers
  stats.py        # Univariate + PCA
  report.py       # Jinja2 HTML report
  logging.py      # Rich logger setup
configs/          # Example configs
templates/        # Jinja2 + static assets
examples/         # Small MS‑DIAL test dataset + sample metadata
tests/            # Pytest for parser, QC, normalization, CLI
```

## Milestones
- M1: Skeleton repo + CLI (`init`, `validate`, `run` stubs), config schema.
- M2: Robust MS‑DIAL parser + metadata validation; golden tests.
- M3: QC filtering (missingness, RSD, blank) + metrics.
- M4: Normalization + imputation options with unit tests.
- M5: Drift correction (QC‑LOESS) + optional ComBat; diagnostics plots.
- M6: Stats (t‑test/ANOVA, FDR) + PCA; CSV exports.
- M7: HTML report with key plots and summaries.
- M8: Cross‑platform packaging; example dataset; quickstart docs.

## Performance & Reliability
- Efficient `pandas` usage; vectorized ops; typed dtypes where possible.
- Deterministic seeds; manifest with config hash and dependency versions.
- Clear, actionable validation errors and warnings (e.g., no QC samples found).

## Open Decisions
- Initial acquisition types: HILIC + C18 (+ Lipidomics), pos/neg both?
- Availability of pooled QC samples for RSD and LOESS?
- Default normalization (PQN vs median vs IS‑based).
- Batch metadata fidelity (Batch and InjectionOrder reliable?).
- Required plots/sections for v1 report.
- Packaging preference: `pipx` vs one‑file exe as primary.

## Next Steps
1) Confirm target datasets (HILIC/C18/Lipidomics, pos/neg).
2) Decide on metadata source: embedded rows vs separate `samples.csv`.
3) Scaffold repo (M1) and implement MS‑DIAL parser (M2).
4) Run `validate` on a provided file to produce a quick summary report (counts, missingness, QC/blank detection).

---
This plan is intended to be living documentation; we’ll update it as we agree on defaults and make implementation choices.
