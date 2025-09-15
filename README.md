# Metabolomics Pipeline (MS-DIAL Post-Processing)

Cross-platform (Windows/macOS) pipeline to parse MS-DIAL alignment outputs, perform QC/filtering/normalization, and export cleaned matrices with a concise report.

## Status
- Project scaffolded with a Typer CLI and a basic MS-DIAL validator.
- See `PLAN.md` for the roadmap and milestones.

## Quickstart
1) Create/activate a Python 3.10+ environment
- Windows (PowerShell):
  - `py -3 -m venv .venv`
  - `.venv\\Scripts\\Activate.ps1`
- macOS/Linux:
  - `python3 -m venv .venv && source .venv/bin/activate`

2) Install in editable mode
- With uv (recommended):
  - `uv venv`
  - `.venv\\Scripts\\Activate.ps1`  (Windows)  or `source .venv/bin/activate` (macOS/Linux)
  - `uv pip install -e .`
- With pip:
  - `pip install -e .`

3) Run CLI
- `metabo --help`
- Validate an MS-DIAL table:
  - `metabo validate --input M2_HILIC_NEG_nistannotated.csv`
- Merge all files (default: WIDE, one feature per row, normalized sample columns):
  - `metabo merge "." --output outputs\\merged.csv`  (pandas required)
  - Leading columns include `isomer_label` (assigned by RT clustering within metabolite + adduct), then `chrom`, and core IDs.

## Commands
- `metabo init` — scaffold a config template and folders.
- `metabo validate` — parse MS-DIAL alignment table(s), summarize samples/features, and check embedded metadata rows.
- `metabo merge` — combine multiple alignment tables:
  - Default: wide format (one feature per row; normalized sample columns across assays).
  - Long format: `--format long` (columns: chrom, annotation_level, alignment_id, rt_min, mz, metabolite_name, adduct, sample_id, intensity).
- `metabo run` — placeholder for the full pipeline (coming milestones M3–M7).
- `metabo diag` — print environment info and versions.

## Inputs
- MS-DIAL alignment tables (`*.csv`/`*.txt`) with the four metadata rows: Class, File type, Injection order, Batch ID.
- Optional: `samples.csv` for richer metadata.

### Sample ID normalization
- During merge, sample column headers are consolidated to a common core by removing assay and polarity tokens and lowercasing. Examples:
  - `M2_Lipids_TV_5+6_pos` → `m2_tv_5+6`
  - `M2_C18_TV_5+6_POS` → `m2_tv_5+6`
  - `M2_Hilic_PE_3+4_neg` → `m2_pe_3+4`

## Processing Logic

**Filtering**
- MS/MS spectrum required with ≥ 3 nonzero fragment ions (parses `MS/MS spectrum`).
- S/N pre-filter: drops rows with `S/N average < 5.0` (if column present).
- Replicate-group QC (per file): groups inferred from normalized sample names (e.g., `m2_tv_1+2`, `m2_tv_3+4`). For each group:
  - `blank_fold_<grp>` = max(replicates) / Blank; NaN if Blank missing or zero.
  - `present_percent_<grp>` = percent of replicates with intensity > 0.
  - `cv_percent_<grp>` = CV among present replicates (std/mean × 100; requires ≥ 2 present).
- Pass flags and gating: `pass_<grp>` is True if `blank_fold ≥ 7`, `present_percent ≥ 60`, `cv_percent ≤ 40`. A feature is kept if any group passes (`pass_any_groups`).

**Merging (Wide Format, default)**
- Reads all alignment tables, applies filters, and appends features (rows) from all sources.
- Normalizes sample columns to a common set across assays (HILIC/C18/Lipidomics), keeping one column per biological sample.
- Leading columns: `isomer_label`, `chrom`, `annotation_level`, `alignment_id`, `rt_min`, `mz`, `metabolite_name`, `adduct`, followed by optional MS‑DIAL metadata, then sample columns, then QC metrics and pass flags.

**Isomer Labeling**
- Within each file, cluster by RT inside (`metabolite_name`, `adduct`) with a 1.0 min window.
- Assign ordered labels by RT: `isomer_1`, `isomer_2`, …
- After merging, suffix the `metabolite_name` with its isomer label (e.g., `Adenine_isomer_1`) so isomers are distinct for downstream steps.

**Deduplication**
- Runs after all files are merged, at the name level (post isomer suffix).
- One row per `metabolite_name`: keep the row with lowest `cv_median_percent` (median of group CVs), breaking ties by highest average intensity across normalized sample columns.
- Other evidence (e.g., scores) is retained as columns but not used to rank by default.

## Dev Notes
- Data files are ignored by `.gitignore` (see rules for `M2_*.csv`, `M2 data/`, and archives).
- Line endings normalized via `.gitattributes`.

## License
TBD.

