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
- Long format (one intensity per row) if needed:
  - `metabo merge "." --output outputs\\merged_long.csv --format long --engine pandas`

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

## Dev Notes
- Data files are ignored by `.gitignore` (see rules for `M2_*.csv`, `M2 data/`, and archives).
- Line endings normalized via `.gitattributes`.

## License
TBD.
