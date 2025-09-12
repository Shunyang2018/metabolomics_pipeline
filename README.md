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

2) Install in editable mode (or use `pipx`)
- `pip install -e .`

3) Run CLI
- `metabo --help`
- Validate an MS-DIAL table:
  - `metabo validate --input M2_HILIC_NEG_nistannotated.csv`
- Merge all files in a folder (HILIC/C18/Lipidomics) to long format:
  - `metabo merge "F:\\Shunyang pipeline" --output outputs\\merged_long.csv`

## Commands
- `metabo init` — scaffold a config template and folders.
- `metabo validate` — parse MS-DIAL alignment table(s), summarize samples/features, and check embedded metadata rows.
- `metabo merge` — combine multiple alignment tables into one long-format CSV, annotating source (chromatography and filename). New columns appear first: `chrom`, `source_file`.
- `metabo run` — placeholder for the full pipeline (coming milestones M3–M7).
- `metabo diag` — print environment info and versions.

## Inputs
- MS-DIAL alignment tables (`*.csv`/`*.txt`) with the four metadata rows: Class, File type, Injection order, Batch ID.
- Optional: `samples.csv` for richer metadata.

## Dev Notes
- Data files are ignored by `.gitignore` (see rules for `M2_*.csv`, `M2 data/`, and archives).
- Line endings normalized via `.gitattributes`.

## License
TBD.
