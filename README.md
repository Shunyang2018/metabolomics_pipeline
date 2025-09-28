# Metabolomics Pipeline (MS-DIAL Post-Processing)

Cross-platform (Windows/macOS) pipeline to parse MS-DIAL alignment outputs, apply basic QC filtering, harmonize sample columns, and export a merged matrix. The merge step also prepares SIRIUS input files for Level 3 unknowns.

## Status
- Typer-based CLI available with validate, merge, and diagnostics.
- SIRIUS integration: exports `.ms` inputs from Level 3 unknowns and can drive SIRIUS runs via `metabo sirius`.
- See `PLAN.md` for the roadmap and milestones.

## Quickstart
1) Create/activate a Python 3.10+ environment
- Windows (PowerShell):
  - `py -3 -m venv .venv`
  - `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`
  - `.venv\\Scripts\\Activate.ps1`
- macOS/Linux:
  - `python3 -m venv .venv && source .venv/bin/activate`

2) Install in editable mode
- With uv (recommended):
  - `uv venv`
  - `.venv\\Scripts\\Activate.ps1` (Windows) or `source .venv/bin/activate` (macOS/Linux)
  - `uv pip install -e .`
- With pip:
  - `pip install -e .`

3) Run CLI
- `metabo --help`
- Validate an MS-DIAL table:
  - `metabo validate --input M2_HILIC_NEG_nistannotated.csv`
- Merge files from a folder (wide format: one feature per row):
  - `metabo merge . --output outputs/merged.csv`
  - Outputs also include SIRIUS inputs: `outputs/sirius_unknown_pos.ms` and `outputs/sirius_unknown_neg.ms` built from Level 3 unknowns (m/z 150–800) after deduplication.

## Commands
- `metabo init`: scaffold a minimal config and folders.
- `metabo validate`: parse MS-DIAL alignment table(s), summarize samples/features, and check embedded metadata rows.
- `metabo merge`: combine multiple alignment tables into a wide matrix; writes `outputs/merged.csv` and SIRIUS `.ms` files for Level 3 unknowns.
- `metabo sirius-init`: write a template YAML at `configs/sirius.yml` for SIRIUS settings (adducts, DBs, tasks).
- `metabo sirius`: run SIRIUS on the exported `.ms` files (formula, fingerprint, structure, CANOPUS).
- `metabo sirius-export`: aggregate SIRIUS identifications and merge back into the table by `feature_id`.
- `metabo classify`: classify Level 1/2 rows via ClassyFire API using unique InChIKeys.
- `metabo classify-check`: probe a ClassyFire endpoint for availability/latency.
- `metabo run`: placeholder for the full pipeline.
- `metabo diag`: print environment info and versions.
- `metabo clean`: remove generated outputs (SIRIUS exports, caches, optional merged tables).

## Inputs
- MS-DIAL alignment tables (`*.csv`/`*.txt`) that include the four metadata rows: Class, File type, Injection order, Batch ID.
- Optional: `samples.csv` with richer sample metadata (not required if embedded rows are sufficient).

### Sample ID normalization
- During merge, sample column headers are normalized to a common core by removing assay and polarity tokens and lowercasing. Examples:
  - `M2_Lipids_TV_5+6_pos` -> `m2_tv_5+6`
  - `M2_C18_TV_5+6_POS` -> `m2_tv_5+6`
  - `M2_Hilic_PE_3+4_neg` -> `m2_pe_3+4`

## Processing

- Filtering: requires an MS/MS spectrum with at least 3 nonzero fragment ions; if `S/N average` is present, drops rows with S/N < 5.0.
- Replicate-group QC (per file): groups inferred from normalized sample names (e.g., `m2_tv_1+2`). For each group:
  - `blank_fold_<grp>` = max intensity across all sample replicates / mean of all blank columns (NaN if blanks missing or zero)
  - `present_percent_<grp>` = percent of replicates with intensity > 0 (default = 60%; configurable via `PRESENT_PERCENT_MIN`)
  - `cv_percent_<grp>` = CV among present replicates (requires at least 2 present; optional filter via `CV_PERCENT_MAX`)
  - A feature passes a group if: blank_fold >= 7 and present_percent >= 60 (default `PRESENT_PERCENT_MIN`). When `CV_PERCENT_MAX` is set, cv_percent must also be <= that value. A row is kept if any group passes (`pass_any_groups`).
- Merging (wide): union all passing rows, normalize sample columns to a common set across assays (HILIC/C18/Lipidomics), keep one column per biological sample.
- Column order starts with identifiers and MS-DIAL metadata such as `feature_id`, `chrom`, `annotation_level`, `Alignment ID`, `Average Rt(min)`, `Average Mz`, `Metabolite name`, `Adduct type`, followed by other MS-DIAL fields, then sample columns, and finally QC metrics.
- Level 3 handling: all Level 3 `Metabolite name` values are set to `Unknown` for clarity.
- Deduplication: name conflicts are resolved within RT/m/z clusters; Level 3 rows are further reduced to representatives and restricted to m/z 150–800 for SIRIUS export.

## SIRIUS Integration

- Merge writes SIRIUS input files automatically:
  - `outputs/sirius_unknown_pos.ms`
  - `outputs/sirius_unknown_neg.ms`
  These are built from Level 3 unknowns using the MS-DIAL `MS1 isotopic spectrum` and `MS/MS spectrum` columns and include metadata like precursor m/z, RT, adduct, and a stable `feature_id`.

- Configure and run SIRIUS:
  - Create a template config: `metabo sirius_init --path configs/sirius.yml`
  - Edit `configs/sirius.yml` to adjust executable path, cores, ppm, adducts, and tasks.
  - Run SIRIUS on exported files (example with defaults):
    - `metabo sirius --config configs/sirius.yml --pos-ms outputs/sirius_unknown_pos.ms --neg-ms outputs/sirius_unknown_neg.ms --out-pos outputs/sirius_pos --out-neg outputs/sirius_neg`
  - Merge results back (adds values, keeps original columns):
    - `metabo sirius-export --pos-dir outputs/sirius_pos --neg-dir outputs/sirius_neg --join-merged outputs/merged.csv`
    - Defaults: reuses `outputs/sirius_identifications.csv` if present (no progress bar). To force a rescan of folders, add `--rescan`.
    - If `outputs/merged_classyfire.csv` exists, it is used automatically as the join target; output defaults to `<join>_ultimate.csv`.

### Installing SIRIUS on Windows and macOS
- Windows: Install SIRIUS to `C:\Program Files\SIRIUS\` or set `SIRIUS_HOME`; the CLI auto-detects common locations and honours `SIRIUS_EXECUTABLE`.
- macOS: Place `sirius.app` under `/Applications` or `~/Applications`; the CLI resolves `<bundle>/Contents/MacOS/sirius` or use `--sirius-exe`/`SIRIUS_HOME` for custom installs.
- Any platform: Keeping `sirius` on `PATH` works out-of-the-box and `metabo sirius` logs the executable it will run.

### Cleaning Generated Artifacts
- `metabo clean --dry-run` previews the cleanup list (SIRIUS exports, caches) without deleting anything.
- `metabo clean` deletes `outputs/sirius_unknown_*.ms`, `outputs/sirius_pos*/`, `outputs/sirius_neg*/`, `outputs/sirius_identifications.csv`, and `outputs/classyfire_cache.*`.
- Add `--prune-tables` to drop derived CSVs (e.g., `merged_classyfire.csv`, `_ultimate.csv`), `--prune-all` to remove the entire `outputs` directory, and `--extra` for custom paths.

## End-to-End Steps

- Merge MS-DIAL alignment outputs to a single table and generate SIRIUS inputs:
  - `metabo merge . --output outputs/merged.csv`

- Classify Level 1/2 compounds (writes a classified copy next to the input):
  - `metabo classify outputs/merged.csv --join-merged`
  - Result: `outputs/merged_classyfire.csv` (adds `cf_*` columns after `SMILES`)
  - Reuse existing classified CSV by default (no extra API calls). Force a refresh with `--force-api`.

- Run SIRIUS on the exported `.ms` files (POS/NEG):
  - `metabo sirius --config configs/sirius.yml --pos-ms outputs/sirius_unknown_pos.ms --neg-ms outputs/sirius_unknown_neg.ms --out-pos outputs/sirius_pos --out-neg outputs/sirius_neg`

- Merge SIRIUS results into the classified table to get an all-in-one “ultimate” table:
  - `metabo sirius-export --join-merged outputs/merged.csv`
  - If `outputs/merged_classyfire.csv` exists, this writes `outputs/merged_classyfire_ultimate.csv`.
  - Add `--rescan` the first time to extract formula/CANOPUS-only evidence for L4.

- Optional: overwrite in place instead of writing derived names
  - Classify in place: `metabo classify outputs/merged.csv --join-merged --output-csv outputs/merged.csv`
  - Then SIRIUS merge defaults to `outputs/merged_ultimate.csv`:
    - `metabo sirius-export --join-merged outputs/merged.csv`

SIRIUS merge behavior
- Keeps all original rows and columns. Only values are updated:
  - `Metabolite name`: if a SIRIUS hit exists, set to `SIRIUS_<candidate_name>`.
  - Adds CANOPUS classes (only new columns): `SIRIUS_canopus_superclass`, `SIRIUS_canopus_class`, `SIRIUS_canopus_most_specific`.
  - Updates `annotation_level` where original is missing or `"3"`:
    - Level `3`: SIRIUS has structure (SMILES and InChIKey).
    - Level `4`: otherwise if SIRIUS has formula or CANOPUS classification.
    - Level `5`: no structure, no formula/classification from SIRIUS.
- Shows a progress bar (rich) and end-of-run stats.
  - When reusing an existing identifications CSV, the progress bar is skipped and the tool derives POS/NEG counts and level distribution from the file.

## Classification (ClassyFire)

- Classify Level 1/2 rows via the Wishart lab ClassyFire API using InChIKeys:
  - Full table with taxonomy columns added (writes `<input>_classyfire.csv` by default):
    - `metabo classify outputs/merged.csv --join-merged`
  - Compact lookup (id + taxonomy only):
    - `metabo classify outputs/merged.csv --results-only true`
- Notes:
  - Default timeout: 15s; adjust pacing with `--sleep-sec`, retries with `--retries`.
  - `--join-merged` adds `cf_*` columns into the table (inserted after `SMILES` when present). Keep your input safe by writing to a new path or accept the derived `<input>_classyfire.csv`.
  - Caching: results are cached by default at `outputs/classyfire_cache.json`. Re-runs reuse cached InChIKeys and only query missing ones. Use `--cache-path` to change location or format (`.json` or `.yaml`).
  - Endpoint: defaults to `http://classyfire.wishartlab.com`. Use `metabo classify-check` to probe availability.

## Dev Notes
- Data files are git-ignored (e.g., `M2_*.csv`, data folders, archives).
- Line endings are normalized via `.gitattributes`.

## License
TBD.
