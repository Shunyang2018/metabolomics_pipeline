# Metabolomics Pipeline (MS-DIAL Post-Processing)

Cross-platform (Windows/macOS) pipeline to parse MS-DIAL alignment outputs, apply basic QC filtering, harmonize sample columns, and export a merged matrix. The merge step also prepares SIRIUS input files for Level 3 unknowns.

## Status
- Typer-based CLI available with validate, merge, SIRIUS integration, and ClassyFire classification.
- SIRIUS integration: exports `.ms` inputs from Level 3 unknowns and can drive SIRIUS runs via `metabo sirius`.
- Configuration: All settings centralized in `src/metabo_pipeline/constants.py` (no config YAML files needed).
- Testing: Comprehensive test suite with unit and integration tests.
- See `PLAN.md` for the roadmap and milestones.

## Quickstart

1) **Install uv** (if not already installed):
   ```bash
   # macOS/Linux
   curl -LsSf https://astral.sh/uv/install.sh | sh

   # Windows
   powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
   ```

2) **Setup and install** (uv handles everything):
   ```bash
   # Create virtual environment and install package
   uv venv
   uv pip install -e .
   ```

3) **Activate environment and run**:
   - Windows (PowerShell): `.venv\Scripts\Activate.ps1`
   - macOS/Linux: `source .venv/bin/activate`
   - Then run: `metabo --help`

**Quick commands:**
- `metabo merge` (uses INPUT_DIR from constants.py, outputs to INPUT_DIR/outputs/)
- `metabo merge /custom/path --output-dir /custom/output` (override paths)
- Outputs include SIRIUS inputs: `outputs/sirius_unknown_pos.ms` and `outputs/sirius_unknown_neg.ms` built from Level 3 unknowns (m/z 150–800)

## Commands
- `metabo run`: **run the complete pipeline** - merge → (classify + sirius in parallel) → final → bioactivity. One command to do it all!
- `metabo merge`: combine multiple alignment tables into a wide matrix; writes `outputs/merged.csv` and SIRIUS `.ms` files for Level 3 unknowns.
- `metabo classify`: classify Level 1/2 rows via ClassyFire API using unique InChIKeys.
- `metabo sirius`: run SIRIUS on the exported `.ms` files (formula, fingerprint, structure, CANOPUS).
- `metabo final`: create the final merged table combining MS-DIAL, ClassyFire, and SIRIUS results.
- `metabo bioactivity`: match features against a bioactivity reference database via InChIKey skeleton; writes `outputs/bioactives.csv`.
- `metabo classify-check`: probe a ClassyFire endpoint for availability/latency.

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
  - Configure SIRIUS path in `src/metabo_pipeline/constants.py` (set `SIRIUS_EXECUTABLE`)
  - Optional: edit `configs/sirius.yml` to adjust adducts, DBs, and tasks.
  - Run SIRIUS on exported files:
    - `metabo sirius` (uses defaults from constants.py)
    - `metabo sirius --test-real-case` (quick test: 4 compounds, ~2-3 min, real SIRIUS execution)
    - `metabo sirius --dry-run` (instant test: creates fake output without running SIRIUS)
  - Merge results back (adds values, keeps original columns):
    - `metabo final`
    - Defaults: reuses `outputs/sirius_identifications.csv` if present (no progress bar). To force a rescan of folders, add `--rescan`.
    - If `outputs/merged_classyfire.csv` exists, it is used automatically as the join target; output defaults to `<join>_final.csv`.

## Bioactivity Matching

- `metabo bioactivity` matches identified/tentatively-identified features against a reference
  bioactivity database (any CSV with an `InChIKey` column, e.g. a natural-product bioactivity
  database) and writes `outputs/bioactives.csv`.
- Matching is done on the **first block of the InChIKey** (the 14-character skeleton hash before
  the first hyphen), so a feature still matches the database even if stereochemistry, protonation
  state, or isotope labeling differ slightly.
- Output shape: one row per `(feature, sample it was detected in, bioactivity hit)` — i.e. MS-DIAL
  feature identity + `sample_type` (the sample column name) + `intensity`, followed by every column
  from the bioactivity database (`Molecule`, `Benefit`, `Endpoint`, `Evidence.Strength`, etc.). A
  feature with several database hits (synonyms, multiple reported benefits) produces one row per
  hit; only samples with nonzero intensity are included.
- Configure the default database path via `BIOACTIVITY_DB_PATH` in `constants.py`, or override with
  `metabo bioactivity --db /path/to/database.csv`.
- Runs automatically as the last step of `metabo run` (skip with `--skip-bioactivity`).

### Installing SIRIUS on Windows and macOS
- Windows: Install SIRIUS to `C:\Program Files\SIRIUS\` or set `SIRIUS_HOME`; the CLI auto-detects common locations and honours `SIRIUS_EXECUTABLE`.
- macOS: Place `sirius.app` under `/Applications` or `~/Applications`; the CLI resolves `<bundle>/Contents/MacOS/sirius` or use `--sirius-exe`/`SIRIUS_HOME` for custom installs.
- Any platform: Keeping `sirius` on `PATH` works out-of-the-box and `metabo sirius` logs the executable it will run.

### SIRIUS Authentication (Required for SIRIUS 6+)

SIRIUS 6+ requires authentication to use the command-line tool. You need a free academic account from https://bio.informatik.uni-jena.de/sirius-register/.

**⚠️ Important: Use GUI Login (Recommended)**

Due to SSL certificate handling issues with CLI login, **login via SIRIUS GUI is required**:

1. **Open SIRIUS application** (double-click `sirius.app` on macOS or run SIRIUS on Windows)

2. **Login with your account:**
   - Look for "Login" or "Account" in the menu
   - Enter your email and password
   - Wait for "Logged in" confirmation

3. **Close SIRIUS** completely

4. **Verify login worked:**
   ```bash
   /path/to/sirius login --show
   ```
   Should show your account info instead of "Not logged in"

5. **Session persists** - you only need to login once; the session is saved and used by CLI commands

**Note:** Login sessions may expire after ~24 hours. If `metabo sirius` fails with "Login ERROR", repeat the GUI login steps above.

**Verify login status:**
```bash
/path/to/sirius login --show
```

**Note:** Your SIRIUS login session is saved in `~/.sirius-6.3/` and persists across runs until it expires (typically 24 hours). If you see "Login ERROR: Please Login to use SIRIUS", simply repeat the GUI login steps above.

## Quick Start: Full Pipeline in One Command

The fastest way to run the entire pipeline:

```bash
# Activate environment (if not already active)
source .venv/bin/activate  # macOS/Linux
# or: .venv\Scripts\Activate.ps1  # Windows

# Ensure SIRIUS is logged in (GUI login - see above)
# Check: /path/to/sirius login --show

# Run everything: merge → (classify + sirius in parallel) → final
metabo run
```

**What `metabo run` does:**
1. 📊 Merges MS-DIAL files into a single table
2. 🔬 Runs ClassyFire (Level 1/2) and SIRIUS (Level 3) **in parallel** to save time
3. 📋 Creates the final comprehensive table with all results

**Options:**
- `metabo run --skip-merge` - Skip merge if you already have `merged.csv`
- `metabo run --skip-classify` - Skip ClassyFire classification
- `metabo run --skip-sirius` - Skip SIRIUS analysis
- `metabo run --input-dir /path` - Custom input directory
- `metabo run --output-dir /path` - Custom output directory

**Output:** `outputs/merged_classyfire_final.csv` - your final table! 🎉

---

## End-to-End Steps (Manual Control)

- Merge MS-DIAL alignment outputs to a single table and generate SIRIUS inputs:
  - `metabo merge` (uses INPUT_DIR from constants.py)

- Classify Level 1/2 compounds (writes a classified copy next to the input):
  - `metabo classify outputs/merged.csv --join-merged`
  - Result: `outputs/merged_classyfire.csv` (adds `cf_*` columns after `SMILES`)
  - Reuse existing classified CSV by default (no extra API calls). Force a refresh with `--force-api`.

- Run SIRIUS on the exported `.ms` files (POS/NEG):
  - `metabo sirius` (uses defaults from constants.py)
  - `metabo sirius --test-real-case` (quick test: 4 compounds, ~2-3 min)
  - `metabo sirius --dry-run` (instant test: fake outputs)
  - `metabo sirius --force` (force re-run even if outputs exist)

- Collect SIRIUS results into CSV for inspection (optional step):
  - `metabo sirius-collect` (exports SIRIUS results from TSV summaries to `sirius_identifications.csv`)
  - `metabo sirius-collect --force` (force re-scan, ignore cached extractions)
  - This allows you to inspect SIRIUS annotations before merging into the final output
  - Requires SIRIUS TSV summary directories (`sirius_pos/`, `sirius_neg/`)

- Create the final merged table with all results (MS-DIAL + ClassyFire + SIRIUS):
  - `metabo final` (uses defaults from constants.py)
  - `metabo final --rescan` (force rescan of SIRIUS folders)
  - If `outputs/merged_classyfire.csv` exists, this writes `outputs/merged_classyfire_final.csv` - your final comprehensive table!

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
  - Caching: results are cached by default at `src/metabo_pipeline/classyfire_cache.json` (package directory with 10K+ pre-cached compounds). Re-runs reuse cached InChIKeys and only query missing ones. Use `--cache-path` to change location or format (`.json` or `.yaml`).
  - Endpoint: defaults to `http://classyfire.wishartlab.com`. Use `metabo classify-check` to probe availability.

## Dev Notes
- Data files are git-ignored (e.g., `M2_*.csv`, data folders, archives).
- Line endings are normalized via `.gitattributes`.
- Formatting is enforced via pre-commit/Black 24.8.0 (`pip install pre-commit`, run `pre-commit install`, and `pre-commit run --all-files` for a full check).

## License
TBD.
