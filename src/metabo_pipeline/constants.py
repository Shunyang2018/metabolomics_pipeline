"""
Centralized constants for pipeline thresholds and defaults.

All QC filtering and processing parameters are defined here.
Adjust these values to customize the pipeline behavior.
"""

# =============================================================================
# MS/MS FILTERING
# =============================================================================
MSMS_MIN_IONS: int = 3
"""Minimum number of nonzero fragment ions required in MS/MS spectrum.
Features with fewer than this many fragments are filtered out."""

# =============================================================================
# SIGNAL-TO-NOISE FILTERING
# =============================================================================
SNR_MIN: float = 5.0
"""Minimum signal-to-noise ratio (if 'S/N average' column is present).
Features with S/N below this threshold are filtered out."""

# =============================================================================
# ISOMER LABELING
# =============================================================================
ISOMER_RT_WINDOW_MIN: float = 1.0
"""RT window (minutes) for clustering potential isomers.
Features with same metabolite name + adduct within this RT window are labeled as isomers."""

# =============================================================================
# REPLICATE-GROUP QC THRESHOLDS
# =============================================================================
BLANK_FOLD_MIN: float = 7.0
"""Minimum blank fold-change threshold.
For each replicate group: max(sample_intensities) / mean(blank_intensities) must be >= this value.
Features that don't pass this threshold in any group are filtered out."""

PRESENT_PERCENT_MIN: float = 60.0
"""Minimum detection percentage within replicate groups.
Percentage of replicates with intensity > 0 must be >= this value.
Features that don't meet this threshold in any group are filtered out."""

CV_PERCENT_MAX: float | None = None
"""Maximum coefficient of variation (CV%) allowed within replicate groups.
Set to None to disable CV filtering.
If set to a value (e.g., 30.0), features with CV > this value are filtered out."""

# =============================================================================
# LEVEL 3 DEDUPLICATION FOR SIRIUS EXPORT
# =============================================================================
DEDUP_RT_WINDOW_MIN: float = 0.6
"""RT window (minutes) for clustering Level 3 unknowns during deduplication.
Features within this RT window are considered potential duplicates."""

DEDUP_MZ_PPM: float = 20.0
"""Mass accuracy window (ppm) for clustering Level 3 unknowns during deduplication.
Features within this m/z tolerance are considered potential duplicates."""

# =============================================================================
# SIRIUS SETTINGS
# =============================================================================
SIRIUS_VERSION: str = "v6"
"""SIRIUS output format version: "v5" (directory structure) or "v6" (TSV summaries).
- v5: SIRIUS 5 creates compound_1/, compound_2/ directories with structure_candidates.tsv
- v6: SIRIUS 6 creates .sirius project files and TSV summary files
Set this to match your SIRIUS installation."""

SIRIUS_EXECUTABLE: str = "/Users/wangs261/Downloads/sirius.app/Contents/MacOS/sirius"
"""Path to SIRIUS executable on your system.
The CLI will auto-detect common locations if this is not found.
macOS: typically /Applications/sirius.app/Contents/MacOS/sirius
Windows: typically C:\\Program Files\\SIRIUS\\sirius.exe
Linux: typically /usr/local/bin/sirius or /opt/sirius/bin/sirius"""

SIRIUS_CORES: int = 8
"""Number of CPU cores for SIRIUS to use during computations."""

SIRIUS_MZDEV_PPM: float = 10.0
"""MS2 mass deviation tolerance (ppm) for SIRIUS analysis."""

# =============================================================================
# INPUT/OUTPUT DIRECTORIES
# =============================================================================
INPUT_DIR: str = "/Users/wangs261/Documents/project/excel_merge/test_csv"
"""Default directory containing MS-DIAL alignment output files (CSV/TXT).
Can be overridden via command line: metabo merge <custom_path>"""

# =============================================================================
# BIOACTIVITY MATCHING
# =============================================================================
BIOACTIVITY_DB_PATH: str = (
    "/Users/ivanablazenovic/Downloads/"
    "natural_product_metabolite_bioactivity_database_may_2026.csv"
)
"""Path to the bioactivity reference database CSV (must include an InChIKey column).
Can be overridden via command line: metabo bioactivity --db <custom_path>"""
