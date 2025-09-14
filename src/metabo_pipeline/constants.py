"""
Centralized constants for pipeline thresholds and defaults.

Adjust here or expose via CLI/config in future iterations.
"""

# MS/MS filtering
MSMS_MIN_IONS: int = 3  # minimum nonzero fragment ions required

# Signal-to-noise filtering
SNR_MIN: float = 5.0  # drop features with S/N average below this (if available)

# Isomer labeling
ISOMER_RT_WINDOW_MIN: float = 1.0  # minutes for RT clustering within metabolite+adduct

# Replicate-group QC thresholds
BLANK_FOLD_MIN: float = 7.0      # max(reps) / blank must be >= this
PRESENT_PERCENT_MIN: float = 60.0  # percent of replicates detected must be >= this
CV_PERCENT_MAX: float = 40.0       # CV among present replicates must be <= this

