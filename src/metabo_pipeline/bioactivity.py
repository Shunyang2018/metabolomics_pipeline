"""Match MS-DIAL/ClassyFire features against a bioactivity reference database.

Matching is done on the first block of the InChIKey (the 14-character
skeleton hash, before the first hyphen), which encodes molecular connectivity
only. Matching on this block (rather than the full InChIKey) still catches
bioactivity hits for a feature whose stereochemistry, protonation state, or
isotope layer differs slightly from the reference entry.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

# MS-DIAL/pipeline columns that are per-feature metadata, never per-sample data.
_KNOWN_NON_SAMPLE_COLS = {
    "feature_id",
    "chrom",
    "annotation_level",
    "mode",
    "Alignment ID",
    "Average Rt(min)",
    "Average Mz",
    "Metabolite name",
    "Adduct type",
    "Post curation result",
    "Fill %",
    "MS/MS assigned",
    "Reference RT",
    "Reference m/z",
    "Formula",
    "Ontology",
    "INCHIKEY",
    "SMILES",
    "Annotation tag (VS1.0)",
    "RT matched",
    "m/z matched",
    "MS/MS matched",
    "Comment",
    "Manually modified for quantification",
    "Manually modified for annotation",
    "Isotope tracking parent ID",
    "Isotope tracking weight number",
    "RT similarity",
    "m/z similarity",
    "Simple dot product",
    "Weighted dot product",
    "Reverse dot product",
    "Matched peaks count",
    "Matched peaks percentage",
    "Total score",
    "S/N average",
    "Spectrum reference file name",
    "MS1 isotopic spectrum",
    "MS/MS spectrum",
    "cf_kingdom",
    "cf_superclass",
    "cf_class",
    "cf_subclass",
    "cf_direct_parent",
}
_METRIC_PREFIXES = ("blank_fold_", "present_percent_", "cv_percent_", "pass_", "SIRIUS_")
_METRIC_EXACT = {"pass_any_groups", "_polarity", "_rt", "_mz", "_sn", "_wd"}


def inchikey_block1(value) -> Optional[str]:
    """Return the first (skeleton) block of an InChIKey, or None if blank/invalid."""
    s = str(value or "").strip().upper()
    if not s or s == "NAN":
        return None
    return s.split("-")[0] or None


def identify_sample_columns(columns: List[str]) -> List[str]:
    """Return the subset of columns that hold per-sample intensities."""
    return [
        c
        for c in columns
        if c not in _KNOWN_NON_SAMPLE_COLS
        and c not in _METRIC_EXACT
        and not str(c).startswith(_METRIC_PREFIXES)
    ]


def match_bioactives(
    merged_csv: Path,
    bioactivity_db_csv: Path,
    output_csv: Path,
    inchikey_col: str = "INCHIKEY",
    db_inchikey_col: str = "InChIKey",
) -> Dict[str, int]:
    """Match msdial features to a bioactivity database on InChIKey skeleton.

    Emits one row per (feature, sample it was detected in, bioactivity hit).
    """
    feats = pd.read_csv(merged_csv, low_memory=False)
    db = pd.read_csv(bioactivity_db_csv, low_memory=False)

    if inchikey_col not in feats.columns:
        raise ValueError(f"'{inchikey_col}' column not found in {merged_csv}")
    if db_inchikey_col not in db.columns:
        raise ValueError(f"'{db_inchikey_col}' column not found in {bioactivity_db_csv}")

    feats = feats.copy()
    db = db.copy()
    feats["_ikey_block1"] = feats[inchikey_col].apply(inchikey_block1)
    db["_ikey_block1"] = db[db_inchikey_col].apply(inchikey_block1)

    matched_feats = feats[feats["_ikey_block1"].notna()].copy()
    matched_db = db[db["_ikey_block1"].notna()].copy()

    sample_cols = identify_sample_columns(list(feats.columns))

    id_cols = [
        c
        for c in (
            "feature_id",
            "chrom",
            "annotation_level",
            "Metabolite name",
            inchikey_col,
            "Formula",
            "Adduct type",
            "Average Rt(min)",
            "Average Mz",
        )
        if c in matched_feats.columns
    ]

    long_df = matched_feats.melt(
        id_vars=id_cols + ["_ikey_block1"],
        value_vars=sample_cols,
        var_name="sample_type",
        value_name="intensity",
    )
    long_df["intensity"] = pd.to_numeric(long_df["intensity"], errors="coerce")
    long_df = long_df[long_df["intensity"].fillna(0) > 0]

    out = long_df.merge(
        matched_db, on="_ikey_block1", how="inner", suffixes=("", "_bioactivity")
    )
    out = out.drop(columns=["_ikey_block1"])

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_csv, index=False)

    return {
        "features_with_inchikey": int(feats[inchikey_col].notna().sum()),
        "unique_matched_skeletons": int(matched_feats["_ikey_block1"].nunique()),
        "bioactivity_db_rows": int(len(matched_db)),
        "output_rows": int(len(out)),
        "unique_features_matched": (
            int(out["feature_id"].nunique()) if "feature_id" in out.columns else 0
        ),
    }
