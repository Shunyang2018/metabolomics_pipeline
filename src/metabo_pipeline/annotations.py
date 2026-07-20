from __future__ import annotations

import pandas as pd


def _to_float(val):
    """Coerce a value to float."""
    return float(val)


def assign_annotation_level_row(row: pd.Series) -> str:
    """Infer an annotation level for a single MS-DIAL row."""
    name = str(row.get("Metabolite name", "")).strip()
    lname = name.lower()
    if (
        (not name)
        or lname.startswith("unknown")
        or lname.startswith("low score")
        or lname.startswith("no ms2")
    ):
        return "3"

    wdot = _to_float(row["Weighted dot product"])
    rdot = _to_float(row["Reverse dot product"])
    mcount = int(round(float(row["Matched peaks count"])))

    # Weak or missing weighted score → Level 3
    if not pd.notna(wdot) or wdot < 0.5:
        return "3"
    # If there are no matched peaks reported by MS-DIAL, treat as Level 3
    if mcount is None or mcount < 1:
        return "3"

    # Level 1 rule (scores assumed on 0–1 scale) — require sufficient matched peaks
    if pd.notna(rdot) and (wdot > 0.75) and (abs(wdot - rdot) < 0.2) and (mcount >= 3):
        return "1"
    # Level 2 requires Reverse dot > 0.5 and matched peaks >= 3
    if pd.notna(rdot) and rdot > 0.5 and mcount >= 3:
        return "2"
    return "3"


def annotate_levels(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate an entire dataframe with MS-DIAL level calls."""
    df = df.copy()
    # Avoid pandas apply edge-case on empty frames (can return empty DataFrame)
    if df.empty:
        df["annotation_level"] = pd.Series([], dtype=str)
        return df
    df["annotation_level"] = df.apply(assign_annotation_level_row, axis=1)
    return df
