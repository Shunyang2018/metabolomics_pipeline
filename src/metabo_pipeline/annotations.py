from __future__ import annotations

import pandas as pd


def assign_annotation_level_row(row: pd.Series) -> str:
    name = str(row.get("Metabolite name", "")).strip()
    if (not name) or name.lower().startswith("unknown"):
        return "3"

    def _to_float(val):
        try:
            return float(val)
        except Exception:
            return float("nan")

    wdot = _to_float(row.get("Weighted dot product"))
    rdot = _to_float(row.get("Reverse dot product"))
    try:
        mcount = int(round(float(row.get("Matched peaks count", float("nan")))))
    except Exception:
        mcount = -1

    # Weak or missing weighted score → Level 3
    if not pd.notna(wdot) or wdot < 0.5:
        return "3"
    # Level 1 rule (scores assumed on 0–1 scale)
    if pd.notna(rdot) and (wdot > 0.75) and (abs(wdot - rdot) < 0.2) and (mcount >= 3):
        return "1"
    # Level 2 now additionally requires Reverse dot > 0.5; otherwise treat as Level 3
    if pd.notna(rdot) and rdot > 0.5:
        return "2"
    return "3"


def annotate_levels(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["annotation_level"] = df.apply(assign_annotation_level_row, axis=1)
    return df
