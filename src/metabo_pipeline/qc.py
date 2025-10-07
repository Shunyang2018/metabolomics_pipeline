from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple


def count_msms_ions(msms: str) -> int:
    """Count nonzero fragment ions encoded in an MS/MS spectrum cell."""
    if not isinstance(msms, str):
        return 0
    s = msms.strip().lower()
    if s in {"", "null", "na", "none"}:
        return 0
    cnt = 0
    for tok in s.split():
        if ":" not in tok:
            continue
        inten = float(tok.split(":", 1)[1])
        if inten > 0:
            cnt += 1
    return cnt


def build_group_cols(sample_cols: List[str]) -> Dict[str, List[str]]:
    """Group replicate sample columns by their normalized identifier."""
    import re

    group_cols: Dict[str, List[str]] = {}
    pat = re.compile(r"^m2_([a-z0-9]+)_(.+)$")
    for c in sample_cols:
        if c == "blank":
            continue
        m = pat.match(c)
        if not m:
            continue
        grp = m.group(1)
        if grp in ("pool", "qc"):
            continue
        group_cols.setdefault(grp, []).append(c)
    return group_cols


def compute_group_metrics(df: pd.DataFrame, group_cols: Dict[str, List[str]], blank_col: str | None) -> pd.DataFrame:
    """Compute blank fold, presence, and CV metrics per replicate group."""
    df = df.copy()

    # Gather blank columns (explicit 'blank' plus any column whose token equals 'blank').
    blank_cols: List[str] = []
    if blank_col is not None and blank_col in df.columns:
        blank_cols.append(blank_col)
    for col in df.columns:
        if col == blank_col:
            continue
        tokens = [tok for tok in str(col).lower().replace('-', '_').split('_') if tok]
        if "blank" in tokens:
            blank_cols.append(col)

    blank_values = None
    if blank_cols:
        blank_values = df[blank_cols].apply(pd.to_numeric, errors="coerce")
        blank_avg = blank_values.mean(axis=1, skipna=True)
        blank_denom = blank_avg.replace(0, np.nan)
    else:
        blank_avg = pd.Series([float('nan')] * len(df), index=df.index)
        blank_denom = blank_avg

    sample_cols_all: List[str] = sorted({c for cols in group_cols.values() for c in cols})
    if sample_cols_all:
        sample_vals = df[sample_cols_all].apply(pd.to_numeric, errors="coerce")
        max_all_samples = sample_vals.max(axis=1, skipna=True)
    else:
        max_all_samples = pd.Series([float('nan')] * len(df), index=df.index)

    for grp, cols_grp in group_cols.items():
        vals = df[cols_grp].apply(pd.to_numeric, errors="coerce")
        nrep = max(1, len(cols_grp))
        present = vals.gt(0)
        present_frac = present.sum(axis=1) / float(nrep)
        present_percent = (present_frac * 100.0).astype(float)

        blank_fold = max_all_samples / blank_denom

        vals_present = vals.where(present)
        mean = vals_present.mean(axis=1, skipna=True)
        std = vals_present.std(axis=1, ddof=1, skipna=True)
        denom_mean = mean.replace(0, np.nan)
        cv_percent = (std / denom_mean) * 100.0
        cv_percent = cv_percent.where(present.sum(axis=1) >= 2)

        df[f"blank_fold_{grp}"] = blank_fold
        df[f"present_percent_{grp}"] = present_percent
        df[f"cv_percent_{grp}"] = cv_percent
    return df


def pass_any_mask(
    df: pd.DataFrame,
    group_cols: Dict[str, List[str]],
    blank_fold_min: float,
    present_min: float,
    cv_max: float | None,
) -> pd.Series:
    """Determine whether each row passes any replicate-group QC criteria."""
    passes = []
    for grp in group_cols.keys():
        bf = df.get(f"blank_fold_{grp}")
        pp = df.get(f"present_percent_{grp}")
        cv = df.get(f"cv_percent_{grp}")

        bf_series = pd.to_numeric(bf, errors="coerce") if bf is not None else pd.Series(np.nan, index=df.index)
        pp_series = pd.to_numeric(pp, errors="coerce") if pp is not None else pd.Series(np.nan, index=df.index)
        mask = (bf_series >= blank_fold_min) & (pp_series >= present_min)
        if cv_max is not None and cv is not None:
            cv_series = pd.to_numeric(cv, errors="coerce")
            mask = mask & (cv_series <= cv_max)
        mask = mask.fillna(False)
        passes.append(mask)
    if not passes:
        return pd.Series([True] * len(df), index=df.index)
    m = passes[0]
    for other in passes[1:]:
        m = m | other
    return m

