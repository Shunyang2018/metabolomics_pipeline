from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple


def count_msms_ions(msms: str) -> int:
    if not isinstance(msms, str):
        return 0
    s = msms.strip().lower()
    if s in {"", "null", "na", "none"}:
        return 0
    cnt = 0
    for tok in s.split():
        if ":" not in tok:
            continue
        try:
            inten = float(tok.split(":", 1)[1])
        except Exception:
            continue
        if inten > 0:
            cnt += 1
    return cnt


def build_group_cols(sample_cols: List[str]) -> Dict[str, List[str]]:
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
    df = df.copy()
    if blank_col is not None and blank_col in df.columns:
        blank_series = pd.to_numeric(df[blank_col], errors="coerce")
    else:
        blank_series = None

    for grp, cols_grp in group_cols.items():
        vals = df[cols_grp].apply(pd.to_numeric, errors="coerce")
        nrep = max(1, len(cols_grp))
        present = vals.gt(0)
        present_frac = present.sum(axis=1) / float(nrep)
        present_percent = (present_frac * 100.0).astype(float)

        max_val = vals.max(axis=1, skipna=True)
        if blank_series is not None:
            denom = blank_series.replace(0, np.nan)
            blank_fold = (max_val / denom)
        else:
            blank_fold = pd.Series([float('nan')]*len(df), index=df.index)

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


def pass_any_mask(df: pd.DataFrame, group_cols: Dict[str, List[str]], blank_fold_min: float, present_min: float, cv_max: float) -> pd.Series:
    passes = []
    for grp in group_cols.keys():
        bf = df.get(f"blank_fold_{grp}")
        pp = df.get(f"present_percent_{grp}")
        cv = df.get(f"cv_percent_{grp}")
        mask = (bf.astype(float) >= blank_fold_min) & (pp.astype(float) >= present_min) & (cv.astype(float) <= cv_max)
        mask = mask.fillna(False)
        passes.append(mask)
    if not passes:
        return pd.Series([True] * len(df), index=df.index)
    m = passes[0]
    for other in passes[1:]:
        m = m | other
    return m

