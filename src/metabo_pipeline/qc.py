from __future__ import annotations

from typing import Dict, List

import numpy as np
import pandas as pd


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
        try:
            inten = float(tok.split(":", 1)[1])
        except ValueError:
            continue
        if inten > 0:
            cnt += 1
    return cnt


# Tokens that mark a column as belonging to the blank baseline rather than a
# real sample under test: explicit blanks, method blanks ('mb'), and any
# column designated a blank-equivalent (e.g. a resuspension/solvent control).
_BLANK_LIKE_TOKENS = {"blank", "mb", "resuspension"}


def build_group_cols(
    sample_cols: List[str], real_sample_tokens: "set[str] | None" = None
) -> Dict[str, List[str]]:
    """Group real-sample columns for blank-fold QC.

    Each distinct normalized sample name is its own group (a trailing number
    is treated as part of the sample's identity, e.g. `crude_1_amide` and
    `crude_2_amide` are different samples, not replicates of one condition —
    grouping never assumes otherwise). Blank-like columns (see
    `_BLANK_LIKE_TOKENS`: blank, method-blank, resuspension/solvent control)
    never form a group — they feed the blank baseline instead.

    `real_sample_tokens`, when given (see `constants.REAL_SAMPLE_TOKENS`), is
    a whitelist: only columns containing at least one of these tokens count
    as a real sample group. Everything else — identification/reference runs,
    QC pool injections, MS1-only injections, or any other column that isn't
    a recognized real-sample or blank-like type — is excluded from QC
    entirely rather than silently forming its own group. When `None`, every
    non-blank-like column is treated as a real sample group (the permissive
    default, for datasets that haven't defined an explicit whitelist).
    """
    import re

    group_cols: Dict[str, List[str]] = {}
    for c in sample_cols:
        tokens = [t for t in re.split(r"[^a-z0-9]+", c.lower()) if t]
        if not tokens:
            continue
        if _BLANK_LIKE_TOKENS.intersection(tokens):
            continue
        if real_sample_tokens is not None and not real_sample_tokens.intersection(tokens):
            continue
        group_cols.setdefault(c.lower(), []).append(c)
    return group_cols


def compute_group_metrics(
    df: pd.DataFrame, group_cols: Dict[str, List[str]], blank_col: str | None
) -> pd.DataFrame:
    """Compute blank fold, presence, and CV metrics per replicate group."""
    df = df.copy()

    # Gather blank columns: explicit 'blank', plus any column tokenized as
    # blank-like (see `_BLANK_LIKE_TOKENS` — blank, method-blank, resuspension/
    # solvent controls).
    blank_cols: List[str] = []
    if blank_col is not None and blank_col in df.columns:
        blank_cols.append(blank_col)
    for col in df.columns:
        if col == blank_col:
            continue
        tokens = [tok for tok in str(col).lower().replace("-", "_").split("_") if tok]
        if _BLANK_LIKE_TOKENS.intersection(tokens):
            blank_cols.append(col)

    blank_values = None
    if blank_cols:
        blank_values = df[blank_cols].apply(pd.to_numeric, errors="coerce")
        blank_avg = blank_values.mean(axis=1, skipna=True)
        blank_denom = blank_avg.replace(0, np.nan)
    else:
        blank_avg = pd.Series([float("nan")] * len(df), index=df.index)
        blank_denom = blank_avg

    sample_cols_all: List[str] = sorted(
        {c for cols in group_cols.values() for c in cols}
    )
    if sample_cols_all:
        sample_vals = df[sample_cols_all].apply(pd.to_numeric, errors="coerce")
        max_all_samples = sample_vals.max(axis=1, skipna=True)
    else:
        max_all_samples = pd.Series([float("nan")] * len(df), index=df.index)

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

        bf_series = (
            pd.to_numeric(bf, errors="coerce")
            if bf is not None
            else pd.Series(np.nan, index=df.index)
        )
        pp_series = (
            pd.to_numeric(pp, errors="coerce")
            if pp is not None
            else pd.Series(np.nan, index=df.index)
        )
        mask = (bf_series >= blank_fold_min) & (pp_series >= present_min)
        if cv_max is not None and cv is not None:
            cv_series = pd.to_numeric(cv, errors="coerce")
            mask = mask & (cv_series <= cv_max)
        mask = mask.fillna(False)
        passes.append(mask)
    if not passes:
        import warnings

        warnings.warn(
            "No replicate groups were recognized for blank-fold QC; "
            "every row will pass unfiltered. Check your sample column naming.",
            stacklevel=2,
        )
        return pd.Series([True] * len(df), index=df.index)
    m = passes[0]
    for other in passes[1:]:
        m = m | other
    return m
