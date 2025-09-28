from __future__ import annotations

import numpy as np
import pandas as pd
from typing import List, Tuple


def l3_representatives(df_l3: pd.DataFrame, rt_window_min: float, mz_ppm: float) -> pd.DataFrame:
    """Return representative L3 rows using global clustering within (chrom, polarity).

    - Window: |ΔRT| <= rt_window_min; |Δm/z| <= mz_ppm (ppm) relative to cluster representative
    - Score: pick higher S/N average, then higher Weighted dot product
    """
    if df_l3.empty:
        return df_l3

    work = df_l3.copy()
    # Prepare helper fields
    if "Adduct type" in work.columns:
        work["_adduct_export"] = work["Adduct type"].astype(str)
    else:
        work["_adduct_export"] = ""
    work["_polarity"] = work["_adduct_export"].map(lambda s: "POS" if "+" in str(s or "") else ("NEG" if "-" in str(s or "") else "UNK"))
    if "mode" in work.columns:
        work.loc[work["_polarity"] == "UNK", "_polarity"] = work.loc[work["_polarity"] == "UNK", "mode"].fillna("UNK")
    work["_rt"] = pd.to_numeric(work.get("Average Rt(min)"), errors="coerce")
    work["_mz"] = pd.to_numeric(work.get("Average Mz"), errors="coerce")
    work["_sn"] = pd.to_numeric(work.get("S/N average"), errors="coerce")
    work["_wd"] = pd.to_numeric(work.get("Weighted dot product"), errors="coerce")

    keep_idx: List[int] = []
    for (chrom, pol), g in work.groupby(["chrom", "_polarity"], dropna=False):
        if g.empty:
            continue
        g = g.sort_values(["_rt", "_mz"])  # deterministic
        # Single-linkage: compare to last member in the current cluster
        cluster_rep: Tuple[int, float, float, float, float] | None = None
        last_rt = np.nan
        last_mz = np.nan
        for idx, row in g.iterrows():
            rt = row["_rt"]; mz = row["_mz"]
            sn = row["_sn"] if pd.notna(row["_sn"]) else -1.0
            wd = row["_wd"] if pd.notna(row["_wd"]) else -1.0
            if cluster_rep is None:
                cluster_rep = (idx, rt, mz, sn, wd)
                last_rt, last_mz = rt, mz
                continue
            _, rep_rt, rep_mz, rep_sn, rep_wd = cluster_rep
            mz_tol = (mz_ppm * 1e-6) * last_mz if pd.notna(last_mz) else np.nan
            if (pd.notna(rt) and pd.notna(last_rt) and abs(rt - last_rt) <= rt_window_min) and \
               (pd.notna(mz) and pd.notna(last_mz) and abs(mz - last_mz) <= mz_tol):
                if (sn, wd) > (rep_sn, rep_wd):
                    cluster_rep = (idx, rt, mz, sn, wd)
                # extend cluster
                last_rt, last_mz = rt, mz
            else:
                keep_idx.append(cluster_rep[0])
                cluster_rep = (idx, rt, mz, sn, wd)
                last_rt, last_mz = rt, mz
        if cluster_rep is not None:
            keep_idx.append(cluster_rep[0])

    reps = work.loc[sorted(set(keep_idx))]
    # Clean helpers
    reps = reps.drop(columns=["_adduct_export", "_polarity", "_rt", "_mz", "_sn", "_wd"], errors="ignore")
    return reps


def dedup_name_conflicts_by_cluster(df: pd.DataFrame, rt_window_min: float, mz_ppm: float) -> pd.DataFrame:
    """Within RT/m/z clusters (per chrom, polarity), keep one row per Metabolite name
    choosing the lowest CV (median of cv_percent_*). If CV is missing, fallback to
    higher S/N average, then higher Weighted dot.
    """
    if df.empty:
        return df

    work = df.copy()
    # Prepare helper fields
    if "Adduct type" in work.columns:
        work["_adduct_export"] = work["Adduct type"].astype(str)
    else:
        work["_adduct_export"] = ""
    work["_polarity"] = work["_adduct_export"].map(lambda s: "POS" if "+" in str(s or "") else ("NEG" if "-" in str(s or "") else "UNK"))
    if "mode" in work.columns:
        work.loc[work["_polarity"] == "UNK", "_polarity"] = work.loc[work["_polarity"] == "UNK", "mode"].fillna("UNK")
    work["_rt"] = pd.to_numeric(work.get("Average Rt(min)"), errors="coerce")
    work["_mz"] = pd.to_numeric(work.get("Average Mz"), errors="coerce")
    work["_sn"] = pd.to_numeric(work.get("S/N average"), errors="coerce")
    work["_wd"] = pd.to_numeric(work.get("Weighted dot product"), errors="coerce")

    # CV median from available cv_percent_* columns
    cv_cols = [c for c in work.columns if str(c).startswith("cv_percent_")]
    if cv_cols:
        work["_cv_med"] = work[cv_cols].astype(float).median(axis=1, skipna=True)
    else:
        work["_cv_med"] = np.nan

    picked_idx: List[int] = []
    for (chrom, pol), g in work.groupby(["chrom", "_polarity"], dropna=False):
        if g.empty:
            continue
        g = g.sort_values(["_rt", "_mz"])  # deterministic
        cluster_members: List[int] = []
        # Single-linkage: compare to last member
        last_rt = np.nan
        last_mz = np.nan
        for idx, row in g.iterrows():
            rt = row["_rt"]; mz = row["_mz"]
            if isinstance(last_rt, float) and np.isnan(last_rt):
                # start cluster
                cluster_members = [idx]
                last_rt, last_mz = rt, mz
                continue
            mz_tol = (mz_ppm * 1e-6) * last_mz if pd.notna(last_mz) else np.nan
            if (pd.notna(rt) and pd.notna(last_rt) and abs(rt - last_rt) <= rt_window_min) and \
               (pd.notna(mz) and pd.notna(last_mz) and abs(mz - last_mz) <= mz_tol):
                cluster_members.append(idx)
            else:
                # finalize previous cluster: pick one per Metabolite name
                sub = work.loc[cluster_members]
                for name, subg in sub.groupby("Metabolite name", dropna=False):
                    subg = subg.copy()
                    subg["_cv_sort"] = subg["_cv_med"].fillna(np.inf)
                    subg = subg.sort_values(["_cv_sort", "_sn", "_wd"], ascending=[True, False, False])
                    picked_idx.append(subg.index[0])
                # start new cluster
                cluster_members = [idx]
                last_rt, last_mz = rt, mz
        # finalize last cluster
        if cluster_members:
            sub = work.loc[cluster_members]
            for name, subg in sub.groupby("Metabolite name", dropna=False):
                subg = subg.copy()
                subg["_cv_sort"] = subg["_cv_med"].fillna(np.inf)
                subg = subg.sort_values(["_cv_sort", "_sn", "_wd"], ascending=[True, False, False])
                picked_idx.append(subg.index[0])

    reps = work.loc[sorted(set(picked_idx))]
    # Clean helpers
    reps = reps.drop(columns=["_adduct_export", "_polarity", "_rt", "_mz", "_sn", "_wd", "_cv_med", "_cv_sort"], errors="ignore")
    return reps
