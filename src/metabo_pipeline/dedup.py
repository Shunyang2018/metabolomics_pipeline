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
        cluster_rep: Tuple[int, float, float, float, float] | None = None
        for idx, row in g.iterrows():
            rt = row["_rt"]; mz = row["_mz"]
            sn = row["_sn"] if pd.notna(row["_sn"]) else -1.0
            wd = row["_wd"] if pd.notna(row["_wd"]) else -1.0
            if cluster_rep is None:
                cluster_rep = (idx, rt, mz, sn, wd)
                continue
            _, rep_rt, rep_mz, rep_sn, rep_wd = cluster_rep
            mz_tol = (mz_ppm * 1e-6) * rep_mz if pd.notna(rep_mz) else np.nan
            if (pd.notna(rt) and pd.notna(rep_rt) and abs(rt - rep_rt) <= rt_window_min) and \
               (pd.notna(mz) and pd.notna(rep_mz) and abs(mz - rep_mz) <= mz_tol):
                if (sn, wd) > (rep_sn, rep_wd):
                    cluster_rep = (idx, rt, mz, sn, wd)
            else:
                keep_idx.append(cluster_rep[0])
                cluster_rep = (idx, rt, mz, sn, wd)
        if cluster_rep is not None:
            keep_idx.append(cluster_rep[0])

    reps = work.loc[sorted(set(keep_idx))]
    # Clean helpers
    reps = reps.drop(columns=["_adduct_export", "_polarity", "_rt", "_mz", "_sn", "_wd"], errors="ignore")
    return reps

