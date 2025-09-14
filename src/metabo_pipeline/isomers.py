from __future__ import annotations

import pandas as pd


def assign_isomer_labels_global(df: pd.DataFrame, rt_window_min: float) -> pd.DataFrame:
    """Assign isomer_label globally by clustering on RT within (Metabolite name, Adduct type).
    Leaves label blank for groups with only one cluster.
    """
    if not all(c in df.columns for c in ("Metabolite name", "Adduct type", "Average Rt(min)")):
        return df

    df = df.copy()
    df["_rt_num"] = pd.to_numeric(df["Average Rt(min)"], errors="coerce")

    def _label_group(g: pd.DataFrame) -> pd.DataFrame:
        g = g.sort_values("_rt_num")
        rts = g["_rt_num"].to_numpy()
        labels = []
        cluster = 1
        prev = None
        for rt in rts:
            if prev is None or pd.isna(prev) or pd.isna(rt) or abs(rt - prev) > rt_window_min:
                cluster = cluster if prev is None else cluster + 1
            labels.append(f"isomer_{cluster}")
            prev = rt
        out = g.copy()
        out["isomer_label"] = labels
        return out

    try:
        df = df.groupby(["Metabolite name", "Adduct type"], dropna=False, group_keys=False).apply(
            _label_group, include_groups=False
        )
    except TypeError:
        df = df.groupby(["Metabolite name", "Adduct type"], dropna=False, group_keys=False).apply(_label_group)

    # Blank singles reliably by using a filled key for grouping
    try:
        name_f = df["Metabolite name"].astype(str).fillna("(NA)")
        adduct_f = df["Adduct type"].astype(str).fillna("(NA)")
        key = name_f + "||" + adduct_f
        sizes = key.map(key.value_counts())
        df.loc[(sizes <= 1) | sizes.isna(), "isomer_label"] = ""
    except Exception:
        pass

    return df.drop(columns=["_rt_num"], errors="ignore")

