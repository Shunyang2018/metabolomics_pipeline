from __future__ import annotations

import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Callable

from .constants import (
    MSMS_MIN_IONS,
    SNR_MIN,
    ISOMER_RT_WINDOW_MIN,
    DEDUP_RT_WINDOW_MIN,
    DEDUP_MZ_PPM,
    BLANK_FOLD_MIN,
    PRESENT_PERCENT_MIN,
    CV_PERCENT_MAX,
)
from .utils import infer_chrom_from_name, infer_mode_from_name, normalize_sample_id_core
from .annotations import annotate_levels
from .qc import count_msms_ions, build_group_cols, compute_group_metrics, pass_any_mask
from .dedup import l3_representatives, dedup_name_conflicts_by_cluster
from .sirius_export import build_ms_entries, write_ms_files


def list_alignment_files(input_dir: Path, recursive: bool = False) -> List[Path]:
    pats = ["*.csv", "*.txt"]
    files: List[Path] = []
    if recursive:
        for p in pats:
            files.extend(input_dir.rglob(p))
    else:
        for p in pats:
            files.extend(input_dir.glob(p))
    return sorted({p.resolve() for p in files})


def merge_folder_to_wide_csv(
    input_dir: Path,
    output_csv: Path,
    recursive: bool = False,
    progress: Optional[Callable[[Dict[str, int]], None]] = None,
) -> Dict[str, int]:
    files = list_alignment_files(input_dir, recursive=recursive)
    totals = {"raw": 0, "after_msms": 0, "after_snr": 0, "after_pass": 0, "ann1": 0, "ann2": 0, "ann3": 0}
    per_file: List[Dict[str, int]] = []

    frames: List[pd.DataFrame] = []
    for p in files:
        # Read with header row at line 5 (0-based index 4)
        try:
            df = pd.read_csv(p, header=4, encoding="utf-8-sig")
        except Exception:
            df = pd.read_csv(p, header=4, encoding="utf-8-sig", engine="python")

        cols = list(df.columns)
        last_meta_col = cols.index("MS/MS spectrum") if "MS/MS spectrum" in cols else max(30, cols.index("Alignment ID") + 1 if "Alignment ID" in cols else 30)
        sample_cols_raw = cols[last_meta_col + 1 :]

        raw_cnt = int(df.shape[0])
        df["_ion_cnt"] = df.get("MS/MS spectrum", "").map(count_msms_ions)
        df = df[df["_ion_cnt"] >= MSMS_MIN_IONS].copy()
        after_msms = int(df.shape[0])
        if "S/N average" in df.columns:
            sn = pd.to_numeric(df["S/N average"], errors="coerce")
            df = df[sn >= SNR_MIN].copy()
        after_snr = int(df.shape[0])

        # Annotation level
        df = annotate_levels(df)

        # Normalize sample column names and keep metadata + samples
        rename_samples = {c: normalize_sample_id_core(c) for c in sample_cols_raw}
        df = df.rename(columns=rename_samples)
        feat_df = df[cols[: last_meta_col + 1]].copy()
        feat_df["annotation_level"] = df["annotation_level"]
        norm_sample_cols = list(rename_samples.values())
        if norm_sample_cols:
            feat_df = pd.concat([feat_df, df[norm_sample_cols]], axis=1)

        # chrom/mode from filename
        feat_df.insert(0, "chrom", infer_chrom_from_name(p.name))
        feat_df.insert(1, "mode", infer_mode_from_name(p.name))

        # Per-group QC metrics and gating
        blank_col = "blank" if "blank" in feat_df.columns else None
        group_cols = build_group_cols(norm_sample_cols)
        feat_df = compute_group_metrics(feat_df, group_cols, blank_col)
        mask = pass_any_mask(feat_df, group_cols, BLANK_FOLD_MIN, PRESENT_PERCENT_MIN, CV_PERCENT_MAX)
        feat_df = feat_df[mask].copy()

        # Per-file stats
        ann = feat_df.get("annotation_level")
        a1 = int((ann == "1").sum()) if ann is not None else 0
        a2 = int((ann == "2").sum()) if ann is not None else 0
        a3 = int((ann == "3").sum()) if ann is not None else 0
        rec = {"file": p.name, "raw": raw_cnt, "after_msms": after_msms, "after_snr": after_snr, "after_pass": int(feat_df.shape[0]), "ann1": a1, "ann2": a2, "ann3": a3}
        per_file.append(rec)
        totals["raw"] += raw_cnt; totals["after_msms"] += after_msms; totals["after_snr"] += after_snr; totals["after_pass"] += int(feat_df.shape[0]); totals["ann1"] += a1; totals["ann2"] += a2; totals["ann3"] += a3
        if progress:
            try:
                progress(rec)
            except Exception:
                pass

        frames.append(feat_df)

    # Union frames and order columns
    if frames:
        merged = pd.concat(frames, ignore_index=True)
    else:
        merged = pd.DataFrame()

    # Isomer labeling removed per request (no isomer_label column)

    # Column ordering: originals first
    id_preferred = [
        "chrom",
        "annotation_level",
        "Alignment ID",
        "Average Rt(min)",
        "Average Mz",
        "Metabolite name",
        "Adduct type",
        "Formula",
        "Ontology",
        "INCHIKEY",
        "SMILES",
        "Annotation tag (VS1.0)",
        "MS/MS assigned",
        "RT matched",
        "m/z matched",
        "MS/MS matched",
        "Weighted dot product",
        "Reverse dot product",
        "Simple dot product",
        "Matched peaks count",
        "Matched peaks percentage",
        "Total score",
        "S/N average",
        "Reference RT",
        "Reference m/z",
        "Spectrum reference file name",
        "Comment",
    ]
    all_cols = set(merged.columns)
    id_order = [c for c in id_preferred if c in all_cols]
    import re
    sample_pat = re.compile(r"^m2_[a-z0-9]+_.+")
    metric_cols = sorted(c for c in all_cols if str(c).startswith(("blank_fold_", "present_percent_", "cv_percent_")))
    pass_cols = sorted(c for c in all_cols if str(c).startswith("pass_") or c in ("pass_any_groups", "mode", "_polarity", "_rt", "_mz", "_sn", "_wd"))
    sample_cols = sorted(c for c in all_cols if c not in set(id_order) and c not in set(metric_cols) and c not in set(pass_cols) and sample_pat.match(str(c)))
    other_cols = sorted(c for c in all_cols if c not in set(id_order + sample_cols + metric_cols + pass_cols))
    out_cols = id_order + other_cols + sample_cols + metric_cols
    if out_cols:
        merged = merged.reindex(columns=out_cols)

    # Deduplicate name conflicts within RT/m/z clusters (all rows)
    merged = dedup_name_conflicts_by_cluster(merged, DEDUP_RT_WINDOW_MIN, DEDUP_MZ_PPM)

    # L3 representatives for merged CSV
    rows_pre = int(merged.shape[0])
    if "annotation_level" in merged.columns and not merged.empty:
        l3 = merged[merged["annotation_level"] == "3"].copy()
        non_l3 = merged[merged["annotation_level"] != "3"].copy()
        if not l3.empty:
            l3_keep = l3_representatives(l3, DEDUP_RT_WINDOW_MIN, DEDUP_MZ_PPM)
            merged = pd.concat([non_l3, l3_keep], axis=0, ignore_index=True)
    rows_post = int(merged.shape[0])

    # SIRIUS export from merged L3 representatives
    l3_all = merged[merged["annotation_level"] == "3"].copy() if "annotation_level" in merged.columns else pd.DataFrame()
    sirius_l3_total = int(l3_all.shape[0]) if not l3_all.empty else 0
    pos_entries, neg_entries = build_ms_entries(l3_all)
    sp, sn = write_ms_files(pos_entries, neg_entries, output_csv.parent)

    # Write CSV
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_csv, index=False)

    return {
        "files": len(files),
        "rows": int(len(merged)),
        "rows_pre_dedup": rows_pre,
        "dedup_dropped": rows_pre - rows_post,
        "totals": totals,
        "per_file": per_file,
        "ann_post": {
            "1": int((merged["annotation_level"] == "1").sum()) if "annotation_level" in merged.columns else 0,
            "2": int((merged["annotation_level"] == "2").sum()) if "annotation_level" in merged.columns else 0,
            "3": int((merged["annotation_level"] == "3").sum()) if "annotation_level" in merged.columns else 0,
        },
        "isomer_post": {
            "metabolites_with_isomers": 0,  # can be recomputed if needed
            "total_isomers": 0,
            "max_isomers_per_metabolite": 0,
        },
        "sirius_pos_count": sp,
        "sirius_neg_count": sn,
        "sirius_l3_total": sirius_l3_total,
    }
