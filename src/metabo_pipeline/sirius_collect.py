from __future__ import annotations

import csv
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import pandas as pd


def _read_text(path: Path) -> str:
    """Return file contents as UTF-8 text."""
    return path.read_text(encoding="utf-8", errors="ignore")


def _parse_feature_id_from_spectrum(spectrum_ms: Path) -> Optional[int]:
    """Extract a feature_id integer from a SIRIUS spectrum file."""
    txt = _read_text(spectrum_ms)
    fid: Optional[int] = None
    for line in txt.splitlines():
        if line.startswith("#feature_id "):
            return int(line.split()[1])
        if line.startswith(">compound "):
            tok = line.split(maxsplit=1)[1].strip()
            if tok.isdigit():
                fid = int(tok)
    return fid


def _parse_info_from_compound(compound_info: Path) -> Dict[str, Optional[str]]:
    """Parse ion metadata from a compound.info file."""
    out: Dict[str, Optional[str]] = {"ionType": None, "ionMass": None}
    txt = _read_text(compound_info)
    for line in txt.splitlines():
        if "\t" not in line:
            continue
        k, v = line.split("\t", 1)
        if k == "ionType":
            out["ionType"] = v.strip()
        elif k == "ionMass":
            out["ionMass"] = v.strip()
    return out


def _extract_canopus(comp_dir: Path) -> Dict[str, Optional[str]]:
    """Gather CANOPUS class annotations if present in a compound directory."""
    out = {"SIRIUS_canopus_superclass": None, "SIRIUS_canopus_class": None, "SIRIUS_canopus_most_specific": None}
    can_dir = comp_dir / "canopus"
    if not can_dir.exists():
        return out
    # Heuristic: find a TSV under canopus with class columns
    for tsv in can_dir.glob("*.tsv"):
        with tsv.open("r", encoding="utf-8", errors="ignore") as f:
            header = f.readline().strip().split("\t")
            line = f.readline().strip().split("\t") if not f.closed else []
        if not header or not line:
            continue
        def _get(colname_substr: str) -> Optional[str]:
            for i, h in enumerate(header):
                if colname_substr.lower() in h.lower():
                    return line[i] if i < len(line) else None
            return None
        out["SIRIUS_canopus_superclass"] = _get("superclass")
        out["SIRIUS_canopus_class"] = _get("class")
        out["SIRIUS_canopus_most_specific"] = _get("most") or _get("specific")
        # If at least one field found, stop
        if any(out.values()):
            return out
    return out


def _pick_top_fingerid(fingerid_dir: Path) -> Optional[Dict[str, Optional[str]]]:
    """Select the top-ranked CSI:FingerID hit from a result directory."""
    if not fingerid_dir.exists():
        return None
    best: Optional[Dict[str, Optional[str]]] = None
    best_score: Optional[float] = None
    for tsv in fingerid_dir.glob("*.tsv"):
        # Read minimally to avoid heavy memory use
        with tsv.open("r", encoding="utf-8", newline="") as f:
            rdr = csv.DictReader(f, delimiter="\t")
            row = next(rdr, None)
            if not row:
                continue
            score_str = row.get("score")
            score = float(score_str) if score_str not in (None, "") else None
            if score is None:
                continue
            if best_score is None or score > best_score:
                best_score = score
                best = {
                    "inchikey2D": row.get("inchikey2D"),
                    "name": row.get("name"),
                    "smiles": row.get("smiles"),
                    "xlogp": row.get("xlogp"),
                    "tanimotoSimilarity": row.get("tanimotoSimilarity"),
                    "rank": row.get("rank"),
                    "score": row.get("score"),
                    "molecularFormula": row.get("molecularFormula"),
                }
    return best


def collect_sirius_results(
    pos_dir: Optional[Path],
    neg_dir: Optional[Path],
    output_csv: Path,
    join_with_merged: Optional[Path] = None,
    progress: Optional[Callable[[int, int, str], None]] = None,
    extract_cache: Optional[Path] = None,
    force_rescan: bool = False,
) -> Dict[str, int]:
    """Aggregate SIRIUS identifications and optionally join them onto the merged table."""
    rows: List[Dict[str, Optional[str]]] = []
    processed = 0
    ann_counts: Optional[Dict[str, int]] = None

    def _scan(root: Path, pol: str, total: int) -> Tuple[int, int]:
        """Iterate over SIRIUS compound folders for a given polarity."""
        nonlocal processed
        n_found = n_kept = 0
        if not root or not root.exists():
            return n_found, n_kept
        for comp_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            spectrum = comp_dir / "spectrum.ms"
            compound_info = comp_dir / "compound.info"
            fingerid_dir = comp_dir / "fingerid"
            if not spectrum.exists():
                continue
            n_found += 1
            fid = _parse_feature_id_from_spectrum(spectrum)
            info = _parse_info_from_compound(compound_info) if compound_info.exists() else {}
            cand = _pick_top_fingerid(fingerid_dir)
            if not cand:
                # Capture formula-only / CANOPUS-only evidence to enable L4 assignment
                cano = _extract_canopus(comp_dir)
                has_canopus = any(v for v in cano.values())
                scores_dir = comp_dir / "scores"
                has_formula = scores_dir.exists() and any(scores_dir.glob("*.info"))
                rows.append(
                    {
                        "feature_id": fid,
                        "_polarity": pol,
                        "_sirius_has_struct": False,
                        "_sirius_has_formula": has_formula,
                        "_sirius_has_canopus": has_canopus,
                        "SIRIUS_canopus_superclass": cano.get("SIRIUS_canopus_superclass"),
                        "SIRIUS_canopus_class": cano.get("SIRIUS_canopus_class"),
                        "SIRIUS_canopus_most_specific": cano.get("SIRIUS_canopus_most_specific"),
                    }
                )
                # advance progress
                nonlocal processed
                processed += 1
                if progress:
                    progress(processed, total, str(comp_dir))
                continue
            n_kept += 1
            cano = _extract_canopus(comp_dir)
            has_struct = bool(cand.get("smiles")) and bool(cand.get("inchikey2D"))
            has_formula = bool(cand.get("molecularFormula"))
            has_canopus = any(v for v in cano.values())
            rows.append(
                {
                    "feature_id": fid,
                    # temporary fields used during join/update
                    "sirius_name": cand.get("name"),
                    "_polarity": pol,
                    "_sirius_has_struct": has_struct,
                    "_sirius_has_formula": has_formula,
                    "_sirius_has_canopus": has_canopus,
                    # keep only the requested CANOPUS class fields as new columns
                    "SIRIUS_canopus_superclass": cano.get("SIRIUS_canopus_superclass"),
                    "SIRIUS_canopus_class": cano.get("SIRIUS_canopus_class"),
                    "SIRIUS_canopus_most_specific": cano.get("SIRIUS_canopus_most_specific"),
                }
            )
            processed += 1
            if progress:
                progress(processed, total, str(comp_dir))
        return n_found, n_kept

    pf = pk = nf = nk = 0
    total_dirs = 0
    if pos_dir and pos_dir.exists():
        total_dirs += sum(1 for _ in pos_dir.iterdir() if _.is_dir())
    if neg_dir and neg_dir.exists():
        total_dirs += sum(1 for _ in neg_dir.iterdir() if _.is_dir())

    # Load cached extraction if available and not forcing rescan
    df: pd.DataFrame
    if extract_cache and (not force_rescan) and extract_cache.exists():
        df = pd.read_csv(extract_cache)
    else:
        # Perform scan
        if pos_dir:
            f, k = _scan(pos_dir, "POS", total_dirs)
            pf, pk = f, k
        if neg_dir:
            f, k = _scan(neg_dir, "NEG", total_dirs)
            nf, nk = f, k
        df = pd.DataFrame(rows)
        # Persist extraction cache for future runs
        if extract_cache:
            extract_cache.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(extract_cache, index=False)

    # Normalize columns from existing identifications (compat with historical headers)
    if not df.empty:
        # Ensure we have a sirius_name column for Metabolite name update
        if "sirius_name" not in df.columns:
            name_col = None
            for c in ("name", "SIRIUS_name"):
                if c in df.columns:
                    name_col = c
                    break
            if name_col:
                df["sirius_name"] = df[name_col]
        # Helper flags for level updates
        if "_sirius_has_struct" not in df.columns:
            smiles_col = None
            ik_col = None
            for c in ("smiles", "SIRIUS_smiles"):
                if c in df.columns:
                    smiles_col = c; break
            for c in ("inchikey2D", "SIRIUS_inchikey2D"):
                if c in df.columns:
                    ik_col = c; break
            if smiles_col and ik_col:
                df["_sirius_has_struct"] = df[smiles_col].astype(str).str.len().gt(0) & df[ik_col].astype(str).str.len().gt(0)
        if "_sirius_has_formula" not in df.columns:
            form_col = None
            for c in ("formula", "molecularFormula", "SIRIUS_formula"):
                if c in df.columns:
                    form_col = c; break
            if form_col:
                df["_sirius_has_formula"] = df[form_col].astype(str).str.len().gt(0)
        if "_sirius_has_canopus" not in df.columns:
            can_cols = [c for c in ("SIRIUS_canopus_superclass","SIRIUS_canopus_class","SIRIUS_canopus_most_specific") if c in df.columns]
            if can_cols:
                mask = None
                for c in can_cols:
                    m = df[c].astype(str).str.len().gt(0)
                    mask = m if mask is None else (mask | m)
                if mask is not None:
                    df["_sirius_has_canopus"] = mask
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    if join_with_merged and join_with_merged.exists():
        merged = pd.read_csv(join_with_merged)
        if "feature_id" in merged.columns:
            # Left-join SIRIUS onto merged to keep all original rows
            out = merged.merge(df, on="feature_id", how="left")
            # Overwrite 'Metabolite name' only where a sirius_name exists
            if "sirius_name" in out.columns and "Metabolite name" in out.columns:
                sn = out["sirius_name"].astype("string")
                mask = sn.notna() & (sn.str.len() > 0)
                out.loc[mask, "Metabolite name"] = "SIRIUS_" + sn[mask]
            # Update annotation_level based on SIRIUS evidence for rows with unknown/level 3
            lvl_col = "annotation_level"
            if lvl_col in out.columns:
                # Ensure string dtype before assigning string labels
                out[lvl_col] = out[lvl_col].astype("string")
                lvl = out[lvl_col]
                # Only update where original is missing or '3'
                tgt = lvl.isna() | (lvl == "3")
                has_struct = out.get("_sirius_has_struct")
                has_formula = out.get("_sirius_has_formula")
                has_canopus = out.get("_sirius_has_canopus")
                # Normalize to pandas BooleanDtype to avoid downcasting warnings
                bs = has_struct.astype("boolean") if has_struct is not None else None
                bf = has_formula.astype("boolean") if has_formula is not None else None
                bc = has_canopus.astype("boolean") if has_canopus is not None else None
                struct_mask = (bs.fillna(False) if bs is not None else False)
                canopus_mask = (bc.fillna(False) if bc is not None else False)
                # Assign levels using explicit boolean masks
                out.loc[tgt & struct_mask, lvl_col] = "3"
                # L4 is CANOPUS-only (no structure but CANOPUS class available)
                out.loc[tgt & ~struct_mask & canopus_mask, lvl_col] = "4"
                # Everything else without structure (including formula-only) becomes L5
                out.loc[tgt & ~struct_mask & ~canopus_mask, lvl_col] = "5"
                # Compute post-merge annotation level counts
                vc = out[lvl_col].astype(str).value_counts()
                ann_counts = {k: int(vc.get(k, 0)) for k in ["1","2","3","4","5"]}
            # Drop helper columns
            out = out.drop(columns=[c for c in ("sirius_name","_sirius_has_struct","_sirius_has_formula","_sirius_has_canopus","_polarity") if c in out.columns], errors="ignore")
            out.to_csv(output_csv, index=False)
        else:
            # No feature_id in merged; just write identifications table
            df.to_csv(output_csv, index=False)
    else:
        # No join requested; write identifications table
        df.to_csv(output_csv, index=False)

    # Basic statistics
    miss_pos = pf - pk
    miss_neg = nf - nk
    # If totals are zero (e.g., reused an existing identifications CSV), derive identified counts from any polarity column
    if (pf == 0 and pk == 0 and nf == 0 and nk == 0) and not df.empty:
        pol_col = None
        for cand in ("_polarity", "polarity", "SIRIUS_polarity"):
            if cand in df.columns:
                pol_col = cand
                break
        if pol_col is not None:
            pos_id = int((df[pol_col].astype(str) == "POS").sum())
            neg_id = int((df[pol_col].astype(str) == "NEG").sum())
            pf = pk = pos_id
            nf = nk = neg_id
            miss_pos = 0
            miss_neg = 0
    return {
        "pos_compounds": pf,
        "pos_identified": pk,
        "neg_compounds": nf,
        "neg_identified": nk,
        "miss_pos": miss_pos,
        "miss_neg": miss_neg,
        "rows": int(df.shape[0]),
        "processed": processed,
        "total_dirs": total_dirs,
        "ann_after": ann_counts or {},
    }
