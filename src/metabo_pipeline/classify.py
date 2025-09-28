from __future__ import annotations

import time
import json
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import requests
import yaml
import pandas as pd


def _fetch_classyfire(
    inchikey: str,
    base_url: str = "http://classyfire.wishartlab.com",
    timeout: float = 15.0,
    retries: int = 1,
    backoff: float = 1.5,
) -> Optional[Dict]:
    """Fetch classification for an InChIKey from ClassyFire.

    Tries common endpoints; returns parsed JSON dict or None.
    """
    inchikey = (inchikey or "").strip()
    if not inchikey:
        return None
    # Try entity endpoint
    urls = [
        f"{base_url}/entities/{inchikey}.json",
        f"{base_url}/entities/{inchikey}",
        f"{base_url}/classification?inchikey={inchikey}",
    ]
    attempt = 0
    while attempt <= max(0, retries):
        for url in urls:
            try:
                r = requests.get(url, timeout=timeout)
                if r.status_code == 200:
                    return r.json()
                if r.status_code in (429, 500, 502, 503, 504):
                    # transient; try next url or next attempt
                    continue
            except Exception:
                # network error; try next url or next attempt
                continue
        attempt += 1
        if attempt <= retries:
            time.sleep(backoff ** (attempt - 1))
    return None


def _extract_taxonomy(cf: Dict) -> Dict[str, Optional[str]]:
    """Extract key taxonomy fields from a ClassyFire response."""
    out = {
        "cf_kingdom": None,
        "cf_superclass": None,
        "cf_class": None,
        "cf_subclass": None,
        "cf_direct_parent": None,
    }
    if not cf:
        return out
    # ClassyFire entities may have 'kingdom', 'superclass', 'class', 'subclass', 'direct_parent' objects with 'name'
    for key, out_key in [
        ("kingdom", "cf_kingdom"),
        ("superclass", "cf_superclass"),
        ("class", "cf_class"),
        ("subclass", "cf_subclass"),
        ("direct_parent", "cf_direct_parent"),
    ]:
        node = cf.get(key)
        if isinstance(node, dict):
            out[out_key] = node.get("name")
        elif isinstance(node, str):
            out[out_key] = node
    # Some endpoints wrap taxonomy under 'classification' or 'taxonomy'
    tax = cf.get("classification") or cf.get("taxonomy")
    if isinstance(tax, dict):
        for key, out_key in [
            ("kingdom", "cf_kingdom"),
            ("superclass", "cf_superclass"),
            ("class", "cf_class"),
            ("subclass", "cf_subclass"),
            ("direct_parent", "cf_direct_parent"),
        ]:
            node = tax.get(key)
            if isinstance(node, dict) and not out[out_key]:
                out[out_key] = node.get("name")
    return out


def _load_cache(path: Optional[str]) -> Dict[str, Dict[str, Optional[str]]]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    try:
        if p.suffix.lower() in (".yml", ".yaml"):
            data = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
        else:
            data = json.loads(p.read_text(encoding="utf-8"))
        # normalize keys to str and values to dict of expected fields
        out: Dict[str, Dict[str, Optional[str]]] = {}
        for k, v in (data or {}).items():
            if isinstance(v, dict):
                out[str(k)] = {
                    "cf_kingdom": v.get("cf_kingdom"),
                    "cf_superclass": v.get("cf_superclass"),
                    "cf_class": v.get("cf_class"),
                    "cf_subclass": v.get("cf_subclass"),
                    "cf_direct_parent": v.get("cf_direct_parent"),
                }
        return out
    except Exception:
        return {}


def _save_cache(path: Optional[str], cache: Dict[str, Dict[str, Optional[str]]]) -> None:
    if not path:
        return
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        if p.suffix.lower() in (".yml", ".yaml"):
            p.write_text(yaml.safe_dump(cache, sort_keys=True, allow_unicode=True), encoding="utf-8")
        else:
            p.write_text(json.dumps(cache, ensure_ascii=False, indent=2, sort_keys=True), encoding="utf-8")
    except Exception:
        # Swallow cache write errors (non-fatal)
        pass


def probe_classyfire(base_url: str, timeout: float = 10.0) -> Dict[str, object]:
    """Probe a ClassyFire-compatible API for availability and latency.

    Tries both entity and classification endpoints with a known InChIKey.
    Returns a dict with statuses and timings.
    """
    key = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"  # acetaminophen (paracetamol)
    urls = [
        f"{base_url.rstrip('/')}/entities/{key}.json",
        f"{base_url.rstrip('/')}/classification?inchikey={key}",
    ]
    out: Dict[str, object] = {"base_url": base_url, "ok": True, "checks": []}
    for u in urls:
        rec: Dict[str, object] = {"url": u, "ok": False, "ms": None, "status": None, "error": None}
        try:
            import time as _t
            t0 = _t.perf_counter()
            r = requests.get(u, timeout=timeout)
            dt = int(( _t.perf_counter() - t0) * 1000)
            rec["ms"] = dt
            rec["status"] = r.status_code
            if r.status_code == 200:
                rec["ok"] = True
            else:
                out["ok"] = False
        except Exception as e:
            rec["error"] = str(e)
            out["ok"] = False
        out["checks"].append(rec)
    return out


def classify_level12_with_classyfire(
    input_csv: str,
    output_csv: str,
    sleep_sec: float = 1.0,
    base_url: str = "http://classyfire.wishartlab.com",
    results_only: bool = False,
    id_column: str = "feature_id",
    cache_path: Optional[str] = "outputs/classyfire_cache.json",
    progress: Optional[Callable[[int, int, str, bool, str], None]] = None,
    offline: bool = False,
    timeout: float = 15.0,
    retries: int = 2,
    backoff: float = 1.5,
) -> Dict[str, int]:
    """Classify Level 1/2 rows via ClassyFire using unique InChIKeys and merge back.

    Returns a summary dict with counts.
    """
    # Derive a friendlier default output if caller used the placeholder path
    try:
        if str(output_csv).replace("\\", "/").endswith("outputs/merged_classified.csv") and input_csv:
            ip = Path(input_csv)
            suffix = ip.suffix or ".csv"
            output_csv = str(ip.with_name(ip.stem + "_classyfire").with_suffix(suffix))
    except Exception:
        pass
    df = pd.read_csv(input_csv)
    if "annotation_level" not in df.columns or "INCHIKEY" not in df.columns:
        cols = ", ".join(df.columns)
        raise ValueError("Input CSV must contain 'annotation_level' and 'INCHIKEY' columns. "
                         f"Found columns: {cols}")

    mask = df["annotation_level"].astype(str).isin(["1", "2"]) & df["INCHIKEY"].astype(str).str.len().gt(0)
    keys = sorted(set(df.loc[mask, "INCHIKEY"].astype(str)))
    # Load cache and prepare results with cached hits
    cache = _load_cache(cache_path)
    results: Dict[str, Dict[str, Optional[str]]] = {}
    hit = miss = 0
    total = len(keys)
    processed = 0
    for ik in keys:
        if ik in cache:
            results[ik] = cache[ik]
            hit += 1
            processed += 1
            if progress:
                try:
                    progress(processed, total, ik, True, "cache")
                except Exception:
                    pass
    to_query = [ik for ik in keys if ik not in results]

    # Only attempt network calls for missing keys
    if offline:
        # In offline mode, do not perform any network calls. Mark remaining keys as skipped.
        skipped = len(to_query)
    else:
        skipped = 0
        for i, ik in enumerate(to_query, 1):
            data = _fetch_classyfire(ik, base_url=base_url, timeout=timeout, retries=retries, backoff=backoff)
            if data is not None:
                tax = _extract_taxonomy(data)
                results[ik] = tax
                cache[ik] = tax
                hit += 1
                status = "hit"
            else:
                results[ik] = {
                    "cf_kingdom": None,
                    "cf_superclass": None,
                    "cf_class": None,
                    "cf_subclass": None,
                    "cf_direct_parent": None,
                }
                # Do not add negative results to cache; try again another time
                miss += 1
                status = "miss"
            processed += 1
            if progress:
                try:
                    progress(processed, total, ik, False, status)
                except Exception:
                    pass
            time.sleep(max(0.0, sleep_sec))

    # Persist updated cache
    _save_cache(cache_path, cache)

    # Merge back
    tax_df = pd.DataFrame.from_dict(results, orient="index").reset_index().rename(columns={"index": "INCHIKEY"})
    out = df.merge(tax_df, on="INCHIKEY", how="left")

    cf_cols = [c for c in ["cf_kingdom", "cf_superclass", "cf_class", "cf_subclass", "cf_direct_parent"] if c in out.columns]

    if results_only:
        # Write only identifier + classification columns
        if id_column not in out.columns:
            # Fallback: try Alignment ID; otherwise error
            if "Alignment ID" in out.columns:
                out[id_column] = out["Alignment ID"]
            else:
                raise ValueError(f"Identifier column '{id_column}' not found in input CSV")
        minimal_cols = [id_column] + cf_cols
        out[minimal_cols].drop_duplicates(subset=[id_column]).to_csv(output_csv, index=False)
        return {"unique_keys": len(keys), "hits": hit, "miss": miss, "rows": int(out.shape[0]), "added_cols": len(cf_cols), "skipped_offline": skipped}
    else:
        # Reorder: insert cf_* columns immediately after 'SMILES' if present,
        # otherwise after 'Metabolite name' if present.
        if cf_cols:
            cols = list(out.columns)
            for c in cf_cols:
                if c in cols:
                    cols.remove(c)
            insert_after = None
            if "SMILES" in cols:
                insert_after = "SMILES"
            elif "Metabolite name" in cols:
                insert_after = "Metabolite name"
            if insert_after is not None:
                insert_pos = cols.index(insert_after) + 1
                new_cols = cols[:insert_pos] + cf_cols + cols[insert_pos:]
                out = out.reindex(columns=new_cols)

        out.to_csv(output_csv, index=False)
        return {"unique_keys": len(keys), "hits": hit, "miss": miss, "rows": int(out.shape[0]), "added_cols": len(cf_cols), "skipped_offline": skipped}
