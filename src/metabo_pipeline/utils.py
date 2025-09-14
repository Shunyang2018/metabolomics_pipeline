from __future__ import annotations

import re
from typing import Iterable, List, Tuple


def infer_chrom_from_name(name: str) -> str:
    n = (name or "").lower()
    if "hilic" in n:
        return "HILIC"
    if "c18" in n:
        return "C18"
    if "lipid" in n:
        return "Lipidomics"
    return "unknown"


def infer_mode_from_name(name: str) -> str:
    n = (name or "").lower()
    if "pos" in n:
        return "POS"
    if "neg" in n:
        return "NEG"
    return "UNK"


_TOKEN_PAT = re.compile(r"(?i)(^|[\W_])(lipidomics|lipids|lipid|hilic|c18|pos|neg|ms1)(?=($|[\W_]))")


def normalize_sample_id_core(name: str) -> str:
    v = str(name or "").strip().lower()
    v = re.sub(r"[\s\-]+", "_", v)
    v = _TOKEN_PAT.sub(lambda m: "_" if m.group(1) else "", v)
    v = re.sub(r"_+", "_", v).strip("_")
    return v


def parse_spectrum(cell: str) -> List[Tuple[float, float]]:
    out: List[Tuple[float, float]] = []
    for tok in str(cell or "").split():
        if ":" not in tok:
            continue
        a, b = tok.split(":", 1)
        try:
            out.append((float(a), float(b)))
        except Exception:
            continue
    return out

