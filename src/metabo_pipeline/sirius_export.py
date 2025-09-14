from __future__ import annotations

from pathlib import Path
from typing import List, Tuple

import pandas as pd

from .utils import parse_spectrum


def build_ms_entries(l3_df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    pos_entries: List[str] = []
    neg_entries: List[str] = []
    for _, r in l3_df.iterrows():
        name = f"{r.get('Metabolite name','Unknown')}|{r.get('Adduct type','')}|{r.get('Average Rt(min)','')}"
        precursor = r.get("Average Mz", "")
        ion = r.get("Adduct type", "")
        ms1 = parse_spectrum(r.get("MS1 isotopic spectrum", ""))
        ms2 = parse_spectrum(r.get("MS/MS spectrum", ""))
        block = [
            f"Name: {name}",
            f"PrecursorMz: {precursor}",
            f"Ionization: {ion}",
        ]
        if ms1:
            block.append("MS1:")
            block.extend([f"{mz} {inten}" for mz, inten in ms1])
        if ms2:
            block.append("MS2:")
            block.extend([f"{mz} {inten}" for mz, inten in ms2])
        block.append("")
        entry = "\n".join(block)
        pol = "POS" if "+" in str(ion or "") else ("NEG" if "-" in str(ion or "") else "UNK")
        if pol == "POS":
            pos_entries.append(entry)
        elif pol == "NEG":
            neg_entries.append(entry)
    return pos_entries, neg_entries


def write_ms_files(pos_entries: List[str], neg_entries: List[str], out_dir: Path) -> Tuple[int, int]:
    out_dir.mkdir(parents=True, exist_ok=True)
    pos_path = out_dir / "sirius_unknown_pos.ms"
    neg_path = out_dir / "sirius_unknown_neg.ms"
    if pos_entries:
        pos_path.write_text("\n".join(pos_entries), encoding="utf-8")
    if neg_entries:
        neg_path.write_text("\n".join(neg_entries), encoding="utf-8")
    return len(pos_entries), len(neg_entries)

