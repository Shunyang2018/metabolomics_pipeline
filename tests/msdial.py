"""Builders for synthetic MS-DIAL alignment tables, shared by the test modules.

An MS-DIAL export has four metadata rows (Class, File type, Injection order,
Batch ID) above the real header row, and sample columns start immediately after
'MS/MS spectrum'.
"""

from __future__ import annotations

# Columns up to and including 'MS/MS spectrum' are metadata; samples follow.
META_COLS = [
    "Alignment ID",
    "Average Rt(min)",
    "Average Mz",
    "Metabolite name",
    "Adduct type",
    "Weighted dot product",
    "Reverse dot product",
    "Matched peaks count",
    "S/N average",
    "MS1 isotopic spectrum",
    "MS/MS spectrum",
]

#: Three fragments, all nonzero, so the default MSMS_MIN_IONS of 3 is met.
GOOD_MSMS = "100.05:200 150.08:800 200.11:400"
GOOD_MS1 = "200.10:1000 201.10:50"

DEFAULT_SAMPLES = ["M2_C18_crude_1_pos", "M2_C18_crude_2_pos", "M2_C18_blank_pos"]


def feature(
    alignment_id=0,
    rt=1.5,
    mz=200.1234,
    name="Glucose",
    adduct="[M+H]+",
    wdot=0.9,
    rdot=0.85,
    peaks=5,
    snr=50.0,
    msms=GOOD_MSMS,
    ms1=GOOD_MS1,
):
    """One MS-DIAL feature row, defaulting to a confident Level 1 match."""
    return [
        str(alignment_id),
        str(rt),
        str(mz),
        name,
        adduct,
        str(wdot),
        str(rdot),
        str(peaks),
        str(snr),
        ms1,
        msms,
    ]


def write_msdial(
    path,
    samples=None,
    features=None,
    sep=",",
    sample_intensity=100000.0,
    blank_intensity=100.0,
):
    """Write a structurally valid MS-DIAL export to ``path``.

    Sample intensities are appended to every feature row. The defaults give a
    blank fold-change of 1000, well above the documented threshold of 7, so
    features survive QC unless a test deliberately makes them fail. Columns
    whose name contains 'blank' are classed as blanks.
    """
    samples = DEFAULT_SAMPLES if samples is None else samples
    features = [feature()] if features is None else features

    lead = [""] * len(META_COLS)
    is_blank = ["blank" in s.lower() for s in samples]
    intensities = [str(blank_intensity if b else sample_intensity) for b in is_blank]
    rows = [
        lead + ["Blank" if b else "Sample" for b in is_blank],
        lead + ["Sample"] * len(samples),
        lead + [str(i + 1) for i in range(len(samples))],
        lead + ["1"] * len(samples),
        META_COLS + samples,
    ]
    rows.extend(list(f) + intensities for f in features)
    path.write_text("\n".join(sep.join(r) for r in rows), encoding="utf-8")
    return path
