"""Tests for building SIRIUS .ms entries from Level 3 unknowns.

Entries are split by polarity, inferred from the adduct sign, because SIRIUS is
run separately per ionization mode.
"""

import pandas as pd

from metabo_pipeline.sirius_export import build_ms_entries


def frame(rows):
    defaults = {
        "feature_id": "F1",
        "Average Mz": 200.1234,
        "Average Rt(min)": 1.5,
        "Adduct type": "[M+H]+",
        "Metabolite name": "Unknown",
        "MS1 isotopic spectrum": "200.1:1000 201.1:50",
        "MS/MS spectrum": "100.05:200 150.08:800",
    }
    return pd.DataFrame([{**defaults, **r} for r in rows])


def test_positive_adduct_goes_to_the_positive_file():
    pos, neg = build_ms_entries(frame([{"Adduct type": "[M+H]+"}]))
    assert len(pos) == 1 and neg == []


def test_negative_adduct_goes_to_the_negative_file():
    pos, neg = build_ms_entries(frame([{"Adduct type": "[M-H]-"}]))
    assert pos == [] and len(neg) == 1


def test_unknown_polarity_is_dropped_from_both():
    pos, neg = build_ms_entries(frame([{"Adduct type": ""}]))
    assert pos == [] and neg == []


def test_entry_contains_the_required_sirius_directives():
    pos, _ = build_ms_entries(frame([{}]))
    entry = pos[0]
    for directive in (">compound", ">parentmass", ">retentiontime", ">ionization", ">ms1"):
        assert directive in entry
    assert ">collision" in entry


def test_precursor_mass_and_rt_are_carried_through():
    pos, _ = build_ms_entries(frame([{"Average Mz": 300.5, "Average Rt(min)": 4.25}]))
    assert ">parentmass\t300.5" in pos[0]
    assert ">retentiontime\t4.25" in pos[0]


def test_feature_id_is_appended_to_the_compound_name():
    """SIRIUS summaries key on the compound name, so it must map back."""
    pos, _ = build_ms_entries(frame([{"feature_id": "F42", "Metabolite name": "Unknown"}]))
    assert ">compound\tUnknown_F42" in pos[0]


def test_feature_id_is_not_appended_twice():
    pos, _ = build_ms_entries(frame([{"feature_id": "F42", "Metabolite name": "Unknown_F42"}]))
    assert pos[0].count("F42") == 2  # once in >compound, once in the #feature_id comment


def test_feature_id_is_recorded_as_a_comment():
    pos, _ = build_ms_entries(frame([{"feature_id": "F7"}]))
    assert "#feature_id F7" in pos[0]


def test_blank_metabolite_name_becomes_unknown():
    """NaN is truthy-safe here: it must not produce the name 'nan'."""
    pos, _ = build_ms_entries(frame([{"Metabolite name": float("nan"), "feature_id": "F1"}]))
    assert ">compound\tUnknown_F1" in pos[0]
    assert "nan" not in pos[0]


def test_empty_metabolite_name_becomes_unknown():
    pos, _ = build_ms_entries(frame([{"Metabolite name": "", "feature_id": "F1"}]))
    assert ">compound\tUnknown_F1" in pos[0]


def test_spectra_peaks_are_written_as_mz_intensity_lines():
    pos, _ = build_ms_entries(
        frame([{"MS1 isotopic spectrum": "200.1:1000", "MS/MS spectrum": "100.05:200"}])
    )
    assert "200.1 1000.0" in pos[0]
    assert "100.05 200.0" in pos[0]


def test_missing_spectra_still_produce_an_entry():
    pos, _ = build_ms_entries(frame([{"MS1 isotopic spectrum": "", "MS/MS spectrum": ""}]))
    assert len(pos) == 1
    assert ">ms1" in pos[0]


def test_empty_frame_produces_no_entries():
    pos, neg = build_ms_entries(frame([]))
    assert pos == [] and neg == []


def test_mixed_polarities_are_split():
    pos, neg = build_ms_entries(
        frame(
            [
                {"Adduct type": "[M+H]+", "feature_id": "F1"},
                {"Adduct type": "[M-H]-", "feature_id": "F2"},
                {"Adduct type": "[M+Na]+", "feature_id": "F3"},
            ]
        )
    )
    assert len(pos) == 2
    assert len(neg) == 1
    assert "F2" in neg[0]
