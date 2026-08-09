"""Tests for Level 3 representative selection.

Level 3 (unknown) features are clustered within (chrom, polarity) by RT and m/z
proximity, and one representative is kept per cluster — the one with the higher
S/N, then the higher weighted dot product. This is what limits how many
unknowns get exported to SIRIUS.
"""

import pandas as pd
import pytest

from metabo_pipeline.dedup import l3_representatives

RT_WINDOW = 0.6
MZ_PPM = 20.0


def frame(rows):
    """Build an L3 frame with the columns l3_representatives reads."""
    defaults = {
        "chrom": "C18",
        "mode": "POS",
        "Average Rt(min)": 1.0,
        "Average Mz": 200.0,
        "S/N average": 10.0,
        "Weighted dot product": 0.5,
    }
    return pd.DataFrame([{**defaults, **r} for r in rows])


def test_empty_frame_returns_empty():
    out = l3_representatives(pd.DataFrame(), RT_WINDOW, MZ_PPM)
    assert out.empty


def test_missing_mode_column_raises():
    df = pd.DataFrame([{"chrom": "C18", "Average Rt(min)": 1.0, "Average Mz": 200.0}])
    with pytest.raises(ValueError, match="mode"):
        l3_representatives(df, RT_WINDOW, MZ_PPM)


def test_single_row_is_its_own_representative():
    out = l3_representatives(frame([{}]), RT_WINDOW, MZ_PPM)
    assert len(out) == 1


def test_near_duplicates_collapse_to_one():
    """Same m/z, RT 0.1 min apart -> one cluster."""
    out = l3_representatives(
        frame([{"Average Rt(min)": 1.0}, {"Average Rt(min)": 1.1}]), RT_WINDOW, MZ_PPM
    )
    assert len(out) == 1


def test_rows_beyond_the_rt_window_are_kept_separately():
    out = l3_representatives(
        frame([{"Average Rt(min)": 1.0}, {"Average Rt(min)": 5.0}]), RT_WINDOW, MZ_PPM
    )
    assert len(out) == 2


def test_rows_beyond_the_mz_tolerance_are_kept_separately():
    """20 ppm of 200 Da is 0.004 Da, so a 1 Da gap is a different feature."""
    out = l3_representatives(
        frame([{"Average Mz": 200.0}, {"Average Mz": 201.0}]), RT_WINDOW, MZ_PPM
    )
    assert len(out) == 2


def test_higher_signal_to_noise_wins_within_a_cluster():
    out = l3_representatives(
        frame(
            [
                {"Average Rt(min)": 1.0, "S/N average": 5.0},
                {"Average Rt(min)": 1.05, "S/N average": 50.0},
            ]
        ),
        RT_WINDOW,
        MZ_PPM,
    )
    assert len(out) == 1
    assert out["S/N average"].iloc[0] == 50.0


def test_weighted_dot_breaks_signal_to_noise_ties():
    out = l3_representatives(
        frame(
            [
                {"Average Rt(min)": 1.0, "S/N average": 10.0, "Weighted dot product": 0.2},
                {"Average Rt(min)": 1.05, "S/N average": 10.0, "Weighted dot product": 0.9},
            ]
        ),
        RT_WINDOW,
        MZ_PPM,
    )
    assert len(out) == 1
    assert out["Weighted dot product"].iloc[0] == 0.9


def test_different_polarities_never_cluster_together():
    out = l3_representatives(frame([{"mode": "POS"}, {"mode": "NEG"}]), RT_WINDOW, MZ_PPM)
    assert len(out) == 2


def test_different_chromatography_never_clusters_together():
    out = l3_representatives(frame([{"chrom": "C18"}, {"chrom": "HILIC"}]), RT_WINDOW, MZ_PPM)
    assert len(out) == 2


def test_unrecognized_mode_is_bucketed_as_unknown():
    out = l3_representatives(frame([{"mode": "weird"}, {"mode": "UNK"}]), RT_WINDOW, MZ_PPM)
    assert len(out) == 1


def test_helper_columns_are_removed_from_the_result():
    out = l3_representatives(frame([{}]), RT_WINDOW, MZ_PPM)
    assert not [c for c in out.columns if c.startswith("_")]


def test_missing_signal_to_noise_does_not_crash():
    out = l3_representatives(
        frame([{"Average Rt(min)": 1.0, "S/N average": None}, {"Average Rt(min)": 1.05}]),
        RT_WINDOW,
        MZ_PPM,
    )
    assert len(out) == 1


def test_missing_rt_or_mz_does_not_cluster():
    """Rows without coordinates cannot be shown to be duplicates."""
    out = l3_representatives(
        frame([{"Average Rt(min)": None}, {"Average Mz": None}]), RT_WINDOW, MZ_PPM
    )
    assert len(out) == 2


def test_input_frame_is_not_mutated():
    df = frame([{}, {}])
    before = list(df.columns)
    l3_representatives(df, RT_WINDOW, MZ_PPM)
    assert list(df.columns) == before


def test_chain_of_near_neighbours_forms_one_cluster():
    """Single-linkage: each row is compared to the previous cluster member."""
    rows = [{"Average Rt(min)": 1.0 + 0.1 * i} for i in range(6)]
    out = l3_representatives(frame(rows), RT_WINDOW, MZ_PPM)
    assert len(out) == 1
