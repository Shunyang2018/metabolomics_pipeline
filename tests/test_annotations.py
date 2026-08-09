"""Tests for MS-DIAL annotation-level assignment.

Level 1 requires a strong weighted score that agrees with the reverse score;
Level 2 a decent reverse score; everything else is Level 3 (unknown).
"""

import pandas as pd
import pytest

from metabo_pipeline.annotations import annotate_levels, assign_annotation_level_row


def row(**kwargs) -> pd.Series:
    """Build an MS-DIAL-shaped row, defaulting to a confident Level 1 match."""
    base = {
        "Metabolite name": "Glucose",
        "Weighted dot product": 0.9,
        "Reverse dot product": 0.85,
        "Matched peaks count": 5,
    }
    base.update(kwargs)
    return pd.Series(base)


class TestLevel3ByName:
    @pytest.mark.parametrize(
        "name",
        [
            "",
            "   ",
            "Unknown",
            "unknown",
            "UNKNOWN",
            "Unknown compound",
            "Low score",
            "low score hit",
            "No MS2",
            "no ms2 spectrum",
        ],
    )
    def test_placeholder_names_are_level_3(self, name):
        assert assign_annotation_level_row(row(**{"Metabolite name": name})) == "3"

    def test_missing_name_is_level_3(self):
        """A blank MS-DIAL cell reads as NaN, which must not pass as a real name."""
        assert assign_annotation_level_row(row(**{"Metabolite name": None})) == "3"
        assert assign_annotation_level_row(row(**{"Metabolite name": float("nan")})) == "3"

    def test_a_real_name_is_not_downgraded(self):
        assert assign_annotation_level_row(row(**{"Metabolite name": "Unknowable acid"})) != "3"


class TestLevel1:
    def test_strong_agreeing_scores_give_level_1(self):
        assert assign_annotation_level_row(row()) == "1"

    def test_weighted_score_must_exceed_075(self):
        assert assign_annotation_level_row(row(**{"Weighted dot product": 0.75})) == "2"

    def test_scores_must_agree_within_02(self):
        """|wdot - rdot| >= 0.2 fails Level 1 but can still make Level 2."""
        assert (
            assign_annotation_level_row(
                row(**{"Weighted dot product": 0.9, "Reverse dot product": 0.6})
            )
            == "2"
        )

    def test_needs_at_least_three_matched_peaks(self):
        assert assign_annotation_level_row(row(**{"Matched peaks count": 2})) == "3"


class TestLevel2:
    def test_moderate_reverse_score_gives_level_2(self):
        assert (
            assign_annotation_level_row(
                row(**{"Weighted dot product": 0.6, "Reverse dot product": 0.7})
            )
            == "2"
        )

    def test_reverse_score_must_exceed_05(self):
        assert (
            assign_annotation_level_row(
                row(**{"Weighted dot product": 0.6, "Reverse dot product": 0.5})
            )
            == "3"
        )

    def test_missing_reverse_score_cannot_reach_level_2(self):
        assert assign_annotation_level_row(row(**{"Reverse dot product": None})) == "3"


class TestLevel3ByScore:
    def test_weak_weighted_score_is_level_3(self):
        assert assign_annotation_level_row(row(**{"Weighted dot product": 0.49})) == "3"

    def test_missing_weighted_score_is_level_3(self):
        assert assign_annotation_level_row(row(**{"Weighted dot product": None})) == "3"

    def test_zero_matched_peaks_is_level_3(self):
        assert assign_annotation_level_row(row(**{"Matched peaks count": 0})) == "3"

    def test_missing_matched_peaks_is_level_3(self):
        assert assign_annotation_level_row(row(**{"Matched peaks count": None})) == "3"

    def test_non_numeric_scores_are_coerced_not_raised(self):
        assert assign_annotation_level_row(row(**{"Weighted dot product": "n/a"})) == "3"

    def test_matched_peaks_given_as_float_string(self):
        assert assign_annotation_level_row(row(**{"Matched peaks count": "5.0"})) == "1"

    def test_completely_empty_row_is_level_3(self):
        assert assign_annotation_level_row(pd.Series(dtype=object)) == "3"


class TestAnnotateLevels:
    def test_adds_a_level_per_row(self):
        df = pd.DataFrame(
            [
                {
                    "Metabolite name": "Glucose",
                    "Weighted dot product": 0.9,
                    "Reverse dot product": 0.85,
                    "Matched peaks count": 5,
                },
                {
                    "Metabolite name": "Unknown",
                    "Weighted dot product": 0.9,
                    "Reverse dot product": 0.85,
                    "Matched peaks count": 5,
                },
                {
                    "Metabolite name": "Citrate",
                    "Weighted dot product": 0.6,
                    "Reverse dot product": 0.7,
                    "Matched peaks count": 4,
                },
            ]
        )
        out = annotate_levels(df)
        assert list(out["annotation_level"]) == ["1", "3", "2"]

    def test_does_not_mutate_the_input(self):
        df = pd.DataFrame([{"Metabolite name": "Glucose"}])
        annotate_levels(df)
        assert "annotation_level" not in df.columns

    def test_empty_frame_gets_the_column_without_raising(self):
        out = annotate_levels(pd.DataFrame(columns=["Metabolite name"]))
        assert "annotation_level" in out.columns
        assert out.empty

    def test_original_columns_are_preserved(self):
        df = pd.DataFrame([{"Metabolite name": "Glucose", "Average Mz": 180.06}])
        out = annotate_levels(df)
        assert "Average Mz" in out.columns
