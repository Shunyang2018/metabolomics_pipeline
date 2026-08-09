"""Tests for replicate-group QC: MS/MS ion counting, grouping, and metrics.

Thresholds under test are the documented defaults: blank fold >= 7, presence
>= 60%, and CV filtering disabled unless a maximum is supplied.
"""

import numpy as np
import pandas as pd
import pytest

from metabo_pipeline.qc import (
    build_group_cols,
    compute_group_metrics,
    count_msms_ions,
    pass_any_mask,
)


class TestCountMsmsIons:
    def test_counts_nonzero_fragments(self):
        assert count_msms_ions("100:50 150:200 200:300") == 3

    def test_zero_intensity_fragments_are_not_counted(self):
        assert count_msms_ions("100:0 150:200 200:0") == 1

    def test_tokens_without_a_colon_are_ignored(self):
        assert count_msms_ions("100:50 junk 150:200") == 2

    def test_non_numeric_intensity_is_ignored(self):
        assert count_msms_ions("100:abc 150:200") == 1

    @pytest.mark.parametrize("value", ["", "   ", "null", "NULL", "na", "NA", "none", "None"])
    def test_null_placeholders_count_zero(self, value):
        assert count_msms_ions(value) == 0

    @pytest.mark.parametrize("value", [None, float("nan"), 123, []])
    def test_non_string_input_counts_zero(self, value):
        assert count_msms_ions(value) == 0

    def test_negative_intensity_is_not_counted(self):
        assert count_msms_ions("100:-5 150:200") == 1


class TestBuildGroupCols:
    def test_each_distinct_sample_is_its_own_group(self):
        groups = build_group_cols(["crude_1_amide", "crude_2_amide"])
        assert len(groups) == 2

    def test_blank_like_columns_never_form_a_group(self):
        groups = build_group_cols(["crude_1", "blank_1", "mb_2", "resuspension_3"])
        assert list(groups) == ["crude_1"]

    def test_whitelist_excludes_unrecognized_columns(self):
        """QC-pool and reference injections must not silently become groups."""
        groups = build_group_cols(
            ["crude_1", "mcx_2", "qcpool_1", "identification_run"],
            real_sample_tokens={"crude", "mcx"},
        )
        assert set(groups) == {"crude_1", "mcx_2"}

    def test_without_a_whitelist_every_non_blank_column_is_a_group(self):
        groups = build_group_cols(["crude_1", "qcpool_1"], real_sample_tokens=None)
        assert set(groups) == {"crude_1", "qcpool_1"}

    def test_group_keys_are_lowercased(self):
        groups = build_group_cols(["Crude_1"])
        assert list(groups) == ["crude_1"]

    def test_columns_with_no_alphanumeric_tokens_are_skipped(self):
        assert build_group_cols(["---", "___"]) == {}

    def test_blank_token_must_be_a_whole_token(self):
        """'blankets' contains 'blank' but is not a blank column."""
        assert "blankets_1" in build_group_cols(["blankets_1"])


class TestComputeGroupMetrics:
    @pytest.fixture
    def frame(self):
        return pd.DataFrame(
            {
                "crude_1": [1000.0, 100.0, 0.0],
                "crude_2": [1200.0, 0.0, 0.0],
                "blank": [100.0, 50.0, 10.0],
            }
        )

    def test_blank_fold_uses_max_sample_over_mean_blank(self, frame):
        groups = build_group_cols(["crude_1", "crude_2"], real_sample_tokens={"crude"})
        out = compute_group_metrics(frame, groups, blank_col="blank")
        # row 0: max(1000, 1200) / 100 = 12
        assert out["blank_fold_crude_1"].iloc[0] == pytest.approx(12.0)

    def test_present_percent_counts_nonzero_replicates(self, frame):
        groups = build_group_cols(["crude_1", "crude_2"], real_sample_tokens={"crude"})
        out = compute_group_metrics(frame, groups, blank_col="blank")
        # each group has a single column here, so presence is 0 or 100
        assert out["present_percent_crude_1"].iloc[0] == pytest.approx(100.0)
        assert out["present_percent_crude_2"].iloc[1] == pytest.approx(0.0)

    def test_zero_blank_yields_nan_rather_than_infinity(self):
        df = pd.DataFrame({"crude_1": [1000.0], "blank": [0.0]})
        groups = build_group_cols(["crude_1"], real_sample_tokens={"crude"})
        out = compute_group_metrics(df, groups, blank_col="blank")
        assert pd.isna(out["blank_fold_crude_1"].iloc[0])

    def test_missing_blank_column_yields_nan_blank_fold(self):
        df = pd.DataFrame({"crude_1": [1000.0]})
        groups = build_group_cols(["crude_1"], real_sample_tokens={"crude"})
        out = compute_group_metrics(df, groups, blank_col=None)
        assert pd.isna(out["blank_fold_crude_1"].iloc[0])

    def test_cv_requires_at_least_two_present_replicates(self):
        df = pd.DataFrame({"crude_1_a": [100.0, 100.0], "crude_1_b": [200.0, 0.0]})
        groups = {"crude_1": ["crude_1_a", "crude_1_b"]}
        out = compute_group_metrics(df, groups, blank_col=None)
        assert not pd.isna(out["cv_percent_crude_1"].iloc[0])
        assert pd.isna(out["cv_percent_crude_1"].iloc[1])

    def test_non_numeric_intensities_are_coerced(self):
        df = pd.DataFrame({"crude_1": ["n/a", "1000"], "blank": [10.0, 10.0]})
        groups = build_group_cols(["crude_1"], real_sample_tokens={"crude"})
        out = compute_group_metrics(df, groups, blank_col="blank")
        assert out["present_percent_crude_1"].iloc[0] == pytest.approx(0.0)

    def test_input_frame_is_not_mutated(self, frame):
        before = list(frame.columns)
        groups = build_group_cols(["crude_1"], real_sample_tokens={"crude"})
        compute_group_metrics(frame, groups, blank_col="blank")
        assert list(frame.columns) == before

    def test_method_blank_columns_join_the_blank_baseline(self):
        """'mb' columns feed the baseline even when blank_col names another column."""
        df = pd.DataFrame({"crude_1": [1000.0], "blank": [100.0], "mb_1": [300.0]})
        groups = build_group_cols(["crude_1"], real_sample_tokens={"crude"})
        out = compute_group_metrics(df, groups, blank_col="blank")
        # mean(100, 300) = 200 -> 1000 / 200 = 5
        assert out["blank_fold_crude_1"].iloc[0] == pytest.approx(5.0)


class TestPassAnyMask:
    def _metrics(self, blank_fold, present_percent, cv=None):
        data = {
            "blank_fold_g": [blank_fold],
            "present_percent_g": [present_percent],
            "cv_percent_g": [cv],
        }
        return pd.DataFrame(data), {"g": ["g_1"]}

    def test_passes_when_both_thresholds_met(self):
        df, groups = self._metrics(10.0, 100.0)
        assert pass_any_mask(df, groups, 7.0, 60.0, None).iloc[0]

    def test_blank_fold_boundary_is_inclusive(self):
        df, groups = self._metrics(7.0, 100.0)
        assert pass_any_mask(df, groups, 7.0, 60.0, None).iloc[0]

    def test_present_percent_boundary_is_inclusive(self):
        df, groups = self._metrics(10.0, 60.0)
        assert pass_any_mask(df, groups, 7.0, 60.0, None).iloc[0]

    def test_fails_below_blank_fold(self):
        df, groups = self._metrics(6.9, 100.0)
        assert not pass_any_mask(df, groups, 7.0, 60.0, None).iloc[0]

    def test_nan_metrics_fail_rather_than_propagate(self):
        df, groups = self._metrics(np.nan, np.nan)
        assert not pass_any_mask(df, groups, 7.0, 60.0, None).iloc[0]

    def test_cv_filter_applies_only_when_a_maximum_is_given(self):
        df, groups = self._metrics(10.0, 100.0, cv=50.0)
        assert pass_any_mask(df, groups, 7.0, 60.0, None).iloc[0]
        assert not pass_any_mask(df, groups, 7.0, 60.0, 30.0).iloc[0]

    def test_a_row_passing_any_group_is_kept(self):
        df = pd.DataFrame(
            {
                "blank_fold_a": [10.0],
                "present_percent_a": [100.0],
                "blank_fold_b": [1.0],
                "present_percent_b": [10.0],
            }
        )
        groups = {"a": ["a_1"], "b": ["b_1"]}
        assert pass_any_mask(df, groups, 7.0, 60.0, None).iloc[0]

    def test_no_recognized_groups_warns_and_keeps_everything(self):
        """Silently filtering nothing would be worse than saying so."""
        df = pd.DataFrame({"x": [1, 2, 3]})
        with pytest.warns(UserWarning, match="No replicate groups"):
            mask = pass_any_mask(df, {}, 7.0, 60.0, None)
        assert mask.all()
        assert len(mask) == 3
