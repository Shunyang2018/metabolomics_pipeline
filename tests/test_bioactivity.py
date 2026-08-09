"""Tests for bioactivity matching on the InChIKey skeleton block.

Matching uses only the first InChIKey block (connectivity), so a feature whose
stereochemistry or protonation differs from the reference still matches.
"""

import pandas as pd
import pytest

from metabo_pipeline.bioactivity import (
    identify_sample_columns,
    inchikey_block1,
    match_bioactives,
)

GLUCOSE = "WQZGKKKJIJFFOK-GASJEMHNSA-N"
CITRATE = "KRKNYBCHXYNGOX-UHFFFAOYSA-N"


class TestInchikeyBlock1:
    def test_returns_the_first_block(self):
        assert inchikey_block1(GLUCOSE) == "WQZGKKKJIJFFOK"

    def test_uppercases_and_strips(self):
        assert inchikey_block1("  wqzgkkkjijffok-gasjemhnsa-n  ") == "WQZGKKKJIJFFOK"

    def test_a_bare_block_passes_through(self):
        assert inchikey_block1("WQZGKKKJIJFFOK") == "WQZGKKKJIJFFOK"

    @pytest.mark.parametrize("value", [None, "", "   ", float("nan"), "nan", "NaN"])
    def test_blank_and_nan_values_return_none(self, value):
        assert inchikey_block1(value) is None

    def test_leading_hyphen_returns_none(self):
        assert inchikey_block1("-ABC-N") is None

    def test_stereochemistry_differences_still_match(self):
        a = "WQZGKKKJIJFFOK-GASJEMHNSA-N"
        b = "WQZGKKKJIJFFOK-VFUOTHLCSA-N"
        assert inchikey_block1(a) == inchikey_block1(b)


class TestIdentifySampleColumns:
    def test_msdial_metadata_columns_are_excluded(self):
        cols = ["feature_id", "Average Mz", "Metabolite name", "INCHIKEY", "crude_1"]
        assert identify_sample_columns(cols) == ["crude_1"]

    def test_qc_metric_columns_are_excluded(self):
        cols = [
            "blank_fold_crude_1",
            "present_percent_crude_1",
            "cv_percent_crude_1",
            "pass_any_groups",
            "crude_1",
        ]
        assert identify_sample_columns(cols) == ["crude_1"]

    def test_sirius_and_classyfire_columns_are_excluded(self):
        cols = ["SIRIUS_formula", "cf_superclass", "cf_class", "crude_1"]
        assert identify_sample_columns(cols) == ["crude_1"]

    def test_internal_helper_columns_are_excluded(self):
        assert identify_sample_columns(["_polarity", "_rt", "_mz", "crude_1"]) == ["crude_1"]

    def test_order_is_preserved(self):
        assert identify_sample_columns(["b_1", "a_1"]) == ["b_1", "a_1"]


class TestMatchBioactives:
    @pytest.fixture
    def files(self, tmp_path):
        merged = tmp_path / "merged.csv"
        db = tmp_path / "db.csv"
        out = tmp_path / "bioactives.csv"
        return merged, db, out

    def _write(self, merged, db, feats, entries):
        pd.DataFrame(feats).to_csv(merged, index=False)
        pd.DataFrame(entries).to_csv(db, index=False)

    def test_matches_on_skeleton_and_emits_one_row_per_sample_hit(self, files):
        merged, db, out = files
        self._write(
            merged,
            db,
            [
                {
                    "feature_id": "F1",
                    "Metabolite name": "Glucose",
                    "INCHIKEY": GLUCOSE,
                    "crude_1": 1000,
                    "crude_2": 2000,
                }
            ],
            [{"InChIKey": GLUCOSE, "activity": "antioxidant"}],
        )
        summary = match_bioactives(merged, db, out)
        result = pd.read_csv(out)
        # detected in two samples -> two rows
        assert len(result) == 2
        assert set(result["sample_type"]) == {"crude_1", "crude_2"}
        assert isinstance(summary, dict)

    def test_zero_intensity_samples_are_not_reported(self, files):
        merged, db, out = files
        self._write(
            merged,
            db,
            [{"feature_id": "F1", "INCHIKEY": GLUCOSE, "crude_1": 1000, "crude_2": 0}],
            [{"InChIKey": GLUCOSE, "activity": "antioxidant"}],
        )
        match_bioactives(merged, db, out)
        result = pd.read_csv(out)
        assert list(result["sample_type"]) == ["crude_1"]

    def test_stereoisomer_in_the_database_still_matches(self, files):
        merged, db, out = files
        self._write(
            merged,
            db,
            [{"feature_id": "F1", "INCHIKEY": "WQZGKKKJIJFFOK-GASJEMHNSA-N", "crude_1": 5}],
            [{"InChIKey": "WQZGKKKJIJFFOK-VFUOTHLCSA-N", "activity": "antioxidant"}],
        )
        match_bioactives(merged, db, out)
        assert len(pd.read_csv(out)) == 1

    def test_unmatched_features_are_absent(self, files):
        merged, db, out = files
        self._write(
            merged,
            db,
            [{"feature_id": "F1", "INCHIKEY": CITRATE, "crude_1": 5}],
            [{"InChIKey": GLUCOSE, "activity": "antioxidant"}],
        )
        match_bioactives(merged, db, out)
        assert pd.read_csv(out).empty

    def test_features_without_an_inchikey_are_skipped(self, files):
        merged, db, out = files
        self._write(
            merged,
            db,
            [{"feature_id": "F1", "INCHIKEY": "", "crude_1": 5}],
            [{"InChIKey": GLUCOSE, "activity": "antioxidant"}],
        )
        match_bioactives(merged, db, out)
        assert pd.read_csv(out).empty

    def test_missing_inchikey_column_raises(self, files):
        merged, db, out = files
        self._write(merged, db, [{"feature_id": "F1", "crude_1": 5}], [{"InChIKey": GLUCOSE}])
        with pytest.raises(ValueError, match="INCHIKEY"):
            match_bioactives(merged, db, out)

    def test_missing_database_column_raises(self, files):
        merged, db, out = files
        self._write(merged, db, [{"INCHIKEY": GLUCOSE, "crude_1": 5}], [{"wrong": GLUCOSE}])
        with pytest.raises(ValueError, match="InChIKey"):
            match_bioactives(merged, db, out)

    def test_custom_column_names_are_honoured(self, files):
        merged, db, out = files
        pd.DataFrame([{"my_key": GLUCOSE, "crude_1": 5}]).to_csv(merged, index=False)
        pd.DataFrame([{"db_key": GLUCOSE, "activity": "x"}]).to_csv(db, index=False)
        match_bioactives(merged, db, out, inchikey_col="my_key", db_inchikey_col="db_key")
        assert len(pd.read_csv(out)) == 1
