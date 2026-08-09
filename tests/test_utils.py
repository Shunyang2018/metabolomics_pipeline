"""Tests for filename inference, sample-name normalization, and spectrum parsing."""

import pytest

from metabo_pipeline.utils import (
    infer_chrom_from_name,
    infer_mode_from_name,
    normalize_sample_id_core,
    parse_spectrum,
)


class TestInferChrom:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("Area_2_HILIC_pos.txt", "HILIC"),
            ("alignment_hilic_neg.csv", "HILIC"),
            ("M2_C18_TV_pos.csv", "C18"),
            ("something_c18.txt", "C18"),
            ("M2_Lipidomics_neg.csv", "Lipidomics"),
            ("lipids_pos.csv", "Lipidomics"),
            ("mystery_file.csv", "unknown"),
        ],
    )
    def test_assay_from_filename(self, name, expected):
        assert infer_chrom_from_name(name) == expected

    def test_hilic_wins_over_c18_when_both_present(self):
        """HILIC is checked first; the function is documented as first-match."""
        assert infer_chrom_from_name("hilic_c18_pos.csv") == "HILIC"

    @pytest.mark.parametrize("value", ["", None])
    def test_blank_input_is_unknown(self, value):
        assert infer_chrom_from_name(value) == "unknown"


class TestInferMode:
    @pytest.mark.parametrize(
        "name,expected",
        [
            ("M2_C18_pos.csv", "POS"),
            ("M2_C18_POS.csv", "POS"),
            ("M2_HILIC_neg.txt", "NEG"),
            ("alignment.csv", "UNK"),
        ],
    )
    def test_polarity_from_filename(self, name, expected):
        assert infer_mode_from_name(name) == expected

    @pytest.mark.parametrize("value", ["", None])
    def test_blank_input_is_unknown(self, value):
        assert infer_mode_from_name(value) == "UNK"


class TestNormalizeSampleId:
    def test_readme_examples(self):
        """The three examples documented in the README must hold."""
        assert normalize_sample_id_core("M2_Lipids_TV_5+6_pos") == "m2_tv_5+6"
        assert normalize_sample_id_core("M2_C18_TV_5+6_POS") == "m2_tv_5+6"
        assert normalize_sample_id_core("M2_Hilic_PE_3+4_neg") == "m2_pe_3+4"

    def test_assay_and_polarity_tokens_are_stripped(self):
        assert normalize_sample_id_core("sample_hilic_neg") == "sample"
        assert normalize_sample_id_core("sample_lipidomics_pos") == "sample"
        assert normalize_sample_id_core("sample_ms1") == "sample"

    def test_columns_differing_only_by_assay_collapse_together(self):
        """This is the property that makes cross-assay merging possible."""
        assert normalize_sample_id_core("M2_C18_crude_1_pos") == normalize_sample_id_core(
            "M2_Hilic_crude_1_neg"
        )

    def test_spaces_and_hyphens_become_underscores(self):
        assert normalize_sample_id_core("M2 crude-1") == "m2_crude_1"

    def test_leading_and_trailing_separators_are_trimmed(self):
        assert normalize_sample_id_core("_pos_crude_") == "crude"

    def test_repeated_separators_collapse(self):
        assert normalize_sample_id_core("a__b___c") == "a_b_c"

    def test_substring_matches_are_not_stripped(self):
        """'position' contains 'pos' but is not a polarity token."""
        assert "position" in normalize_sample_id_core("sample_position")

    @pytest.mark.parametrize("value", ["", None])
    def test_blank_input(self, value):
        assert normalize_sample_id_core(value) == ""


class TestParseSpectrum:
    def test_parses_mz_intensity_pairs(self):
        assert parse_spectrum("100.5:200 150.25:300.5") == [
            (100.5, 200.0),
            (150.25, 300.5),
        ]

    def test_tokens_without_a_colon_are_skipped(self):
        assert parse_spectrum("100:200 garbage 150:300") == [(100.0, 200.0), (150.0, 300.0)]

    def test_non_numeric_tokens_are_skipped(self):
        assert parse_spectrum("abc:def 150:300") == [(150.0, 300.0)]

    def test_intensity_containing_extra_colons_is_kept_whole(self):
        """split(':', 1) means only the first colon separates m/z from intensity."""
        assert parse_spectrum("100:200:300") == []

    @pytest.mark.parametrize("value", ["", None, "   "])
    def test_blank_input_is_empty(self, value):
        assert parse_spectrum(value) == []
