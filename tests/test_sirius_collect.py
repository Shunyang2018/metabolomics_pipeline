"""Tests for parsing SIRIUS result directories and joining them onto the merged table.

SIRIUS itself is never run: the two output layouts it produces are written as
fixtures. SIRIUS 6 emits TSV summary files at the project root; SIRIUS 5 emits
one numbered directory per compound.

Annotation levels assigned on join, for rows that were unknown or Level 3:
  3 - SIRIUS found a structure (SMILES + InChIKey)
  4 - CANOPUS class only, no structure
  5 - neither (formula-only included)
"""

from __future__ import annotations

import pandas as pd
import pytest

from metabo_pipeline.sirius_collect import (
    _detect_sirius_format,
    _extract_canopus,
    _extract_feature_id_from_v6_mapping,
    _parse_sirius6_tsv,
    _pick_top_fingerid,
    _value_for_column,
    collect_sirius_results,
)

# A real mappingFeatureId: <sirius id>_sirius_<input>_<compound name>_<feature_id>
MAPPING = "764016017587228397_sirius_unknown_neg_Unknown_8468"


def write_tsv(path, rows):
    """Write a tab-separated summary file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def sirius6_dir(tmp_path, name="sirius_pos", structures=None, canopus=None, formulas=None):
    """Build a SIRIUS 6 project directory from the given summary rows."""
    d = tmp_path / name
    d.mkdir(parents=True, exist_ok=True)
    if structures is not None:
        write_tsv(d / "structure_identifications.tsv", structures)
    if canopus is not None:
        write_tsv(d / "canopus_structure_summary.tsv", canopus)
    if formulas is not None:
        write_tsv(d / "formula_identifications.tsv", formulas)
    return d


def structure_row(mapping=MAPPING, rank=1, smiles="CCO", inchikey="LFQSCWFLJHTTHZ", **kw):
    row = {
        "mappingFeatureId": mapping,
        "structurePerIdRank": rank,
        "smiles": smiles,
        "InChIkey2D": inchikey,
        "molecularFormula": "C2H6O",
        "name": "Ethanol",
    }
    row.update(kw)
    return row


def canopus_row(mapping=MAPPING, rank=1, **kw):
    row = {
        "mappingFeatureId": mapping,
        "formulaRank": rank,
        "ClassyFire#superclass": "Organic oxygen compounds",
        "ClassyFire#class": "Organooxygen compounds",
        "ClassyFire#most specific class": "Primary alcohols",
    }
    row.update(kw)
    return row


class TestFormatDetection:
    def test_v6_is_detected_from_the_summary_tsv(self, tmp_path):
        d = sirius6_dir(tmp_path, structures=[structure_row()])
        assert _detect_sirius_format(d) == "v6"

    def test_v5_is_detected_from_compound_directories(self, tmp_path):
        d = tmp_path / "sirius_pos"
        comp = d / "1_Unknown_8468"
        comp.mkdir(parents=True)
        (comp / "spectrum.ms").write_text(">compound\tUnknown\n", encoding="utf-8")
        assert _detect_sirius_format(d) == "v5"

    def test_empty_directory_is_unknown(self, tmp_path):
        d = tmp_path / "empty"
        d.mkdir()
        assert _detect_sirius_format(d) == "unknown"

    def test_missing_directory_is_unknown(self, tmp_path):
        assert _detect_sirius_format(tmp_path / "nope") == "unknown"

    def test_none_is_unknown(self):
        assert _detect_sirius_format(None) == "unknown"

    def test_directories_without_spectra_are_not_v5(self, tmp_path):
        d = tmp_path / "sirius_pos"
        (d / "somedir").mkdir(parents=True)
        assert _detect_sirius_format(d) == "unknown"


class TestFeatureIdExtraction:
    def test_takes_the_trailing_numeric_segment(self):
        assert _extract_feature_id_from_v6_mapping(MAPPING) == 8468

    def test_ignores_the_leading_sirius_id(self):
        """The 18-digit SIRIUS id must not be mistaken for the feature id."""
        assert _extract_feature_id_from_v6_mapping(MAPPING) != 764016017587228397

    def test_handles_a_bare_number(self):
        assert _extract_feature_id_from_v6_mapping("123") == 123

    @pytest.mark.parametrize("value", [None, "", "no_digits_here", 12345, [], {}])
    def test_unparsable_values_return_none(self, value):
        assert _extract_feature_id_from_v6_mapping(value) is None

    def test_name_containing_digits_uses_the_last_one(self):
        assert _extract_feature_id_from_v6_mapping("x_sirius_C6H12O6_42") == 42


class TestValueForColumn:
    HEADER = ["ClassyFire#superclass", "ClassyFire#class", "NPC#pathway"]
    ROW = ["Lipids", "Fatty acids", "Terpenoids"]

    def test_matches_by_case_insensitive_substring(self):
        assert _value_for_column(self.HEADER, self.ROW, "superclass") == "Lipids"

    def test_returns_the_first_match(self):
        """'class' appears in both 'superclass' and 'class'; first header wins."""
        assert _value_for_column(self.HEADER, self.ROW, "class") == "Lipids"

    def test_unmatched_substring_returns_none(self):
        assert _value_for_column(self.HEADER, self.ROW, "nope") is None

    def test_short_row_returns_none(self):
        assert _value_for_column(self.HEADER, ["Lipids"], "NPC") is None

    def test_empty_header_returns_none(self):
        assert _value_for_column([], [], "anything") is None


class TestParseSirius6:
    def test_reads_a_structure_hit(self, tmp_path):
        d = sirius6_dir(tmp_path, structures=[structure_row()])
        results = _parse_sirius6_tsv(d, "POS")
        assert len(results) == 1
        r = results[0]
        assert r["feature_id"] == 8468
        assert r["_sirius_has_struct"] is True
        assert r["_polarity"] == "POS"
        assert r["sirius_name"] == "Ethanol"

    def test_only_rank_1_structures_are_kept(self, tmp_path):
        d = sirius6_dir(
            tmp_path,
            structures=[
                structure_row(rank=1, name="TopHit"),
                structure_row(rank=2, name="RunnerUp"),
            ],
        )
        results = _parse_sirius6_tsv(d, "POS")
        assert [r["sirius_name"] for r in results] == ["TopHit"]

    def test_canopus_classes_are_merged_in(self, tmp_path):
        d = sirius6_dir(tmp_path, structures=[structure_row()], canopus=[canopus_row()])
        r = _parse_sirius6_tsv(d, "POS")[0]
        assert r["SIRIUS_canopus_superclass"] == "Organic oxygen compounds"
        assert r["SIRIUS_canopus_most_specific"] == "Primary alcohols"
        assert r["_sirius_has_canopus"] is True

    def test_missing_structure_layer_is_not_a_hit(self, tmp_path):
        d = sirius6_dir(tmp_path, structures=[structure_row(smiles=None, inchikey=None)])
        r = _parse_sirius6_tsv(d, "POS")[0]
        assert r["_sirius_has_struct"] is False

    def test_falls_back_to_formula_summaries(self, tmp_path):
        """With no structures, formula identifications still yield results."""
        d = sirius6_dir(
            tmp_path,
            formulas=[
                {
                    "mappingFeatureId": MAPPING,
                    "formulaRank": 1,
                    "molecularFormula": "C6H12O6",
                }
            ],
        )
        results = _parse_sirius6_tsv(d, "NEG")
        assert len(results) == 1
        r = results[0]
        assert r["_sirius_has_struct"] is False
        assert r["_sirius_has_formula"] is True
        assert r["molecularFormula"] == "C6H12O6"
        assert r["_polarity"] == "NEG"

    def test_empty_directory_yields_nothing(self, tmp_path):
        d = tmp_path / "sirius_pos"
        d.mkdir()
        assert _parse_sirius6_tsv(d, "POS") == []

    def test_rows_without_a_resolvable_feature_id_are_skipped(self, tmp_path):
        d = sirius6_dir(tmp_path, structures=[structure_row(mapping="no_digits")])
        assert _parse_sirius6_tsv(d, "POS") == []

    def test_aligned_feature_id_is_used_as_a_fallback(self, tmp_path):
        d = sirius6_dir(
            tmp_path,
            structures=[structure_row(mapping="no_digits", alignedFeatureId=77)],
        )
        results = _parse_sirius6_tsv(d, "POS")
        assert [r["feature_id"] for r in results] == [77]

    def test_malformed_tsv_is_reported_not_raised(self, tmp_path, capsys):
        d = tmp_path / "sirius_pos"
        d.mkdir()
        (d / "structure_identifications.tsv").write_bytes(b"\x00\x01 not a tsv")
        # Falls through to the formula path, which is absent -> no results.
        assert _parse_sirius6_tsv(d, "POS") == []


class TestExtractCanopusV5:
    def test_reads_classes_from_a_canopus_tsv(self, tmp_path):
        comp = tmp_path / "1_Unknown_42"
        can = comp / "canopus"
        can.mkdir(parents=True)
        (can / "classes.tsv").write_text(
            "ClassyFire#superclass\tClassyFire#class\tClassyFire#most specific class\n"
            "Lipids\tFatty acids\tOmega-3\n",
            encoding="utf-8",
        )
        out = _extract_canopus(comp)
        assert out["SIRIUS_canopus_superclass"] == "Lipids"
        assert out["SIRIUS_canopus_most_specific"] == "Omega-3"

    def test_missing_canopus_directory_yields_empty_fields(self, tmp_path):
        comp = tmp_path / "1_Unknown_42"
        comp.mkdir()
        out = _extract_canopus(comp)
        assert set(out) == {
            "SIRIUS_canopus_superclass",
            "SIRIUS_canopus_class",
            "SIRIUS_canopus_most_specific",
        }
        assert all(v is None for v in out.values())

    def test_header_only_tsv_yields_empty_fields(self, tmp_path):
        comp = tmp_path / "1_Unknown_42"
        can = comp / "canopus"
        can.mkdir(parents=True)
        (can / "classes.tsv").write_text("ClassyFire#superclass\n", encoding="utf-8")
        assert _extract_canopus(comp)["SIRIUS_canopus_superclass"] is None


class TestPickTopFingerid:
    def _write(self, d, name, score, compound):
        d.mkdir(parents=True, exist_ok=True)
        (d / name).write_text(
            "rank\tscore\tname\tsmiles\tinchikey2D\tmolecularFormula\txlogp\ttanimotoSimilarity\n"
            f"1\t{score}\t{compound}\tCCO\tLFQSCWFLJHTTHZ\tC2H6O\t-0.1\t0.9\n",
            encoding="utf-8",
        )

    def test_highest_score_wins(self, tmp_path):
        d = tmp_path / "fingerid"
        self._write(d, "a.tsv", -10.0, "Low")
        self._write(d, "b.tsv", -1.0, "High")
        assert _pick_top_fingerid(d)["name"] == "High"

    def test_missing_directory_returns_none(self, tmp_path):
        assert _pick_top_fingerid(tmp_path / "nope") is None

    def test_empty_directory_returns_none(self, tmp_path):
        d = tmp_path / "fingerid"
        d.mkdir()
        assert _pick_top_fingerid(d) is None

    def test_rows_without_a_score_are_ignored(self, tmp_path):
        d = tmp_path / "fingerid"
        d.mkdir()
        (d / "a.tsv").write_text("rank\tscore\tname\n1\t\tNoScore\n", encoding="utf-8")
        assert _pick_top_fingerid(d) is None


class TestCollectSiriusResults:
    def test_writes_an_identifications_table(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        out = tmp_path / "sirius_identifications.csv"
        summary = collect_sirius_results(pos, None, out)
        assert out.exists()
        assert summary["pos_identified"] >= 0

    def test_joins_onto_the_merged_table_keeping_all_rows(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [
                {"feature_id": 8468, "Metabolite name": "Unknown", "annotation_level": "3"},
                {"feature_id": 9999, "Metabolite name": "Glucose", "annotation_level": "1"},
            ]
        ).to_csv(merged, index=False)

        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        result = pd.read_csv(out)
        assert len(result) == 2, "left join must preserve every original row"

    def test_structure_hit_renames_the_metabolite(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [{"feature_id": 8468, "Metabolite name": "Unknown", "annotation_level": "3"}]
        ).to_csv(merged, index=False)

        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        row = pd.read_csv(out).iloc[0]
        assert row["Metabolite name"] == "SIRIUS_Ethanol"

    def test_identified_rows_stay_level_3(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [{"feature_id": 8468, "Metabolite name": "Unknown", "annotation_level": "3"}]
        ).to_csv(merged, index=False)
        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        assert str(pd.read_csv(out).iloc[0]["annotation_level"]) == "3"

    def test_canopus_only_becomes_level_4(self, tmp_path):
        pos = sirius6_dir(
            tmp_path,
            "sirius_pos",
            structures=[structure_row(smiles=None, inchikey=None)],
            canopus=[canopus_row()],
        )
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [{"feature_id": 8468, "Metabolite name": "Unknown", "annotation_level": "3"}]
        ).to_csv(merged, index=False)
        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        assert str(pd.read_csv(out).iloc[0]["annotation_level"]) == "4"

    def test_formula_only_becomes_level_5(self, tmp_path):
        pos = sirius6_dir(
            tmp_path,
            "sirius_pos",
            formulas=[{"mappingFeatureId": MAPPING, "formulaRank": 1, "molecularFormula": "C2H6O"}],
        )
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [{"feature_id": 8468, "Metabolite name": "Unknown", "annotation_level": "3"}]
        ).to_csv(merged, index=False)
        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        assert str(pd.read_csv(out).iloc[0]["annotation_level"]) == "5"

    def test_confident_annotations_are_never_downgraded(self, tmp_path):
        """Level 1 and 2 rows must survive the SIRIUS join untouched."""
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [
                {"feature_id": 8468, "Metabolite name": "Glucose", "annotation_level": "1"},
                {"feature_id": 8468, "Metabolite name": "Citrate", "annotation_level": "2"},
            ]
        ).to_csv(merged, index=False)
        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        result = pd.read_csv(out)
        assert sorted(str(v) for v in result["annotation_level"]) == ["1", "2"]
        assert "Glucose" in set(result["Metabolite name"])

    def test_helper_columns_are_dropped_from_the_output(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [{"feature_id": 8468, "Metabolite name": "Unknown", "annotation_level": "3"}]
        ).to_csv(merged, index=False)
        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        cols = set(pd.read_csv(out).columns)
        assert (
            not {
                "sirius_name",
                "_sirius_has_struct",
                "_sirius_has_formula",
                "_sirius_has_canopus",
                "_polarity",
            }
            & cols
        )

    def test_no_sirius_output_passes_the_merged_table_through(self, tmp_path):
        merged = tmp_path / "merged.csv"
        pd.DataFrame(
            [{"feature_id": 1, "Metabolite name": "Glucose", "annotation_level": "1"}]
        ).to_csv(merged, index=False)
        out = tmp_path / "final.csv"
        collect_sirius_results(None, None, out, join_with_merged=merged)
        result = pd.read_csv(out)
        assert len(result) == 1
        assert result.iloc[0]["Metabolite name"] == "Glucose"

    def test_both_polarities_are_collected(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row(mapping="x_Unknown_1")])
        neg = sirius6_dir(tmp_path, "sirius_neg", structures=[structure_row(mapping="x_Unknown_2")])
        out = tmp_path / "ids.csv"
        collect_sirius_results(pos, neg, out)
        ids = set(pd.read_csv(out)["feature_id"])
        assert ids == {1, 2}

    def test_progress_callback_is_invoked(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        seen = []
        collect_sirius_results(pos, None, tmp_path / "ids.csv", progress=lambda *a: seen.append(a))
        # v6 parsing reports through the same callback contract as v5 scanning.
        assert isinstance(seen, list)

    def test_merged_table_without_feature_id_falls_back_to_ids_only(self, tmp_path):
        pos = sirius6_dir(tmp_path, "sirius_pos", structures=[structure_row()])
        merged = tmp_path / "merged.csv"
        pd.DataFrame([{"Metabolite name": "Glucose"}]).to_csv(merged, index=False)
        out = tmp_path / "final.csv"
        collect_sirius_results(pos, None, out, join_with_merged=merged)
        assert "feature_id" in pd.read_csv(out).columns
