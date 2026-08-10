"""End-to-end tests for the merge step.

Exercises the documented pipeline: read MS-DIAL tables (comma- or
tab-delimited), filter on MS/MS ion count and S/N, compute replicate-group QC,
harmonize sample columns across assays, and write SIRIUS inputs for Level 3
unknowns.
"""

import pandas as pd
import pytest
from msdial import feature, write_msdial

from metabo_pipeline.merge import list_alignment_files, merge_folder_to_wide_csv


@pytest.fixture
def input_dir(tmp_path):
    d = tmp_path / "in"
    d.mkdir()
    return d


class TestListAlignmentFiles:
    def test_finds_csv_and_txt(self, input_dir):
        (input_dir / "a.csv").write_text("x", encoding="utf-8")
        (input_dir / "b.txt").write_text("x", encoding="utf-8")
        (input_dir / "c.xlsx").write_text("x", encoding="utf-8")
        found = [p.name for p in list_alignment_files(input_dir)]
        assert found == ["a.csv", "b.txt"]

    def test_non_recursive_by_default(self, input_dir):
        sub = input_dir / "sub"
        sub.mkdir()
        (sub / "deep.csv").write_text("x", encoding="utf-8")
        assert list_alignment_files(input_dir) == []
        assert len(list_alignment_files(input_dir, recursive=True)) == 1


class TestMerge:
    def _run(self, input_dir, tmp_path):
        out = tmp_path / "outputs" / "merged.csv"
        out.parent.mkdir(exist_ok=True)
        summary = merge_folder_to_wide_csv(input_dir, out)
        return summary, out

    def test_merges_a_single_file(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_crude_2_pos", "M2_C18_blank_pos"],
            [feature(0), feature(1, rt=5.0, mz=300.2, name="Citrate")],
        )
        summary, out = self._run(input_dir, tmp_path)
        assert summary["files"] == 1
        assert out.exists()
        assert not pd.read_csv(out).empty

    def test_features_with_too_few_msms_ions_are_dropped(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [
                feature(0, name="Keeper"),
                feature(1, name="TooFewIons", msms="100.05:200"),
            ],
        )
        _, out = self._run(input_dir, tmp_path)
        names = set(pd.read_csv(out)["Metabolite name"])
        assert "Keeper" in names
        assert "TooFewIons" not in names

    def test_low_signal_to_noise_features_are_dropped(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [
                feature(0, name="Keeper", snr=50.0),
                feature(1, name="Noisy", snr=1.0, rt=5.0, mz=300.0),
            ],
        )
        _, out = self._run(input_dir, tmp_path)
        names = set(pd.read_csv(out)["Metabolite name"])
        assert "Keeper" in names
        assert "Noisy" not in names

    def test_annotation_level_and_provenance_columns_are_added(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0)],
        )
        _, out = self._run(input_dir, tmp_path)
        cols = pd.read_csv(out).columns
        for expected in ("feature_id", "chrom", "mode", "annotation_level"):
            assert expected in cols

    def test_chrom_and_mode_come_from_the_filename(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_Hilic_neg.csv",
            ["M2_Hilic_crude_1_neg", "M2_Hilic_blank_neg"],
            [feature(0)],
        )
        _, out = self._run(input_dir, tmp_path)
        df = pd.read_csv(out)
        assert set(df["chrom"]) == {"HILIC"}
        assert set(df["mode"]) == {"NEG"}

    def test_level_3_names_are_replaced_with_unknown(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0, name="Low score hit", wdot=0.2, rdot=0.1, peaks=0)],
        )
        _, out = self._run(input_dir, tmp_path)
        df = pd.read_csv(out)
        l3 = df[df["annotation_level"].astype(str) == "3"]
        assert not l3.empty, "expected the weak-scoring feature to be called Level 3"
        assert set(l3["Metabolite name"]) == {"Unknown"}

    def test_sample_columns_harmonize_across_assays(self, input_dir, tmp_path):
        """The same biological sample run on C18 and HILIC yields one column."""
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0)],
        )
        write_msdial(
            input_dir / "M2_Hilic_neg.csv",
            ["M2_Hilic_crude_1_neg", "M2_Hilic_blank_neg"],
            [feature(0, rt=2.0, mz=250.0, name="Citrate")],
        )
        summary, out = self._run(input_dir, tmp_path)
        cols = list(pd.read_csv(out).columns)
        assert summary["files"] == 2
        assert "m2_crude_1" in cols
        assert "M2_C18_crude_1_pos" not in cols

    def test_tab_delimited_txt_exports_are_read(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.txt",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0)],
            sep="\t",
        )
        summary, out = self._run(input_dir, tmp_path)
        assert summary["files"] == 1
        assert not pd.read_csv(out).empty

    def test_non_msdial_files_are_skipped_not_fatal(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0)],
        )
        (input_dir / "notes.csv").write_text("a,b\n1,2\n3,4\n5,6\n7,8\n9,10\n", encoding="utf-8")
        summary, out = self._run(input_dir, tmp_path)
        assert summary["files"] == 1
        assert out.exists()

    def test_empty_input_directory_does_not_crash(self, input_dir, tmp_path):
        summary, _ = self._run(input_dir, tmp_path)
        assert summary["files"] == 0

    def test_sirius_inputs_are_written_for_unknowns(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0, name="Unknown", wdot=0.2, rdot=0.1, peaks=0)],
        )
        _, out = self._run(input_dir, tmp_path)
        # Written next to the merged output, per the README.
        assert (out.parent / "sirius_unknown_pos.ms").exists()

    def test_summary_reports_per_stage_counts(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0), feature(1, msms="100:1", rt=5.0, mz=300.0)],
        )
        summary, _ = self._run(input_dir, tmp_path)
        totals = summary.get("totals", summary)
        assert totals["raw"] == 2
        assert totals["after_msms"] == 1

    def test_progress_callback_is_invoked(self, input_dir, tmp_path):
        write_msdial(
            input_dir / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_blank_pos"],
            [feature(0)],
        )
        seen = []
        out = tmp_path / "outputs" / "merged.csv"
        out.parent.mkdir(exist_ok=True)
        merge_folder_to_wide_csv(input_dir, out, progress=seen.append)
        assert seen
