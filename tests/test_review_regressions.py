"""Regression cases found while preparing the public reproducible example."""

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest
from msdial import write_msdial
from typer.testing import CliRunner

from metabo_pipeline.bioactivity import identify_sample_columns, match_bioactives
from metabo_pipeline.classify import classify_level12_with_classyfire
from metabo_pipeline.cli import app
from metabo_pipeline.cli.run_cmd import _create_final_merge
from metabo_pipeline.io_msdial import summarize_alignment_table
from metabo_pipeline.merge import merge_folder_to_wide_csv
from metabo_pipeline.reporting import write_report
from metabo_pipeline.sirius_export import write_ms_files

runner = CliRunner()


def test_unknown_sample_vocabulary_does_not_bypass_qc(tmp_path):
    write_msdial(tmp_path / "C18_pos.csv", samples=["case_pos", "blank_pos"])
    with pytest.raises(ValueError, match="No real sample groups"):
        merge_folder_to_wide_csv(tmp_path, tmp_path / "outputs" / "merged.csv")


def test_explicit_sample_token_enables_other_datasets(tmp_path):
    write_msdial(tmp_path / "C18_pos.csv", samples=["case_pos", "blank_pos"])
    out = tmp_path / "outputs" / "merged.csv"
    result = runner.invoke(app, ["merge", str(tmp_path), "-o", str(out), "--sample-token", "case"])
    assert result.exit_code == 0, result.output
    assert len(pd.read_csv(out)) == 1


def test_colliding_sample_names_are_rejected(tmp_path):
    write_msdial(tmp_path / "C18_pos.csv", samples=["case_C18_pos", "case_HILIC_neg", "blank"])
    with pytest.raises(ValueError, match="Ambiguous normalized"):
        merge_folder_to_wide_csv(tmp_path, tmp_path / "outputs" / "merged.csv")


def test_recursive_merge_does_not_read_its_own_outputs(tmp_path):
    write_msdial(tmp_path / "C18_pos.csv")
    out = tmp_path / "outputs" / "merged.csv"
    merge_folder_to_wide_csv(tmp_path, out, recursive=True)
    again = merge_folder_to_wide_csv(tmp_path, out, recursive=True)
    assert again["files_found"] == 1
    assert again["rows"] == 1
    assert pd.read_csv(out)["source_file"].tolist() == ["C18_pos.csv"]


def test_empty_polarity_clears_a_previous_sirius_export(tmp_path):
    write_ms_files([">compound old"], [">compound old"], tmp_path)
    write_ms_files([], [], tmp_path)
    assert (tmp_path / "sirius_unknown_pos.ms").read_text() == ""
    assert (tmp_path / "sirius_unknown_neg.ms").read_text() == ""


def test_tab_delimited_validation_matches_merge_contract(tmp_path):
    path = write_msdial(tmp_path / "C18_pos.txt", sep="\t")
    assert summarize_alignment_table(path).n_features == 1
    result = runner.invoke(app, ["validate", str(path)])
    assert result.exit_code == 0, result.output
    assert json.loads(result.output)["n_features"] == 1


@pytest.mark.parametrize("keys", [[None], ["AAAAAAAAAAAAAA-BBBBBBBBBB-C"], []])
def test_offline_classification_with_no_cache_hits_preserves_rows(tmp_path, keys, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("offline classification contacted the network")

    monkeypatch.setattr("requests.get", forbidden)
    src, dst = tmp_path / "in.csv", tmp_path / "nested" / "out.csv"
    pd.DataFrame(
        {"feature_id": range(len(keys)), "annotation_level": ["1"] * len(keys), "INCHIKEY": keys}
    ).to_csv(src, index=False)
    stats = classify_level12_with_classyfire(
        str(src), str(dst), offline=True, cache_path=str(tmp_path / "cache.json")
    )
    assert len(pd.read_csv(dst)) == len(keys)
    assert stats["hits"] == 0
    assert stats["unique_keys"] == sum(k is not None for k in keys)


def test_reference_matches_exclude_controls_and_new_metadata():
    cols = [
        "source_file",
        "new_numeric_metadata",
        "blank",
        "mb_1",
        "qc_pool",
        "case_1",
        "present_percent_case_1",
    ]
    assert identify_sample_columns(cols) == ["case_1"]


def test_match_summary_counts_only_actual_reference_hits(tmp_path):
    feats, db, out = (tmp_path / name for name in ("feats.csv", "db.csv", "out.csv"))
    pd.DataFrame(
        [
            {"feature_id": 1, "INCHIKEY": "AAAAAAAAAAAAAA-X", "case_1": 10},
            {"feature_id": 2, "INCHIKEY": "BBBBBBBBBBBBBB-X", "case_1": 10},
        ]
    ).to_csv(feats, index=False)
    pd.DataFrame([{"InChIKey": "AAAAAAAAAAAAAA-Y"}]).to_csv(db, index=False)
    stats = match_bioactives(feats, db, out)
    assert stats["unique_matched_skeletons"] == 1


def test_skipped_stages_do_not_reuse_old_classifications_or_predictions(tmp_path):
    src = tmp_path / "merged.csv"
    pd.DataFrame([{"feature_id": 1, "Metabolite name": "Current", "annotation_level": 3}]).to_csv(
        src, index=False
    )
    pd.DataFrame([{"feature_id": 1, "Metabolite name": "Stale"}]).to_csv(
        tmp_path / "merged_classyfire.csv", index=False
    )
    (tmp_path / "sirius_identifications.csv").write_text("stale", encoding="utf-8")
    out = _create_final_merge(str(src), str(tmp_path), use_classified=False, use_sirius=False)
    assert pd.read_csv(out)["Metabolite name"].tolist() == ["Current"]


def test_run_failure_does_not_announce_completion(tmp_path, monkeypatch):
    import metabo_pipeline.cli.run_cmd as run_cmd

    source = tmp_path / "input"
    source.mkdir()
    write_msdial(source / "C18_pos.csv")
    monkeypatch.setattr(run_cmd, "_run_parallel_tasks", lambda *args: {"sirius": (False, "failed")})
    result = runner.invoke(
        app,
        [
            "run",
            "--input-dir",
            str(source),
            "--output-dir",
            str(tmp_path / "out"),
            "--skip-bioactivity",
        ],
    )
    assert result.exit_code == 1
    assert "PIPELINE COMPLETE" not in result.output
    assert not (tmp_path / "out" / "merged_classyfire_final.csv").exists()


def test_full_offline_run_produces_current_final_and_reference_matches(tmp_path):
    source = Path(__file__).resolve().parents[1] / "examples" / "synthetic"
    out = tmp_path / "run"
    result = runner.invoke(
        app,
        [
            "run",
            "--input-dir",
            str(source / "alignments"),
            "--output-dir",
            str(out),
            "--skip-classify",
            "--skip-sirius",
            "--bioactivity-db",
            str(source / "reference.csv"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert len(pd.read_csv(out / "merged_final.csv")) == 5
    assert set(pd.read_csv(out / "bioactives.csv")["sample_type"]) == {
        "demo_crude_1",
        "demo_crude_2",
    }
    # A second unqualified run cannot mix with these existing artifacts.
    again = runner.invoke(
        app,
        [
            "run",
            "--input-dir",
            str(source / "alignments"),
            "--output-dir",
            str(out),
            "--skip-classify",
            "--skip-sirius",
        ],
    )
    assert again.exit_code == 2


def test_public_demo_is_offline_reproducible_and_escaped(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("the public demo must not contact services")

    monkeypatch.setattr("requests.get", forbidden)
    script = Path(__file__).resolve().parents[1] / "examples" / "run_demo.py"
    spec = importlib.util.spec_from_file_location("public_demo", script)
    demo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo)
    manifest = demo.run_demo(tmp_path / "demo")
    assert manifest["merge"]["totals"] == {
        "raw": 8,
        "after_msms": 7,
        "after_snr": 6,
        "after_pass": 5,
        "ann1": 2,
        "ann2": 1,
        "ann3": 2,
    }
    assert manifest["merge"]["rows"] == 5
    assert manifest["merge"]["sirius_l3_total"] == 2
    assert manifest["bioactivity"]["unique_features_matched"] == 2
    assert manifest["bioactivity"]["output_rows"] == 4
    assert all(len(entry["sha256"]) == 64 for entry in manifest["inputs"])
    assert demo.run_demo(tmp_path / "second") == manifest
    write_report(
        pd.DataFrame({"Metabolite name": ["<script>alert(1)</script>"]}),
        manifest,
        tmp_path / "safe.html",
        synthetic=True,
    )
    html = (tmp_path / "safe.html").read_text(encoding="utf-8")
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


@pytest.mark.parametrize("failure", ["positive_process", "summary_export"])
def test_sirius_partial_failure_returns_nonzero(tmp_path, monkeypatch, failure):
    import metabo_pipeline.cli.sirius_cmd as sirius_cmd

    monkeypatch.setattr(sirius_cmd, "_check_sirius_login", lambda *args: None)
    monkeypatch.setattr(
        sirius_cmd, "export_tsv_summaries", lambda *args, **kw: failure != "summary_export"
    )
    for mode in ("pos", "neg"):
        (tmp_path / f"sirius_unknown_{mode}.ms").write_text(
            ">compound Unknown_1\n>parentmass 200\n", encoding="utf-8"
        )

    class FakeProcess:
        def __init__(self, command, **kwargs):
            project = Path(command[command.index("-o") + 1])
            project.write_bytes(b"synthetic test project" * 100)
            self.returncode = 1 if failure == "positive_process" and "pos" in project.name else 0
            self.stdout = iter(["Finished\n"])
            self.pid = 1

        def wait(self):
            return self.returncode

    monkeypatch.setattr(sirius_cmd.subprocess, "Popen", FakeProcess)
    result = runner.invoke(
        app, ["sirius", "--output-dir", str(tmp_path), "--sirius-exe", "mock-sirius"]
    )
    assert result.exit_code == 1, result.output
    assert "SIRIUS runs finished" not in result.output


def test_short_unrelated_csv_does_not_break_a_valid_merge(tmp_path):
    write_msdial(tmp_path / "C18_pos.csv")
    (tmp_path / "notes.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    stats = merge_folder_to_wide_csv(tmp_path, tmp_path / "outputs" / "merged.csv")
    assert stats["files"] == 1
    assert stats["files_skipped"] == ["notes.csv"]


def test_merge_output_can_be_in_the_parent_of_input(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    write_msdial(inputs / "C18_pos.csv")
    stats = merge_folder_to_wide_csv(inputs, tmp_path / "merged.csv")
    assert stats["rows"] == 1


def test_exact_duplicate_headers_are_rejected_before_pandas_mangles_them(tmp_path):
    write_msdial(tmp_path / "C18_pos.csv", samples=["crude_1", "crude_1", "blank"])
    with pytest.raises(ValueError, match="Duplicate raw sample names"):
        merge_folder_to_wide_csv(tmp_path, tmp_path / "outputs" / "merged.csv")
