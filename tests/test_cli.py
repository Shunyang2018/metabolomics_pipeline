"""Tests for the `metabo` command line interface.

These cover the Typer wiring, argument/default resolution and exit codes. No
external service is contacted: SIRIUS and the ClassyFire API are never reached,
and every command is given explicit paths so the machine-specific defaults in
constants.py are not consulted.
"""

from __future__ import annotations

import pandas as pd
import pytest
from msdial import feature, write_msdial
from typer.testing import CliRunner

from metabo_pipeline.cli import app

runner = CliRunner()

GLUCOSE = "WQZGKKKJIJFFOK-GASJEMHNSA-N"

ALL_COMMANDS = [
    "run",
    "merge",
    "classify",
    "classify-check",
    "sirius",
    "sirius-collect",
    "final",
    "bioactivity",
]


class TestAppWiring:
    def test_top_level_help(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "Metabolomics pipeline" in result.output

    def test_every_command_is_registered(self):
        result = runner.invoke(app, ["--help"])
        for name in ALL_COMMANDS:
            assert name in result.output, f"{name} missing from --help"

    @pytest.mark.parametrize("command", ALL_COMMANDS)
    def test_command_help_renders(self, command):
        """A broken Typer signature only shows up when --help is rendered."""
        result = runner.invoke(app, [command, "--help"])
        assert result.exit_code == 0, result.output

    def test_unknown_command_is_rejected(self):
        result = runner.invoke(app, ["no-such-command"])
        assert result.exit_code != 0


class TestMergeCommand:
    @pytest.fixture
    def in_dir(self, tmp_path):
        d = tmp_path / "in"
        d.mkdir()
        write_msdial(
            d / "M2_C18_pos.csv",
            ["M2_C18_crude_1_pos", "M2_C18_crude_2_pos", "M2_C18_blank_pos"],
            [feature(0), feature(1, rt=5.0, mz=300.2, name="Citrate")],
        )
        return d

    def test_merges_and_writes_the_requested_output(self, in_dir, tmp_path):
        out = tmp_path / "out" / "merged.csv"
        result = runner.invoke(app, ["merge", str(in_dir), "--output", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert not pd.read_csv(out).empty

    def test_reports_rows_written(self, in_dir, tmp_path):
        out = tmp_path / "out" / "merged.csv"
        result = runner.invoke(app, ["merge", str(in_dir), "--output", str(out)])
        assert "Rows written" in result.output

    def test_short_output_flag(self, in_dir, tmp_path):
        out = tmp_path / "out" / "merged.csv"
        result = runner.invoke(app, ["merge", str(in_dir), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_missing_input_directory_exits_2(self, tmp_path):
        result = runner.invoke(app, ["merge", str(tmp_path / "nope")])
        assert result.exit_code == 2
        assert "not found" in result.output

    def test_a_file_given_instead_of_a_directory_exits_2(self, tmp_path):
        f = tmp_path / "a.csv"
        f.write_text("x", encoding="utf-8")
        result = runner.invoke(app, ["merge", str(f)])
        assert result.exit_code == 2

    def test_recursive_flag_descends(self, tmp_path):
        d = tmp_path / "in"
        (d / "sub").mkdir(parents=True)
        write_msdial(d / "sub" / "M2_C18_pos.csv")
        out = tmp_path / "out" / "merged.csv"

        shallow = runner.invoke(app, ["merge", str(d), "--output", str(out)])
        assert "Merged 0 files" in shallow.output

        deep = runner.invoke(app, ["merge", str(d), "--output", str(out), "--recursive"])
        assert deep.exit_code == 0, deep.output
        assert "Merged 1 files" in deep.output

    def test_skipped_files_are_not_counted_as_merged(self, in_dir, tmp_path):
        """The count printed to the user must exclude non-MS-DIAL tables."""
        (in_dir / "notes.csv").write_text("a,b\n1,2\n3,4\n5,6\n7,8\n9,10\n", encoding="utf-8")
        out = tmp_path / "out" / "merged.csv"
        result = runner.invoke(app, ["merge", str(in_dir), "--output", str(out)])
        assert "Merged 1 files" in result.output


class TestBioactivityCommand:
    @pytest.fixture
    def merged_csv(self, tmp_path):
        path = tmp_path / "outputs" / "merged.csv"
        path.parent.mkdir(parents=True)
        pd.DataFrame([{"feature_id": "F1", "INCHIKEY": GLUCOSE, "crude_1": 1000}]).to_csv(
            path, index=False
        )
        return path

    @pytest.fixture
    def db_csv(self, tmp_path):
        path = tmp_path / "db.csv"
        pd.DataFrame([{"InChIKey": GLUCOSE, "Benefit": "antioxidant"}]).to_csv(path, index=False)
        return path

    def test_matches_and_writes_output(self, merged_csv, db_csv, tmp_path):
        out = tmp_path / "bioactives.csv"
        result = runner.invoke(
            app,
            ["bioactivity", str(merged_csv), "--db", str(db_csv), "--output", str(out)],
        )
        assert result.exit_code == 0, result.output
        assert len(pd.read_csv(out)) == 1

    def test_reports_match_counts(self, merged_csv, db_csv, tmp_path):
        out = tmp_path / "bioactives.csv"
        result = runner.invoke(
            app,
            ["bioactivity", str(merged_csv), "--db", str(db_csv), "--output", str(out)],
        )
        assert "Features with InChIKey" in result.output

    def test_missing_database_exits_1(self, merged_csv, tmp_path):
        result = runner.invoke(
            app, ["bioactivity", str(merged_csv), "--db", str(tmp_path / "nope.csv")]
        )
        assert result.exit_code == 1
        assert "not found" in result.output

    def test_discovers_merged_csv_in_the_output_dir(self, merged_csv, db_csv, tmp_path):
        """With no positional argument it falls back to OUTPUT_DIR discovery."""
        out = tmp_path / "bioactives.csv"
        result = runner.invoke(
            app,
            [
                "bioactivity",
                "--output-dir",
                str(merged_csv.parent),
                "--db",
                str(db_csv),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_prefers_the_most_processed_table(self, tmp_path, db_csv):
        """merged_classyfire_final.csv wins over merged.csv when both exist."""
        out_dir = tmp_path / "outputs"
        out_dir.mkdir()
        for name, fid in [
            ("merged.csv", "PLAIN"),
            ("merged_classyfire_final.csv", "FINAL"),
        ]:
            pd.DataFrame([{"feature_id": fid, "INCHIKEY": GLUCOSE, "crude_1": 1000}]).to_csv(
                out_dir / name, index=False
            )

        out = tmp_path / "bioactives.csv"
        result = runner.invoke(
            app,
            [
                "bioactivity",
                "--output-dir",
                str(out_dir),
                "--db",
                str(db_csv),
                "--output",
                str(out),
            ],
        )
        assert result.exit_code == 0, result.output
        assert set(pd.read_csv(out)["feature_id"]) == {"FINAL"}

    def test_no_merged_csv_anywhere_exits_1(self, tmp_path, db_csv):
        empty = tmp_path / "empty"
        empty.mkdir()
        result = runner.invoke(
            app, ["bioactivity", "--output-dir", str(empty), "--db", str(db_csv)]
        )
        assert result.exit_code == 1
        assert "No merged CSV found" in result.output


class TestClassifyCheckCommand:
    def test_reports_a_reachable_endpoint(self, monkeypatch):
        import requests

        class Resp:
            status_code = 200
            text = "{}"

            def raise_for_status(self):
                pass

        monkeypatch.setattr(requests, "get", lambda *a, **k: Resp())
        result = runner.invoke(app, ["classify-check"])
        assert result.exit_code == 0, result.output

    def test_unreachable_endpoint_does_not_traceback(self, monkeypatch):
        import requests

        def boom(*a, **k):
            raise requests.ConnectionError("unreachable")

        monkeypatch.setattr(requests, "get", boom)
        result = runner.invoke(app, ["classify-check"])
        # It may exit non-zero, but must not surface an unhandled exception.
        assert not isinstance(result.exception, requests.ConnectionError)


class TestClassifyCommand:
    def test_missing_input_file_is_reported(self, tmp_path):
        result = runner.invoke(app, ["classify", str(tmp_path / "nope.csv")])
        assert result.exit_code != 0

    def test_input_without_required_columns_is_reported(self, tmp_path):
        bad = tmp_path / "bad.csv"
        pd.DataFrame([{"something": 1}]).to_csv(bad, index=False)
        result = runner.invoke(app, ["classify", str(bad), "--offline"])
        assert result.exit_code != 0

    def test_offline_mode_needs_no_network(self, tmp_path, monkeypatch):
        """--offline must serve from cache only; a network call would be a bug."""
        import requests

        def forbidden(*a, **k):
            raise AssertionError("--offline must not contact the network")

        monkeypatch.setattr(requests, "get", forbidden)
        monkeypatch.setattr(requests, "post", forbidden)

        src = tmp_path / "merged.csv"
        pd.DataFrame([{"feature_id": "F1", "annotation_level": "1", "INCHIKEY": GLUCOSE}]).to_csv(
            src, index=False
        )
        cache = tmp_path / "cache.json"
        cache.write_text("{}", encoding="utf-8")

        result = runner.invoke(
            app,
            [
                "classify",
                str(src),
                "--offline",
                "--cache-path",
                str(cache),
                "--output",
                str(tmp_path / "out.csv"),
            ],
        )
        assert not isinstance(result.exception, AssertionError), result.exception


class TestFinalCommand:
    def test_no_inputs_warns_that_the_table_will_be_empty(self, tmp_path):
        """Succeeding silently here makes an upstream failure look like success.

        The step still exits 0 on purpose: a dataset with no Level 3 unknowns
        legitimately has no SIRIUS output to join.
        """
        result = runner.invoke(app, ["final", "--output-dir", str(tmp_path)])
        assert result.exit_code == 0, result.output
        assert "will be empty" in result.output

    def test_joins_an_existing_merged_table(self, tmp_path):
        pd.DataFrame(
            [{"feature_id": "F1", "Metabolite name": "Glucose", "annotation_level": "1"}]
        ).to_csv(tmp_path / "merged.csv", index=False)
        result = runner.invoke(app, ["final", "--output-dir", str(tmp_path)])
        assert result.exit_code == 0, result.output
        assert "will be empty" not in result.output
        assert (tmp_path / "merged_final.csv").exists()


class TestSiriusCommand:
    @pytest.fixture
    def ms_dir(self, tmp_path):
        out_dir = tmp_path / "outputs"
        out_dir.mkdir()
        (out_dir / "sirius_unknown_pos.ms").write_text(
            ">compound\tUnknown_F1\n>parentmass\t200.1\n\n"
            ">compound\tUnknown_F2\n>parentmass\t300.1\n",
            encoding="utf-8",
        )
        return out_dir

    def test_dry_run_needs_no_sirius_install(self, ms_dir):
        """--dry-run is documented as producing output without running SIRIUS.

        It previously ignored the flag entirely and shelled out to a login
        check, dying with FileNotFoundError wherever SIRIUS was not installed.
        """
        result = runner.invoke(app, ["sirius", "--dry-run", "--output-dir", str(ms_dir)])
        assert result.exit_code == 0, result.output
        assert not isinstance(result.exception, FileNotFoundError)

    def test_dry_run_skips_the_login_check(self, ms_dir):
        result = runner.invoke(app, ["sirius", "--dry-run", "--output-dir", str(ms_dir)])
        assert "skipping SIRIUS login check" in result.output

    def test_dry_run_counts_the_compounds_it_would_process(self, ms_dir):
        result = runner.invoke(app, ["sirius", "--dry-run", "--output-dir", str(ms_dir)])
        assert "2 compounds" in result.output

    def test_dry_run_creates_the_output_directory(self, ms_dir):
        runner.invoke(app, ["sirius", "--dry-run", "--output-dir", str(ms_dir)])
        assert (ms_dir / "sirius_pos.sirius").exists()

    def test_dry_run_reports_missing_ms_files_without_failing(self, ms_dir):
        """Only POS exists in the fixture, so NEG must be reported and skipped."""
        result = runner.invoke(app, ["sirius", "--dry-run", "--output-dir", str(ms_dir)])
        assert "Missing .ms file" in result.output
        assert result.exit_code == 0

    def test_missing_executable_exits_cleanly_without_dry_run(self, ms_dir):
        """A real run with no SIRIUS installed must not raise FileNotFoundError."""
        result = runner.invoke(
            app,
            ["sirius", "--output-dir", str(ms_dir), "--sirius-exe", str(ms_dir / "nope")],
        )
        assert not isinstance(result.exception, FileNotFoundError)
        assert result.exit_code == 2
        assert "not found" in result.output


class TestSiriusCollectCommand:
    def test_missing_result_directory_is_reported(self, tmp_path):
        result = runner.invoke(app, ["sirius-collect", str(tmp_path / "nope")])
        assert result.exit_code != 0
        assert result.exception is None or isinstance(result.exception, SystemExit)
