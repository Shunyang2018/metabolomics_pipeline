"""Tests for configuration defaults.

Paths used to be hardcoded to two developers' home directories
(`/Users/wangs261/...`, `/Users/ivanablazenovic/...`), so a fresh clone of this
"cross-platform" pipeline had defaults pointing at a specific Mac. These tests
pin the portable behaviour.
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path

import pytest

from metabo_pipeline import constants


def reload_constants(monkeypatch, **env):
    """Re-import constants with the given environment applied."""
    for key in ("METABO_INPUT_DIR", "METABO_BIOACTIVITY_DB", "SIRIUS_EXECUTABLE"):
        monkeypatch.delenv(key, raising=False)
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    return importlib.reload(constants)


@pytest.fixture(autouse=True)
def _restore_constants():
    """Leave the module as the rest of the suite expects to find it."""
    yield
    importlib.reload(constants)


class TestNoMachineSpecificDefaults:
    @pytest.mark.parametrize("name", ["INPUT_DIR", "BIOACTIVITY_DB_PATH", "SIRIUS_EXECUTABLE"])
    def test_default_is_not_an_absolute_user_path(self, monkeypatch, name):
        mod = reload_constants(monkeypatch)
        value = getattr(mod, name)
        assert "/Users/" not in value
        assert "C:\\Users" not in value
        assert not Path(value).is_absolute() or value == "", (
            f"{name} defaults to the absolute path {value!r}, which cannot be "
            "correct on another machine"
        )

    def test_no_developer_usernames_are_committed(self, monkeypatch):
        mod = reload_constants(monkeypatch)
        blob = " ".join(
            str(getattr(mod, n)) for n in ("INPUT_DIR", "BIOACTIVITY_DB_PATH", "SIRIUS_EXECUTABLE")
        )
        for leaked in ("wangs261", "ivanablazenovic"):
            assert leaked not in blob


class TestEnvironmentOverrides:
    def test_input_dir(self, monkeypatch):
        mod = reload_constants(monkeypatch, METABO_INPUT_DIR="/data/msdial")
        assert mod.INPUT_DIR == "/data/msdial"

    def test_bioactivity_db(self, monkeypatch):
        mod = reload_constants(monkeypatch, METABO_BIOACTIVITY_DB="/db/bio.csv")
        assert mod.BIOACTIVITY_DB_PATH == "/db/bio.csv"

    def test_sirius_executable(self, monkeypatch):
        mod = reload_constants(monkeypatch, SIRIUS_EXECUTABLE="/opt/sirius/bin/sirius")
        assert mod.SIRIUS_EXECUTABLE == "/opt/sirius/bin/sirius"

    def test_defaults_apply_when_unset(self, monkeypatch):
        mod = reload_constants(monkeypatch)
        assert mod.INPUT_DIR == "data"
        assert mod.BIOACTIVITY_DB_PATH == "data/bioactivity_database.csv"
        assert mod.SIRIUS_EXECUTABLE == ""


class TestEmptyExecutableTriggersAutodetection:
    def test_guess_falls_back_to_platform_search(self, monkeypatch):
        """An empty SIRIUS_EXECUTABLE must mean 'go look', not 'use empty path'."""
        from metabo_pipeline.sirius_utils import guess_sirius_executable

        monkeypatch.delenv("SIRIUS_EXECUTABLE", raising=False)
        monkeypatch.delenv("SIRIUS_EXE", raising=False)
        resolved, note = guess_sirius_executable("")
        assert resolved
        assert "sirius" in resolved.lower()

    def test_explicit_path_still_wins(self, tmp_path, monkeypatch):
        from metabo_pipeline.sirius_utils import guess_sirius_executable

        exe = tmp_path / "sirius"
        exe.write_text("#!/bin/sh\n", encoding="utf-8")
        exe.chmod(0o755)
        monkeypatch.delenv("SIRIUS_EXECUTABLE", raising=False)
        resolved, note = guess_sirius_executable(str(exe))
        assert resolved == str(exe)
        assert note == "cli-option"

    def test_env_var_is_consulted(self, tmp_path, monkeypatch):
        from metabo_pipeline.sirius_utils import guess_sirius_executable

        exe = tmp_path / "sirius"
        exe.write_text("#!/bin/sh\n", encoding="utf-8")
        exe.chmod(0o755)
        monkeypatch.setenv("SIRIUS_EXECUTABLE", str(exe))
        resolved, note = guess_sirius_executable(None)
        assert resolved == str(exe)
        assert note == "env:SIRIUS_EXECUTABLE"


class TestThresholdsRemainDocumented:
    """The README quotes these defaults; drift would make it wrong."""

    def test_documented_qc_thresholds(self):
        assert constants.MSMS_MIN_IONS == 3
        assert constants.SNR_MIN == 5.0
        assert constants.BLANK_FOLD_MIN == 7.0
        assert constants.PRESENT_PERCENT_MIN == 60.0

    def test_cv_filtering_is_off_by_default(self):
        assert constants.CV_PERCENT_MAX is None

    def test_documented_dedup_window(self):
        assert constants.DEDUP_RT_WINDOW_MIN == 0.6
        assert constants.DEDUP_MZ_PPM == 20.0

    def test_env_does_not_leak_into_thresholds(self, monkeypatch):
        monkeypatch.setenv("METABO_INPUT_DIR", "/elsewhere")
        mod = importlib.reload(constants)
        assert mod.BLANK_FOLD_MIN == 7.0
        assert os.environ["METABO_INPUT_DIR"] == "/elsewhere"
