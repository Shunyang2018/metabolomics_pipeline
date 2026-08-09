"""Tests for reading MS-DIAL alignment-table metadata.

An MS-DIAL export carries four metadata rows (Class, File type, Injection
order, Batch ID) above the real header row, and sample columns begin after
'MS/MS spectrum'.
"""

import pytest

from metabo_pipeline.io_msdial import summarize_alignment_table

FIXED_COLS = [
    "Alignment ID",
    "Average Rt(min)",
    "Average Mz",
    "Metabolite name",
    "MS/MS spectrum",
]
SAMPLES = ["M2_C18_crude_1_pos", "M2_C18_crude_2_pos", "M2_C18_blank_pos"]


def write_table(path, samples=SAMPLES, n_features=3, injection_order=None, pad=0):
    """Write a minimal but structurally valid MS-DIAL alignment table."""
    blank_lead = [""] * (len(FIXED_COLS) - pad)
    order = injection_order or [str(i + 1) for i in range(len(samples))]
    rows = [
        blank_lead + ["Sample"] * (len(samples) - 1) + ["Blank"],
        blank_lead + ["Sample"] * (len(samples) - 1) + ["Blank"],
        blank_lead + order,
        blank_lead + ["1"] * len(samples),
        FIXED_COLS + samples,
    ]
    for i in range(n_features):
        rows.append([str(i), "1.5", "180.06", "Glucose", "100:200"] + ["1000"] * len(samples))
    path.write_text("\n".join(",".join(r) for r in rows), encoding="utf-8")
    return path


def test_reads_sample_names_after_msms_column(tmp_path):
    summary = summarize_alignment_table(write_table(tmp_path / "a.csv"))
    assert summary.samples == SAMPLES


def test_counts_feature_rows_excluding_headers(tmp_path):
    summary = summarize_alignment_table(write_table(tmp_path / "a.csv", n_features=7))
    assert summary.n_features == 7


def test_maps_metadata_per_sample(tmp_path):
    summary = summarize_alignment_table(write_table(tmp_path / "a.csv"))
    assert summary.metadata.classes["M2_C18_blank_pos"] == "Blank"
    assert summary.metadata.classes["M2_C18_crude_1_pos"] == "Sample"
    assert summary.metadata.injection_order["M2_C18_crude_1_pos"] == 1
    assert summary.metadata.batch_id["M2_C18_crude_1_pos"] == "1"


def test_records_the_source_path(tmp_path):
    path = write_table(tmp_path / "a.csv")
    assert summarize_alignment_table(path).path == path


def test_too_few_header_rows_raises_a_clear_error(tmp_path):
    path = tmp_path / "short.csv"
    path.write_text("a,b\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Not enough header rows"):
        summarize_alignment_table(path)


def test_empty_file_raises(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="Not enough header rows"):
        summarize_alignment_table(path)


def test_table_with_no_feature_rows_reports_zero(tmp_path):
    summary = summarize_alignment_table(write_table(tmp_path / "a.csv", n_features=0))
    assert summary.n_features == 0
    assert summary.samples == SAMPLES


def test_blank_injection_order_does_not_crash(tmp_path):
    """A blank cell in the injection-order row previously raised ValueError."""
    path = write_table(tmp_path / "a.csv", injection_order=["1", "", "3"])
    summary = summarize_alignment_table(path)
    order = summary.metadata.injection_order
    assert order["M2_C18_crude_1_pos"] == 1
    assert "M2_C18_crude_2_pos" not in order


def test_non_numeric_injection_order_does_not_crash(tmp_path):
    path = write_table(tmp_path / "a.csv", injection_order=["1", "n/a", "3"])
    summary = summarize_alignment_table(path)
    assert "M2_C18_crude_2_pos" not in summary.metadata.injection_order


def test_utf8_bom_is_tolerated(tmp_path):
    """MS-DIAL on Windows writes a BOM; the first column name must still match."""
    path = write_table(tmp_path / "bom.csv")
    raw = path.read_bytes()
    path.write_bytes(b"\xef\xbb\xbf" + raw)
    assert summarize_alignment_table(path).samples == SAMPLES


def test_metadata_row_shorter_than_header_is_tolerated(tmp_path):
    """Truncated metadata rows yield fewer entries rather than raising."""
    path = write_table(tmp_path / "a.csv")
    lines = path.read_text(encoding="utf-8").splitlines()
    # drop the final cell from the Class row
    lines[0] = ",".join(lines[0].split(",")[:-1])
    path.write_text("\n".join(lines), encoding="utf-8")
    summary = summarize_alignment_table(path)
    assert len(summary.metadata.classes) == len(SAMPLES) - 1


def test_missing_msms_column_names_the_problem(tmp_path):
    """Sample columns are located relative to 'MS/MS spectrum'."""
    path = write_table(tmp_path / "a.csv")
    lines = path.read_text(encoding="utf-8").splitlines()
    lines[4] = lines[4].replace("MS/MS spectrum", "Something Else")
    path.write_text("\n".join(lines), encoding="utf-8")
    with pytest.raises(ValueError, match="MS/MS spectrum"):
        summarize_alignment_table(path)
