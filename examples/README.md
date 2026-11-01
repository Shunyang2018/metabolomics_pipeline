# Offline example

From the repository root after `uv sync --locked --extra dev`:

```bash
uv run python examples/run_demo.py --output-dir outputs/demo
```

Open `outputs/demo/report.html`. No HTTP requests or external processes are used.
Use a fresh output directory for each run; the script refuses to overwrite a run.

## Data and expected results

`synthetic/alignments/` contains two small CSVs modeled on the MS-DIAL export
layout: C18 positive and HILIC negative. There are two sample columns and one
blank per file. Sample identities harmonize across the two methods.

All spectra, intensities, names, InChIKey-shaped identifiers, taxonomy labels
and reference records are invented for testing. They do not describe real
molecules or bioactivity. The cache is a synthetic stand-in for the shape of
ClassyFire results, not a downloaded or scientifically verified taxonomy.

| Check | Expected |
| --- | --- |
| Input features | 8 |
| Rejected for sparse MS/MS | 1 |
| Rejected for low S/N | 1 |
| Rejected for blank background | 1 |
| Final features | 5 |
| Pipeline categories 1 / 2 / 3 | 2 / 1 / 2 |
| SIRIUS POS / NEG inputs | 1 / 1 |
| Features with illustrative reference matches | 2 |
| Feature/sample/reference output rows | 4 |

The demo exports unknown spectra but does not run or fabricate SIRIUS results.
`manifest.json` records every fixture's SHA-256 hash, library versions,
thresholds and stage counts. Repeating with the same code, fixtures and
environment produces the same manifest (no timestamps or absolute paths).

## Artifacts to inspect

- `merged.csv`: QC-filtered matrix, `source_file`, original `Alignment ID`,
  harmonized sample intensities and QC metrics.
- `annotated.csv`: the same features plus illustrative cached taxonomy.
- `sirius_unknown_pos.ms` and `sirius_unknown_neg.ms`: actual export format.
- `reference_matches.csv`: only detected, QC-evaluated sample columns; controls
  are excluded. Connectivity matches are review candidates, not activity findings.
- `manifest.json`: reproduction record.
- `report.html`: self-contained feature-retention summary and preview.

Fixture files and this example code are covered by the repository MIT license.
