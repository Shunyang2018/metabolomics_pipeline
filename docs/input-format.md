# Input contract

Inputs are **MS-DIAL alignment exports**, after peak picking and alignment.
Instrument raw files, mzML, vendor formats and unaligned peak lists are not read.

## Layout

1. Class metadata row.
2. File-type metadata row.
3. Injection-order metadata row.
4. Batch-ID metadata row.
5. Column names.
6. Feature rows.

CSV uses commas; TXT uses tabs. Text is UTF-8 (a BOM is accepted).
`metabo validate PATH` inspects this layout and reports sample metadata; it is
a structural inspection, not a scientific validation of spectra or scores.

The merge requires `Alignment ID` and `MS/MS spectrum`. For meaningful annotation
and export, provide `Average Rt(min)`, `Average Mz`, `Metabolite name`,
`Adduct type`, `Weighted dot product`, `Reverse dot product`, `Matched peaks count`,
`MS1 isotopic spectrum`, and optionally `S/N average`, `INCHIKEY` and `SMILES`.
Spectra are whitespace-separated `mz:intensity` pairs. Scores are assumed to be
on a **0–1** scale; exports on another scale must be converted before use.

Every column after `MS/MS spectrum` is treated as an intensity column, so no
additional feature metadata may follow it. Filenames must contain `pos` or `neg`
(case-insensitive); C18/HILIC/Lipid tokens determine the chromatographic method.
Do not mix CSVs from unrelated schemas in the input folder. Short/non-MS-DIAL
files are skipped; a folder without usable alignment tables raises an error.

## Sample names and groups

Normalization lowercases names, converts spaces/hyphens to underscores, and
removes assay/polarity tokens (`c18`, `hilic`, `lipid`, `lipids`, `lipidomics`,
`pos`, `neg`, `ms1`). Thus `demo_C18_crude_1_pos` and
`demo_HILIC_crude_1_neg` align to `demo_crude_1` across separate tables.

Within one table, names must remain unique after normalization. The pipeline
rejects collisions rather than summing injections or guessing which column to use.
A numeric suffix is part of sample identity: `crude_1` and `crude_2` are separate
groups, not inferred biological replicates. The low-level QC API accepts explicit
replicate groups, but the CLI does not yet read a replicate-design sheet.

The default real-sample whitelist is the token set `{crude, mcx}`. Adapt it with
repeated `--sample-token` options. Blank-like tokens (`blank`, `mb`,
`resuspension`) feed the baseline and never form sample groups. Other columns,
including pooled QC and reference injections, are ignored for QC unless explicitly
included by your sample vocabulary. Zero recognized groups raises an error.

The four metadata rows are available to `validate`; QC sample roles currently
come from names, not those rows. Check vocabulary and blank assignments before
interpreting results from a new study.

## Provenance

`source_file` is relative to the input directory. Together with `Alignment ID`,
it maps an output row back to its source. `feature_id` is sequential within a
run and can change when inputs or filters change; it is not a permanent molecular
identifier. Keep SIRIUS outputs with the run that produced their exported IDs.
