# Changelog

## Unreleased

- Add an offline synthetic example, input hashes, an HTML report, and reproducible
  expected outputs checked by tests and CI.
- Explain input formats, sample roles, annotation-score categories, algorithmic
  limitations, and external-data provenance.
- Add a `validate` command and `python -m metabo_pipeline` entry point.
- Reject unrecognized sample vocabulary and within-file normalized-name collisions;
  support repeated `--sample-token` options and retain `source_file` provenance.
- Exclude a merge's own output directory during recursive input discovery and
  report empty/unusable input folders as errors.
- Replace stale SIRIUS exports even for empty polarities; fail on incomplete
  required polarity runs or TSV export failures.
- Prevent skipped stages from reusing previous annotations, propagate task failures,
  and require a fresh output directory unless explicitly resuming a merge.
- Handle zero-hit offline taxonomy joins and missing keys, keep taxonomy caches
  outside the package, and remove the bundled historical taxonomy seed data.
- Exclude controls/metadata from reference matching and correct the matched-key
  summary to count actual matches.
- Support tab-delimited validation and honor explicitly configured SIRIUS paths.
