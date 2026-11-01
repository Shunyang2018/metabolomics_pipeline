# External integrations

The public demo is self-contained. These adapters are optional and must be
validated with the particular installed software and data used in a study.

## ClassyFire

`metabo classify` queries an InChIKey entity endpoint or uses an explicit JSON
cache. `--offline` makes no HTTP requests; missing keys remain unclassified.
`--cache-path` selects the run's cache. Without it, the cache is under
`LOCALAPPDATA/metabo_pipeline` on Windows or
`XDG_CACHE_HOME/metabo_pipeline` (fallback `~/.cache/metabo_pipeline`).
The package directory is never used as a writable default cache.

Real taxonomy records are not bundled. The example cache is explicitly synthetic.
The default external endpoint is `http://classyfire.wishartlab.com`; availability,
terms, data coverage and connectivity are outside this repository's control.
`metabo classify-check` probes the endpoint. Use local cache snapshots when
reproducibility matters and record their source and permitted use.

## SIRIUS

Install SIRIUS separately and follow its documentation for authentication and
access to the requested tasks. Set `SIRIUS_EXECUTABLE`, put it on PATH, or pass
`--sirius-exe`. A configured executable is authoritative; a typo does not silently
fall back to another installation.

The adapter launches the configured tasks from `configs/sirius.yml` and collects
TSV summaries. `sirius.cores`, `sirius.executable` and `tasks` are consumed by the
current driver. Historical `adducts`, `common_config` and `mzdev_ppm` settings
were not forwarded by this adapter; the sample configuration omits them so it
does not imply those settings are effective. Set instrument-specific options in
your SIRIUS workflow and verify them in its logs. The retained `--mzdev-ppm`
compatibility option currently warns rather than claiming to apply a setting.

`--dry-run` creates placeholder project directories and checks command plumbing;
it produces no predictions and is not a scientific validation. Use a separate
directory from real results. The legacy `--test-real-case` option expects local
`test_data/test_pos.ms` and `test_data/test_neg.ms`; that dataset is not distributed.

`metabo run` stops if an annotation task fails. Fresh runs require an empty output
directory. `--skip-merge` is an explicit resume of that run, not permission to
reuse results across changed inputs. The standalone `final`/`sirius-collect`
commands can reuse existing files: ensure they belong to the same feature IDs.

Unit tests exercise output parsing and mock failure paths. No claim is made that
the current external adapters pass live tests against every SIRIUS version or
that an upstream service was reachable during this review.
