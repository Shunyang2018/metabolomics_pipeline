# Data and dependency provenance

- Source code: MIT, with authorship preserved in Git history and [LICENSE](../LICENSE).
- `examples/synthetic/*`: generated specifically for this repository; MIT. All
  signals, spectra, taxonomy and reference associations are illustrative.
- No instrument raw files, employer datasets, commercial reports or real
  bioactivity database are included in the public example.
- An earlier revision contained a real ClassyFire seed cache. The showcase
  revision removes it from the current source tree and package: downloaded
  records are now stored in a user cache, with no dependency on historical data.
  Git history still records earlier revisions; this is not a history rewrite.
- SIRIUS is an external installation, not included or relicensed here.
- ClassyFire responses and any user-supplied bioactivity database retain their
  respective source terms. Record provenance and permitted redistribution before
  adding such data to an example.
- Python dependencies and exact resolved versions are recorded in `pyproject.toml`
  and `uv.lock`; their licenses remain those of the individual projects.

This repository demonstrates downstream MS-DIAL processing. It should be described
on its own implementation and evidence, without transferring performance or
experimental-validation claims from a separate proprietary pipeline.
