# Contributing

Open an issue with the command, expected behavior, actual behavior, package/Python
versions and a minimal synthetic input. Do not attach private or restricted data.
For changes, branch from `main` and keep scientific behavior changes explicit.

```bash
uv sync --locked --extra dev
uv run ruff check .
uv run ruff format --check .
uv run pytest --cov=metabo_pipeline --cov-report=term-missing
uv run python examples/run_demo.py --output-dir outputs/contribution-demo
```

Add a regression test for a bug or changed scientific rule. Prefer deterministic,
offline fixtures. External-service tests must be opt-in and must not require
credentials in CI. If dependency declarations change, refresh and commit `uv.lock`.

Explain which output changes and why. Update methods/input documentation when
thresholds, sample roles, annotation rules or formats change. Preserve source
provenance and avoid presenting model/reference matches as validated identities.
Keep implementation, review and test evidence visible in the pull request; use
the same standards for human-written and agent-assisted contributions.

Before contributing a dataset, state its origin, license, scientific purpose and
expected results. Synthetic examples must be labeled, including any invented
taxonomy or reference entries. Original authorship is preserved in Git history.
