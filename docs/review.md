# Repository review and showcase preparation

Reviewed baseline: `be846f24d4fc16f64fc05977e8826a277ea3923f`.

The repository already separated scientific logic from CLI commands and had
285 tests. The main gaps were reproducible public inputs, precise scientific
interpretation, and failure/stale-result behavior. The baseline run on this
Windows/Python 3.12 environment passed 284 tests; one missing-executable test
instead reached an installed fallback SIRIUS and failed in the sandbox. An
additional direct CLI help check exposed a legacy-console encoding error.

## Findings addressed

| Priority | Finding and effect | Change |
| --- | --- | --- |
| High | Unrecognized sample vocabulary bypassed blank/presence QC | The merge rejects zero recognized real-sample groups; `--sample-token` adapts naming |
| High | Normalization could give two injections the same sample column | Reject raw and normalized within-file collisions |
| High | Old SIRIUS input/prediction/classification files could contaminate a rerun | Replace empty exports, require fresh run output directories, and avoid skipped-stage caches |
| High | Orchestration could print success after an external task failed | Propagate exceptions/failures and incomplete required-polarity/summary-export outcomes |
| High | Score-based category 1 could be read as standard-confirmed identification | Explicitly document project categories and scientific validation limits |
| Medium | Offline classification with no cache hits could fail its join | Preserve a typed empty taxonomy schema and missing-key handling |
| Medium | Reference output could include controls/metadata and overstated matched-key counts | Use QC-recorded sample membership, exclude controls, and count actual matched keys |
| Medium | Recursive input discovery could ingest its own outputs | Exclude the output subtree when nested inside input and exclude the exact output file |
| Medium | TXT validation disagreed with TXT merge support | Respect tab delimiters in both metadata inspection and row counts |
| Medium | Package data cache needed writable source paths and undocumented seed data | Move caching to a user/explicit path; remove the historical seed from the current tree |
| Medium | Explicit executable typos could select another local SIRIUS installation | Honor configured paths and check missing executables before launch |
| Medium | CLI output could crash under legacy Windows encodings | Use a plain ASCII command summary and stream-safe logging |
| Medium | Example/documentation files were absent from the source distribution | Add an explicit source manifest; build both distribution formats |

## Reproducible demonstration

The public fixture has eight synthetic input features across two methods and
polarities. It exercises three independent rejection reasons and retains five
features, exports two unknowns and produces four illustrative reference rows.
The script records input hashes, effective thresholds, package versions and
stage counts, plus source/script hashes. Regression tests run it without HTTP access and compare repeated
manifests. HTML output is escaped and contains no remote assets or scripts.

The README now leads with the problem, runnable example and expected results.
Input/method/integration documentation, citation metadata, contribution guidance,
data provenance and an unreleased changelog support external review.

## Remaining scientific and integration limitations

- No redistributable experimental validation set, measured identification
  precision/recall, or production-scale benchmark is provided. No such claim is made.
- Default CLI groups are individual named samples, not a replicate-design sheet.
  Shared maximum-sample/mean-blank fold metrics and disabled CV filtering are
  preserved and documented rather than silently changing the study's QC policy.
- Unknown-feature deduplication remains a sequential neighbor heuristic that can
  merge chains or miss nonadjacent neighbors. It needs independent evaluation
  before interpreting representatives as unique molecules.
- The standalone reference matcher uses a metadata blacklist when no QC sample
  membership is available. User-supplied tables need role review.
- Live SIRIUS/HTTP integration is not tested in the offline suite. The minimal
  configuration lists only settings the current driver actually consumes.
- `--skip-merge` and standalone result-collection commands require the user to
  keep results with their original input/feature IDs. This is not full automatic
  cache-invalidation or a workflow engine.

These are explicit research-development opportunities, not hidden evidence of
validated accuracy. A suitable next scientific milestone is a versioned public
experimental dataset with reference labels, followed by a documented comparison
of QC and representative-selection strategies.

## Verification of this revision

- Windows, Python 3.12.14: **306 tests passed**, including the offline demonstration
  and failure-path regressions. Overall measured coverage was approximately 70%;
  live external execution is not covered.
- Ruff lint and format checks passed; Git whitespace checks passed.
- Locked dependency resolution and both source/wheel builds succeeded.
- The built wheel was extracted and imported outside the checkout; module CLI,
  a merge under a cp1252 piped console, and the example from the source distribution
  succeeded using those wheel modules. Expected counts were five retained
  features and four illustrative reference rows.
- Local Markdown targets and the presence of examples, documentation and citation
  metadata in the source distribution were checked. The historical taxonomy seed
  was absent from the wheel.
- HTML escaping and deterministic output were tested. Interactive visual preview
  was not performed because the connected browser disallowed local `file:` URLs.
- Multi-platform GitHub Actions results are separate from these local results;
  the workflow is configured for all six OS/Python combinations.
