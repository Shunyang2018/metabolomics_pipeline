# Methods and interpretation

This is a research workflow with explicit heuristics. Thresholds need assessment
against the instrument, assay and sample design of each study. Synthetic test
results are not a performance benchmark on experimental metabolomics data.

## Filtering and grouping

The merge keeps spectra with at least three positive-intensity MS/MS fragments.
When `S/N average` exists, its value must be at least 5. Missing/non-numeric S/N
values fail that filter. There is no S/N filter if the column is absent.

For a row, blank fold is **the maximum over all recognized real-sample columns
divided by the mean over all blank-like columns** in that input file. This shared
fold value is attached to every group; it is not a per-group numerator.
Missing or zero blank means produce NaN and fail the blank-fold threshold (7).

Presence is the percentage of positive-intensity columns in each explicit group
(default minimum 60%). CV uses sample standard deviation (`ddof=1`) among positive
intensities and requires at least two present replicates. CV filtering is disabled
by default. The CLI treats each unique normalized sample as its own group, so
presence is usually 0/100 and CV undefined; these values must not be presented as
biological replicate statistics. A row survives if any recognized group passes.

## Pipeline annotation categories

The historical `annotation_level` column is retained for compatibility.
**Its values are project-specific score categories, not evidence of
standard-confirmed identity or automatic assignment to a community confidence standard.**
In particular, category 1 does not establish a match to an authentic standard.

| Category before SIRIUS | Rule |
| --- | --- |
| 1 | Non-placeholder name, weighted dot >0.75, absolute difference from reverse dot <0.2, at least 3 matched peaks |
| 2 | Non-placeholder name, weighted dot >=0.5, reverse dot >0.5, at least 3 matched peaks, and not category 1 |
| 3 | Missing/placeholder name, weak/missing scores, or insufficient peak evidence |

Category 3 names become `Unknown`. When SIRIUS results are joined, only unknown
or missing categories are eligible for updates: 3 for a predicted structure,
4 for CANOPUS class evidence without a structure, and 5 otherwise (including a
formula-only result). Without SIRIUS results, categories pass through unchanged.
These software rules are neither calibrated probabilities nor measured false
discovery rates. Category 3 consequently has different pre/post-SIRIUS meanings;
inspect the stage and evidence rather than treating the integer as a universal label.

## Merge, deduplication and export

Tables are stacked as features and aligned by normalized intensity-column names.
This does not sum abundances across methods, calibrate response factors, correct
batch effects, or guarantee that rows from different assays are unique molecules.

Unknowns are sorted by RT and m/z within chromatographic-method/polarity buckets.
Sequential neighboring rows are clustered within 0.6 minutes and 20 ppm, and the
highest S/N (then weighted dot score) supplies a representative. This is a heuristic
scan, not exhaustive multidimensional clustering. Chains can span more than the
nominal window; isomers can be collapsed. Other nearby pairs may be missed when
an intervening row breaks the sorted chain. Retained known-name rows are not deduplicated.

All retained unknown representatives stay in the matrix. Only those with precursor
m/z 150–800 are exported for SIRIUS; this is a project setting, not a universal
SIRIUS accuracy or applicability guarantee. The current report counts retained
features, not unique compounds. Review discarded representatives before using
the table for quantitative inference.

## Reference matching

Matching uses the first InChIKey block as a connectivity key. This deliberately
does not preserve full stereochemical identity. A reference association is not
proof of the exact compound, concentration, causal activity, or efficacy in a new
biological context. One feature can produce several reference rows, and each
detected sample adds rows; output rows are not independent experimental hits.

For pipeline-generated matrices, reference matching uses the QC-recorded sample
columns and excludes blank/QC controls. For standalone tables without QC columns,
column roles are inferred by a metadata blacklist: inspect or subset those tables
to avoid treating unfamiliar numeric metadata as sample intensities.

## Validation still needed

Use a redistributable experimental dataset with documented reference standards,
sample design and ground truth to assess sensitivity, precision and reproducibility.
Such validation, exhaustive clustering, a replicate-design sheet, and live
version-pinned external-service integration tests remain future work. No new
scientific accuracy or speed claim is made by the public synthetic demonstration.
