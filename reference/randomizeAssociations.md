# Test association-network co-occurrences against a null model

This function tests the null hypothesis of temporally independent space
use (i.e., each animal occurs independently of the other) using Monte
Carlo permutation tests. The algorithm permutes entries within each
column, keeping the total number of detections of each individual and
the relative occurrence frequencies across receivers unchanged.
Permutations can be constrained by given factors (e.g., diel phase) and
are performed only within the monitoring period of each animal (interval
between its release and the date of its last detection). Empirical
p-values are calculated by comparing the observed spatiotemporal
overlaps against their null distribution. The p-values can be one-tailed
or two-tailed based on the specified alternative hypothesis (a
continuity correction is applied to avoid p-values of exactly 0 or 1).

## Usage

``` r
randomizeAssociations(
  data,
  overlaps,
  constraint.by = NULL,
  iterations = 1000,
  alternative = c("two.sided"),
  conf.level = 0.95,
  p.adjust.method = "fdr",
  cores = 1,
  random.seed = NULL,
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- data:

  A data frame containing binned detections in the wide format (time bin
  x individual matrix, with values corresponding to the receiver/station
  with the highest number of detections), as returned by
  [`createWideTable`](https://miguelgandra.github.io/moby/reference/createWideTable.md).

- overlaps:

  Data frame containing paiwise overlaps, as returned by
  [`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md).

- constraint.by:

  Optional. Variable(s) to constrain permutations, e.g. to account for
  potential diel or seasonal trends in animal occurrences. If supplied,
  permutations across time bins will be restricted to the same
  combination of levels of these variables.

- iterations:

  Number of Monte Carlo iterations (simulated datasets). Defaults to
  1000.

- alternative:

  Character string specifying the alternative hypothesis. Must be one of
  "two.sided", "less", or "greater". "two.sided" tests for deviation in
  either direction, "less" tests if the observed value is significantly
  less than the null distribution, and "greater" tests if the observed
  value is significantly greater than the null distribution. Defaults to
  "two.sided".

- conf.level:

  Confidence level for the test. Defaults to 0.95 (95%).

- p.adjust.method:

  Method used to correct the per-dyad p-values for multiple comparisons,
  passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Because
  every dyad in the network is tested simultaneously, some correction is
  recommended to control false positives. Defaults to `"fdr"`
  (Benjamini-Hochberg false discovery rate); use `"none"` to retain raw
  p-values, or any other method accepted by `p.adjust` (e.g.
  `"bonferroni"`, `"holm"`). The dyad-level `Association` label is based
  on the adjusted p-value, while both the raw (`P-value`) and adjusted
  (`Adjusted p-value`) values are reported.

- cores:

  Number of CPU cores to use for the computations. Defaults to 1, which
  means no parallel computing (single core). If set to a value greater
  than 1, the function will use parallel computing to speed up
  calculations. Run
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to check the number of available cores.

- random.seed:

  Optional. Set the seed for a reproducible randomization. See
  [`set.seed`](https://rdrr.io/r/base/Random.html).

- verbose:

  Logical; print progress (a permutation progress bar). Defaults to
  `getOption("moby.verbose", TRUE)`.

## Value

A list containing:

- **summary**:

  A summary table with overall results.

- **pairwise_results**:

  A data frame similar to that returned by the
  [`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md)
  function, but with additional columns indicating the mean of the null
  distribution, empirical p-values, and association type for each animal
  dyad.

- **randomized_overlaps**:

  A numeric matrix of simulated pairwise overlaps, with one row per
  animal dyad (the row names give the dyad, e.g. `"A-B"`) and one column
  per iteration. If a subset variable was specified in the supplied
  overlaps, this is instead a named list of such matrices, one per level
  of the subset variable.

The summary table includes the following columns:

- Type: The type of comparison. This column is included only when
  `id.groups` are identified in the supplied overlaps.

- Subset: This column (e.g., a temporal subset) is included only when a
  subset variable was used to calculate the supplied overlaps.

- N dyads: The number of dyads (pairs of individuals).

- Mean interval (d): The mean monitoring period (in days).

- Mean overlap (%): The mean observed overlap percentage with standard
  deviation.

- Mean null distr (%): The mean overlap percentage from the null
  distribution with standard deviation.

- P-value: The estimated p-value(s), indicating the significance of the
  observed overlaps.

- Association: The association type (e.g., positive, negative,
  non-significant) based on the p-value.

- Pairs Non Sig: The number of pairs with non-significant.

- Pairs \> Random: The number of pairs with observed overlap greater
  than the null distribution.

- Pairs \< Random: The number of pairs with observed overlap less than
  the null distribution.

## References

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. Ecology, 81(9), 2606-2621.

Farine, D. R. (2017). A guide to null models for animal social network
analysis. Methods Ecol Evol 8: 1309-1320.

## See also

[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
[`plotAssociations`](https://miguelgandra.github.io/moby/reference/plotAssociations.md)

## Examples

``` r
# \donttest{
data(rays)
wide <- createWideTable(rays, value.col = "station")
#> Warning: - 'id.col' converted to factor.
#> Warning: 3 (ID, time-bin) combination(s) had multiple differing values; the first was kept. Aggregate upstream (e.g. calculateCOAs) to control this.
#> Tied (ID, time-bin) instances (first value kept):
#>                  timebin  ID                ties
#> 1898 2023-06-14 11:00:00 D03 ST01 (1) | ST06 (1)
#> 2163 2023-04-24 05:00:00 D04 ST01 (4) | ST05 (4)
#> 5302 2023-06-25 16:00:00 R04 ST03 (1) | ST05 (1)
assoc <- calculateAssociations(wide)
#> Calculating overlap - complete monitoring duration
#> Total execution time: 0.03 secs
# test the observed co-occurrences against a permutation null model
# (iterations kept low here for speed; use the default 1000 in practice)
rand <- randomizeAssociations(wide, assoc, iterations = 100, random.seed = 1)
#> Total execution time: 0.14 secs
rand$summary
#>   Type N dyads Mean interval (d) Mean overlap (%) Mean null distr (%) P-value
#> 1  All      28                77      0.47 ± 0.81         0.42 ± 0.09   0.515
#>       Association Pairs Non Sig Pairs > Random Pairs < Random
#> 1 non-significant            28              0              0
# }
```
