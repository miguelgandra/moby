# Test movement-network transitions against a null model

Tests whether the directed transitions in a movement network (from
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md))
occur more (or less) often than expected under random movement, using a
within-individual permutation null model. For each individual, the
time-ordered sequence of visited locations is randomly permuted and the
transition network recomputed; repeating this many times builds a null
distribution for each edge. Empirical, multiple-comparison-corrected
p-values then identify transitions (and the overall network) that are
stronger or weaker than expected.

The null permutes the order in which each individual visited its set of
locations, preserving that set (and hence each animal's site fidelity)
while randomising connectivity. It therefore tests the non-randomness of
*movement order / connectivity*, not of site use itself. Note that
movement-network null models are an active area of research; this is one
defensible, transparent choice rather than a universally agreed
standard.

## Usage

``` r
randomizeTransitions(
  network,
  iterations = 1000,
  alternative = c("two.sided", "greater", "less"),
  p.adjust.method = "fdr",
  conf.level = 0.95,
  random.seed = NULL,
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- network:

  A `mobyNetwork` object of type `"movement"`, from
  [`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md).

- iterations:

  Number of permutations. Defaults to 1000.

- alternative:

  Alternative hypothesis: `"two.sided"` (default), `"greater"` or
  `"less"`.

- p.adjust.method:

  Multiple-comparison correction across edges, passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Defaults to
  `"fdr"`; use `"none"` for raw p-values.

- conf.level:

  Confidence level used for the significance label. Defaults to 0.95.

- random.seed:

  Optional integer for reproducibility.

- verbose:

  Logical; print progress (a permutation progress bar). Defaults to
  `getOption("moby.verbose", TRUE)`.

## Value

A list with

- edges:

  A data frame of per-transition results: `group`, `from`, `to`,
  `n_movements` (observed), `mean_null`, `p_value`, `p_adjusted`, and
  `association` (`"more"`, `"less"` or `"non-significant"`, based on the
  adjusted p-value).

- network:

  A data frame with a network-level test per group: the observed number
  of realised transitions versus the null, with its empirical p-value.

## See also

[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md),
[`randomizeAssociations`](https://miguelgandra.github.io/moby/reference/randomizeAssociations.md)

## Examples

``` r
# \donttest{
data(rays)
trans <- calculateTransitions(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> - Segmenting residence events with max.gap = 48 hours; tune/justify per system (Inf = split on location change only).
# test each transition against a within-individual permutation null
# (iterations kept low here for speed; use the default 1000 in practice)
rand <- randomizeTransitions(trans, iterations = 100, random.seed = 1)
#> Running 100 permutations across 2 group(s)...
head(rand$edges)
#>                       group from   to n_movements mean_null   p_value
#> Raja clavata.1 Raja clavata ST01 ST02           1      1.39 1.0000000
#> Raja clavata.2 Raja clavata ST01 ST03           5      2.64 0.1188119
#> Raja clavata.3 Raja clavata ST01 ST04           2      2.44 1.0000000
#> Raja clavata.4 Raja clavata ST01 ST05           1      1.16 1.0000000
#> Raja clavata.5 Raja clavata ST01 ST06           1      1.25 1.0000000
#> Raja clavata.6 Raja clavata ST02 ST01           2      1.19 0.7326733
#>                p_adjusted     association
#> Raja clavata.1          1 non-significant
#> Raja clavata.2          1 non-significant
#> Raja clavata.3          1 non-significant
#> Raja clavata.4          1 non-significant
#> Raja clavata.5          1 non-significant
#> Raja clavata.6          1 non-significant
rand$network
#>                                 group n_transitions mean_null    p_value
#> Raja clavata             Raja clavata           114     89.98 0.01980198
#> Dasyatis pastinaca Dasyatis pastinaca           100     79.68 0.01980198
# }
```
