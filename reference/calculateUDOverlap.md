# Pairwise overlap between utilization distributions

Computes pairwise overlap between the individual utilization
distributions (UDs) produced by
[`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md),
returning a tidy table of dyadic overlap values (one row per pair of
animals). It is the home-range counterpart to
[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md)
(which measures spatiotemporal co-occurrence): here the question is how
much two animals' *space use* overlaps, irrespective of timing.

The overlap computation is delegated to the estimator that produced the
UDs, so the result is method-consistent: for autocorrelated KDE (AKDE,
the
[`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
default) it uses [`overlap`](https://rdrr.io/pkg/ctmm/man/overlap.html),
which returns the Bhattacharyya coefficient **with confidence
intervals**; for classic kernel density (`method = "kde"`) it uses
[`kerneloverlaphr`](https://rdrr.io/pkg/adehabitatHR/man/kerneloverlap.html),
which offers a wider set of indices but no CIs.

## Usage

``` r
calculateUDOverlap(
  uds,
  index = "BA",
  contour = 95,
  conf.level = 0.95,
  id.groups = NULL,
  verbose = TRUE
)
```

## Arguments

- uds:

  The output of
  [`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
  (a list with a `$ud` element). The estimation method (AKDE vs KDE) is
  detected automatically.

- index:

  Overlap index to compute (case-insensitive). One of `"BA"` (default),
  and for KDE-estimated UDs additionally `"UDOI"`, `"HR"`, `"PHR"`,
  `"VI"`, `"HD"`.

- contour:

  KDE only: the home-range isopleth percentage the overlap is restricted
  to (a value in (0, 100\]). Defaults to 95. Ignored for AKDE (which
  integrates the full UD).

- conf.level:

  AKDE only: the confidence level for the reported overlap CIs. Defaults
  to 0.95. Ignored for KDE (which provides no CIs).

- id.groups:

  Optional named list of ID groups. When supplied, each pair is
  annotated with its `group1`/`group2` membership and a `pair_type` of
  `"within"` or `"between"`.

- verbose:

  Logical; print progress messages. Defaults to TRUE.

## Value

A tidy data frame with one row per pair of individuals (self-pairs
excluded) — unordered pairs for symmetric indices, ordered pairs for the
directional indices `HR`/`PHR`:

- id1, id2:

  the two animals in the pair.

- :

  the overlap value (column named after `index`, e.g. `BA` or `UDOI`).
  For `HR`/`PHR` this is the directional value for `(id1, id2)`.

- \_lower, \_upper:

  (AKDE only) the confidence-interval bounds.

- group1, group2, pair_type:

  (when `id.groups` is supplied or `uds` was grouped) the group
  membership of each animal and whether the pair is within or between
  groups.

When a single unit produces pairs, the estimator's own overlap matrix is
attached as attribute `"matrix"` (symmetric except for the directional
`HR`/`PHR` indices, which stay asymmetric); `"method"`, `"index"`, and
`"contour"`/`"conf.level"` record how it was computed.

## Details

Available `index` values depend on how the UDs were estimated:

- **AKDE** (`ctmm`): only `"BA"` (Bhattacharyya coefficient), reported
  with a lower/upper confidence interval at the requested `conf.level`.
  This autocorrelation-aware estimate is the recommended choice for
  tracking data.

- **KDE** (`adehabitatHR`): `"BA"` (Bhattacharyya), `"UDOI"`
  (utilization distribution overlap index), `"HR"` (home-range area
  overlap), `"PHR"` (probability of finding animal j in animal i's
  range), `"VI"` (volume of intersection) or `"HD"` (Hellinger
  distance). See Fieberg & Kochanny (2005) for definitions.

`"BA"` is the default because it is the one index available under *both*
methods, so the default works whether the UDs are AKDE or KDE. All
indices are bounded in `[0, 1]` except `UDOI`, which can exceed 1 when
two ranges overlap and are both non-uniform. Overlap is computed within
each unit (e.g. species / `id.groups` block) that `uds` was estimated
over, never across units (whose UDs live on separate grids); pairs are
formed among individuals of the same unit only.

**Symmetric vs directional indices.** Most indices (`"BA"`, `"UDOI"`,
`"VI"`, `"HD"`) are symmetric — the overlap of A with B equals that of B
with A — and are reported once per unordered pair. `"HR"` and `"PHR"`
are *directional*: `HR[i, j]` is the proportion of animal `i`'s range
covered by animal `j`, and `PHR[i, j]` the probability of finding animal
`j` inside animal `i`'s range, so `[i, j] != [j, i]`. For these two
indices the result therefore has one row per *ordered* pair (both
`id1->id2` and `id2->id1`), the value in each row being the index
evaluated for `(id1, id2)`.

## References

Fieberg, J. & Kochanny, C. O. (2005). Quantifying home-range overlap:
the importance of the utilization distribution. Journal of Wildlife
Management, 69(4), 1346-1359.

Winner, K., Noonan, M. J., Fleming, C. H., et al. (2018). Statistical
inference for home range overlap. Methods in Ecology and Evolution,
9(7), 1679-1691.

## See also

[`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md),
[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
[`plotMaps`](https://miguelgandra.github.io/moby/reference/plotMaps.md)

## Examples

``` r
# \donttest{
if (requireNamespace("adehabitatHR", quietly = TRUE)) {
  data(rays)
  # estimate UDs for a few animals (kde keeps the example fast), then measure pairwise overlap
  sub <- rays[rays$ID %in% head(levels(factor(rays$ID)), 3), ]
  uds <- calculateUDs(sub, method = "kde", bandwidth = 500, verbose = FALSE)
  calculateUDOverlap(uds, index = "UDOI")
}
#> Warning: - 'id.col' converted to factor.
#> Warning: - Some of the ID(s) in id.groups don't match the IDs in the data.
#> - Computed UDOI overlap for 3 pair(s) (kde).
#>   id1 id2      UDOI
#> 1 D01 D02 0.9402257
#> 2 D01 D03 0.8684309
#> 3 D02 D03 0.9359016
# }
```
