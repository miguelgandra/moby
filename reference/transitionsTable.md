# Summarise a movement network as a transitions table

Produces a publication-ready summary table from a movement network (the
output of
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md)).
Each row is a directed transition between two locations, with the number
of movements, the number (and percentage) of distinct individuals
performing it, and the mean transit duration. When per-animal metadata
is supplied, numeric variables are summarised as mean +/- error and
categorical variables as level counts, per transition type.

This is the formatting counterpart to
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md)
(which holds the numeric network), mirroring the
[`calculateResidency`](https://miguelgandra.github.io/moby/reference/calculateResidency.md)
/
[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md)
split. It is purely a table: network visualisation is handled by
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the network
object, and temporal distributions of transition timing are available
from the network's `transition_records` attribute.

## Usage

``` r
transitionsTable(network, id.metadata = NULL, error.stat = "se")
```

## Arguments

- network:

  A `mobyNetwork` object of type `"movement"`, from
  [`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md).

- id.metadata:

  Optional data frame of per-animal metadata. Must contain an animal-ID
  column matching the network's `id.col`. Numeric columns are summarised
  as mean +/- error and categorical columns as counts, per transition
  type.

- error.stat:

  Error statistic for numeric metadata summaries: `"se"` (standard
  error, default) or `"sd"` (standard deviation).

## Value

A data frame with one row per directed transition (and group-label rows
when the network was built with `id.groups`).

## See also

[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md),
[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md)

## Examples

``` r
data(rays)
trans <- calculateTransitions(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> - Segmenting residence events with max.gap = 48 hours; tune/justify per system (Inf = split on location change only).
# publication-ready summary of the directed transitions
transitionsTable(trans)
#>                  Type Movements Individuals Mean duration (h)
#> 1        Raja clavata                                        
#> 2       ST01 --> ST02         1     1 (25%)              58.4
#> 3       ST01 --> ST03         5     3 (75%)       30.7 ± 12.6
#> 4       ST01 --> ST04         2     2 (50%)       96.2 ± 48.0
#> 5       ST01 --> ST05         1     1 (25%)              32.5
#> 6       ST01 --> ST06         1     1 (25%)              32.8
#> 7       ST02 --> ST01         2     1 (25%)       60.5 ± 11.9
#> 8       ST02 --> ST03         4     2 (50%)       73.0 ± 34.0
#> 9       ST02 --> ST04        10     3 (75%)       23.1 ± 11.1
#> 10      ST02 --> ST05         3     2 (50%)        12.0 ± 8.3
#> 11      ST02 --> ST06         3     2 (50%)       18.1 ± 11.6
#> 12      ST03 --> ST01         3     2 (50%)       94.9 ± 44.0
#> 13      ST03 --> ST02         9     3 (75%)        28.7 ± 9.6
#> 14      ST03 --> ST04         8    4 (100%)        30.2 ± 8.0
#> 15      ST03 --> ST05         2     1 (25%)       52.8 ± 34.7
#> 16      ST03 --> ST06         4     3 (75%)        34.5 ± 9.5
#> 17      ST04 --> ST01         2     2 (50%)         6.5 ± 2.7
#> 18      ST04 --> ST02         7     2 (50%)       47.4 ± 30.2
#> 19      ST04 --> ST03         7     3 (75%)      103.3 ± 18.1
#> 20      ST04 --> ST05         5     3 (75%)       48.5 ± 28.5
#> 21      ST04 --> ST06         7     2 (50%)       45.2 ± 19.5
#> 22      ST05 --> ST01         1     1 (25%)              51.9
#> 23      ST05 --> ST02         3     2 (50%)        36.3 ± 8.8
#> 24      ST05 --> ST03         5    4 (100%)       37.1 ± 27.9
#> 25      ST05 --> ST04         2     2 (50%)        22.7 ± 0.1
#> 26      ST05 --> ST06         1     1 (25%)              36.4
#> 27      ST06 --> ST01         1     1 (25%)               7.0
#> 28      ST06 --> ST02         3     3 (75%)       48.1 ± 17.8
#> 29      ST06 --> ST03         5    4 (100%)       69.0 ± 17.4
#> 30      ST06 --> ST04         6     1 (25%)        17.6 ± 7.6
#> 31      ST06 --> ST05         1     1 (25%)             150.2
#> 32 Dasyatis pastinaca                                        
#> 33      ST01 --> ST02         3     2 (50%)        55.3 ± 3.2
#> 34      ST01 --> ST03         1     1 (25%)              98.1
#> 35      ST01 --> ST04         1     1 (25%)             144.6
#> 36      ST01 --> ST05         7     3 (75%)       25.2 ± 12.7
#> 37      ST01 --> ST06         6     3 (75%)       42.4 ± 16.1
#> 38      ST02 --> ST01         3     2 (50%)       69.5 ± 35.5
#> 39      ST02 --> ST03         1     1 (25%)              16.6
#> 40      ST02 --> ST04         1     1 (25%)              46.4
#> 41      ST02 --> ST05         3     2 (50%)       26.5 ± 26.3
#> 42      ST02 --> ST06         4     2 (50%)      101.5 ± 20.2
#> 43      ST03 --> ST02         3     2 (50%)       39.1 ± 22.9
#> 44      ST03 --> ST05         7     3 (75%)       67.2 ± 21.7
#> 45      ST03 --> ST06         5     2 (50%)        18.7 ± 7.4
#> 46      ST04 --> ST01         1     1 (25%)               3.9
#> 47      ST04 --> ST02         1     1 (25%)              41.9
#> 48      ST04 --> ST03         1     1 (25%)              19.7
#> 49      ST04 --> ST05         2     1 (25%)       122.9 ± 3.5
#> 50      ST04 --> ST06         1     1 (25%)              28.0
#> 51      ST05 --> ST01         8     3 (75%)       81.0 ± 28.8
#> 52      ST05 --> ST02         2     1 (25%)         0.4 ± 0.3
#> 53      ST05 --> ST03         4     3 (75%)       87.3 ± 18.6
#> 54      ST05 --> ST04         3     3 (75%)       63.1 ± 15.8
#> 55      ST05 --> ST06         8     3 (75%)       36.4 ± 15.1
#> 56      ST06 --> ST01         6     3 (75%)       48.5 ± 21.8
#> 57      ST06 --> ST02         3     1 (25%)      137.4 ± 56.0
#> 58      ST06 --> ST03         6     3 (75%)       70.4 ± 32.2
#> 59      ST06 --> ST04         2     2 (50%)        52.3 ± 2.3
#> 60      ST06 --> ST05         7    4 (100%)       96.5 ± 34.1
```
