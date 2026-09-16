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
transitionsTable(
  network,
  id.metadata = NULL,
  error.stat = "se",
  verbose = getOption("moby.verbose", TRUE)
)
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

- verbose:

  Logical; print a summary of the operation. Defaults to
  `getOption("moby.verbose", TRUE)`.

## Value

A `mobyTable`: a TYPED data frame (counts stay integer, indices numeric,
dates POSIXct), one row per directed transition, so the result can be
computed on directly. Presentation - fixed precision, the display-only
`mean +/- error` row, group headings - is applied by
[`format`](https://miguelgandra.github.io/moby/reference/format.mobyTable.md)
and
[`print`](https://miguelgandra.github.io/moby/reference/print.mobyTable.md).
Export the rendered version with
`write.csv(format(x), file, row.names = FALSE)`.

Columns use stable snake_case names, which is what `format(decimals=)`
and `format(group.by=)` are keyed on; the publication headers live in
`format(style = "report")`.

- `transition` (`"A --> B"`), and a `group` factor when the network was
  built with `id.groups`

- `n_movements`, `n_individuals`, `pct_individuals`

- `mean_duration`, `error_duration` (hours, or days when transits are
  long - the unit is recorded on the table and named in the header)

- `mean_<var>` for each numeric `id.metadata` column, and the raw column
  for each categorical one

The count and the share of individuals are separate columns, as are the
mean duration and its error: written as one string (`"3 (75%)"`,
`"30.7 +/- 12.6"`) neither could be computed on.

## See also

[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md),
[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md)

## Examples

``` r
data(rays)
trans <- calculateTransitions(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> ── calculateTransitions() ────────────────────────────────────────────── moby ──
#> 
#> ℹ Building a directed movement network between locations
#> • Input: 1,643 detections · 8 individuals · 6 nodes (station)
#> 
#> → Method
#>   • max.gap  48 hours (a longer absence starts a new visit; tune per system)
# publication-ready summary of the directed transitions
transitionsTable(trans)
#> ── transitionsTable() ────────────────────────────────────────────────── moby ──
#> 
#> ℹ Summarising directed transitions between sites
#> • Input: 6 sites · 58 transitions
#> 
#> → Method
#>   • error  standard error (se)
#> <mobyTable: transitions> 58 rows (2 groups; one mean ± se row per group)
#> 
#> ─ Raja clavata
#>     transition n_movements n_individuals pct_individuals mean_duration
#>  ST01 --> ST02           1             1              25          58.4
#>  ST01 --> ST03           5             3              75          30.7
#>  ST01 --> ST04           2             2              50          96.2
#>  ST01 --> ST05           1             1              25          32.5
#>  ST01 --> ST06           1             1              25          32.8
#>  ST02 --> ST01           2             1              25          60.5
#>  ST02 --> ST03           4             2              50          73.0
#>  ST02 --> ST04          10             3              75          23.1
#>  ST02 --> ST05           3             2              50          12.0
#>  ST02 --> ST06           3             2              50          18.1
#>  ST03 --> ST01           3             2              50          94.9
#>  ST03 --> ST02           9             3              75          28.7
#>  ST03 --> ST04           8             4             100          30.2
#>  ST03 --> ST05           2             1              25          52.8
#>  ST03 --> ST06           4             3              75          34.5
#>  ST04 --> ST01           2             2              50           6.5
#>  ST04 --> ST02           7             2              50          47.4
#>  ST04 --> ST03           7             3              75         103.3
#>  ST04 --> ST05           5             3              75          48.5
#>  ST04 --> ST06           7             2              50          45.2
#>  ST05 --> ST01           1             1              25          51.9
#>  ST05 --> ST02           3             2              50          36.3
#>  ST05 --> ST03           5             4             100          37.1
#>  ST05 --> ST04           2             2              50          22.7
#>  ST05 --> ST06           1             1              25          36.4
#>  ST06 --> ST01           1             1              25           7.0
#>  ST06 --> ST02           3             3              75          48.1
#>  ST06 --> ST03           5             4             100          69.0
#>  ST06 --> ST04           6             1              25          17.6
#>  ST06 --> ST05           1             1              25         150.2
#>      mean ± se       4 ± 0         2 ± 0               -    46.8 ± 5.8
#>  error_duration
#>               -
#>            12.6
#>            48.0
#>               -
#>               -
#>            11.9
#>            34.0
#>            11.1
#>             8.3
#>            11.6
#>            44.0
#>             9.6
#>             8.0
#>            34.7
#>             9.5
#>             2.7
#>            30.2
#>            18.1
#>            28.5
#>            19.5
#>               -
#>             8.8
#>            27.9
#>             0.1
#>               -
#>               -
#>            17.8
#>            17.4
#>             7.6
#>               -
#>      18.3 ± 2.7
#> 
#> ─ Dasyatis pastinaca
#>     transition n_movements n_individuals pct_individuals mean_duration
#>  ST01 --> ST02           3             2              50          55.3
#>  ST01 --> ST03           1             1              25          98.1
#>  ST01 --> ST04           1             1              25         144.6
#>  ST01 --> ST05           7             3              75          25.2
#>  ST01 --> ST06           6             3              75          42.4
#>  ST02 --> ST01           3             2              50          69.5
#>  ST02 --> ST03           1             1              25          16.6
#>  ST02 --> ST04           1             1              25          46.4
#>  ST02 --> ST05           3             2              50          26.5
#>  ST02 --> ST06           4             2              50         101.5
#>  ST03 --> ST02           3             2              50          39.1
#>  ST03 --> ST05           7             3              75          67.2
#>  ST03 --> ST06           5             2              50          18.7
#>  ST04 --> ST01           1             1              25           3.9
#>  ST04 --> ST02           1             1              25          41.9
#>  ST04 --> ST03           1             1              25          19.7
#>  ST04 --> ST05           2             1              25         122.9
#>  ST04 --> ST06           1             1              25          28.0
#>  ST05 --> ST01           8             3              75          81.0
#>  ST05 --> ST02           2             1              25           0.4
#>  ST05 --> ST03           4             3              75          87.3
#>  ST05 --> ST04           3             3              75          63.1
#>  ST05 --> ST06           8             3              75          36.4
#>  ST06 --> ST01           6             3              75          48.5
#>  ST06 --> ST02           3             1              25         137.4
#>  ST06 --> ST03           6             3              75          70.4
#>  ST06 --> ST04           2             2              50          52.3
#>  ST06 --> ST05           7             4             100          96.5
#>      mean ± se       4 ± 0         2 ± 0               -    58.6 ± 7.3
#>  error_duration
#>             3.2
#>               -
#>               -
#>            12.7
#>            16.1
#>            35.5
#>               -
#>               -
#>            26.3
#>            20.2
#>            22.9
#>            21.7
#>             7.4
#>               -
#>               -
#>               -
#>             3.5
#>               -
#>            28.8
#>             0.3
#>            18.6
#>            15.8
#>            15.1
#>            21.8
#>            56.0
#>            32.2
#>             2.3
#>            34.1
#>      19.7 ± 3.1
```
