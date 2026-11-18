# Create movement stats table

Creates a publication-ready table of per-animal movement metrics (total
distance travelled, rate of movement, linearity index and home-range
areas), with a summary mean +/- SE row. This is a formatter: the
underlying numeric values are computed by
[`calculateROM`](https://miguelgandra.github.io/moby/reference/calculateROM.md)
(total distance and rate of movement) and
[`calculateLinearityIndex`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md)
(movement directness); use those functions directly when you need the
raw values rather than a formatted table.

## Usage

``` r
movementTable(
  data,
  uds,
  id.col = NULL,
  timebin.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  dist.col = "dist_m",
  id.groups = NULL,
  land.shape = NULL,
  epsg.code = NULL,
  discard.missing = TRUE,
  verbose = getOption("moby.verbose", TRUE),
  ...
)
```

## Arguments

- data:

  A data frame containing binned animal detections and distances
  traveled, as returned by
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).

- uds:

  Output of
  [`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md).

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- timebin.col:

  Name of the column containing time bins (in POSIXct format). Defaults
  to `"timebin"`.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- dist.col:

  Name of the column containing distance values (in meters). Defaults to
  'dist_m'.

- id.groups:

  Optional. A list containing ID groups, used to visually aggregate
  animals belonging to the same class (e.g. different species).

- land.shape:

  Optional. A projected shape file containing coastlines, used (when
  supplied) to compute net displacements along the shortest in-water
  path for the linearity index.

- epsg.code:

  Coordinate reference system used to project positions (class 'CRS').
  If not supplied, CRS is assumed to be the same as in land.shape.

- discard.missing:

  If true, only individuals with detections are included.

- verbose:

  Logical; print a summary of the operation. Defaults to
  `getOption("moby.verbose", TRUE)`.

- ...:

  Additional arguments passed to
  [`calculateLinearityIndex`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md)
  (and onwards to
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)),
  used to calculate distances between the first and last recorded
  detections for each individual (e.g., `grid.resolution`,
  `mov.directions` and `cores`).

## Value

A `mobyTable`: a TYPED data frame (counts stay integer, indices numeric,
dates POSIXct), one row per individual, so the result can be computed on
directly. Presentation - fixed precision, the display-only
`mean +/- error` row, group headings - is applied by
[`format`](https://miguelgandra.github.io/moby/reference/format.mobyTable.md)
and
[`print`](https://miguelgandra.github.io/moby/reference/print.mobyTable.md).
Export the rendered version with
`write.csv(format(x), file, row.names = FALSE)`.

Columns use stable snake_case names, which is what `format(decimals=)`
and `format(group.by=)` are keyed on; the publication headers live in
`format(style = "report")`.

- your `id.col`, and a `group` factor when `id.groups` names more than
  one group

- `distance_km`, `rom`, `rom_max`, `linearity_index`

- the home-range columns from
  [`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)

`rom`/`rom_max` are in m/h unless every group is fast enough to warrant
km/h, in which case they are scaled and the unit is recorded on the
table - so the column NAME never changes with the data, and
[`format()`](https://rdrr.io/r/base/format.html) names the unit in the
header.

## See also

[`calculateROM`](https://miguelgandra.github.io/moby/reference/calculateROM.md),
[`calculateLinearityIndex`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md)

## Examples

``` r
# \donttest{
data(rays)

# build per-time-bin tracks with stepwise distances
coas <- calculateCOAs(rays)
#> Warning: - 'id.col' converted to factor.
#> ── calculateCOAs() ───────────────────────────────────────────────────── moby ──
#> 
#> ℹ Estimating centres of activity per individual and time bin
#> • Input: 1,643 detections · 8 individuals
#> 
#> ✔ 794 positions estimated across 8 individuals
tracks <- calculateStepDistances(coas, verbose = FALSE)

if (requireNamespace("adehabitatHR", quietly = TRUE)) {
  # home-range areas (a coarse estimation grid keeps this example fast)
  grid <- terra::rast(terra::ext(-9.05, -8.95, 38.43, 38.48),
                      ncol = 60, nrow = 60, crs = "EPSG:4326")
  terra::values(grid) <- 0
  grid <- terra::project(grid, "EPSG:32629")
  kud <- calculateUDs(coas, method = "kde", bandwidth = 500,
                       spatial.grid = grid)

  # publication-ready movement metrics table (one row per animal + mean +/- SE)
  movementTable(tracks, uds = kud)
}
#> ── calculateUDs() ────────────────────────────────────────────────────── moby ──
#> 
#> ℹ Estimating utilization distributions
#> • Input: 794 positions · 8 individuals
#> 
#> → Method
#>   • estimator  kernel density (KDE)
#>   • bandwidth  500 m
#>   • contours   50% · 95%
#>   • grouping   id.groups
#> 
#> ✔ 8 utilization distributions estimated
#> ⏱ runtime: 1.5s
#> ── movementTable() ───────────────────────────────────────────────────── moby ──
#> 
#> ℹ Summarising distance, rate of movement and space use per individual
#> • Input: 794 positions · 8 individuals
#> 
#> → Method
#>   • net displacement  straight-line (great-circle)
#> 
#> ℹ Irregular time-bin widths detected · distances interpolated to a common interval
#> <mobyTable: movement> 8 rows (2 groups; one mean ± se row per group)
#> 
#> ─ Raja clavata
#>         ID distance_km        rom        rom_max linearity_index         N COAs
#>        R01        62.2       31.3         2019.6            0.02         130.00
#>        R02        42.1       21.5         2106.6            0.08          85.00
#>        R03        71.6       43.3         1799.8            0.00         106.00
#>        R04        51.6       25.8         1467.6            0.03         126.00
#>  mean ± se  56.8 ± 6.4 30.5 ± 4.7 1848.4 ± 142.4     0.03 ± 0.02 111.75 ± 10.35
#>  UD 50% (Km2) UD 95% (Km2)
#>          3.64        16.94
#>          3.69        15.91
#>          4.65        18.18
#>          3.61        16.28
#>             -            -
#> 
#> ─ Dasyatis pastinaca
#>         ID distance_km        rom        rom_max linearity_index        N COAs
#>        D01        60.2       30.9         2275.3            0.02        116.00
#>        D02        46.7       22.5          273.3            0.00         72.00
#>        D03        45.5       23.2          979.4            0.00         83.00
#>        D04        54.0       28.3         2097.9            0.04         76.00
#>  mean ± se  51.6 ± 3.4 26.2 ± 2.0 1406.5 ± 474.3     0.02 ± 0.01 86.75 ± 10.01
#>  UD 50% (Km2) UD 95% (Km2)
#>          4.53        17.65
#>          3.57        17.42
#>          3.54        17.19
#>          3.90        17.57
#>             -            -
# }
```
