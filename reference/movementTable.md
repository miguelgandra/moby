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
  ud.results,
  land.shape = NULL,
  epsg.code = NULL,
  id.groups = NULL,
  id.col = NULL,
  timebin.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  dist.col = "dist_m",
  discard.missing = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame containing binned animal detections and distances
  traveled, as returned by
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).

- ud.results:

  Output of
  [`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md).

- land.shape:

  Optional. A projected shape file containing coastlines, used (when
  supplied) to compute net displacements along the shortest in-water
  path for the linearity index.

- epsg.code:

  Coordinate reference system used to project positions (class 'CRS').
  If not supplied, CRS is assumed to be the same as in land.shape.

- id.groups:

  Optional. A list containing ID groups, used to visually aggregate
  animals belonging to the same class (e.g. different species).

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

- discard.missing:

  If true, only individuals with detections are included.

- ...:

  Additional arguments passed to
  [`calculateLinearityIndex`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md)
  (and onwards to
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)),
  used to calculate distances between the first and last recorded
  detections for each individual (e.g., `grid.resolution`,
  `mov.directions` and `cores`).

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
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: datetime.
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
  movementTable(tracks, ud.results = kud)
}
#> Grouping data by: id.groups
#> Estimating kernel utilization distributions [Dasyatis pastinaca]...
#> Calculating 50% contours...
#> Calculating 95% contours...
#> Estimating kernel utilization distributions [Raja clavata]...
#> Calculating 50% contours...
#> Calculating 95% contours...
#> Total execution time: 0.59 secs
#> Interpolating distances
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#>                    ID Distance (km)  ROM (m/h)  Max ROM (m/h)          LI
#> 1        Raja clavata                                                    
#> 2                 R01          62.2       31.3         2019.6        0.02
#> 3                 R02          42.1       21.5         2106.6        0.08
#> 4                 R03          71.6       43.3         1799.8        0.00
#> 5                 R04          51.6       25.8         1467.6        0.03
#> 6                mean    56.9 ± 6.4 30.5 ± 4.7 1848.4 ± 142.4 0.03 ± 0.02
#> 7  Dasyatis pastinaca                                                    
#> 8                 D01          60.2       30.9         2275.3        0.02
#> 9                 D02          46.7       22.5          273.3        0.00
#> 10                D03          45.5       23.2          979.4        0.00
#> 11                D04          54.0       28.3         2097.9        0.04
#> 12               mean    51.6 ± 3.4 26.2 ± 2.0 1406.5 ± 474.3 0.01 ± 0.01
#>      N COAs UD 50% (Km2) UD 95% (Km2)
#> 1                                    
#> 2       130         3.64        16.94
#> 3        85         3.69        15.91
#> 4       106         4.65        18.18
#> 5       126         3.61        16.28
#> 6  112 ± 10  3.90 ± 0.25 16.83 ± 0.50
#> 7                                    
#> 8       116         4.53        17.65
#> 9        72         3.57        17.42
#> 10       83         3.54        17.19
#> 11       76         3.90        17.57
#> 12  87 ± 10  3.88 ± 0.23 17.46 ± 0.10
# }
```
