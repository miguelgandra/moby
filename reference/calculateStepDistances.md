# Estimate step distances and reconstruct movement tracks

This function calculates the distance travelled between consecutive
positions (step distances) for multiple individuals, reconstructing the
underlying movement tracks (shortest in-water paths). If no land
shapefile is provided, the function defaults to calculating linear
(great-circle) distances, which are If no land shapefile is provided,
the function defaults to calculating linear (great-circle) distances,
which are faster to compute. When a land shapefile is provided, the
function can utilize parallel computing to expedite least-cost path
estimation, depending on the number of CPU cores specified. However,
depending on the chosen grid resolution and the spatial extent of the
provided positions, this process may take a long time to run, even with
parallel computing enabled.

## Usage

``` r
calculateStepDistances(
  data,
  id.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  land.shape = NULL,
  epsg.code = NULL,
  grid.resolution = 100,
  mov.directions = 16,
  cores = 1,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame with animal positions, containing longitude and latitude
  values.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- land.shape:

  Optional. A shapefile containing coastlines or landmasses. It can be
  supplied as an 'sf' object or as an object of class
  'SpatialPolygonsDataFrame' or 'SpatialPolygons'. If the provided
  object is not of class 'sf', the function will attempt to convert it
  to an 'sf' object for compatibility with subsequent spatial
  operations. If not supplied, the function will calculate linear
  (great-circle) distances instead of shortest in-water distances.
  Linear distances are much faster to compute but may be less realistic
  or accurate in cases where positions are separated by landmasses or
  other spatial obstacles.

- epsg.code:

  The EPSG code (integer) representing the coordinate reference system
  (CRS) to be used for projecting the positions. If not specified, the
  function will attempt to use the CRS from the provided land.shape (if
  available).

- grid.resolution:

  The grid cell size (in meters) used to estimate shortest paths. A
  higher resolution leads to more precise path calculations but may
  increase computation time.

- mov.directions:

  Size of the movement neighbourhood used when building the least-cost
  graph: `4` (rook), `8` (rook + bishop), or `16` (adds knight moves,
  the default; smoother least-cost paths).

- cores:

  Number of CPU cores to use for the computations. Defaults to 1, which
  means no parallel computing (single core). If set to a value greater
  than 1, the function will use parallel computing to speed up
  calculations. This parameter is only relevant for estimating
  least-cost paths, specifically when a `land.shape` is provided. Run
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to check the number of available cores.

- verbose:

  Logical. Should the function output process information and display a
  progress bar? Defaults to TRUE.

## Value

The input data (a
[`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
when one was supplied, otherwise a data frame) with an added `dist_m`
column giving the distance from each position to the next consecutive
position of the same individual (in metres; `NA` for each individual's
last position). The movement trajectories (one spatial-line geometry per
individual) are attached as the `"trajectories"` attribute and retrieved
with
[`getTrajectories`](https://miguelgandra.github.io/moby/reference/getTrajectories.md).
The distance-enriched output pipes directly into
[`calculateROM`](https://miguelgandra.github.io/moby/reference/calculateROM.md)
and
[`calculateLinearityIndex`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md).

## Note

Least-cost in-water paths are computed on a **terra**-rasterised cost
surface routed with **igraph** (Dijkstra); great-circle distances use
**geosphere**. The graph is built once and reused across a call's
segments (and across calls when the cache is supplied).

## See also

[`getTrajectories`](https://miguelgandra.github.io/moby/reference/getTrajectories.md),
[`calculateROM`](https://miguelgandra.github.io/moby/reference/calculateROM.md),
[`calculateLinearityIndex`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md),
[`plotMaps`](https://miguelgandra.github.io/moby/reference/plotMaps.md)

## Examples

``` r
data(rays)

# great-circle step distances between consecutive positions
# (no land shape supplied: fast linear paths)
rays_dist <- calculateStepDistances(rays)
#> Warning: - 'id.col' converted to factor.
#>   |                                                                              |                                                                      |   0%Calculating linear paths between consecutive positions...
#>   |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> Total execution time: 0.07 secs 
head(rays_dist$dist_m)
#> [1]    0.000    0.000    0.000    0.000 1639.958    0.000

# \donttest{
# shortest in-water paths around a coastline (least-cost routing)
land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(-8.99, 38.44), c(-8.96, 38.44), c(-8.96, 38.47),
  c(-8.99, 38.47), c(-8.99, 38.44)))), crs = 4326))
rays_lc <- calculateStepDistances(rays[1:60, ], land.shape = land,
                                  grid.resolution = 200)
#> Warning: - 'id.col' converted to factor.
#> Warning: - 27 detection(s) at 3 location(s) fall on the supplied land shape. A
#> detection's position is the position of the receiver that logged it, so this
#> usually means the deployment metadata places a receiver on land. Audit the
#> receiver log with checkDeployments() and correct the coordinates at source; for
#> a genuine near-shore position that a coarse coastline overlaps,
#> correctPositions() can relocate points to the nearest marine cell.
#> Building least-cost graph (200m grid | 16 directions)
#>   |                                                                              |                                                                      |   0%Calculating least-cost paths between consecutive positions...
#>   |                                                                              |======================================================================| 100%
#> Total execution time: 0.31 secs 
# }
```
