# Calculate distance to nearest shore/land feature and related movement metrics.

This function calculates the distance to the nearest coastline or land
feature for each animal position. Instead of computing pairwise
distances between all points and coastline vertices, it uses a
rasterised land mask and a pre-computed distance-to-land raster,
providing major gains in speed and memory efficiency for large datasets.
If step lengths between consecutive positions are supplied, the function
also quantifies movement direction relative to land (inshore, offshore,
alongshore).

## Usage

``` r
calculateLandDists(
  data,
  land.shape,
  epsg.code = NULL,
  mov.threshold = 0.5,
  id.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  dist.col = NULL,
  raster.res = 100
)
```

## Arguments

- data:

  A data frame containing animal positions.

- land.shape:

  An sf object representing landmass.

- epsg.code:

  Integer. The projected EPSG code (must use meters). If not supplied,
  it is taken from the `mobyData` metadata (see
  [`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md));
  a projected (metric) CRS is required.

- mov.threshold:

  Numeric (0-1). Proportion of movement perpendicular to shore required
  for inshore/offshore classification.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- dist.col:

  Optional. Column name for step length (meters).

- raster.res:

  Numeric. Resolution of the distance raster in meters.

## Value

A data frame with added spatial and movement columns.

## Examples

``` r
# \donttest{
data(rays)
# a small land polygon spanning the study area
land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(-9.05, 38.49), c(-8.90, 38.49), c(-8.90, 38.43),
  c(-9.05, 38.49)))), crs = 4326))
rays_land <- calculateLandDists(rays, land.shape = land)
#> Calculating distances to nearest land (raster method)...
#> Done! Total execution time: 2.68 secs 
head(rays_land$land_dist)
#> [1] 1298.88 1298.88 1298.88 1298.88 1298.88 2383.32
# }
```
