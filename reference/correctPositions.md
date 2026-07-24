# Relocate points on land to the nearest marine cell.

This function relocates positions on land to the nearest marine cell
using either a coastline shapefile or a raster containing land surfaces
or bathymetry values.

## Usage

``` r
correctPositions(
  data,
  lon.col = NULL,
  lat.col = NULL,
  spatial.layer,
  epsg.code = NULL,
  raster.type = "land",
  depth.threshold = 0,
  max.distance.km = 50,
  plot = FALSE,
  cores = 1,
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- data:

  A data frame with animal positions, containing longitude and latitude
  values.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- spatial.layer:

  A spatial object used to determine land and marine areas. This can be
  one of the following:

  - A `SpatRaster` (terra; a `RasterLayer` is also accepted) with binary
    values or bathymetric data.

  - An `sf` or `SpatialPolygons` object representing coastlines or
    landmasses.

- epsg.code:

  Optional integer EPSG code of a **projected** (metre-based) coordinate
  reference system, used when projecting coordinates or computing
  distances/areas.

- raster.type:

  A character string indicating the type of raster when `spatial.layer`
  is a `Raster`. This parameter can be either 'land', 'water', or
  'bathy'.

  - If 'land', the raster is assumed to contain NA values in water
    surfaces.

  - If 'water', the raster is assumed to contain NA values in land
    surfaces.

  - If 'bathy', the raster is assumed to range between 0 (land) to
    maximum depth (positive or negative values).

- depth.threshold:

  A numeric value. If set, shortest in-water paths will be calculated
  only for cells with depths \>= this threshold. Only takes effect when
  raster is of type 'bathy'. Defaults to 0 (water vs land).

- max.distance.km:

  A numeric value specifying the maximum distance (in kilometers) to
  consider when relocating points. This parameter limits the search
  radius for the nearest marine cell, ensuring that only cells within
  the specified distance are evaluated. Points that are further than
  this distance from the nearest water will have their coordinates set
  to NA.

- plot:

  A logical value indicating whether to generate a plot of the
  relocation process. If set to TRUE, the function will plot the
  `spatial.layer` along with the original and updated positions. If set
  to FALSE, no plot will be generated.

- cores:

  Number of CPU cores to use for the computations. Defaults to 1, which
  means no parallel computing (single core). If set to a value greater
  than 1, the function will use parallel computing to speed up
  calculations. Run
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to check the number of available cores.

- verbose:

  Logical; print progress and a relocation summary to the console.
  Defaults to `getOption("moby.verbose", TRUE)`.

## Value

A list with two elements:

- data:

  The original data frame with updated positions for the points that
  were relocated from land.

- summary:

  A summary of the input data and changes made to the positions.

## Examples

``` r
# \donttest{
data(rays)
# a small land polygon overlapping part of the study area
land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
  c(-8.99, 38.44), c(-8.96, 38.44), c(-8.96, 38.47),
  c(-8.99, 38.47), c(-8.99, 38.44)))), crs = 4326))
corrected <- correctPositions(rays[1:50, ], spatial.layer = land)
#> Relocating positions on land to the nearest marine cell
#> Points relocated: 23
#> Mean relocation distance: 133 m (0 m — 436 m)
#> Total execution time: 0.38 secs
#> Warning: Coordinates were initially in a geographic CRS. They have been projected for spatial processing and converted back to geographic coordinates.
attr(corrected, "points.relocated")
#> [1] 23
head(corrected$summary)
#>   index original.lon original.lat new.lon new.lat distance_m
#> 1     6       -8.990       38.442   -8.99  38.442        0.0
#> 2     7       -8.990       38.442   -8.99  38.442        0.0
#> 3     8       -8.990       38.442   -8.99  38.442        0.0
#> 4    23       -8.985       38.466   -8.99  38.466      436.2
#> 5    24       -8.985       38.466   -8.99  38.466      436.2
#> 6    25       -8.985       38.466   -8.99  38.466      436.2
# }
```
