# Calculate animals' kernel utilization areas

This function estimates animal utilization distributions (UDs) and
home-range areas. By default it uses **autocorrelated kernel density
estimation (AKDE)** via the `ctmm` package (`method = "akde"`), which
fits a continuous-time movement model to account for the serial
autocorrelation inherent in tracking data and returns area estimates
with confidence intervals. A classic (IID) fixed-bandwidth kernel
density estimator (`method = "kde"`, via
[`kernelUD`](https://rdrr.io/pkg/adehabitatHR/man/kernelUD.html)) is
also available for speed or backward compatibility. Utilization areas
are returned for the specified isopleths (default 50% and 95%); for
display, isopleth polygons can be clipped to land, and UDs can be
estimated independently for different groups or time periods.

## Usage

``` r
calculateUDs(
  data,
  bandwidth = NULL,
  method = c("akde", "kde"),
  spatial.grid = NULL,
  subset = NULL,
  id.groups = NULL,
  land.shape = NULL,
  id.col = NULL,
  time.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  epsg.code = NULL,
  contour.percent = c(50, 95),
  model.selection = c("fit", "select"),
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- data:

  A data frame with animal positions (e.g. COAs), containing
  longitude/latitude (or projected x/y) columns and, for
  `method = "akde"`, a POSIXct time column.

- bandwidth:

  Numeric smoothing parameter (h) for `method = "kde"` only; ignored for
  AKDE (which estimates smoothing from the fitted movement model).
  Larger values produce smoother, broader distributions; smaller values
  give more localized, potentially fragmented estimates.

- method:

  Estimation method: `"akde"` (default; autocorrelated KDE via `ctmm`,
  with confidence intervals) or `"kde"` (classic fixed-bandwidth KDE via
  `adehabitatHR`).

- spatial.grid:

  Optional. A `Raster` or `SpatialPixels` object representing the grid
  over which the animal kernel utilization distributions (UDs) will be
  estimated (see the `grid` argument in
  [`kernelUD`](https://rdrr.io/pkg/adehabitatHR/man/kernelUD.html)). If
  set to `NULL`, the function will automatically generate an appropriate
  grid based on the spatial extent of the supplied animal's positions.

- subset:

  Optional. A variable used to subset the data, allowing UDs to be
  calculated independently for each level of this variable. This should
  be the name of a column in the provided dataset. If left `NULL`, UDs
  are calculated for the whole monitoring period.

- id.groups:

  Optional. A named list where each element represents a group (e.g.,
  species, sex, or age class), containing a vector of IDs for that
  group. If supplied, UDs will be calculated independently for each
  group. The names of the list correspond to the group labels.

- land.shape:

  Optional. A shapefile containing coastlines or landmasses, provided
  either as an 'sf' object or as a
  'SpatialPolygonsDataFrame'/'SpatialPolygons' object. If the input is
  not in 'sf' format, the function will automatically convert it to 'sf'
  to ensure compatibility with subsequent spatial operations. Used to
  clip and exclude any portions of the estimated areas that overlap with
  landmasses.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- time.col:

  Name of the POSIXct time column used by `method = "akde"` to model
  temporal autocorrelation. If `NULL`, the function uses the `mobyData`
  time-bin/date-time column and, failing that, a `"timebin"` or
  `"datetime"` column present in the data.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- epsg.code:

  Optional integer EPSG code of a **projected** (metre-based) coordinate
  reference system, used when projecting coordinates or computing
  distances/areas.

- contour.percent:

  Numeric vector. The percentages for which isopleths (contour areas)
  are calculated. Defaults to 50% and 95%, representing core and total
  areas of utilization.

- model.selection:

  For `method = "akde"`: `"fit"` (default) fits a single movement model
  from an automated guess (faster); `"select"` runs
  [`ctmm::ctmm.select`](https://rdrr.io/pkg/ctmm/man/ctmm.fit.html) to
  choose among candidate models (more thorough, slower).

- verbose:

  Logical. If TRUE, the function will print detailed processing
  information. Defaults to `getOption("moby.verbose", TRUE)`.

## Value

A list containing:

- ud:

  The estimated utilization distributions: an `adehabitatHR` "estUDm"
  object for `method = "kde"`, or a named list of `ctmm` UD objects for
  `method = "akde"`.

- summary_table:

  A data frame summarizing the (point-estimate) area for each contour
  (isopleth) per individual, in km2. Suitable for
  [`movementTable`](https://miguelgandra.github.io/moby/reference/movementTable.md).

- `K50`, `K95`, ...:

  One 'sf' object per requested `contour.percent`, named
  `paste0("K", contour.percent)`, holding the isopleth polygons. For
  `method = "akde"` each individual contributes three polygons tagged by
  a `ci` column (`"low"`, `"est"`, `"high"`) - the confidence envelope
  of the isopleth, which
  [`plotMaps`](https://miguelgandra.github.io/moby/reference/plotMaps.md)
  can draw. If a subset/group was provided, a column distinguishes the
  groups.

- area_estimates:

  (`method = "akde"` only) A tidy data frame of area estimates with
  lower/upper confidence limits (in km2) and the effective sample size
  (DOF) for each individual and contour.

The results list also contains multiple attributes to store relevant
metadata, such as function options and processing details. These
attributes might be useful for tracking parameters and ensuring
reproducibility of the analysis.

## Details

This function estimates kernel utilization distributions (UDs) based on
animal location data, allowing for quick analyses of space-use patterns.
For a comprehensive overview of other home-range estimation methods,
check out Kraft et al. (2023) (full reference below in the References
section). The function also includes options for handling landmasses and
for grouping data by subsets or groups for independent analysis.

**Land clipping**: Land clipping is applied post-hoc, after kernel
density estimation. If you need to account for physical barriers like
land during UD estimation, consider alternative methods (e.g. dynamic
Brownian Bridge Movement Models as provided in the `RSP` package; Niella
et al. 2020).

**Bandwidth (h)**: The smoothing factor, or bandwidth (h), is a critical
parameter in kernel utilization distribution (UD) analysis, representing
the standard deviation of the kernel. It defines the extent to which a
location can influence the home range estimation and the overall density
estimate. The choice of bandwidth significantly impacts the results of
the analysis. The bandwidth can either be fixed (using a single value
for all data points) or variable (adapting based on point density). The
`method = "kde"` pathway uses a fixed bandwidth; the default
`method = "akde"` instead estimates smoothing from a fitted
continuous-time movement model and is generally preferred for
autocorrelated tracking data (Fleming et al. 2015).

- A larger bandwidth increases the influence of more distant data points
  on the home range estimation, leading to a wider utilization
  distribution (UD) and a larger overall home range size. This increased
  smoothing can help mitigate sampling errors but may obscure finer
  details, retaining only the most prominent features of the spatial
  data.

- Conversely, a smaller bandwidth allows for greater detail at smaller
  spatial scales but might result in more fragmented home range
  estimates.

### Bandwidth Selection Methods

The optimal bandwidth is not universally defined and depends on the
specific context and dataset. Several methods are available for
selecting the bandwidth, including:

- **Reference Bandwidth (Href):** A constant based on data variance,
  assuming normally distributed data. It provides consistent smoothing
  but may overestimate home range size in sparse datasets.

- **LSCV (Least-Squares Cross-Validation):** Minimizes error between
  observed and predicted data, but may undersmooth in sparse datasets.

- **Ad-Hoc Choice:** Researchers manually select a bandwidth based on
  prior knowledge or exploratory data analysis.

- **Direct Plug-In:** Estimates bandwidth by solving equations based on
  the second derivative of the kernel density. Computationally intensive
  but can yield robust results, especially in complex datasets.

In some cases, bandwidth is adjusted based on the detection range in the
array, particularly when spatial estimates produce smaller, disconnected
isopleths. Adjusting the bandwidth can help form continuous areas that
better reflect biological relevance.

## References

Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., &
Abecasis, D. (2023). Residency and space use estimation methods based on
passive acoustic telemetry data. Movement Ecology, 11(1), 12.
https://doi.org/10.1186/s40462-023-00349-y

Niella, Y., Flávio, H., Smoothey, A. F., Aarestrup, K., Taylor, M. D.,
Peddemors, V. M., & Harcourt, R. (2020). Refined Shortest Paths (RSP):
Incorporation of topography in space use estimation from node-based
telemetry data. Methods in Ecology and Evolution, 11(12), 1733-1742.
https://doi.org/10.1111/2041-210X.13484

Worton, B. J. (1989). Kernel methods for estimating the utilization
distribution in home-range studies. Ecology, 70(1), 164-168.
https://doi.org/10.2307/1938423

Fleming, C. H., Fagan, W. F., Mueller, T., Olson, K. A., Leimgruber, P.,
& Calabrese, J. M. (2015). Rigorous home range estimation with movement
data: a new autocorrelated kernel density estimator. Ecology, 96(5),
1182-1188. https://doi.org/10.1890/14-2010.1

## See also

[`kernelUD`](https://rdrr.io/pkg/adehabitatHR/man/kernelUD.html),
[`akde`](https://rdrr.io/pkg/ctmm/man/akde.html)

## Examples

``` r
# \donttest{
data(rays)
if (requireNamespace("adehabitatHR", quietly = TRUE)) {
  # a coarse estimation grid keeps this example fast
  grid <- terra::rast(terra::ext(-9.05, -8.95, 38.43, 38.48),
                      ncol = 60, nrow = 60, crs = "EPSG:4326")
  terra::values(grid) <- 0
  grid <- terra::project(grid, "EPSG:32629")
  kud <- calculateUDs(rays, method = "kde", bandwidth = 500,
                       spatial.grid = grid)
  kud$summary_table
}
#> Warning: - 'id.col' converted to factor.
#> Grouping data by: id.groups
#> Estimating kernel utilization distributions [Dasyatis pastinaca]...
#> Calculating 50% contours...
#> Calculating 95% contours...
#> Estimating kernel utilization distributions [Raja clavata]...
#> Calculating 50% contours...
#> Calculating 95% contours...
#> Total execution time: 0.63 secs
#>                group  ID N COAs UD 50% (Km2) UD 95% (Km2)
#> 1 Dasyatis pastinaca D01    249         4.61        17.56
#> 2 Dasyatis pastinaca D02    154         3.42        16.55
#> 3 Dasyatis pastinaca D03    160         3.48        17.26
#> 4 Dasyatis pastinaca D04    169         3.77        17.29
#> 5       Raja clavata R01    283         3.67        16.87
#> 6       Raja clavata R02    160         3.60        15.68
#> 7       Raja clavata R03    207         4.59        18.11
#> 8       Raja clavata R04    261         3.63        16.44
# }
```
