# Create a moby telemetry dataset

Wraps a detection data frame into a `mobyData` object: a lightweight
`data.frame` subclass that carries dataset-level metadata (column
mapping, coordinate reference system, tagging dates, ID groups and an
optional land layer) as attributes. Downstream `moby` functions read
this metadata automatically, so column names, EPSG codes and tagging
dates do not have to be re-specified on every call.

Because a `mobyData` object **is** a `data.frame`, it can be used
anywhere a data frame is expected (including with base subsetting,
`dplyr`, etc.), and plain data frames remain fully supported by all
`moby` functions. This replaces the former global
`setDefaults()`/`getDefaults()` mechanism with reproducible,
object-carried metadata.

Calling `as_moby()` on an existing `mobyData` object updates only the
supplied fields and inherits the rest, making it easy to add metadata
(e.g. `tagging.dates`) after construction.

## Usage

``` r
as_moby(
  data,
  id.col = .mobyDefaults[["id.col"]],
  datetime.col = .mobyDefaults[["datetime.col"]],
  timebin.col = .mobyDefaults[["timebin.col"]],
  station.col = .mobyDefaults[["station.col"]],
  lon.col = .mobyDefaults[["lon.col"]],
  lat.col = .mobyDefaults[["lat.col"]],
  epsg.code = NULL,
  tagging.dates = NULL,
  nominal.delay = NULL,
  id.groups = NULL,
  land.shape = NULL
)
```

## Arguments

- data:

  A data frame of detections (one row per detection, or per binned
  record).

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- datetime.col:

  Name of the column containing date-times in POSIXct format. Defaults
  to `"datetime"`.

- timebin.col:

  Name of the column containing time bins (in POSIXct format). Defaults
  to `"timebin"`.

- station.col:

  Name of the column containing station/receiver IDs. Defaults to
  `"station"`.

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

- tagging.dates:

  Optional POSIXct vector of tagging/release dates. Either a single
  value (applied to all individuals) or a named vector whose names match
  the animal IDs.

- nominal.delay:

  Optional transmitter nominal (mean) delay, in seconds. Either a single
  value (applied to all individuals) or a named numeric vector whose
  names match the animal IDs (for arrays mixing tag families, e.g. 60 s
  and 120 s tags). Stored in the metadata and read automatically by
  [`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
  to scale its short-interval (min_lag) false-detection filter. Usually
  populated for you by
  [`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
  when the tag table carries a delay column.

- id.groups:

  Optional named list grouping IDs (e.g. by species, sex or life stage),
  used by many functions to compute metrics or draw plots independently
  per group.

- land.shape:

  Optional `sf` (or `SpatialPolygons*`) object representing landmasses,
  used by spatial functions (e.g. in-water distances, UD land clipping).

## Value

A `mobyData` object (a `data.frame` with a `"moby"` metadata attribute).

## See also

[`mobyMeta`](https://miguelgandra.github.io/moby/reference/mobyMeta.md),
[`is_moby`](https://miguelgandra.github.io/moby/reference/is_moby.md)

## Examples

``` r
df <- data.frame(
  ID = c("A", "A", "B"),
  datetime = as.POSIXct(c("2023-01-01 00:00", "2023-01-01 01:00", "2023-01-01 00:00"),
                        tz = "UTC"),
  lon = c(-8.1, -8.2, -8.0),
  lat = c(37.0, 37.1, 37.0),
  station = c("R1", "R2", "R1")
)
md <- as_moby(df, tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: timebin.
md
#> <mobyData> 3 records x 5 columns
#>   individuals: 2  (id.col = 'ID')
#>   period: 2023-01-01 to 2023-01-01 01:00:00 (tz = UTC)
#>   columns: datetime.col='datetime', timebin.col='timebin', station.col='station', lon.col='lon', lat.col='lat'
#>   metadata: tagging.dates (1)

# add or update metadata later
md <- as_moby(md, id.groups = list(grp1 = c("A", "B")))
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: timebin.
```
