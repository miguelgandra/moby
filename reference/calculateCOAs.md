# Calculate Centers of Activity (COAs)

This function calculates position estimates (Centers of Activity, COAs)
from presence data collected at multiple receivers. It uses weighted
means of longitude and latitude based on detection counts during
specified time bins. Additionally, it aggregates all remaining columns
dynamically. For numeric columns, the mean is calculated, while for
character or factor columns, unique values are concatenated and
separated by "\|".

## Usage

``` r
calculateCOAs(
  data,
  id.col = NULL,
  timebin.col = NULL,
  station.col = NULL,
  lon.col = NULL,
  lat.col = NULL
)
```

## Arguments

- data:

  A data frame containing animal detections and including a time bin
  column (as specified by the `timebin.col` argument). Time bins can be
  created using the
  [`getTimeBins`](https://miguelgandra.github.io/moby/reference/getTimeBins.md)
  function.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

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

## Value

A data frame of per-(ID, time bin) positions with columns:

- The unique identifier (`id.col`) for each individual.

- The time bin (`timebin.col`) for which the COA was calculated.

- Mean longitude and latitude (`lon.col` and `lat.col`) for each ID and
  time bin.

- The number of detections (`detections`) for each ID and time bin.

- The number of unique stations visited (`stations`) for each ID and
  time bin.

- For numeric columns: Mean values for each ID and time bin.

- For character or factor columns: Concatenated unique values (separated
  by "\|") for each ID and time bin.

When the input is a
[`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md),
the result is also a `mobyData`: the input's metadata (CRS/`epsg.code`,
tagging dates, `id.groups`, land layer) is carried forward - so
downstream spatial analyses such as
[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
and
[`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
inherit it automatically. (The raw datetime column does not survive
aggregation, so `datetime.col` is not retained.)

## Examples

``` r
data(rays)

# calculate Centers of Activity (one position per ID and time bin)
coas <- calculateCOAs(rays)
#> Warning: - 'id.col' converted to factor.
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: datetime.
head(coas)
#> <mobyData> 6 records x 10 columns
#>   individuals: 3  (id.col = 'ID')
#>   columns: datetime.col='datetime', timebin.col='timebin', station.col='station', lon.col='lon', lat.col='lat'
#>   metadata: epsg=32629, tagging.dates (8), id.groups (2)
```
