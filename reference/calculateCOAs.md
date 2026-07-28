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
  lat.col = NULL,
  verbose = getOption("moby.verbose", TRUE)
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

- verbose:

  Logical; print a summary of the operation. Defaults to
  `getOption("moby.verbose", TRUE)`.

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
#> ── calculateCOAs() ───────────────────────────────────────────────────── moby ──
#> 
#> ℹ Estimating centres of activity per individual and time bin
#> • Input: 1,643 detections · 8 individuals
#> 
#> ✔ 794 positions estimated across 8 individuals
head(coas)
#>    ID             timebin    lon    lat detections stations station
#> 1 D02 2023-04-02 17:00:00 -9.008 38.464          3        1    ST02
#> 2 D02 2023-04-02 19:00:00 -9.008 38.464          3        1    ST02
#> 3 D02 2023-04-02 20:00:00 -9.008 38.464          4        1    ST02
#> 4 R04 2023-04-03 09:00:00 -8.985 38.466          1        1    ST04
#> 5 R04 2023-04-03 10:00:00 -8.985 38.466          1        1    ST04
#> 6 D04 2023-04-05 21:00:00 -8.996 38.456          1        1    ST03
#>              species  receiver    transmitter
#> 1 Dasyatis pastinaca VR2W-1002 A69-1602-30006
#> 2 Dasyatis pastinaca VR2W-1002 A69-1602-30006
#> 3 Dasyatis pastinaca VR2W-1002 A69-1602-30006
#> 4       Raja clavata VR2W-1004 A69-1602-30004
#> 5       Raja clavata VR2W-1004 A69-1602-30004
#> 6 Dasyatis pastinaca VR2W-1003 A69-1602-30008
```
