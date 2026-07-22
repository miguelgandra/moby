# Retrieve diel phase boundary times (hours)

Returns a table containing average sunrise, sunset and twilight times
for a given period (study duration). '

## Usage

``` r
getSunTimes(coords, start.time, end.time, by = "%m", solar.depth = 18)
```

## Arguments

- coords:

  A SpatialPoints, matrix or numeric object containing longitude and
  latitude coordinates (in that order) at which to estimate sunrise and
  sunset times.

- start.time:

  A POSIXct object containing the earliest date of the monitoring
  period.

- end.time:

  A POSIXct object containing the latest date of the monitoring period.

- by:

  Date-time format (as defined by
  [`strptime`](https://rdrr.io/r/base/strptime.html)) containing the
  time frame used to average sunrise, sunset and twilight times.
  Defaults to month ("%m").

- solar.depth:

  Angle of the sun below the horizon (in degrees). Passed the solarDep
  argument in
  [`crepuscule`](https://rdrr.io/pkg/suntools/man/crepuscule.html)
  function. Defaults to 18 (astronomical twilight).

## Value

Data frame with diel phase' boundary times (in hours)

## Examples

``` r
# average diel-phase boundary times (hours) off SW Portugal, early May
getSunTimes(c(-9, 38.4),
            as.POSIXct("2023-05-01", tz = "UTC"),
            as.POSIXct("2023-05-07", tz = "UTC"))
#>   interval    dawns sunrises  sunsets    dusks
#> 1       05 3.916368 5.602511 19.50174 21.19206
```
