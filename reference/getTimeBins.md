# Assign time bins

This function assigns time bins to a dataset based on a specified time
interval, enabling the grouping of detection data into fixed time
periods. Time binning is particularly valuable in telemetry studies,
where detection events are often aggregated over regular intervals
(e.g., hourly or every 30 minutes).

The function determines time bins by rounding datetime values to the
nearest specified interval using one of three methods: `"floor"` (round
down), `"ceiling"` (round up), or `"round"` (round to the nearest). It
essentially serves as a wrapper for the `round_date()`, `floor_date()`,
and `ceiling_date()` functions from the `lubridate` package (Grolemund &
Wickham, 2011).

## Usage

``` r
getTimeBins(datetimes, interval = "30 mins", rounding.method = "floor")
```

## Arguments

- datetimes:

  A POSIXct vector containing datetimes.

- interval:

  A character string specifying the time unit or multiple of a unit to
  round to. Valid base units are `"second"`, `"minute"`, `"hour"`,
  `"day"`, `"week"`, `"month"`, `"bimonth"`, `"quarter"`, `"season"`,
  `"halfyear"`, and `"year"`. Arbitrary unique English abbreviations as
  in the
  [`lubridate::period()`](https://lubridate.tidyverse.org/reference/period.html)
  constructor are also allowed. Defaults to `"30 mins"`.

- rounding.method:

  The method for assigning time bins. Options are `"floor"` (default),
  `"ceiling"`, or `"round"`.

## Value

A vector of datetime values in POSIXct format representing the assigned
time bins, rounded according to the specified `interval` and
`rounding.method`.

## References

Grolemund, G., & Wickham, H. (2011). Dates and times made easy with
lubridate. Journal of statistical software, 40, 1-25.

## See also

[`round_date`](https://lubridate.tidyverse.org/reference/round_date.html)

## Examples

``` r
# Sample dataset with datetime column
data <- data.frame(
  id = 1:6,
  datetime = as.POSIXct(c(
    "2024-11-15 08:23:45",
    "2024-11-15 08:45:00",
    "2024-11-15 09:05:30",
    "2024-11-15 09:20:00",
    "2024-11-15 09:59:59",
    "2024-11-15 10:10:00"
  ))
)

# Assign time bins using a 30-minute interval
data$timebin <- getTimeBins(datetimes = data$datetime, interval = "30 mins")
```
