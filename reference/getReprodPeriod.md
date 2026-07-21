# Assign reproductive status

Assign spawning vs resting periods based on given start and end
interval.

## Usage

``` r
getReprodPeriod(
  datetimes,
  spawning.start,
  spawning.end,
  format = "%m",
  tz = "UTC"
)
```

## Arguments

- datetimes:

  A POSIXct object containing the respective datetimes or time-bins.

- spawning.start:

  Start of the spawning season. Can be supplied as a POSIXct object or
  as month (e.g., "08" for August).

- spawning.end:

  End of the spawning season. Can be supplied as a POSIXct object or as
  month (e.g., "09" for September).

- format:

  Character string giving the datetimes-time format of the supplied
  spawning.start and spawning.end. See
  [`strptime`](https://rdrr.io/r/base/strptime.html) details for further
  info on available formats. Defaults to "%m" (month).

- tz:

  A character string specifying the time zone to be used for the
  conversion. Defaults to "UTC".

## Value

A factor indicating the reproductive state (resting vs spawning).

## Examples

``` r
# Using integer month format (spawning period between June and September)
datetimes <- as.POSIXct("2024-05-30")
getReprodPeriod(datetimes, spawning.start="05", spawning.end="09", format="%m")
#> [1] spawning
#> Levels: resting spawning

# Using abbreviated month format (spawning period between June and September)
datetimes <- as.POSIXct("2024-05-30")
getReprodPeriod(datetimes, spawning.start="Jun", spawning.end="Sep", format="%b")
#> [1] resting
#> Levels: resting spawning

# Using day/month format (spawning period between 15th August and 31th August)
datetimes <- as.POSIXct("2024-09-01")
getReprodPeriod(datetimes, "15/08", "31/08", format="%d/%m")
#> [1] resting
#> Levels: resting spawning

# Using POSIXct objects
datetimes <- as.POSIXct("2024-05-30")
spawning.start <- as.POSIXct("2024-04-01")
spawning.end <- as.POSIXct("2024-09-30")
getReprodPeriod(datetimes, spawning.start, spawning.end, format="%Y-%m-%d")
#> [1] spawning
#> Levels: resting spawning
```
