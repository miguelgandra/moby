# Introduction to moby

``` r

library(moby)
#> ======================================================
#>  moby v1.0.0
#> ------------------------------------------------------
#>   - filter & quality-control acoustic detections
#>   - estimate residency, home range & movement networks
#> ------------------------------------------------------
#>   dive in:          help(package = "moby")
#>   chart the depths: https://github.com/miguelgandra/moby
#>   cite your voyage: citation("moby")
#> ======================================================
```

## Overview

`moby` is an R package for the analysis and visualization of passive
acoustic telemetry data, with a focus on marine environments. It
provides a unified set of tools spanning data cleaning, temporal and
spatial classification, movement and home-range analysis, residency
metrics, and spatiotemporal overlap / social-network analysis.

Most functions accept user-defined column names (e.g. `id.col`,
`datetime.col`, `station.col`). When several functions are applied to
the same dataset, you can wrap the data once with
[`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md):
the resulting `mobyData` object carries the column mapping (and,
optionally, the coordinate reference system, tagging dates and ID
groups) as metadata, so these do not need to be repeated on every call.
A `mobyData` object is still an ordinary `data.frame`, and plain data
frames remain fully supported.

``` r

df <- data.frame(
  ID = c("A", "A", "B"),
  datetime = as.POSIXct(c("2023-06-01 00:00", "2023-06-01 01:00", "2023-06-01 00:00"),
                        tz = "UTC"),
  station = c("R1", "R2", "R1")
)
md <- as_moby(df)
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: timebin, lon, lat.
md
#> <mobyData> 3 records x 3 columns
#>   individuals: 2  (id.col = 'ID')
#>   period: 2023-06-01 to 2023-06-01 01:00:00 (tz = UTC)
#>   columns: datetime.col='datetime', timebin.col='timebin', station.col='station', lon.col='lon', lat.col='lat'
```

## A small example dataset

The functions below operate on ordinary data frames, making it easy to
integrate `moby` into existing workflows. Here we build a tiny synthetic
detection dataset.

``` r

set.seed(1)
detections <- data.frame(
  ID = factor(rep(c("A", "B"), each = 24)),
  datetime = rep(seq(as.POSIXct("2023-06-01 00:00:00", tz = "UTC"),
                     by = "1 hour", length.out = 24), 2),
  station = sample(c("R1", "R2", "R3"), 48, replace = TRUE)
)
head(detections)
#>   ID            datetime station
#> 1  A 2023-06-01 00:00:00      R1
#> 2  A 2023-06-01 01:00:00      R3
#> 3  A 2023-06-01 02:00:00      R1
#> 4  A 2023-06-01 03:00:00      R2
#> 5  A 2023-06-01 04:00:00      R1
#> 6  A 2023-06-01 05:00:00      R3
```

## Temporal classification

[`getTimeBins()`](https://miguelgandra.github.io/moby/reference/getTimeBins.md)
rounds timestamps to a regular interval, while
[`getSeason()`](https://miguelgandra.github.io/moby/reference/getSeason.md)
and
[`getDielPhase()`](https://miguelgandra.github.io/moby/reference/getDielPhase.md)
assign seasonal and diel classes.

``` r

# assign 2-hour time bins
detections$timebin <- getTimeBins(detections$datetime, interval = "2 hours")

# meteorological season (northern hemisphere)
detections$season <- getSeason(detections$datetime, hemisphere = "Northern")

table(detections$season)
#> 
#> spring summer autumn winter 
#>      0     48      0      0
```

## Reshaping to wide format

[`createWideTable()`](https://miguelgandra.github.io/moby/reference/createWideTable.md)
converts the long-format detections into a time-bin x individual matrix,
the structure used by the overlap and visualization functions.

``` r

wide <- createWideTable(detections, value.col = "detections", verbose = FALSE)
#> Warning: - No 'detections' column found, assuming one detection per row
dim(wide)
#> [1] 12  3
```

## Where to go next

See
[`help(package = "moby")`](https://miguelgandra.github.io/moby/reference)
for the complete list of functions, including
[`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md),
[`calculateStepDistances()`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md),
[`calculateUDs()`](https://miguelgandra.github.io/moby/reference/calculateUDs.md),
[`calculateAssociations()`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
and the family of `plot*()` visualization tools.
