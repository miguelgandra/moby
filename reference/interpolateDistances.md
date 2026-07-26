# Distance interpolation

Function to interpolate (average) distances across all timebins. If an
animal goes undetected for large periods of time, large distance spikes
may occur after running the
[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
function. Hence, to avoid large biases when averaging distances per hour
or day, this function 'dilutes' traveled distances by the total number
of time-bins passed between each two consecutive positions.

## Usage

``` r
interpolateDistances(
  data,
  id.col = NULL,
  timebin.col = NULL,
  dist.col = "dist_m",
  keep.intermediate = FALSE
)
```

## Arguments

- data:

  Output
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
  function: a data frame with animal IDS, positions and distances.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- timebin.col:

  Name of the column containing time bins (in POSIXct format). Defaults
  to `"timebin"`.

- dist.col:

  Name of the column containing distances. Defaults to 'dist_m', the
  output of
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).

- keep.intermediate:

  Boolean indicating if intermediate distances (assigned to time-bins
  without detections) should be kept or discarded. Defaults to false.

## Value

Original data frame plus missing time-bins (with diluted distances)

## Examples

``` r
data(rays)

# build per-time-bin tracks with stepwise distances
coas <- calculateCOAs(rays)
#> Warning: - 'id.col' converted to factor.
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: datetime.
tracks <- calculateStepDistances(coas, verbose = FALSE)

# dilute distances across time-bins with no detections
interp <- interpolateDistances(tracks)
#> Interpolating distances
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
head(interp)
#>               timebin  ID    lon    lat detections stations station
#> 1 2023-04-08 15:00:00 D01 -8.996 38.456          3        1    ST03
#> 2 2023-04-08 16:00:00 D01 -8.996 38.456          2        1    ST03
#> 3 2023-04-08 17:00:00 D01     NA     NA          0        0    <NA>
#> 4 2023-04-08 18:00:00 D01     NA     NA          0        0    <NA>
#> 5 2023-04-08 19:00:00 D01     NA     NA          0        0    <NA>
#> 6 2023-04-08 20:00:00 D01     NA     NA          0        0    <NA>
#>              species  receiver    transmitter   dist_m
#> 1 Dasyatis pastinaca VR2W-1003 A69-1602-30005  0.00000
#> 2 Dasyatis pastinaca VR2W-1003 A69-1602-30005 81.99792
#> 3               <NA>      <NA> A69-1602-30005       NA
#> 4               <NA>      <NA> A69-1602-30005       NA
#> 5               <NA>      <NA> A69-1602-30005       NA
#> 6               <NA>      <NA> A69-1602-30005       NA
```
