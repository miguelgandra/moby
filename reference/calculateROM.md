# Calculate rates of movement (and total distance travelled)

Summarises each animal's stepwise distances (as returned by
[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md))
into a tidy, fully numeric table — one row per individual — suitable for
plotting and downstream statistical analysis. This is the numeric core
used internally by
[`movementTable`](https://miguelgandra.github.io/moby/reference/movementTable.md)
(which formats these values, joins home-range areas and adds summary
rows for publication); use `calculateROM()` directly when you need the
raw values rather than a formatted table.

Two related quantities are derived from the per-step distance column:

- **Total distance** — the sum of consecutive in-water (or linear) step
  distances, i.e. the total path length travelled by each individual.

- **Rate of movement (ROM)** — the per-step distance rescaled to an
  hourly rate using the (modal) time-bin interval. Both the mean and
  maximum ROM are reported, in metres per hour.

When the detection series contains irregular time-bin intervals (more
than two distinct bin widths), distances are first reconciled to a
common interval with
[`interpolateDistances`](https://miguelgandra.github.io/moby/reference/interpolateDistances.md)
so that rates are comparable across individuals.

## Usage

``` r
calculateROM(data, id.col = NULL, timebin.col = NULL, dist.col = "dist_m")
```

## Arguments

- data:

  A data frame (or
  [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md))
  of binned detections with stepwise distances, as returned by
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- timebin.col:

  Name of the column containing time bins (in POSIXct format). Defaults
  to `"timebin"`.

- dist.col:

  Name of the column containing the (stepwise) distance values, in
  metres. Defaults to `"dist_m"`.

## Value

A data frame with one row per individual containing: the ID column,
`n_steps` (number of movement steps, i.e. non-missing distances),
`total_distance_m` (total path length in metres), and `mean_rom` /
`max_rom` (mean and maximum rate of movement in metres per hour). The
(modal) time-bin interval used to convert distances to rates, in
minutes, is stored in the `"interval"` attribute.

## See also

[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md),
[`calculateLinearityIndex`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md),
[`movementTable`](https://miguelgandra.github.io/moby/reference/movementTable.md),
[`interpolateDistances`](https://miguelgandra.github.io/moby/reference/interpolateDistances.md)

## Examples

``` r
data(rays)

# build per-time-bin tracks with stepwise distances
coas <- calculateCOAs(rays)
#> Warning: - 'id.col' converted to factor.
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: datetime.
tracks <- calculateStepDistances(coas, verbose = FALSE)

# summarise total distance travelled and rate of movement per individual
calculateROM(tracks)
#> Interpolating distances
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#>    ID n_steps total_distance_m mean_rom   max_rom
#> 1 D01    1947         60186.62 30.91249 2275.3137
#> 2 D02    2077         46676.77 22.47317  273.3264
#> 3 D03    1967         45540.16 23.15209  979.4320
#> 4 D04    1912         54023.53 28.25499 2097.8850
#> 5 R01    1987         62193.49 31.30020 2019.5640
#> 6 R02    1953         42052.62 21.53232 2106.6194
#> 7 R03    1651         71563.60 43.34561 1799.7578
#> 8 R04    2002         51580.68 25.76458 1467.6392
```
