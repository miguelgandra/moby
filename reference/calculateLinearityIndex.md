# Calculate the linearity index of individual movements

Computes a movement **linearity (directness) index** for each animal,
returning a tidy, fully numeric table (one row per individual). The
index is the ratio between the net displacement (the shortest in-water
path between an individual's first and last positions) and the total
path length actually travelled (the sum of consecutive step distances):
\$\$LI = net\\ displacement / total\\ distance.\$\$ Values approach 1
for highly directional (near straight-line) movements and approach 0 for
convoluted, back-and-forth or resident movements.

Net displacements are obtained with
[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md),
so when a `land.shape` is supplied they follow the shortest in-water
route (consistent with the total distance); otherwise great-circle
distances are used. This is the numeric core behind the `LI` column of
[`movementTable`](https://miguelgandra.github.io/moby/reference/movementTable.md).

## Usage

``` r
calculateLinearityIndex(
  data,
  land.shape = NULL,
  epsg.code = NULL,
  id.col = NULL,
  timebin.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  dist.col = "dist_m",
  ...
)
```

## Arguments

- data:

  A data frame (or
  [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md))
  of binned detections with stepwise distances, as returned by
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).

- land.shape:

  Optional. A projected coastline/landmass layer (`sf` or convertible).
  When supplied, net displacements follow the shortest in-water path;
  otherwise linear (great-circle) distances are used.

- epsg.code:

  Optional. Coordinate reference system used to project positions. If
  not supplied, the CRS is taken from `land.shape` (when available).

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- timebin.col:

  Name of the column containing time bins (in POSIXct format). Defaults
  to `"timebin"`.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- dist.col:

  Name of the column containing the (stepwise) distance values, in
  metres. Defaults to `"dist_m"`.

- ...:

  Additional arguments passed to
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
  (e.g. `grid.resolution`, `mov.directions`, `cores`), used when
  computing net displacements.

## Value

A data frame with one row per individual containing: the ID column,
`net_distance_m` (shortest distance between the first and last
positions, in metres), `total_distance_m` (total path length, in metres)
and `linearity_index` (their ratio, between 0 and 1).

## See also

[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md),
[`calculateROM`](https://miguelgandra.github.io/moby/reference/calculateROM.md),
[`movementTable`](https://miguelgandra.github.io/moby/reference/movementTable.md)

## Examples

``` r
data(rays)

# build per-time-bin tracks with stepwise distances
coas <- calculateCOAs(rays)
#> Warning: - 'id.col' converted to factor.
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: datetime.
tracks <- calculateStepDistances(coas, verbose = FALSE)

# movement directness per individual (net displacement / total path length)
calculateLinearityIndex(tracks)
#>    ID net_distance_m total_distance_m linearity_index
#> 1 D01       1467.639         60186.62      0.02438481
#> 2 D02          0.000         46676.77      0.00000000
#> 3 D03          0.000         45540.16      0.00000000
#> 4 D04       2141.456         54023.53      0.03963931
#> 5 R01       1373.174         62193.49      0.02207906
#> 6 R02       3332.627         42052.62      0.07924898
#> 7 R03          0.000         71563.60      0.00000000
#> 8 R04       1467.639         51580.68      0.02845327
```
