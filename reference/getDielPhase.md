# Estimate diel phase

Function to retrieve diel phase (e.g. day/night) for given datetimes and
coordinates.  

The number of retrieved levels can be set by the 'phases' argument:  

- phases=2: day \| night  

- phases=3: day \| crepuscule \| night  

- phases=4: dawn \| day \| dusk \| night  

Crepuscular periods can be defined based on different solar elevation
angles:  

- solar.depth=6: civil twilight  

- solar.depth=12: nautical twilight  

- solar.depth=18: astronomical twilight  

## Usage

``` r
getDielPhase(datetimes, coords, phases = 2, solar.depth = 18)
```

## Arguments

- datetimes:

  A POSIXct object containing the respective datetimes or time-bins.

- coords:

  A SpatialPoints, matrix, or data frame object containing geographic
  (unprojected) longitude and latitude coordinates (in that order) for
  which to estimate sunrise and sunset times. If a single point or a
  matrix/data frame with one row is provided, the same coordinates will
  be used for all calculations.

- phases:

  Integer indicating the number of diel phases to return (2, 3, or 4).

- solar.depth:

  Numeric value indicating the angle of the sun below the horizon (in
  degrees). Passed to the
  [`crepuscule`](https://rdrr.io/pkg/suntools/man/crepuscule.html)
  function.

## Value

A factor indicating the diel phase.

## Examples

``` r
datetimes <- as.POSIXct("2024-05-30 12:00:00", tz = "UTC")
coords <- matrix(c(-7.997, 37.008), ncol = 2)
getDielPhase(datetimes, coords, phases = 4, solar.depth = 12)
#> [1] day
#> Levels: dawn day dusk night
```
