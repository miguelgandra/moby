# Return season boundaries for background shading.

Returns a table containing season boundaries/limits, used to add
background polygons in
[`plotAbacus`](https://miguelgandra.github.io/moby/reference/plotAbacus.md).

## Usage

``` r
shadeSeasons(
  start.time,
  end.time,
  interval,
  color.pal = c("white", "grey96", "grey83", "grey90"),
  hemisphere = "Northern"
)
```

## Arguments

- start.time:

  A POSIXct object containing the earliest date used in the plot region.

- end.time:

  A POSIXct object containing the latest date used in the plot region.

- interval:

  Time-bins interval (in minutes).

- color.pal:

  Vector of 4 colors, one for each season (in the following order:
  winter, spring, summer and autumn).

- hemisphere:

  Earth hemisphere for which to calculate seasons.

## See also

[`getSeason`](https://miguelgandra.github.io/moby/reference/getSeason.md)

## Examples

``` r
# season boundaries over one year (daily resolution), used for background shading
shadeSeasons(as.POSIXct("2023-01-01", tz = "UTC"),
             as.POSIXct("2023-12-31", tz = "UTC"),
             interval = 1440)
#>   season               start                 end  color
#> 1 winter 2023-01-01 00:00:00 2023-02-28 12:00:00  white
#> 2 spring 2023-02-28 12:00:00 2023-05-31 12:00:00 grey96
#> 3 summer 2023-05-31 12:00:00 2023-08-31 12:00:00 grey83
#> 4 autumn 2023-08-31 12:00:00 2023-11-30 12:00:00 grey90
#> 5 winter 2023-11-30 12:00:00 2023-12-31 00:00:00  white
```
