# Retrieve mobyData metadata

Returns the metadata list stored on a `mobyData` object (column mapping,
EPSG code, tagging dates, ID groups and land layer), or `NULL` for a
plain data frame.

## Usage

``` r
mobyMeta(x)
```

## Arguments

- x:

  A `mobyData` object (or any object).

## Value

A named list of metadata, or `NULL`.

## See also

[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md)

## Examples

``` r
data(rays)
meta <- mobyMeta(rays)
names(meta)
#>  [1] "id.col"        "datetime.col"  "timebin.col"   "station.col"  
#>  [5] "lon.col"       "lat.col"       "epsg.code"     "tagging.dates"
#>  [9] "id.groups"     "land.shape"   
meta$epsg.code
#> [1] 32629
meta$id.groups
#> $`Raja clavata`
#> [1] "R01" "R02" "R03" "R04"
#> 
#> $`Dasyatis pastinaca`
#> [1] "D01" "D02" "D03" "D04"
#> 
```
