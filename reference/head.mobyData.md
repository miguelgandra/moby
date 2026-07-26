# Return the first or last rows of a `mobyData`

[`head()`](https://rdrr.io/r/utils/head.html) and
[`tail()`](https://rdrr.io/r/utils/head.html) behave as they do for any
data frame: they return the first (or last) `n` rows of the detection
table for inspection. They return a plain `data.frame`, so the rows
print as rows - the study metadata summary is what
[`print()`](https://rdrr.io/r/base/print.html) on the `mobyData` itself
is for.

To subset a dataset and *keep* it a `mobyData` (metadata and all), use
`[` instead: `x[1:100, ]`.

## Usage

``` r
# S3 method for class 'mobyData'
head(x, n = 6L, ...)

# S3 method for class 'mobyData'
tail(x, n = 6L, ...)
```

## Arguments

- x:

  A
  [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
  object.

- n:

  Number of rows. Defaults to 6, as for
  [`utils::head`](https://rdrr.io/r/utils/head.html).

- ...:

  Further arguments passed to
  [`utils::head`](https://rdrr.io/r/utils/head.html) /
  [`utils::tail`](https://rdrr.io/r/utils/head.html).

## Value

A `data.frame` of the first (or last) `n` rows.

## See also

[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md),
[`mobyMeta`](https://miguelgandra.github.io/moby/reference/mobyMeta.md)

## Examples

``` r
data(rays)
head(rays, 3)      # rows, like any data frame
#>    ID            datetime station    lon    lat            species  receiver
#> 1 D01 2023-04-08 15:05:10    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#> 2 D01 2023-04-08 15:14:59    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#> 3 D01 2023-04-08 15:52:10    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#>      transmitter             timebin
#> 1 A69-1602-30005 2023-04-08 15:00:00
#> 2 A69-1602-30005 2023-04-08 15:00:00
#> 3 A69-1602-30005 2023-04-08 15:00:00
dim(rays[1:100, ]) # '[' keeps the mobyData class and its metadata
#> [1] 100   9
```
