# Demote a mobyData to a plain data frame

Drops the `mobyData` class **and** its metadata, returning an ordinary
data frame. Without this method
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) would
remove the class but keep the metadata attribute, leaving an object that
[`is_moby()`](https://miguelgandra.github.io/moby/reference/is_moby.md)
calls plain while moby's column resolution still reads its stored roles.

## Usage

``` r
# S3 method for class 'mobyData'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `mobyData` object.

- ...:

  Ignored.

## Value

A plain `data.frame`, carrying no moby metadata.

## See also

[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md),
[`is_moby`](https://miguelgandra.github.io/moby/reference/is_moby.md),
[`mobyMeta`](https://miguelgandra.github.io/moby/reference/mobyMeta.md)

## Examples

``` r
data(rays)
plain <- as.data.frame(rays)
is_moby(plain)         # FALSE
#> [1] FALSE
mobyMeta(plain)        # NULL - the metadata is gone, not merely hidden
#> NULL
```
