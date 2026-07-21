# Check whether an object is a mobyData dataset

Check whether an object is a mobyData dataset

## Usage

``` r
is_moby(x)
```

## Arguments

- x:

  An object.

## Value

A logical value.

## See also

[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md)

## Examples

``` r
data(rays)
is_moby(rays)                # TRUE
#> [1] TRUE
is_moby(as.data.frame(rays)) # FALSE (plain data frame)
#> [1] FALSE
```
