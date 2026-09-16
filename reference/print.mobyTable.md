# Print a moby summary table

Renders the formatted table (see
[`format`](https://miguelgandra.github.io/moby/reference/format.mobyTable.md))
beneath a one-line banner, with grouped tables broken by a blank line
and a group heading. The object itself stays typed - this affects only
what is shown.

## Usage

``` r
# S3 method for class 'mobyTable'
print(x, ...)
```

## Arguments

- x:

  A `mobyTable`.

- ...:

  Passed to
  [`format`](https://miguelgandra.github.io/moby/reference/format.mobyTable.md),
  so the display can be tuned in place, e.g.
  `print(x, style = "report", decimals = c(distance_km = 2))`.

## Value

`x`, invisibly.
