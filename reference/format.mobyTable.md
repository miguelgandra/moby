# Format a moby summary table for display or export

Renders a `mobyTable` as a character data frame: fixed per-metric
precision, missing values as `"-"`, and (where a group has more than one
row) a display-only `mean +/- error` row. This is the export route -
`write.csv(format(x), "table.csv", row.names = FALSE)` - and the same
rendering [`print()`](https://rdrr.io/r/base/print.html) shows.

The returned table stays rectangular. Group headings and the blank lines
between groups belong to the console rendering, not to the exported
object, so a CSV never carries an empty record.

## Usage

``` r
# S3 method for class 'mobyTable'
format(
  x,
  style = c("internal", "report", "concise"),
  symbols = c("ascii", "unicode"),
  decimals = NULL,
  group.by = NULL,
  include.summary.row = TRUE,
  datetime.format = "%d/%m/%Y",
  ...
)
```

## Arguments

- x:

  A `mobyTable` (from
  [`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md),
  [`movementTable`](https://miguelgandra.github.io/moby/reference/movementTable.md)
  or
  [`transitionsTable`](https://miguelgandra.github.io/moby/reference/transitionsTable.md)).

- style:

  Column-name style. `"internal"` (default) keeps the snake_case names
  the API guarantees; `"report"` uses publication-ready headers
  (`n_detections` -\> "N Detect"); `"concise"` uses the same headers
  abbreviated for narrow tables. Only the names differ.

- symbols:

  Whether the rendered table may use typographic symbols: `"ascii"`
  (default) writes `+/-`, `"unicode"` the plus-minus sign. ASCII is the
  default because this table is usually written to a file, and a
  spreadsheet opening a UTF-8 CSV with no byte-order mark guesses the
  encoding. [`print()`](https://rdrr.io/r/base/print.html) picks the
  right one for the terminal on its own.

- decimals:

  Optional per-column override of the display precision, as a named
  numeric vector of decimal places - `c(distance_km = 2)`. Merged OVER
  the built-in precision, so naming one column leaves the rest
  untouched, and both the values and the `mean +/- error` row follow it.
  Named by the INTERNAL column names, which do not change with `style`.

- group.by:

  Column to group the rendering by. Defaults to the table's own `group`
  column when it has one (the `id.groups` split), or to no grouping.
  Each group gets its own `mean +/- error` row. Pass `FALSE` to render
  ungrouped.

- include.summary.row:

  Logical; append the display-only `mean +/- error` row(s). Default
  `TRUE`. A group of ONE row never gets one: the "mean" would just
  restate that row.

- datetime.format:

  `strftime` format for date/time columns. Default `"%d/%m/%Y"`.

- ...:

  Unused.

## Value

A character `data.frame`.

## See also

[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md),
[`movementTable`](https://miguelgandra.github.io/moby/reference/movementTable.md),
[`transitionsTable`](https://miguelgandra.github.io/moby/reference/transitionsTable.md)

## Examples

``` r
data(rays)
tbl <- summaryTable(rays, last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"))
#> Warning: - 'id.col' converted to factor.
#> ── summaryTable() ────────────────────────────────────────────────────── moby ──
#> 
#> ℹ Summarising monitoring and residency metrics per individual
#> • Input: 1,643 detections · 8 individuals
#> 
#> → Method
#>   • residency index  IR1 · IR2 · IR2/IR1
#>   • start point      release date
#>   • error            standard deviation (sd)
#> Warning: - No 'detections' column found, assuming one detection per row.
# the object itself is typed - you can compute on it
mean(tbl$n_detections)
#> [1] 205.375
# ...and format() renders the version you export
head(format(tbl, style = "report"))
#>            ID Tagging date Last detection   N Detect N Receiv
#> 1         R01   03/04/2023     29/06/2023        283        6
#> 2         R02   08/04/2023     28/06/2023        160        5
#> 3         R03   10/04/2023     29/06/2023        207        6
#> 4         R04   01/04/2023     25/06/2023        261        6
#> 5 mean +/- sd            -              - 228 +/- 55  6 +/- 0
#> 6         D01   07/04/2023     28/06/2023        249        6
#>   Monitoring duration (d) Detection span (d) N days detected           IR1
#> 1                     272                 88              38          0.43
#> 2                     267                 82              23          0.28
#> 3                     265                 81              30          0.37
#> 4                     274                 86              36          0.42
#> 5               270 +/- 4           84 +/- 3        32 +/- 7 0.38 +/- 0.07
#> 6                     268                 83              31          0.37
#>             IR2       IR2/IR1              Group
#> 1          0.14          0.32       Raja clavata
#> 2          0.09          0.31       Raja clavata
#> 3          0.11          0.31       Raja clavata
#> 4          0.13          0.31       Raja clavata
#> 5 0.12 +/- 0.02 0.31 +/- 0.01       Raja clavata
#> 6          0.12          0.31 Dasyatis pastinaca
```
