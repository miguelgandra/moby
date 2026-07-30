# Import and harmonise tag / animal metadata

Reads a tag/animal metadata table and harmonises it into a consistent
schema (`transmitter`, `ID`, `tagging_date`, plus biometrics such as
`species`, `sex`, `length`, and the transmitter `nominal_delay` when the
source provides one). Use together with
[`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
to attach animal IDs, tagging dates and nominal delays to a detection
dataset.

## Usage

``` r
importTags(
  x,
  source = c("vue", "glatos", "otn", "etn", "generic"),
  tz = "UTC",
  col.map = NULL,
  datetime.format = NULL,
  keep.extra = TRUE,
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- x:

  A path to a `.csv`/`.xlsx` tag-metadata file, or a data frame (e.g.
  from `etn::get_tags()` / `etn::get_animals()`).

- source:

  One of `"vue"`, `"glatos"`, `"otn"`, `"etn"` or `"generic"`.

- tz:

  Time zone used to parse the tagging date. Defaults to `"UTC"`.

- col.map:

  Optional named list mapping canonical tag fields to the column name(s)
  in `x`, merged over (and extending) the `source` preset. See
  [moby_import_schema](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
  for the full list of tag fields and examples.

- datetime.format:

  Optional explicit `strptime` format for the tagging-date column. When
  `NULL` (default) the layout is inferred, but only when the answer is
  unambiguous: a column such as `"07/06/2013"`, where day-first and
  month-first both fit and disagree, raises an error asking for this
  argument rather than silently choosing one. A column no layout can
  read is an error too, so an unreadable column can never travel on as
  silent `NA`s.

  Timezone handling depends on the input type. TEXT carries no zone, so
  `tz` is applied as the zone while parsing. A column that is ALREADY
  `POSIXct` (as when a data frame is passed in, e.g. from an API) is an
  absolute instant chosen by the caller: it is never reinterpreted, and
  `tz` changes only how it is displayed.

- keep.extra:

  Logical; retain unmapped source columns. Defaults to `TRUE` so that
  additional biometric fields are preserved.

- verbose:

  Logical; print a summary of the operation. Defaults to
  `getOption("moby.verbose", TRUE)`.

## Value

A data frame with at least `transmitter` and (when available) `ID`,
`tagging_date` (POSIXct), `nominal_delay` (seconds) and biometric
columns.

## Details

When the source specifies a delay RANGE (`min_delay`/`max_delay`) rather
than a nominal delay - as many tag-specification exports do -
`nominal_delay` is derived as their midpoint. The nominal delay is what
[`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
uses to scale its short-interval (min_lag) false-detection filter, so
carrying it through here means that filter can be enabled automatically.

## See also

[`moby_import_schema`](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
for the canonical field list;
[`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md),
[`importDetections`](https://miguelgandra.github.io/moby/reference/importDetections.md)

## Examples

``` r
# harmonise a tag-metadata table (here the bundled 'rays_tags' data frame)
tags <- importTags(rays_tags, source = "generic",
                   col.map = list(ID = "ID", transmitter = "transmitter",
                                  tagging_date = "tagging_date", species = "species"))
#> ── importTags() ──────────────────────────────────────────────────────── moby ──
#> 
#> ℹ Reading and harmonising tag and animal metadata
#> • Input: data frame
#> 
#> → Method
#>   • source    generic (user col.map)
#>   • timezone  UTC
#> 
#> ✔ 8 tags imported
#> ℹ Fields resolved: ID, transmitter, tagging_date, species
head(tags)
#>    ID    transmitter tagging_date            species tagging_station
#> 1 R01 A69-1602-30001   2023-04-03       Raja clavata            ST02
#> 2 R02 A69-1602-30002   2023-04-08       Raja clavata            ST03
#> 3 R03 A69-1602-30003   2023-04-10       Raja clavata            ST06
#> 4 R04 A69-1602-30004   2023-04-01       Raja clavata            ST03
#> 5 D01 A69-1602-30005   2023-04-07 Dasyatis pastinaca            ST03
#> 6 D02 A69-1602-30006   2023-04-02 Dasyatis pastinaca            ST02
```
