# Assign animal IDs (and tagging dates) to detections from tag metadata

Joins tag metadata (see
[`importTags`](https://miguelgandra.github.io/moby/reference/importTags.md))
to a detection dataset on the transmitter code, assigning each detection
an animal `ID`. Matching is tolerant: it first tries the full
transmitter string and then falls back to the trailing numeric code (so
`"A69-1602-111"` matches a tag stored as `"111"`). When the tag table
carries tagging dates, these are attached to the returned
[`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
object's metadata (and can flow automatically into functions such as
[`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
and
[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md)).
Optional biometric columns can be joined in as well.

## Usage

``` r
assignAnimalIDs(
  detections,
  tags,
  transmitter.col = "transmitter",
  id.col = NULL,
  keep.cols = NULL,
  set.tagging.dates = TRUE,
  set.nominal.delay = TRUE
)
```

## Arguments

- detections:

  A detection dataset (`mobyData` or data frame) with a transmitter
  column.

- tags:

  A harmonised tag table from
  [`importTags`](https://miguelgandra.github.io/moby/reference/importTags.md)
  (or a data frame with at least a `transmitter` column).

- transmitter.col:

  Name of the transmitter column in `detections`. Defaults to
  `"transmitter"`.

- id.col:

  Name of the animal-ID column to (re)create in `detections`. Resolved
  from the `mobyData` metadata or `"ID"` when `NULL`.

- keep.cols:

  Optional character vector of additional `tags` columns (e.g. `"sex"`,
  `"length"`, `"species"`) to join into the detections.

- set.tagging.dates:

  Logical; if `TRUE` (default) and `tags` has a `tagging_date` column,
  attach per-individual tagging dates to the returned object's metadata.

- set.nominal.delay:

  Logical; if `TRUE` (default) and `tags` has a `nominal_delay` column
  (see
  [`importTags`](https://miguelgandra.github.io/moby/reference/importTags.md)),
  attach per-individual transmitter nominal delays (seconds) to the
  returned object's metadata.
  [`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
  reads these automatically to scale its short-interval (min_lag)
  false-detection filter, so arrays mixing tag families (e.g. 60 s and
  120 s tags) are handled per animal.

## Value

A [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
object with the `ID` column assigned (and, optionally, tagging dates and
biometric columns attached). Detections whose transmitter is absent from
`tags` keep `NA` IDs (with a warning).

## See also

[`importTags`](https://miguelgandra.github.io/moby/reference/importTags.md),
[`importDetections`](https://miguelgandra.github.io/moby/reference/importDetections.md),
[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md)

## Examples

``` r
# join tag metadata to detections to assign animal IDs (and tagging dates)
tags <- importTags(rays_tags, source = "generic",
                   col.map = list(ID = "ID", transmitter = "transmitter",
                                  tagging_date = "tagging_date"))
# detections carrying the tagged rays' transmitters
det <- rays_detections[rays_detections$transmitter %in% rays_tags$transmitter, ]
det <- assignAnimalIDs(det, tags)
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: timebin.
levels(det$ID)
#> [1] "D01" "D02" "D03" "D04" "R01" "R02" "R03" "R04"
```
