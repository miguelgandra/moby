# Assign animal IDs (and tagging dates) to detections from tag metadata

Resolves the transmitter codes in a detection dataset to the animals
that carried them, using a tag table (see
[`importTags`](https://miguelgandra.github.io/moby/reference/importTags.md)).
This is more than a column join: several transmitter codes can map to
the same animal, and the function also derives the per-animal metadata -
tagging dates and transmitter nominal delays - that later steps such as
[`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
and
[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md)
read automatically. Optional biometric columns can be joined in as well.

## Usage

``` r
assignAnimalIDs(
  detections,
  tags,
  id.col = NULL,
  transmitter.col = "transmitter",
  keep.cols = NULL,
  set.tagging.dates = TRUE,
  set.nominal.delay = TRUE,
  verbose = getOption("moby.verbose", TRUE)
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

- id.col:

  Name of the animal-ID column to (re)create in `detections`. Resolved
  from the `mobyData` metadata or `"ID"` when `NULL`.

- transmitter.col:

  Name of the transmitter column in `detections`. Defaults to
  `"transmitter"`.

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

- verbose:

  Logical; print a summary of the operation. Defaults to
  `getOption("moby.verbose", TRUE)`.

## Value

A [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
object with the `ID` column assigned (and, optionally, tagging dates and
biometric columns attached). Detections whose transmitter is absent from
`tags` keep `NA` IDs (with a warning).

## Details

**Why several codes can mean one animal.** Most often this is because a
single physical transmitter emits on more than one code space: sensor
tags (for example depth-sensing Vemco V16P transmitters) report their
identity and their sensor data under different code spaces, so one tag
appears in the detection file under two or more transmitter codes.
Animals genuinely carrying several transmitters at once are uncommon;
the other routine case is an animal that was recaptured and re-tagged,
so its history spans consecutive tags. List each code as its own row of
`tags`, mapped to the same animal `ID`, and all of them resolve to that
animal.

**Matching.** The full transmitter string is matched first. A
trailing-numeric fallback then covers tag tables that store the bare
code (so `"A69-1602-111"` matches a tag stored as `"111"`). The fallback
is deliberately conservative: a detection whose code space differs from
the candidate tag's is left unmatched rather than assigned, because
`"A69-9999-111"` and `"A69-1602-111"` are different transmitters and
guessing would attribute another project's animal to yours. Such codes
are reported, so you can add them to `tags` if they are further code
spaces of your own tags.

**Where it sits in the pipeline.** This step is independent of
[`matchDeployments`](https://miguelgandra.github.io/moby/reference/matchDeployments.md):
one resolves transmitters to animals, the other resolves receivers to
deployment windows, and they key on different columns. Either order
gives the same result.

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
#> ℹ Fields resolved: ID, transmitter, tagging_date
# detections carrying the tagged rays' transmitters
det <- rays_detections[rays_detections$transmitter %in% rays_tags$transmitter, ]
det <- assignAnimalIDs(det, tags)
#> ── assignAnimalIDs() ─────────────────────────────────────────────────── moby ──
#> 
#> ℹ Joining tag metadata to replace placeholder IDs with real animal IDs
#> • Input: 1,643 detections · 8 transmitters
#> 
#> ✔ 1,643 detections matched to 8 animals
#> ℹ Metadata attached: tagging.dates (8 individuals)
levels(det$ID)
#> [1] "D01" "D02" "D03" "D04" "R01" "R02" "R03" "R04"
```
