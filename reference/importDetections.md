# Import and harmonise acoustic detection data

Reads acoustic-telemetry detections from a range of common sources and
harmonises them into a single consistent schema - moby's canonical
column names, with date-times parsed - ready to be turned into a
[`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
by
[`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
/ [`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md).
Supported sources are Innovasea/VEMCO VUE exports, Innovasea VDAT/Fathom
`DET.csv` files, and detection extracts from the GLATOS, OTN and ETN
(`etn` package) systems. A `generic` mode plus a user-supplied `col.map`
handle non-standard layouts.

## Usage

``` r
importDetections(
  x,
  source = c("vue", "vdat", "glatos", "otn", "etn", "generic"),
  tz = "UTC",
  col.map = NULL,
  datetime.format = NULL,
  keep.extra = FALSE
)
```

## Arguments

- x:

  A path to a `.csv` (or `.xlsx`) detection file, or a data frame
  already loaded in R (e.g. the output of
  `etn::get_acoustic_detections()` or
  `glatos::read_glatos_detections()`).

- source:

  One of `"vue"`, `"vdat"`, `"glatos"`, `"otn"`, `"etn"` or `"generic"`.
  For `"generic"`, supply `col.map`.

- tz:

  Time zone used to parse date-times. Defaults to `"UTC"` (the
  convention for GLATOS/OTN/ETN); set explicitly for VUE/VDAT exports
  recorded in another zone.

- col.map:

  Optional named list mapping canonical fields to the column name(s) in
  `x`, merged over (and overriding) the chosen `source` preset. See
  [moby_import_schema](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
  for the full list of detection fields, which are required, and worked
  examples.

- datetime.format:

  Optional explicit `strptime` format for the datetime column; if
  `NULL`, common layouts are auto-detected.

- keep.extra:

  Logical; retain source columns that were not mapped to a canonical
  field. Defaults to `FALSE`.

## Value

A data frame with harmonised columns (`ID`, `datetime`, `transmitter`,
`receiver`, `station`, `lon`, `lat`, ...), sorted by animal and time.
When the source has no animal identifier, `ID` is initialised from
`transmitter` (assign true animal IDs afterwards with
[`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)).

This is a *harmonised table*, not yet a
[`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md):
importing reshapes columns, it does not attach study metadata. Turn it
into a `mobyData` with
[`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
(animal IDs, tagging dates, nominal delays) and/or
[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md)
(`id.groups`, CRS, land layer) - see the ‘Which function do I use?’
section of
[`moby_import_schema`](https://miguelgandra.github.io/moby/reference/moby_import_schema.md).
Because the output already uses moby's canonical column names, those
functions need no explicit column arguments.

## See also

[`moby_import_schema`](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
for the canonical field list and how importing relates to
[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md);
[`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md),
[`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md),
[`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md),
[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md)

## Examples

``` r
# harmonise a generic-layout detection CSV via an explicit column map
csv <- system.file("extdata", "rays_detections.csv", package = "moby")
det <- importDetections(csv, source = "generic",
                        col.map = list(ID = "animal_id", datetime = "timestamp",
                                       station = "station_name", lon = "deploy_longitude",
                                       lat = "deploy_latitude", receiver = "receiver_id",
                                       transmitter = "transmitter"))
head(det)
#>    ID            datetime    transmitter  receiver station    lon    lat
#> 1 D01 2023-04-08 15:05:10 A69-1602-30005 VR2W-1003    ST03 -8.996 38.456
#> 2 D01 2023-04-08 15:14:59 A69-1602-30005 VR2W-1003    ST03 -8.996 38.456
#> 3 D01 2023-04-08 15:52:10 A69-1602-30005 VR2W-1003    ST03 -8.996 38.456
#> 4 D01 2023-04-08 16:51:29 A69-1602-30005 VR2W-1003    ST03 -8.996 38.456
#> 5 D01 2023-04-08 16:51:56 A69-1602-30005 VR2W-1003    ST03 -8.996 38.456
#> 6 D01 2023-04-09 12:54:09 A69-1602-30005 VR2W-1006    ST06 -8.990 38.442
```
