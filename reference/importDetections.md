# Import and harmonise acoustic detection data

Reads acoustic-telemetry detections from a range of common sources and
harmonises them into a single consistent schema, returning a
[`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
object ready for the rest of the `moby` workflow. Supported sources are
Innovasea/VEMCO VUE exports, Innovasea VDAT/Fathom `DET.csv` files, and
detection extracts from the GLATOS, OTN and ETN (`etn` package) systems.
A `generic` mode plus a user-supplied `col.map` handle non-standard
layouts.

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

  Optional named list mapping canonical fields (`datetime`,
  `transmitter`, `receiver`, `station`, `lon`, `lat`, `ID`,
  `sensor_value`, `sensor_unit`, ...) to the column name(s) in `x`.
  Merged over (and overriding) the chosen `source` preset.

- datetime.format:

  Optional explicit `strptime` format for the datetime column; if
  `NULL`, common layouts are auto-detected.

- keep.extra:

  Logical; retain source columns that were not mapped to a canonical
  field. Defaults to `FALSE`.

## Value

A [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
object with harmonised columns (`ID`, `datetime`, `transmitter`,
`receiver`, `station`, `lon`, `lat`, ...). When the source has no animal
identifier, `ID` is initialised from `transmitter` (assign true animal
IDs later by joining tag metadata).

## See also

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
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: timebin.
head(det)
#> <mobyData> 6 records x 7 columns
#>   individuals: 1  (id.col = 'ID')
#>   period: 2023-04-08 15:05:10 to 2023-04-09 12:54:09 (tz = UTC)
#>   columns: datetime.col='datetime', timebin.col='timebin', station.col='station', lon.col='lon', lat.col='lat'
```
