# Generate summary table for tagged animals

This function generates a comprehensive summary table for tagged
animals, including tagging and last detection dates, along with various
monitoring and residency metrics. It also allows the inclusion of
additional sensor data (e.g., depth, temperature) and provides an option
to summarize these readings by calculating the mean, minimum, and
maximum values for each specified sensor column. Additionally, the
function can incorporate additional metadata and automatically computes
overall means and error metrics.

## Usage

``` r
summaryTable(
  data,
  id.metadata = NULL,
  id.col = NULL,
  datetime.col = NULL,
  station.col = NULL,
  id.groups = NULL,
  tagging.dates = NULL,
  tag.durations = NULL,
  sensor.cols = NULL,
  sensor.titles = NULL,
  residency.index = c("IR1", "IR2", "IR2/IR1"),
  start.point = "release",
  last.monitoring.date = NULL,
  residency.by = NULL,
  error.stat = "sd",
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- data:

  A data frame containing animal detections. Each row should represent
  an individual detection event, unless a 'detections' column is
  included to indicate the number of detections for each row.

- id.metadata:

  A data frame containing metadata about the tagged animals, such as
  their length, sex, or transmitter type. All columns in this data frame
  will be summarized and included in the final table. If there are
  multiple rows per animal, variables will be collapsed before merging
  with other statistics.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- datetime.col:

  Name of the column containing date-times in POSIXct format. Defaults
  to `"datetime"`.

- station.col:

  Name of the column containing station/receiver IDs. Defaults to
  `"station"`.

- id.groups:

  Optional. A list where each element is a group of IDs, used for
  visually aggregating animals belonging to the same class (e.g.,
  different species or life stages). If supplied, averages will be
  calculated independently for each group.

- tagging.dates:

  Optional POSIXct vector of tagging/release dates. Either a single
  value (applied to all individuals) or a named vector whose names match
  the animal IDs.

- tag.durations:

  Optional. A numeric vector containing the estimated battery duration
  of the deployed tags (in days). This parameter must be either:

  - A single numeric value, which will be applied to all unique animal
    IDs; or

  - A named numeric vector, where the names correspond to the animal IDs
    in the `id.col` column. If multiple tag durations are provided, the
    vector must include all IDs and will be reordered to align with the
    levels of `id.col`.

- sensor.cols:

  Optional. A character vector specifying column names in `data` that
  contain additional sensor readings (e.g., depth, temperature). For
  each column specified, the mean, minimum, and maximum values will be
  calculated and incorporated into the summary.

- sensor.titles:

  Optional. Titles for the sensor readings in `sensor.cols`. If not
  provided, defaults to the column names.

- residency.index:

  A character string specifying the type of residency index to
  calculate. Options include:

  - "IR1": Residency Index 1, calculated as the number of days the
    animal was detected (Dd) divided by the detection interval (Di),
    i.e., the number of days between release/first detection and last
    detection (days at liberty). This represents a maximum residency
    value, considering only the period for which the animal was known to
    be alive and the tag operational.

  - "IR2": Residency Index 2, calculated as the number of days the
    animal was detected (Dd) divided by the study interval (Dt), i.e.,
    the total number of days between release/first detection and last
    data download or tag expiration date. This approach provides a
    minimum residency value, assuming the animal was alive and
    detectable throughout the study period.

  - "IWR": Weighted Residency Index, which corresponds to the IR2 index
    weighted by the ratio between the detection interval (Di, the number
    of days between the first and last detection) and the study interval
    (Dt, the total monitoring period). This accounts for the number of
    days detected and the spread of detections within the monitoring
    period, providing a measure of residency that balances the frequency
    of detections with their temporal distribution.

  - "IR2/IR1": The ratio of IR2 to IR1, providing a measure of the gap
    between the last detection and the end of the monitoring period.

  The choice of index can affect the interpretation of residency
  patterns, so it's important to select the one(s) that best fits the
  study objectives. Further information on residency estimation methods
  can be found in Kraft et al. (2023) and Appert et al. (2023) - see
  below in the references section.

- start.point:

  A character string specifying the starting point for calculating days
  at liberty. Options include:

  - `"release"`: The release day is used as the starting point for
    calculating days at liberty.

  - `"first.detection"`: The first detection after tagging is used as
    the starting point.

  Defaults to `"release"`.

- last.monitoring.date:

  Optional. A POSIXct object or a named vector of POSIXct objects
  specifying the last timestamp when data could be retrieved, typically
  corresponding to the last data download date or the final day
  receivers were operational. If a single value is provided, it will be
  applied to all individuals. If a named vector is provided, the names
  should correspond to individual IDs, allowing for unique timestamps
  per individual. When `tag.durations` are also supplied, the total
  monitoring duration for each individual will be estimated based on the
  shortest of the two values: the tag expiration date or the last
  monitoring day.

- residency.by:

  Optional. Variable used to calculate partial residencies (e.g. array
  or habitat). Defaults to NULL.

- error.stat:

  The statistic to use for variability/error calculation, either 'sd'
  (standard deviation) or 'se' (standard error). Defaults to 'sd'.

- verbose:

  Logical; print a summary of the operation. Defaults to
  `getOption("moby.verbose", TRUE)`.

## Value

A `mobyTable`: a TYPED data frame (counts stay integer, indices numeric,
dates POSIXct), one row per tagged animal, so the result can be computed
on directly. Presentation - fixed precision, the display-only
`mean +/- error` row, group headings - is applied by
[`format`](https://miguelgandra.github.io/moby/reference/format.mobyTable.md)
and
[`print`](https://miguelgandra.github.io/moby/reference/print.mobyTable.md).
Export the rendered version with
`write.csv(format(x), file, row.names = FALSE)`.

Columns use stable snake_case names, which is what `format(decimals=)`
and `format(group.by=)` are keyed on; the publication headers live in
`format(style = "report")`.

- `ID` (or your `id.col`), and a `group` factor when `id.groups` names
  more than one group

- `tagging_date`, `last_detection` (POSIXct)

- `n_detections`, `n_receivers`, `monitoring_duration_d`,
  `detection_span_d`, `n_days_detected`

- one column per requested `residency.index` (`IR1`, `IR2`, `IR2/IR1`,
  ...), plus the `residency.by` variants when supplied

- `<sensor>_mean`, `<sensor>_min`, `<sensor>_max` for each of
  `sensor.cols`

- any columns carried over from `id.metadata`

## See also

[`calculateResidency`](https://miguelgandra.github.io/moby/reference/calculateResidency.md)
for the underlying numeric residency indices (useful for plotting or
downstream statistical analysis).

## Examples

``` r
data(rays)
# per-animal summary with residency indices, grouped by species
summaryTable(rays,
             last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"),
             id.groups = mobyMeta(rays)$id.groups)
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
#> <mobyTable: summary> 8 rows (2 groups; one mean ± sd row per group)
#> 
#> ─ Raja clavata
#>         ID tagging_date last_detection n_detections n_receivers
#>        R01   03/04/2023     29/06/2023          283           6
#>        R02   08/04/2023     28/06/2023          160           5
#>        R03   10/04/2023     29/06/2023          207           6
#>        R04   01/04/2023     25/06/2023          261           6
#>  mean ± sd            -              -     228 ± 55       6 ± 0
#>  monitoring_duration_d detection_span_d n_days_detected         IR1         IR2
#>                    272               88              38        0.43        0.14
#>                    267               82              23        0.28        0.09
#>                    265               81              30        0.37        0.11
#>                    274               86              36        0.42        0.13
#>                270 ± 4           84 ± 3          32 ± 7 0.38 ± 0.07 0.12 ± 0.02
#>      IR2/IR1
#>         0.32
#>         0.31
#>         0.31
#>         0.31
#>  0.31 ± 0.01
#> 
#> ─ Dasyatis pastinaca
#>         ID tagging_date last_detection n_detections n_receivers
#>        D01   07/04/2023     28/06/2023          249           6
#>        D02   02/04/2023     28/06/2023          154           6
#>        D03   05/04/2023     30/06/2023          160           6
#>        D04   02/04/2023     24/06/2023          169           6
#>  mean ± sd            -              -     183 ± 44       6 ± 0
#>  monitoring_duration_d detection_span_d n_days_detected         IR1         IR2
#>                    268               83              31        0.37        0.12
#>                    273               88              25        0.28        0.09
#>                    270               87              24        0.28        0.09
#>                    273               84              30        0.36        0.11
#>                271 ± 2           86 ± 2          28 ± 4 0.32 ± 0.05 0.10 ± 0.01
#>      IR2/IR1
#>         0.31
#>         0.32
#>         0.32
#>         0.31
#>  0.32 ± 0.01
```
