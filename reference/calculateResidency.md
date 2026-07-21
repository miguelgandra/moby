# Calculate residency indices

Computes individual residency indices and the temporal building blocks
they are derived from, returning a tidy, fully numeric table (one row
per animal) suitable for plotting and downstream statistical analysis.
This is the numeric core used internally by
[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md)
(which formats these values for publication); use `calculateResidency()`
directly when you need the raw values rather than a formatted table.

Three indices, widely used in acoustic telemetry studies (Kraft et al.
2023; Appert et al. 2023), are available:

- **IR1** = Dd / Di, the proportion of days detected over the detection
  span (days at liberty); a maximum residency, considering only the
  period the animal was known to be alive and the tag operational.

- **IR2** = Dd / Dt, the proportion of days detected over the full study
  interval (release to last data download or tag expiration); a minimum
  residency.

- **IWR** = (Dd / Dt) \* (Di / Dt), the IR2 index weighted by the ratio
  of the detection interval to the study interval.

- **IR2/IR1**, the ratio of IR2 to IR1 (gap between last detection and
  end of monitoring).

where Dd = number of days with detections, Di = detection span (days at
liberty, first/release to last detection, inclusive), and Dt = study
interval (release to monitoring end).

## Usage

``` r
calculateResidency(
  data,
  tagging.dates = NULL,
  tag.durations = NULL,
  id.col = NULL,
  datetime.col = NULL,
  last.monitoring.date = NULL,
  residency.index = c("IR1", "IR2", "IWR"),
  start.point = "release",
  residency.by = NULL,
  cap = TRUE
)
```

## Arguments

- data:

  A data frame (or
  [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md))
  of detections.

- tagging.dates:

  A POSIXct vector of tagging/release dates (single value or named by
  ID). Inherited from the `mobyData` metadata when available.

- tag.durations:

  Optional numeric vector of tag battery durations (in days), used (with
  `last.monitoring.date`) to define the study interval Dt. Required
  (together with, or instead of, `last.monitoring.date`) when `IR2` or
  `IWR` are requested.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- datetime.col:

  Name of the column containing date-times in POSIXct format. Defaults
  to `"datetime"`.

- last.monitoring.date:

  Optional POSIXct value or named vector giving the last date data could
  be retrieved (last download / receivers operational). When both this
  and `tag.durations` are supplied, the shorter of the two defines the
  monitoring end per individual.

- residency.index:

  Character vector of indices to compute; any of `"IR1"`, `"IR2"`,
  `"IWR"`, `"IR2/IR1"`. Defaults to `c("IR1", "IR2", "IWR")`.

- start.point:

  Starting point for the detection span: `"release"` (default) or
  `"first.detection"`.

- residency.by:

  Optional column name used to additionally compute partial (spatially
  structured) residencies, one column per level (e.g. per habitat or
  array).

- cap:

  Logical; cap index values at 1 (their theoretical maximum). Defaults
  to TRUE. Set to FALSE to retain raw values (useful for diagnosing edge
  effects).

## Value

A data frame with one row per individual containing: the ID column,
`tagging_date`, `first_detection`, `last_detection`, `monitoring_end`
(POSIXct); `days_detected` (Dd), `detection_span` (Di) and
`monitoring_duration` (Dt) in days; one numeric column per requested
index; and, if `residency.by` is set, additional `"<index> <level>"`
partial-residency columns.

## References

Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., &
Abecasis, D. (2023). Residency and space use estimation methods based on
passive acoustic telemetry data. Movement Ecology, 11(1), 12.
https://doi.org/10.1186/s40462-023-00349-y

Appert, C., Udyawer, V., Simpfendorfer, C. A., et al. (2023). Use,
misuse, and ambiguity of indices of residence in acoustic telemetry
studies. Marine Ecology Progress Series, 714, 27-44.

## See also

[`summaryTable`](https://miguelgandra.github.io/moby/reference/summaryTable.md)

## Examples

``` r
data(rays)
# residency indices, using a fixed last-monitoring date to define the study interval (Dt)
res <- calculateResidency(rays,
         last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"))
#> Warning: - 'id.col' converted to factor.
head(res)
#>    ID tagging_date     first_detection      last_detection monitoring_end
#> 1 D01   2023-04-07 2023-04-08 15:05:10 2023-06-28 18:33:40     2023-12-31
#> 2 D02   2023-04-02 2023-04-02 17:11:13 2023-06-28 06:55:50     2023-12-31
#> 3 D03   2023-04-05 2023-04-09 11:28:11 2023-06-30 10:10:52     2023-12-31
#> 4 D04   2023-04-02 2023-04-05 21:58:36 2023-06-24 13:59:54     2023-12-31
#> 5 R01   2023-04-03 2023-04-07 20:38:30 2023-06-29 15:38:57     2023-12-31
#> 6 R02   2023-04-08 2023-04-08 07:08:58 2023-06-28 16:57:07     2023-12-31
#>   days_detected detection_span monitoring_duration       IR1        IR2
#> 1            31             83                 268 0.3734940 0.11567164
#> 2            25             88                 273 0.2840909 0.09157509
#> 3            24             87                 270 0.2758621 0.08888889
#> 4            30             84                 273 0.3571429 0.10989011
#> 5            38             88                 272 0.4318182 0.13970588
#> 6            23             82                 267 0.2804878 0.08614232
#>          IWR
#> 1 0.03582368
#> 2 0.02951871
#> 3 0.02864198
#> 4 0.03381234
#> 5 0.04519896
#> 6 0.02645569
```
