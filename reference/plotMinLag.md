# Diagnose the min_lag false-detection threshold

Plots the empirical distribution of `min_lag` - each detection's time
gap to the nearest other detection of the same transmitter on the same
receiver - normalised by the transmitter nominal delay, and reports the
proportion of detections that the short-interval false-detection filter
of
[`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
would flag at several candidate thresholds. It is a quick check on
whether the default `min.lag.factor = 30` is appropriate for a given
dataset, rather than assuming it.

## Usage

``` r
plotMinLag(
  data,
  id.col = NULL,
  datetime.col = NULL,
  station.col = NULL,
  nominal.delay = NULL,
  factors = c(20, 30, 50),
  remove.duplicates = TRUE,
  bar.color = "#0072B2",
  background.color = "grey96",
  breaks = 40,
  main = NULL,
  cex = 1,
  file = NULL,
  width = NULL,
  height = NULL,
  res = 300
)
```

## Arguments

- data:

  A `mobyData` object or a data frame of detections.

- id.col:

  Name of the column with animal IDs. Resolved from the metadata /
  canonical default ("ID") when not supplied.

- datetime.col:

  Name of the column with detection timestamps (`POSIXct`). Resolved
  from the metadata / canonical default ("datetime") when not supplied.

- station.col:

  Name of the receiver/station column. Resolved from the metadata /
  canonical default ("station") when not supplied.

- nominal.delay:

  Transmitter nominal (mean) delay, in seconds. A single value applied
  to all individuals, or a vector named by `id.col` for mixed tag
  families. Read from the `mobyData` metadata (`nominal.delay`) when not
  supplied. Required: the diagnostic normalises `min_lag` by it.

- factors:

  Reference threshold multipliers (of the nominal delay) to overlay on
  the plot and tabulate. Defaults to `c(20, 30, 50)`; the middle value
  is the
  [`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
  default.

- remove.duplicates:

  Logical; drop exact-duplicate records (same ID, timestamp and station)
  before computing `min_lag`, matching the stage-0 behaviour of
  [`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md).
  Defaults to TRUE.

- bar.color:

  Fill colour for the histogram. Defaults to a colourblind-safe blue.

- background.color:

  Panel background colour. Defaults to "grey96"; `NA` draws none.

- breaks:

  Number of histogram bins on the log10 axis. Defaults to 40.

- main:

  Optional plot title.

- cex:

  Global expansion factor for labels. Defaults to 1.

- file:

  Optional output file. If `NULL` (the default), the figure is drawn on
  the current graphics device - the usual interactive behaviour. If a
  file path is supplied, moby opens a graphics device chosen from the
  file extension (`.pdf`, `.svg`, `.png`, `.jpg`/`.jpeg`, `.tif`/`.tiff`
  or `.bmp`), draws the figure to it, and closes the device
  automatically (also if an error occurs). For multi-page or batch
  workflows, keep `file = NULL` and manage the device yourself (e.g.
  `pdf(...); plot...(); dev.off()`).

- width, height:

  Output size in inches. Used *only* when `file` is supplied. If `NULL`
  (default), a size is derived from the figure's structure (e.g. the
  number of individuals or panels). These defaults are sensible starting
  points, *not* guarantees: dense or unusual figures may still need you
  to set `width`/`height` explicitly. The two can be set independently.

- res:

  Resolution in pixels per inch, for raster formats only (`.png`,
  `.jpg`, `.tif`, `.bmp`); ignored for vector formats (`.pdf`, `.svg`).
  Used only when `file` is supplied. Defaults to 300.

## Value

Invisibly, a data frame with one row per `factors` value: the multiplier
(`factor`), the threshold in seconds (`threshold_s`; `NA` when
`nominal.delay` varies across individuals), and the number (`n_flagged`)
and percentage (`pct_flagged`) of detections flagged - i.e. with
`min_lag` above the threshold, or a lone decode - at that multiplier.

## Details

A genuinely present tag produces bursts of closely spaced detections
(small `min_lag`), whereas an isolated spurious decode has no nearby
companion (large `min_lag`). The two typically appear as separate modes
on the log axis, and a defensible threshold sits in the valley between
them. Read the plot rather than a single number: if the proportion
flagged changes little across the reference multipliers, the choice is
robust; if it changes sharply, the threshold matters and should be
justified for that dataset.

`min_lag` is computed with the same internal helper used by
[`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md),
so - on the same input - the flagged proportions reported here match
what the filter would remove. The diagnostic applies only the optional
duplicate-removal step (not the pre-tagging / cut-off steps), so run it
on the detections you intend to filter.

## See also

[`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md)

## Examples

``` r
# assess whether the default 30x min_lag threshold suits the data (120 s nominal delay)
flagged <- plotMinLag(rays, nominal.delay = 120)
#> Warning: - 'id.col' converted to factor.

flagged
#>   factor threshold_s n_flagged pct_flagged
#> 1     20        2400       206       12.54
#> 2     30        3600        87        5.30
#> 3     50        6000        26        1.58
```
