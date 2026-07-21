# Chronogram (hour x date heatmap)

Draws a chronogram: a two-dimensional hour-of-day by date plot of a
chosen metric (number of detections, individuals or co-occurring
animals) per time bin. Detections can be shown as size-coded points
(`style = "points"`) or colour-coded raster cells (`style = "raster"`),
with the diel cycle (dawn / sunrise / sunset / dusk boundaries, and
optional day/night/season shading) overlaid to reveal temporal patterns
of habitat use. When `split.by` is supplied, a separate panel is drawn
for each group (e.g. species).

## Usage

``` r
plotChronogram(
  data,
  id.col = NULL,
  timebin.col = NULL,
  station.col = NULL,
  metric = "detections",
  split.by = NULL,
  style = "points",
  color.by = NULL,
  color.pal = NULL,
  sunriset.coords = NULL,
  solar.depth = 18,
  diel.lines = 4,
  shade = FALSE,
  date.interval = "auto",
  date.format = NULL,
  date.start = 1,
  pt.cex = c(0.5, 3),
  alpha = 0.7,
  background.color = "white",
  grid = TRUE,
  grid.color = "white",
  highlight.isolated = FALSE,
  shared.scale = FALSE,
  shared.dates = TRUE,
  main = NULL,
  legend = TRUE,
  legend.cols = NULL,
  cex = 1,
  ncol = NULL,
  file = NULL,
  width = NULL,
  height = NULL,
  res = 300,
  ...
)
```

## Arguments

- data:

  A data frame of detections including a time-bin column
  (`timebin.col`).

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- timebin.col:

  Name of the column containing time bins (in POSIXct format). Defaults
  to `"timebin"`.

- station.col:

  Name of the column containing station/receiver IDs. Defaults to
  `"station"`.

- metric:

  The metric to plot: one of `"detections"` (the default),
  `"individuals"` or `"co-occurrences"` (see Details for the
  co-occurrence definition and its caveats).

- split.by:

  Optional column name; when supplied, a separate panel is drawn for
  each of its levels (e.g. species, life stage, habitat).

- style:

  `"points"` (metric encoded by point size) or `"raster"` (metric
  encoded by colour).

- color.by:

  For `style = "points"`, a column used to colour the points (factor or
  numeric).

- color.pal:

  Colours (or a palette function) for the points/raster. If NULL, a
  colourblind-safe default is used (Okabe-Ito for categorical, viridis
  for continuous/raster).

- sunriset.coords:

  Longitude/latitude (a length-2 numeric, matrix or `SpatialPoints`) at
  which to compute diel-phase times. Required only when `diel.lines > 0`
  or `shade = "diel"`.

- solar.depth:

  Sun angle below the horizon (degrees) defining twilight. Defaults to
  18.

- diel.lines:

  Number of diel-boundary lines to draw: 0, 2 (sunrise/sunset) or 4
  (dawn, sunrise, sunset, dusk). Defaults to 4.

- shade:

  Background shading for `style = "points"`. One of `"diel"`
  (day/night/twilight), `"season"`, FALSE (none, default), or a data
  frame of custom periods with columns `start` and `end` (POSIXct or
  coercible) and, optionally, `label` and `color` (same structure as
  [`plotAbacus`](https://miguelgandra.github.io/moby/reference/plotAbacus.md)).

- date.interval:

  Date-axis labels: `"auto"` (default) for span-appropriate "pretty"
  breaks, or a numeric value selecting every n-th formatted date (manual
  mode).

- date.format:

  Date format for the labels
  ([`strftime`](https://rdrr.io/r/base/strptime.html)). If NULL
  (default), chosen automatically in `"auto"` mode and defaults to
  `"%d/%b"` in manual mode.

- date.start:

  In manual mode, the index of the first labelled date; in `"auto"`
  mode, a phase offset re-anchoring the labels within each period.
  Defaults to 1.

- pt.cex:

  Length-2 numeric giving the min and max point sizes
  (`style = "points"`). Defaults to c(0.5, 3).

- alpha:

  Opacity of the points, from 0 (transparent) to 1 (opaque). Defaults to
  0.7.

- background.color:

  Background colour (`style = "points"`, when `shade` is not used).

- grid:

  Logical; draw faint hour/date guide lines. Defaults to TRUE.

- grid.color:

  Colour of the grid lines. Defaults to "white".

- highlight.isolated:

  Logical; for `style = "points"` with `color.by`, plot isolated colour
  levels on top so they are not hidden behind denser clusters. Defaults
  to FALSE.

- shared.scale:

  Use a common metric (colour/size) scale across panels. Defaults to
  FALSE.

- shared.dates:

  Use a common date range across panels. Defaults to TRUE.

- main:

  Panel title(s). If NULL (default), titles are generated automatically
  (the group name when `split.by` is used, otherwise the metric); a
  character vector overrides them (recycled across panels); FALSE omits
  titles.

- legend:

  Logical; draw the legend(s). Defaults to TRUE.

- legend.cols:

  Number of columns for the `color.by` legend. If NULL, chosen
  automatically.

- cex:

  Global expansion factor for all plot text. Defaults to 1.

- ncol:

  Number of columns in the panel layout. If NULL, panels are stacked
  (one column).

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

- ...:

  Further arguments passed to
  [`image`](https://rdrr.io/r/graphics/image.html) (raster) or
  [`points`](https://rdrr.io/r/graphics/points.html) (points).

## Value

Called for its side effect: it generates chronogram panel(s).

## Details

Date-axis labels are chosen automatically to suit the temporal span
(`date.interval = "auto"`); a numeric `date.interval` reverts to manual
mode, and `date.start` shifts the labels (a phase offset in auto mode).
For `style = "points"`, `color.by` colours points by any variable while
point size encodes the metric; `shade` shades the diel cycle or seasons
in the background.

The `"co-occurrences"` metric counts the number of distinct individuals
sharing the same station within the same time bin (i.e. the local group
size, restricted to cells with two or more animals). It is therefore an
exploratory descriptor whose magnitude is sensitive to the time-bin
width and the spatial resolution of `station.col` (coarser bins/stations
inflate apparent co-occurrence) and to detection effort; it is not
corrected for the number of animals available. For quantitative
co-occurrence/association analysis (with permutation-based null models)
use
[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md).

## Note

The right margin is sized from the measured legend width and legends are
positioned in inch units (via
[`grconvertX`](https://rdrr.io/r/graphics/convertXY.html) /
[`grconvertY`](https://rdrr.io/r/graphics/convertXY.html)), so the
layout stays consistent across datasets and graphics devices.

## Examples

``` r
# chronogram of detections with the diel cycle overlaid
plotChronogram(rays, sunriset.coords = c(-9, 38.4))
#> Warning: - 'id.col' converted to factor.
#> 
#> Chronogram
#> ------------------------------------------------------
#>   Individuals:     8
#>   Detections:      1,643
#>   Period:          2023-04-02 to 2023-06-30 (88 d)
#>   Time bin:        60 min
#>   Metric:          detections
#>   Style:           points
#>   Diel:            4 lines
#>   Date axis:       auto
#> ------------------------------------------------------


# without the diel overlay (no coordinates required)
plotChronogram(rays, diel.lines = 0)
#> Warning: - 'id.col' converted to factor.
#> 
#> Chronogram
#> ------------------------------------------------------
#>   Individuals:     8
#>   Detections:      1,643
#>   Period:          2023-04-02 to 2023-06-30 (88 d)
#>   Time bin:        60 min
#>   Metric:          detections
#>   Style:           points
#>   Diel:            off
#>   Date axis:       auto
#> ------------------------------------------------------
```
