# Abacus plot

Produces an abacus plot: a timeline of colour-coded detections for each
tagged individual over the monitoring period. Animals can be grouped
(e.g. by species, sex or age) via `id.groups`, detections can be
coloured by any variable (e.g. receiver/site) through `color.by`, and
per-animal release and estimated tag-expiry dates are marked with
dedicated symbols. Optional background shading (annual seasons or custom
periods) and a labelled top time-band provide temporal context.

## Usage

``` r
plotAbacus(
  data,
  id.col = NULL,
  datetime.col = NULL,
  id.groups = NULL,
  tagging.dates = NULL,
  tag.durations = NULL,
  color.by = NULL,
  color.pal = NULL,
  discard.missing = FALSE,
  shade = TRUE,
  top.band = "%Y",
  legend = TRUE,
  main = NULL,
  bin = NULL,
  scale.by.count = FALSE,
  date.interval = "auto",
  date.format = NULL,
  date.start = 1,
  pch = 16,
  pt.cex = 1,
  alpha = 1,
  event.pch = c(release = 8, expiry = 4),
  highlight.isolated = TRUE,
  background.color = "grey96",
  legend.cols = NULL,
  cex = 1,
  file = NULL,
  width = NULL,
  height = NULL,
  res = 300,
  ...
)
```

## Arguments

- data:

  A data frame containing animal detections.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- datetime.col:

  Name of the column containing date-times in POSIXct format. Defaults
  to `"datetime"`.

- id.groups:

  Optional named list of ID groups used to visually aggregate animals
  belonging to the same class (e.g. species, sex or age). Each element
  is a vector of IDs in that group.

- tagging.dates:

  A POSIXct vector of tagging/release dates (single value or named by
  ID). Inherited from the `mobyData` metadata when available.

- tag.durations:

  Optional numeric vector of estimated tag battery durations (in days),
  used to mark the estimated tag-expiry date of each animal. Either a
  single value (applied to all IDs) or a named vector matching the IDs
  in `id.col`.

- color.by:

  Name of the column used to colour-code detections (e.g. station,
  habitat or a temporal class). If NULL, all detections are drawn in a
  single colour.

- color.pal:

  Colours used to plot detections, one per `color.by` level. May be a
  plain vector (matched to the levels by position) or a **named** vector
  (matched by name, which is safer when the palette order may not follow
  the factor-level order).

- discard.missing:

  Logical. If TRUE, animals without detections are omitted. Defaults to
  FALSE.

- shade:

  Controls the background shading. One of: `TRUE` (default) for
  annual-season shading (see
  [`shadeSeasons`](https://miguelgandra.github.io/moby/reference/shadeSeasons.md));
  `FALSE` for a plain `background.color`; or a data frame of custom
  periods with columns `start` and `end` (POSIXct or coercible) and,
  optionally, `label` (used for the legend and to assign a colour per
  category) and `color`.

- top.band:

  A date format ([`strftime`](https://rdrr.io/r/base/strptime.html)) for
  the labelled time band drawn above the panel (e.g. `"%Y"` for years).
  Set to FALSE to omit it. Defaults to `"%Y"`.

- legend:

  Logical. If TRUE (default), a legend (shading periods, release/expiry
  symbols and `color.by` levels) is drawn in the right margin.

- main:

  Optional plot title.

- bin:

  Optional aggregation of detections before plotting, useful for large
  datasets (reduces overplotting, file size and render time). One of:
  `NULL` (default, one point per detection); a unit string (`"hour"`,
  `"day"`, `"week"`, `"month"`); a number of minutes; or `"auto"` (a bin
  width chosen from the temporal span). Detections are aggregated by
  `id.col`, `color.by` (if set) and time bin, so one point represents
  all detections of that group within the interval.

- scale.by.count:

  Logical. When `bin` is set, if TRUE the point *area* scales with the
  number of detections in each bin (`pt.cex` is the reference size); if
  FALSE (default) all points share a constant size. Ignored when
  `bin = NULL`.

- date.interval:

  Controls the x-axis date labels. Either `"auto"` (default), which
  generates "pretty" breaks appropriate to the temporal span, or a
  numeric value selecting every n-th formatted date (manual mode; used
  together with `date.start` and `date.format`).

- date.format:

  A date format ([`strftime`](https://rdrr.io/r/base/strptime.html)) for
  the x-axis labels. If NULL (default), it is chosen automatically in
  `"auto"` mode and falls back to `"%b"` in manual mode.

- date.start:

  Integer controlling label placement. In manual mode (numeric
  `date.interval`), the index of the first displayed date. In `"auto"`
  mode, a phase/offset that re-anchors the automatically chosen labels
  within each period without changing the interval - e.g. for yearly
  labels it selects the month (1 = January ... 12 = December); for
  monthly labels, the day of month. Defaults to 1.

- pch:

  Plotting symbol for detections. Defaults to 16 (filled circle).

- pt.cex:

  Expansion factor for the detection points, or `"auto"` to scale the
  points to the available row height (keeping density consistent across
  datasets and devices). Defaults to 1.

- alpha:

  Opacity of the detection points, from 0 (transparent) to 1 (opaque).
  Defaults to 1.

- event.pch:

  Length-2 vector of symbols for the release and tag-expiry markers,
  ideally named `c(release = , expiry = )`. Defaults to
  `c(release = 8, expiry = 4)`.

- highlight.isolated:

  Logical. If TRUE and `color.by` is set, isolated detections are
  brought forward (plotted on top) so they are not hidden behind denser
  point clusters. Defaults to TRUE.

- background.color:

  Background colour used when `shade = FALSE` (or for uncovered gaps
  with custom periods). Defaults to "grey96".

- legend.cols:

  Integer number of columns for the `color.by` legend. If NULL
  (default), it is chosen automatically from the estimated legend height
  versus the device height.

- cex:

  Global expansion factor applied to all plot text (axis labels, titles,
  legend and the top band). Detection point size is controlled
  separately via `pt.cex`. Defaults to 1.

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

  Further graphical parameters passed to
  [`points`](https://rdrr.io/r/graphics/points.html) when drawing the
  detections (e.g. `lwd`, `bg`).

## Value

Called for its side effect: it generates an abacus plot.

## Details

Rows are laid out in the order of `id.groups` (and, within a single
group, in the factor-level order of `id.col`); to reorder animals,
reorder the `id.col` factor levels or the IDs within `id.groups`.

Date-axis labels are, by default, chosen automatically to suit the
temporal span of the data (from a few days to several years). Supplying
a numeric `date.interval` switches to manual mode, where every n-th
formatted date (per `date.format`, starting at `date.start`) is
labelled. When `highlight.isolated = TRUE` and a `color.by` variable is
set, isolated detections are drawn on top of denser point clusters so
they remain visible.

## Note

The layout is device-stable: decorative elements (legend, background
shading, the top band and the inter-legend spacing) are sized in
physical (inch) units and positioned with
[`grconvertX`](https://rdrr.io/r/graphics/convertXY.html) /
[`grconvertY`](https://rdrr.io/r/graphics/convertXY.html), while the
right margin is sized from the measured legend width. Only the data
panel stretches with the device, so the appearance stays consistent
across datasets and graphics devices of differing dimensions.

## Examples

``` r
# abacus plot of the bundled ray detections
# (release dates are read from the mobyData metadata)
plotAbacus(rays)
#> Warning: - 'id.col' converted to factor.
#> 
#> Abacus plot
#> ------------------------------------------------------
#>   Individuals:     8 (2 groups)
#>   Detections:      1,643
#>   Period:          2023-04-01 to 2023-06-30 (90 d)
#>   Points:          cex 1.00
#>   Date axis:       auto -> weekly ("%d %b")
#>   Top band:        "%Y"
#>   Shading:         seasonal
#>   Legend:          1 col
#> ------------------------------------------------------


# colour the detections by receiver station
plotAbacus(rays, color.by = "station")
#> Warning: - 'id.col' converted to factor.
#> Warning: - 'color.by' variable converted to factor.
#> 
#> Abacus plot
#> ------------------------------------------------------
#>   Individuals:     8 (2 groups)
#>   Detections:      1,643
#>   Period:          2023-04-01 to 2023-06-30 (90 d)
#>   Colour:          station (6 levels)
#>   Points:          cex 1.00
#>   Date axis:       auto -> weekly ("%d %b")
#>   Top band:        "%Y"
#>   Shading:         seasonal
#>   Legend:          1 col
#> ------------------------------------------------------
```
