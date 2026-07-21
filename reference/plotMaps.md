# Plot home-range maps with movement trajectories

Draws one map per individual showing its kernel utilisation distribution
(UD home range) as graded utilisation isopleths, optionally over a
coastline and a background raster (e.g. bathymetry or temperature),
together with inferred movement tracks and time-ordered detection
points. Panels can be grouped by `id.groups`.

## Usage

``` r
plotMaps(
  data,
  uds = NULL,
  animal.tracks = NULL,
  id.groups = NULL,
  id.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  land.shape = NULL,
  coastline = getOption("moby.coastline", TRUE),
  epsg.code = NULL,
  land.color = "gray50",
  background.layer = NULL,
  background.color = "#F3F7F7",
  background.pal = NULL,
  background.title = NULL,
  discard.missing = TRUE,
  color.pal = NULL,
  ud.contour = 0.95,
  ud.uncertainty = c("band", "outline", "none"),
  tracks.color = "black",
  tracks.lty = 2,
  tracks.lwd = 0.2,
  plot.detections = TRUE,
  pts.color = c("#009E73", "white", "#D55E00"),
  pts.cex = c(1, 0.6, 1),
  ud.legend = TRUE,
  title.color = "black",
  title.pos = "topleft",
  title.inset = c(-0.08, 0),
  scale.km = NULL,
  scale.color = "black",
  scale.pos = "bottomright",
  scale.inset = 0.05,
  extent.factor = 1.1,
  main = NULL,
  cex = 1,
  ncol = NULL,
  file = NULL,
  width = NULL,
  height = NULL,
  res = 300
)
```

## Arguments

- data:

  A `mobyData` or data frame of animal positions (centres of activity)
  with lon/lat.

- uds:

  Utilisation distributions from
  [`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md) -
  either KDE (`estUD`/`estUDm`, drawn as graded volume isopleths) or
  AKDE (drawn from the stored isopleth contours, with a confidence
  envelope; see `ud.uncertainty`). Pass the whole
  [`calculateUDs()`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
  output; the method is detected automatically.

- animal.tracks:

  The distance-enriched output of
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
  (tracks are read from its `"trajectories"` attribute), or a
  pre-extracted list of per-individual linestring geometries.

- id.groups:

  Optional named list of ID groups (one block of panels each).

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- land.shape:

  Optional `sf` coastline/landmass polygon drawn over the background.

- coastline:

  When no `land.shape` is supplied, draw a default coastline for the map
  extent. `TRUE` auto-picks a scale; a string
  (`"small"`/`"medium"`/`"large"`) forces an `rnaturalearth` scale;
  `FALSE` draws no coastline. Defaults to
  `getOption("moby.coastline", TRUE)` (set that option to switch the
  fallback off package-wide). Requires the (Suggested) `rnaturalearth`
  (with `rnaturalearthdata`) or `maps` package; if neither is installed,
  no coastline is drawn.

- epsg.code:

  Projected CRS (numeric EPSG code or `crs` object, metre units) for the
  map. If NULL, inferred from `land.shape`, then from `uds`.

- land.color:

  Colour of land areas. Defaults to "gray50".

- background.layer:

  Optional projected raster displayed behind the home range.

- background.color:

  Solid background colour when no `background.layer` is supplied.
  Defaults to "#F3F7F7".

- background.pal:

  Colour palette for the `background.layer` raster. If NULL, a muted
  bathymetry palette is used.

- background.title:

  Legend title for the background layer. Defaults to "Layer values".

- discard.missing:

  Logical; drop individuals with fewer than 5 detections (UD needs \>=
  5). Defaults to TRUE.

- color.pal:

  Colour palette for the home-range isopleths. If NULL, viridis is used.

- ud.contour:

  Outer utilisation isopleth to display, as a proportion in (0, 1\].
  Cells outside this contour are made transparent. e.g. 0.95 shows the
  95% home range. Defaults to 0.95.

- ud.uncertainty:

  For AKDE UDs only: how to draw the confidence envelope of the outer
  isopleth. `"band"` (default) shades the region between the low and
  high contours; `"outline"` draws the low and high contours as dashed
  outlines; `"none"` draws only the estimate. Ignored for KDE (which has
  no CI).

- tracks.color, tracks.lty, tracks.lwd:

  Colour, line type and width of the movement trajectories.

- plot.detections:

  Logical; overlay detection points. Defaults to TRUE.

- pts.color:

  Colour(s) for detection points: one colour, or three for the first /
  intermediate / last detection. Defaults to c("#009E73", "white",
  "#D55E00").

- pts.cex:

  Size(s) for detection points: one value, or three for first /
  intermediate / last. Defaults to c(1, 0.6, 1).

- ud.legend:

  Logical; draw the home-range isopleth legend. Defaults to TRUE.

- title.color:

  Colour of the per-panel ID label. Defaults to "black".

- title.pos:

  Keyword position of the ID label (e.g. "topleft"). Defaults to
  "topleft".

- title.inset:

  Inset of the ID label from the plot edges (single value or x/y
  vector). Defaults to c(-0.08, 0).

- scale.km:

  Scale-bar length in kilometres. If NULL (default), ~20% of the map
  width.

- scale.color:

  Colour of the scale-bar labels. Defaults to "black".

- scale.pos:

  Keyword position of the scale bar. Defaults to "bottomright".

- scale.inset:

  Inset of the scale bar from the plot edges. Defaults to 0.05.

- extent.factor:

  Bounding-box expansion factor around the positions. Defaults to 1.1.

- main:

  Optional overall title above the panel grid.

- cex:

  Global expansion factor for all plot text (ID labels, legends, scale
  bar). Defaults to 1.

- ncol:

  Number of panel columns. If NULL, set from the number of individuals.

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

Invisibly, a tidy data frame with one row per mapped individual: the
isopleth level, the home-range area (km^2, when coordinates are metric)
and the number of detections.

## Details

The home range is drawn as true **utilisation isopleths**: the UD is
converted to a volume distribution with
[`getvolumeUD`](https://rdrr.io/pkg/adehabitatHR/man/kernelUD.html), so
`ud.contour = 0.95` really masks everything outside the 95% home range,
and the colour ramp and legend encode the isopleth level (core =
brightest). All layers (coordinates, coastline, background raster) are
reconciled to one projected CRS via the shared spatial pipeline before
plotting.

## See also

[`calculateUDs`](https://miguelgandra.github.io/moby/reference/calculateUDs.md),
[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md),
[`plotMovements`](https://miguelgandra.github.io/moby/reference/plotMovements.md)

## Examples

``` r
if (FALSE) { # \dontrun{
pdf("./home-range-maps.pdf", width = 12, height = 20)
plotMaps(data = animal_coas, uds = kud_output, animal.tracks = tracks,
         land.shape = coastline, id.col = "ID", lon.col = "lon", lat.col = "lat")
dev.off()
} # }
```
