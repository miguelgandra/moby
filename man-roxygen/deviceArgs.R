#' @param file Optional output file. If `NULL` (the default), the figure is drawn on the current
#' graphics device - the usual interactive behaviour. If a file path is supplied, moby opens a
#' graphics device chosen from the file extension (`.pdf`, `.svg`, `.png`, `.jpg`/`.jpeg`,
#' `.tif`/`.tiff` or `.bmp`), draws the figure to it, and closes the device automatically (also if
#' an error occurs). For multi-page or batch workflows, keep `file = NULL` and manage the device
#' yourself (e.g. `pdf(...); plot...(); dev.off()`).
#' @param width,height Output size in inches. Used *only* when `file` is supplied. If `NULL`
#' (default), a size is derived from the figure's structure (e.g. the number of individuals or
#' panels). These defaults are sensible starting points, *not* guarantees: dense or unusual figures
#' may still need you to set `width`/`height` explicitly. The two can be set independently.
#' @param res Resolution in pixels per inch, for raster formats only (`.png`, `.jpg`, `.tif`,
#' `.bmp`); ignored for vector formats (`.pdf`, `.svg`). Used only when `file` is supplied. Defaults
#' to 300.
