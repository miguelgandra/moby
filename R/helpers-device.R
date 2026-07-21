##################################################################################################
## Internal helpers: graphics-device handling for the `file` export argument #####################
## Shared machinery behind the `file`/`width`/`height`/`res` arguments common to every moby plot. #
##################################################################################################


## The shared `file`/`width`/`height`/`res` argument documentation lives in the roxygen template
## `man-roxygen/deviceArgs.R`; every plotting function pulls it in with `@template deviceArgs`, so
## the wording stays identical package-wide.


##################################################################################################
## Device-safe par() save/restore ################################################################

#' Capture par() without opening a device
#'
#' @description `par(no.readonly = TRUE)` silently opens the default device (e.g. `Rplots.pdf`) when
#' no device is open - a footgun in non-interactive/batch sessions. `.savePar()` returns the current
#' par only if a device is already open, and `NULL` otherwise, so a plotting function can capture par
#' for later restoration without ever creating a stray device. Pair with .restorePar().
#' @return A par list, or `NULL` if no device was open.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.savePar <- function() if(is.null(grDevices::dev.list())) NULL else graphics::par(no.readonly = TRUE)

#' Restore par() captured by .savePar()
#'
#' @description Restores the full par list (matching base R's usual `par(par(no.readonly=TRUE))`
#' round-trip), wrapped in `suppressWarnings` so the benign "calling par(new=TRUE) with no plot"
#' warning - which can fire when restoring after an error, before any plot was drawn - is not
#' surfaced to the user.
#' @param p A par list from .savePar(), or `NULL` (a no-op).
#' @return Invisibly `NULL`.
#' @keywords internal
#' @noRd
.restorePar <- function(p) { if(!is.null(p)) suppressWarnings(graphics::par(p)); invisible(NULL) }


##################################################################################################
## Resolve an auto/user dimension ################################################################

#' Resolve one output dimension (user value, or a clamped structural default)
#'
#' @description Returns the user-supplied dimension unchanged, or - when `NULL` - a value derived
#' from a simple `base + slope * n` rule clamped to `[lo, hi]`. Also reports whether the structural
#' default was clamped at the upper bound, which the caller uses to flag a possibly-crowded figure.
#' @param user The user-supplied dimension (inches), or `NULL` to auto-size.
#' @param base,slope,n The linear rule: `base + slope * n`, where `n` is a structural count.
#' @param lo,hi Lower/upper clamps (inches).
#' @return A list with `value` (inches) and `clamped` (TRUE if the structural default hit `hi`).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.autoDim <- function(user, base, slope, n, lo, hi){
  if(!is.null(user)) return(list(value = user, clamped = FALSE))
  raw <- base + slope * n
  list(value = max(lo, min(hi, raw)), clamped = raw > hi + 1e-9)
}


##################################################################################################
## Open a file-based graphics device from a path #################################################

#' Open a graphics device inferred from a file extension
#'
#' @description Opens the graphics device matching the extension of `file`, sized `width` x `height`
#' inches (raster formats at `res` ppi). The caller is responsible for closing it (typically via
#' `on.exit(grDevices::dev.off())`).
#' @param file Output file path (extension selects the format).
#' @param width,height Device size in inches (already resolved; non-NULL).
#' @param res Resolution (ppi) for raster formats.
#' @return Invisibly `NULL`; called for the side effect of opening a device.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.openDevice <- function(file, width, height, res = 300){
  ext <- tolower(tools::file_ext(file))
  if(!nzchar(ext))
    stop("Cannot determine the output format: 'file' has no extension (e.g. use \"plot.pdf\").", call. = FALSE)
  raster <- function(fun) fun(file, width = width, height = height, units = "in", res = res)
  before <- as.integer(grDevices::dev.cur())
  switch(ext,
    pdf  = grDevices::pdf(file, width = width, height = height),
    svg  = grDevices::svg(file, width = width, height = height),
    png  = raster(grDevices::png),
    jpg  = ,
    jpeg = raster(grDevices::jpeg),
    tif  = ,
    tiff = raster(grDevices::tiff),
    bmp  = raster(grDevices::bmp),
    stop(sprintf("Unsupported output format '.%s'. Supported: pdf, svg, png, jpg/jpeg, tif/tiff, bmp.", ext),
         call. = FALSE))
  # A device can fail to start WITHOUT signalling an R error - notably cairo-backed devices (svg, and
  # png/jpeg/tiff when type = "cairo") whose surface cannot be created, e.g. on a headless machine
  # with no usable fonts. The device then shuts itself straight back down, leaving the caller's
  # device current. Callers register on.exit(dev.off()) immediately after this returns, so failing to
  # detect it would silently draw into the user's active device and then CLOSE it. Fail loudly instead.
  if(as.integer(grDevices::dev.cur()) == before)
    stop(sprintf(paste0("Could not open a '%s' graphics device; the file was not written. This format ",
                        "may be unsupported by this R build (see capabilities()). Try a different ",
                        "extension, e.g. \"%s.png\"."),
                 ext, tools::file_path_sans_ext(basename(file))), call. = FALSE)
  invisible(NULL)
}


##################################################################################################
## Set up file output for a plotting function ####################################################

#' Resolve auto-sizing, warn on crowding, and open the output device
#'
#' @description Convenience wrapper used at the top of every plotting function's file-output block.
#' Given the user's `width`/`height` (possibly `NULL`) and the per-function structural sizing rule,
#' it resolves the dimensions, emits a single advisory `message()` if an auto-sized dimension was
#' clamped (a possibly-crowded figure), and opens the device. Returns nothing; the caller registers
#' `on.exit(grDevices::dev.off(), add = TRUE)` to guarantee the device is closed.
#' @param file,width,height,res As passed by the user.
#' @param w.rule,h.rule Lists `list(base=, slope=, n=, lo=, hi=)` giving the width/height rules.
#' @param crowd.unit Optional label (e.g. "individuals", "panels") naming the structural driver; when
#' supplied and an auto dimension is clamped, an advisory message is emitted.
#' @return Invisibly `NULL`.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.beginFileOutput <- function(file, width, height, res, w.rule, h.rule, crowd.unit = NULL){
  wd <- do.call(.autoDim, c(list(user = width),  w.rule))
  ht <- do.call(.autoDim, c(list(user = height), h.rule))
  if(!is.null(crowd.unit) && (wd$clamped || ht$clamped)){
    n <- if(ht$clamped) h.rule$n else w.rule$n
    message(sprintf(paste0("moby: %d %s exceed the comfortable default figure size; the saved figure may be ",
                           "crowded. Set 'width'/'height' explicitly for a larger canvas."), n, crowd.unit))
  }
  .openDevice(file, wd$value, ht$value, res)
  invisible(NULL)
}

##################################################################################################
##################################################################################################
##################################################################################################
