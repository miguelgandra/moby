# Shared `file`/`width`/`height`/`res` export argument, across plotting functions.

fo_data <- function() {
  set.seed(1); n <- 2000
  t0 <- as.POSIXct("2023-01-01", tz = "UTC"); dt <- t0 + sort(runif(n, 0, 120 * 86400))
  hr <- as.numeric(format(dt, "%H"))
  data.frame(
    ID = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    datetime = dt,
    timebin = as.POSIXct(round(as.numeric(dt) / 3600) * 3600, origin = "1970-01-01", tz = "UTC"),
    depth = 15 - 8 * cos(2 * pi * hr / 24) + rnorm(n, 0, 2),
    station = factor(sample(c("S1", "S2", "S3"), n, replace = TRUE)),
    species = factor(sample(c("sp1", "sp2"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
}
fo_coords <- c(-8.9, 37.0)

# each entry: a one-arg function taking the output `file` path
fo_calls <- function(d) list(
  abacus       = function(f) plotAbacus(d, id.col = "ID", datetime.col = "datetime",
                                        tagging.dates = rep(min(d$datetime), 3), file = f),
  actograms    = function(f) plotActograms(d, id.col = "ID", datetime.col = "datetime",
                                           tagging.dates = rep(min(d$datetime), 3), coords = fo_coords, file = f),
  chronogram   = function(f) plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                                            coords = fo_coords, file = f),
  contours     = function(f) plotContours(d, variables = "depth", id.col = "ID", datetime.col = "datetime",
                                          coords = fo_coords, file = f),
  stationStats = function(f) plotStationStats(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                                              type = c("detections", "individuals"), file = f)
)

test_that("file = <path> writes a figure and never leaks or litters a device", {
  d <- fo_data()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)  # a device is open (interactive-like)
  before <- length(dev.list())
  owd <- setwd(tempdir()); on.exit(setwd(owd), add = TRUE)
  for(nm in names(fo_calls(d))){
    f <- tempfile(fileext = ".pdf")
    expect_no_error(suppressWarnings(suppressMessages(fo_calls(d)[[nm]](f))))
    expect_true(file.exists(f) && file.info(f)$size > 0, info = nm)
    expect_equal(length(dev.list()), before, info = paste(nm, "leaked a device"))
  }
  expect_false(file.exists(file.path(tempdir(), "Rplots.pdf")))
})

test_that("file = NULL leaves the active device untouched (draws to it)", {
  d <- fo_data()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  cur <- dev.cur()
  expect_no_error(suppressWarnings(suppressMessages(
    plotStationStats(d, id.col = "ID", timebin.col = "timebin", station.col = "station"))))
  expect_equal(dev.cur(), cur)                      # same device still current
})

# Can this platform's bare svg() device actually write a file? `capabilities("cairo")` is not a
# reliable answer: it reports what R was compiled with, not whether a cairo surface can be created
# at runtime. On a headless machine with no usable fonts the surface fails, and the device shuts
# itself straight back down WITHOUT signalling an R error.
#
# Hence the dev.cur() bookkeeping: if no new device appeared, calling dev.off() here would close the
# CALLER's device. Probing through grDevices only keeps the verdict independent of moby.
fo_svgWorks <- function() {
  f <- tempfile(fileext = ".svg")
  before <- as.integer(grDevices::dev.cur())
  on.exit({                                     # safety net: never leak, never close someone else's
    if(as.integer(grDevices::dev.cur()) != before && as.integer(grDevices::dev.cur()) != 1L)
      grDevices::dev.off()
  }, add = TRUE)
  isTRUE(tryCatch({
    grDevices::svg(f, width = 4, height = 3)
    if(as.integer(grDevices::dev.cur()) == before) FALSE   # never opened; nothing of ours to close
    else {
      plot(1:10)
      grDevices::dev.off()
      file.exists(f) && file.info(f)$size > 0
    }
  }, error = function(e) FALSE))
}

test_that("output format is inferred from the extension; bad extensions error", {
  d <- fo_data()
  # Decide up front, while no device of ours is open, so a broken backend cannot disturb the loop.
  exts <- c(".pdf", ".png", ".jpeg", ".tiff", ".svg")
  if(!fo_svgWorks()) exts <- setdiff(exts, ".svg")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  for(ext in exts){
    f <- tempfile(fileext = ext)
    expect_no_error(suppressWarnings(suppressMessages(
      plotAbacus(d, id.col = "ID", datetime.col = "datetime", tagging.dates = rep(min(d$datetime), 3),
                 file = f, width = 8, height = 5, res = 100))))
    expect_true(file.exists(f) && file.info(f)$size > 0,
                info = sprintf("%s (exists=%s, size=%s, cairo=%s)", ext, file.exists(f),
                               if(file.exists(f)) file.info(f)$size else NA, capabilities("cairo")))
  }
  expect_error(
    plotAbacus(d, id.col = "ID", datetime.col = "datetime", tagging.dates = rep(min(d$datetime), 3),
               file = tempfile(fileext = ".xyz")),
    "Unsupported output format")
  expect_error(
    plotAbacus(d, id.col = "ID", datetime.col = "datetime", tagging.dates = rep(min(d$datetime), 3),
               file = tempfile()),   # no extension
    "no extension")
})

test_that(".autoDim resolves user values and clamps structural defaults", {
  expect_equal(.autoDim(12, base = 2, slope = 1, n = 100, lo = 3, hi = 20)$value, 12)   # user value wins
  expect_false(.autoDim(12, base = 2, slope = 1, n = 100, lo = 3, hi = 20)$clamped)
  expect_equal(.autoDim(NULL, base = 2, slope = 0.5, n = 4, lo = 3, hi = 20)$value, 4)  # 2 + 0.5*4
  expect_true(.autoDim(NULL, base = 2, slope = 1, n = 100, lo = 3, hi = 20)$clamped)    # 102 -> clamp 20
  expect_equal(.autoDim(NULL, base = 2, slope = 1, n = 100, lo = 3, hi = 20)$value, 20)
})

test_that(".savePar does not open a device when none is open", {
  # in a fresh context with no device, .savePar must return NULL without creating one
  if(!is.null(dev.list())) for(i in rev(dev.list())) dev.off()
  expect_null(.savePar())
  expect_null(dev.list())
})

test_that(".openDevice errors instead of hijacking the caller's device when the device never opens", {
  # Regression: a cairo-backed device whose surface cannot be created shuts itself back down without
  # signalling an R error. .openDevice used to return successfully in that case, so the caller's
  # on.exit(dev.off()) drew into - and then closed - the USER's active device.
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  before <- dev.cur()
  local_mocked_bindings(svg = function(...) invisible(NULL), .package = "grDevices")
  expect_error(moby:::.openDevice(tempfile(fileext = ".svg"), width = 4, height = 3),
               "Could not open")
  expect_equal(dev.cur(), before)          # caller's device untouched
})
