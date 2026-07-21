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
                                           tagging.dates = rep(min(d$datetime), 3), sunriset.coords = fo_coords, file = f),
  chronogram   = function(f) plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                                            sunriset.coords = fo_coords, file = f),
  contours     = function(f) plotContours(d, variables = "depth", id.col = "ID", datetime.col = "datetime",
                                          sunriset.coords = fo_coords, file = f),
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

# Does this platform's bare svg() device actually produce a file? `capabilities("cairo")` is not a
# reliable answer: a cairo surface that fails to initialise (e.g. no usable fonts on a headless CI
# runner) discards its output silently, so svg() opens and closes without error yet writes nothing.
# Probing with grDevices alone keeps the check independent of moby.
fo_svgWorks <- function() {
  f <- tempfile(fileext = ".svg")
  ok <- tryCatch({
    grDevices::svg(f, width = 4, height = 3)
    plot(1:10)
    grDevices::dev.off()
    file.exists(f) && file.info(f)$size > 0
  }, error = function(e) FALSE)
  isTRUE(ok)
}

test_that("output format is inferred from the extension; bad extensions error", {
  d <- fo_data()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  for(ext in c(".pdf", ".png", ".jpeg", ".tiff", ".svg")){
    f <- tempfile(fileext = ext)
    expect_no_error(suppressWarnings(suppressMessages(
      plotAbacus(d, id.col = "ID", datetime.col = "datetime", tagging.dates = rep(min(d$datetime), 3),
                 file = f, width = 8, height = 5, res = 100))))
    written <- file.exists(f) && file.info(f)$size > 0
    # If moby produced no SVG, decide whether that is moby's fault or the platform's: re-run the
    # same format through grDevices directly. Only excuse the failure if the bare device fails too.
    if(ext == ".svg" && !written && !fo_svgWorks()) next
    expect_true(written, info = sprintf("%s (exists=%s, size=%s, cairo=%s)", ext, file.exists(f),
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
