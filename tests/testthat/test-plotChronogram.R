chrono_dataset <- function() {
  set.seed(1)
  n <- 400
  t0 <- as.POSIXct("2023-04-01", tz = "UTC")
  dt <- t0 + sort(runif(n, 0, 60 * 86400))
  data.frame(
    ID = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    timebin = as.POSIXct(round(as.numeric(dt) / 3600) * 3600, origin = "1970-01-01", tz = "UTC"),
    station = factor(sample(c("S1", "S2", "S3"), n, replace = TRUE)),
    species = factor(sample(c("sp1", "sp2"), n, replace = TRUE)),
    depth = runif(n, 1, 30),
    stringsAsFactors = FALSE
  )
}

coords <- c(-8.9, 37.0)

test_that("plotChronogram runs across styles, metrics and overlays", {
  d <- chrono_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  cc <- function(...) plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                                     coords = coords, ...)
  expect_no_error(suppressWarnings(suppressMessages(cc())))
  expect_no_error(suppressWarnings(suppressMessages(cc(style = "raster"))))
  expect_no_error(suppressWarnings(suppressMessages(cc(color.by = "station", shade = "diel"))))
  expect_no_error(suppressWarnings(suppressMessages(cc(shade = "season"))))
  expect_no_error(suppressWarnings(suppressMessages(cc(metric = "co-occurrences"))))
  expect_no_error(suppressWarnings(suppressMessages(cc(metric = "individuals"))))
  expect_no_error(suppressWarnings(suppressMessages(cc(split.by = "species", color.by = "station"))))
  expect_no_error(suppressWarnings(suppressMessages(cc(color.by = "depth"))))          # numeric color.by
  expect_no_error(suppressWarnings(suppressMessages(cc(color.by = "station", legend = FALSE))))
})

test_that("plotChronogram does not require coords when diel is off", {
  d <- chrono_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(
    plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station", diel.lines = 0))))
  # but it errors clearly when diel lines are requested without coords
  expect_error(
    plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station", diel.lines = 4),
    "coords")
})

test_that("plotChronogram validates inputs and returns invisibly", {
  d <- chrono_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  expect_error(plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                              diel.lines = 0, metric = "bogus"), "metric")
  expect_error(plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                              diel.lines = 0, style = "x"), "style")
  ret <- suppressWarnings(suppressMessages(
    plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station", diel.lines = 0)))
  expect_null(ret)
})
