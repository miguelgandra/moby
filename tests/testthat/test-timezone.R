# Timezone-regression tests. The tz cluster (day/hour derived in UTC instead of the data's own zone)
# is PROVABLY LATENT on the UTC-only `rays` data, so these use a non-UTC fixture (nonutc_detections()).

test_that(".dataTZ reports the data's own timezone, not UTC", {
  d <- nonutc_detections()
  expect_identical(moby:::.dataTZ(d$datetime), "Etc/GMT-10")
  expect_identical(moby:::.dataTZ(as.POSIXct("2023-01-01", tz = "UTC")), "UTC")
})

test_that("filterDetections min.days buckets days in the DATA timezone (not UTC)", {
  d <- nonutc_detections()   # A: 1 local day (2 UTC days); B: 2 local days
  tags <- stats::setNames(rep(as.POSIXct("2023-06-01", tz = "Etc/GMT-10"), 2), c("A", "B"))
  res <- suppressWarnings(suppressMessages(
    filterDetections(d, tagging.dates = tags, id.col = "ID", datetime.col = "datetime",
                     lon.col = "lon", lat.col = "lat",
                     min.days = 2)))
  kept <- unique(as.character(res$data$ID))
  # A spans 1 local day -> dropped by min.days=2; B spans 2 local days -> kept.
  # Under the old UTC-bucketing bug A would count 2 UTC days and be wrongly kept.
  expect_false("A" %in% kept)
  expect_true("B" %in% kept)
})

test_that("filterDetections re-wraps its output as mobyData (metadata preserved)", {
  d <- nonutc_detections()
  tags <- stats::setNames(rep(as.POSIXct("2023-06-01", tz = "Etc/GMT-10"), 2), c("A", "B"))
  md <- as_moby(d, id.col = "ID", datetime.col = "datetime", lon.col = "lon", lat.col = "lat",
                tagging.dates = tags)
  res <- suppressWarnings(suppressMessages(
    filterDetections(md, min.days = 0)))
  expect_true(is_moby(res$data))
  expect_equal(mobyMeta(res$data)$datetime.col, "datetime")
})

test_that("chronogram / abacus render on non-UTC data without error", {
  # a longer non-UTC series (3 animals over ~40 days) so the plots have enough span for their defaults
  tz <- "Etc/GMT-10"
  times <- as.POSIXct("2023-06-01 08:00", tz = tz) + sort(sample(0:(40 * 86400), 90))
  d <- data.frame(ID = factor(sample(c("A", "B", "C"), 90, replace = TRUE)), datetime = times,
                  station = sample(c("S1", "S2"), 90, replace = TRUE), stringsAsFactors = FALSE)
  d$timebin <- getTimeBins(d$datetime, "1 hour")
  tags <- stats::setNames(rep(as.POSIXct("2023-06-01", tz = tz), 3), c("A", "B", "C"))
  tmp <- tempfile(fileext = ".pdf"); grDevices::pdf(tmp); on.exit({ grDevices::dev.off(); unlink(tmp) })
  expect_no_error(suppressWarnings(suppressMessages(
    plotChronogram(d, id.col = "ID", timebin.col = "timebin", station.col = "station", diel.lines = 0))))
  expect_no_error(suppressWarnings(suppressMessages(
    plotAbacus(d, id.col = "ID", datetime.col = "datetime", station.col = "station",
               tagging.dates = tags, shade = FALSE))))
})
