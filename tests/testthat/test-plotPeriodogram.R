# Synthetic binned detections with an injected rhythm of known period (hours)
rhythm_data <- function(ids = c("A", "B"), pers = c(24, 12), days = 40, amp = 3, base = 3, seed = 1) {
  set.seed(seed)
  do.call(rbind, Map(function(id, per) {
    t <- seq(as.POSIXct("2020-01-01 00:00", tz = "UTC"), by = "1 hour", length.out = days * 24)
    lambda <- pmax(0, base + amp * sin(2 * pi * as.numeric(t) / 3600 / per) + rnorm(length(t), 0, 0.5))
    data.frame(ID = id, timebin = t, detections = rpois(length(t), lambda), stringsAsFactors = FALSE)
  }, ids, pers))
}

test_that(".buildRhythmSeries represents gaps instead of deleting them (the core FFT bug)", {
  # 10 days hourly, then remove a 2-day block in the middle -> the gap must survive as filled bins,
  # so the series length still spans the whole observed range (not collapsed to contiguous samples).
  t <- seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "1 hour", length.out = 10 * 24)
  d <- data.frame(ID = "A", timebin = t, detections = 1L)
  d <- d[!(d$timebin >= as.POSIXct("2020-01-05", tz = "UTC") &
           d$timebin <  as.POSIXct("2020-01-07", tz = "UTC")), ]
  built <- .buildRhythmSeries(d, "ID", "timebin", "detections", gap.handling = "zero", min.days = NULL)
  expect_equal(built$dt, 60)                                  # 60-min sampling inferred
  # observed span is first->last bin = 10 days of hourly bins minus none at the ends -> ~240 samples
  expect_equal(length(built$series$A$values), length(t))     # gap represented, length preserved
  expect_true(any(built$series$A$values != 0))               # gap-filled zeros coexist with signal
})

test_that("plotPeriodogram recovers the injected period and returns a tidy table (fft)", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressMessages(plotPeriodogram(rhythm_data(), id.col = "ID", timebin.col = "timebin", min.days = 5))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("id", "dominant_period", "power", "fap", "band") %in% names(res)))
  expect_equal(round(res$dominant_period[res$id == "A"]), 24) # diel rhythm recovered
  expect_equal(round(res$dominant_period[res$id == "B"]), 12) # tidal rhythm recovered
  expect_equal(res$band, c("diel", "tidal"))
})

test_that("method = 'lomb' runs when available and returns a false-alarm probability", {
  skip_if_not_installed("lomb")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressMessages(plotPeriodogram(rhythm_data(), id.col = "ID", timebin.col = "timebin",
                                          method = "lomb", min.days = 5))
  expect_equal(round(res$dominant_period[res$id == "A"]), 24)
  expect_true(all(is.finite(res$fap)))                        # analytic FAP present for lomb
  expect_true(all(res$fap < 0.001))                           # strong injected signal
})

test_that("options run without error: presences, detrend, shared.scale, id.groups, bands off", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- rhythm_data(ids = c("A", "B", "C", "D"), pers = c(24, 12, 24, 8))
  pp <- function(...) suppressMessages(plotPeriodogram(d, id.col = "ID", timebin.col = "timebin", min.days = 5, ...))
  expect_no_error(pp(type = "presences"))
  expect_no_error(pp(detrend = "linear"))
  expect_no_error(pp(detrend = "diff", shared.scale = FALSE))
  expect_no_error(pp(highlight.bands = FALSE, period.range = c(4, 30)))
  expect_no_error(pp(id.groups = list(bay = c("A", "B"), reef = c("C", "D"))))
})

test_that("plotPeriodogram validates inputs and writes files without leaking a device", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- rhythm_data()
  expect_error(plotPeriodogram(d[, c("ID", "timebin")], id.col = "ID", timebin.col = "timebin"), "detections")
  expect_error(plotPeriodogram(d, id.col = "ID", timebin.col = "timebin", period.range = 5), "length-2")
  before <- length(dev.list()); f <- tempfile(fileext = ".png")
  expect_no_error(suppressMessages(plotPeriodogram(d, id.col = "ID", timebin.col = "timebin", min.days = 5, file = f)))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)                    # no device leaked
})
