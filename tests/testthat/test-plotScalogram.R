# Synthetic binned detections with an injected rhythm of known period (hours)
scal_data <- function(ids = c("A", "B"), pers = c(24, 12), days = 45, amp = 3, base = 3, seed = 1) {
  set.seed(seed)
  do.call(rbind, Map(function(id, per) {
    t <- seq(as.POSIXct("2020-01-01 00:00", tz = "UTC"), by = "1 hour", length.out = days * 24)
    lambda <- pmax(0, base + amp * sin(2 * pi * as.numeric(t) / 3600 / per) + rnorm(length(t), 0, 0.6))
    data.frame(ID = id, timebin = t, detections = rpois(length(t), lambda), stringsAsFactors = FALSE)
  }, ids, pers))
}

test_that("plotScalogram computes a CWT, recovers the dominant period, and returns tidy output", {
  skip_if_not_installed("wavScalogram")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(
    plotScalogram(scal_data(), variable = "detections", id.col = "ID", timebin.col = "timebin")))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("id", "dominant_period", "power") %in% names(res)))
  # global spectrum peaks near the injected period (tolerant: wavelet period resolution is coarse)
  expect_lt(abs(res$dominant_period[res$id == "A"] - 24), 3)
  expect_lt(abs(res$dominant_period[res$id == "B"] - 12), 2)
  expect_type(attr(res, "spectra"), "list")
  expect_type(attr(res, "periods"), "list")
})

test_that("plotScalogram options run: log/linear/quantile scaling, detrend, mask.coi, shared.scale", {
  skip_if_not_installed("wavScalogram")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- scal_data(ids = c("A", "B", "C"), pers = c(24, 12, 24))
  ps <- function(...) suppressWarnings(suppressMessages(
    plotScalogram(d, variable = "detections", id.col = "ID", timebin.col = "timebin", ...)))
  expect_no_error(ps(power.scaling = "linear"))
  expect_no_error(ps(power.scaling = "quantile"))
  expect_no_error(ps(power.scaling = "sqrt", detrend = "linear"))
  expect_no_error(ps(mask.coi = TRUE, shared.scale = TRUE))
  expect_no_error(ps(id.groups = list(g1 = c("A", "B"), g2 = "C"), ncol = 2))
  expect_no_error(ps(upper.quant = 0.99, time.unit = "hours"))
})

test_that("plotScalogram validates inputs and writes files without leaking a device", {
  skip_if_not_installed("wavScalogram")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- scal_data()
  d$chr <- as.character(d$detections)
  expect_error(plotScalogram(d, variable = "chr", id.col = "ID", timebin.col = "timebin"), "numeric")
  expect_error(plotScalogram(d, variable = "detections", id.col = "ID", timebin.col = "timebin",
                             period.range = 5), "length-2")
  before <- length(dev.list()); f <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(
    plotScalogram(d, variable = "detections", id.col = "ID", timebin.col = "timebin", file = f))))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)                    # no device leaked
})
