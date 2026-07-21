# plotMinLag(): the min_lag threshold diagnostic. The key guarantee is that it computes min_lag from
# the SAME internal helper as filterDetections(), so its flagged counts match what the filter removes.

# dense bursts (small min_lag -> kept) + scattered isolated decodes (large min_lag -> flagged)
mk_minlag_data <- function(seed = 1) {
  set.seed(seed)
  t0 <- as.POSIXct("2023-01-01", tz = "UTC")
  one <- function(id) rbind(
    data.frame(ID = id, datetime = t0 + cumsum(runif(40, 80, 160)), station = "R1"),   # burst at R1
    data.frame(ID = id, datetime = t0 + sort(runif(4, 5 * 86400, 60 * 86400)), station = "R2"))  # isolated at R2
  d <- do.call(rbind, lapply(c("A", "B", "C"), one))
  d$ID <- factor(d$ID)
  # lon/lat are unused by min_lag but filterDetections() validates them (for the consistency test)
  d$lon <- -8; d$lat <- 37
  d[order(d$ID, d$datetime), ]
}

test_that("plotMinLag flagged counts equal filterDetections' min_lag removals (shared helper)", {
  d <- mk_minlag_data()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(plotMinLag(d, nominal.delay = 120, factors = 30)))
  # same data, no pre-tagging / cut-off removal -> the filter's min_lag stage sees the same detections
  filt <- suppressWarnings(suppressMessages(filterDetections(
    d, tagging.dates = as.POSIXct("2022-01-01", tz = "UTC"), nominal.delay = 120,
    min.lag.factor = 30, verbose = FALSE)))
  n_false <- sum(grepl("false detection", filt$data_discarded$reason))
  expect_gt(n_false, 0)                                       # the case actually exercises the filter
  expect_equal(res$n_flagged[res$factor == 30], n_false)
})

test_that("plotMinLag returns a monotonic flagged table and requires nominal.delay", {
  d <- mk_minlag_data(2)
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(plotMinLag(d, nominal.delay = 120, factors = c(20, 30, 50))))
  expect_s3_class(res, "data.frame")
  expect_named(res, c("factor", "threshold_s", "n_flagged", "pct_flagged"))
  expect_true(all(diff(res$n_flagged) <= 0))                 # higher threshold -> no more flagged
  expect_equal(res$threshold_s, c(20, 30, 50) * 120)          # single nominal -> concrete thresholds
  expect_error(suppressWarnings(plotMinLag(d)), "nominal.delay")
})

test_that("plotMinLag writes to a file and reports NA thresholds for mixed tag families", {
  d <- mk_minlag_data(3)
  f <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(
    plotMinLag(d, nominal.delay = c(A = 60, B = 120, C = 120), file = f))))
  expect_true(file.exists(f))
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(plotMinLag(d, nominal.delay = c(A = 60, B = 120, C = 120))))
  expect_true(all(is.na(res$threshold_s)))                    # nominal delay varies across IDs
})
