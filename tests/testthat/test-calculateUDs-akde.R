# AKDE relies on the (suggested) ctmm package and model fitting; keep small and guarded.

make_track <- function(id, n = 80, seed = 1) {
  set.seed(seed)
  t <- as.POSIXct("2023-01-01", tz = "UTC") + (1:n) * 3600
  x <- 5e5 + cumsum(rnorm(n, 0, 80))
  y <- 4.1e6 + cumsum(rnorm(n, 0, 80))
  data.frame(ID = id, timebin = t, lon = x, lat = y, stringsAsFactors = FALSE)
}

test_that("AKDE (default method) returns areas with confidence intervals", {
  skip_if_not_installed("ctmm")
  d <- rbind(make_track("A", seed = 1), make_track("B", seed = 2))
  d$ID <- factor(d$ID)
  md <- as_moby(d, timebin.col = "timebin", lon.col = "lon", lat.col = "lat", epsg.code = 32629)

  res <- suppressWarnings(calculateUDs(md, verbose = FALSE))
  expect_equal(attr(res, "method"), "akde")
  expect_true(all(c("ud", "summary_table", "area_estimates", "K50", "K95") %in% names(res)))
  expect_true(all(c("area_est", "area_low", "area_high", "DOF_area") %in% colnames(res$area_estimates)))
  expect_gt(nrow(res$area_estimates), 0)
  # confidence intervals must bracket the estimate
  expect_true(all(res$area_estimates$area_low <= res$area_estimates$area_est, na.rm = TRUE))
  expect_true(all(res$area_estimates$area_high >= res$area_estimates$area_est, na.rm = TRUE))
  # summary table stays movementTable-compatible (numeric, coercible area columns)
  expect_false(any(is.na(suppressWarnings(as.numeric(res$summary_table[["UD 95% (Km2)"]])))))
  expect_s3_class(res$K95, "sf")
  # the isopleth contours keep the confidence envelope (low/est/high), for plotMaps to draw
  expect_true("ci" %in% names(res$K95))
  expect_setequal(unique(res$K95$ci), c("low", "est", "high"))
  expect_equal(nrow(res$K95), 2L * 3L)                       # 2 individuals x (low, est, high)
})

test_that("AKDE skips individuals with too few locations", {
  skip_if_not_installed("ctmm")
  d <- rbind(make_track("A", n = 80, seed = 1), make_track("B", n = 3, seed = 9))
  d$ID <- factor(d$ID)
  md <- as_moby(d, timebin.col = "timebin", lon.col = "lon", lat.col = "lat", epsg.code = 32629)
  res <- suppressWarnings(calculateUDs(md, verbose = FALSE))
  expect_false("B" %in% res$area_estimates$ID)  # B (3 locations) dropped
})

test_that("method='kde' still works and requires a bandwidth", {
  skip_if_not_installed("adehabitatHR")
  d <- rbind(make_track("A", seed = 1), make_track("B", seed = 2))
  d$ID <- factor(d$ID)
  md <- as_moby(d, timebin.col = "timebin", lon.col = "lon", lat.col = "lat", epsg.code = 32629)

  expect_error(calculateUDs(md, method = "kde"), "bandwidth")
  res <- suppressWarnings(calculateUDs(md, method = "kde", bandwidth = 200, verbose = FALSE))
  expect_true("summary_table" %in% names(res))
})
