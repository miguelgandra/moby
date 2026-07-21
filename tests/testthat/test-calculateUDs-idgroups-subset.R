# Regression: calculateUDs() must accept id.groups AND subset together. The internal group column is
# created as "id.groups" but was read as "id.group" (singular), so interaction() saw NULL and the call
# died with "replacement has 0 rows". Note that a mobyData carrying id.groups in its metadata hits this
# even without an explicit id.groups argument.

test_that("calculateUDs runs with id.groups and subset combined (id.group column-name regression)", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(42)
  mk <- function(id, cx, cy, m, k = 40) data.frame(
    ID = id, timebin = as.POSIXct("2021-01-01", tz = "UTC") + seq_len(k) * 1800,
    x = rnorm(k, cx, 300), y = rnorm(k, cy, 300), month = m)
  coas <- rbind(mk("A", 0, 0, "Jan"),    mk("A", 0, 0, "Feb"),
                mk("B", 5000, 5000, "Jan"), mk("B", 5000, 5000, "Feb"),
                mk("C", 0, 5000, "Jan"),  mk("C", 0, 5000, "Feb"),
                mk("D", 5000, 0, "Jan"),  mk("D", 5000, 0, "Feb"))
  coas$ID <- factor(coas$ID)
  grid <- terra::rast(terra::ext(-3000, 8000, -3000, 8000), res = 250, crs = "EPSG:32629")

  expect_no_error(
    suppressWarnings(suppressMessages(calculateUDs(
      coas, method = "kde", bandwidth = 800, spatial.grid = grid, subset = "month",
      id.groups = list(g1 = c("A", "B"), g2 = c("C", "D")),
      lon.col = "x", lat.col = "y", epsg.code = 32629, verbose = FALSE)))
  )
})
