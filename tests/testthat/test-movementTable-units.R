# Regression: movementTable() must combine id.groups that move at different speeds. The rate-of-
# movement unit (m/h vs km/h) was chosen per group and baked into the column name, so a fast group and
# a slow group produced "ROM (km/h)" vs "ROM (m/h)" and the final rbind() died with
# "names do not match previous names". The unit is now decided once for the whole table.

test_that("movementTable combines fast and slow id.groups (per-group rom_units rbind regression)", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7)
  mkmove <- function(id, step, k = 40) {
    ang <- cumsum(rnorm(k, 0, 0.3))
    data.frame(ID = id, timebin = as.POSIXct("2021-01-01", tz = "UTC") + seq_len(k) * 1800,
               x = cumsum(step * cos(ang)), y = cumsum(step * sin(ang)))
  }
  # fast movers ~4000 m/h -> km/h; slow movers ~75 m/h -> m/h (under the old per-group logic)
  mv <- rbind(mkmove("F1", 2000), mkmove("F2", 2100), mkmove("S1", 40), mkmove("S2", 35))
  mv$ID <- factor(mv$ID)

  tracks <- suppressWarnings(suppressMessages(
    calculateStepDistances(mv, lon.col = "x", lat.col = "y", epsg.code = 32629, verbose = FALSE)))
  rng <- range(c(tracks$x, tracks$y), na.rm = TRUE)
  grid <- terra::rast(terra::ext(rng[1] - 10000, rng[2] + 10000, rng[1] - 10000, rng[2] + 10000),
                      res = 1000, crs = "EPSG:32629")
  uds <- suppressWarnings(suppressMessages(
    calculateUDs(tracks, method = "kde", bandwidth = 3000, spatial.grid = grid,
                 lon.col = "x", lat.col = "y", epsg.code = 32629, verbose = FALSE)))

  mt <- suppressWarnings(suppressMessages(movementTable(
    tracks, uds, id.groups = list(fast = c("F1", "F2"), slow = c("S1", "S2")),
    epsg.code = 32629, lon.col = "x", lat.col = "y")))

  expect_s3_class(mt, "data.frame")
  # a single, shared rate-of-movement unit is used across the combined table
  expect_length(grep("^ROM \\(", colnames(mt)), 1)
  expect_length(grep("^Max ROM \\(", colnames(mt)), 1)
})
