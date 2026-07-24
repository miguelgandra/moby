test_that("calculateCOAs carries mobyData metadata forward so the pipeline stays metadata-driven", {
  coas <- suppressWarnings(suppressMessages(calculateCOAs(rays)))
  expect_true(is_moby(coas))
  expect_equal(mobyMeta(coas)$epsg.code, mobyMeta(rays)$epsg.code)     # CRS carried forward
  expect_false(is.null(mobyMeta(coas)$id.groups))                     # id.groups carried
  expect_false(is.null(mobyMeta(coas)$tagging.dates))                 # tagging dates carried
  # the payoff: downstream inherits the CRS -> calculateUDs needs no explicit epsg.code
  expect_no_error(suppressWarnings(suppressMessages(calculateUDs(coas, method = "kde", bandwidth = 250))))
})

test_that("calculateCOAs on a plain data frame returns a plain data frame", {
  plain <- as.data.frame(rays); attr(plain, "moby") <- NULL
  res <- suppressWarnings(suppressMessages(calculateCOAs(
    plain, id.col = "ID", timebin.col = "timebin", lon.col = "lon", lat.col = "lat", station.col = "station")))
  expect_false(is_moby(res))
  expect_true(all(c("ID", "timebin", "lon", "lat", "detections", "stations") %in% names(res)))
})

test_that("calculateCOAs tolerates an all-NA extra column (e.g. an empty sensor field)", {
  # position-only tags carry an empty sensor column; the dynamic per-column aggregation must not drop
  # every time bin (regression: aggregate()'s default na.omit -> "no rows to aggregate")
  plain <- as.data.frame(rays); attr(plain, "moby") <- NULL
  plain$sensor_value <- NA_real_          # all-NA numeric
  plain$sensor_unit  <- NA_character_     # all-NA character
  res <- suppressWarnings(suppressMessages(calculateCOAs(
    plain, id.col = "ID", timebin.col = "timebin", lon.col = "lon", lat.col = "lat", station.col = "station")))
  expect_gt(nrow(res), 0)
  expect_true(all(c("ID", "timebin", "lon", "lat", "detections", "sensor_value") %in% names(res)))
  expect_true(all(is.nan(res$sensor_value)))   # all-NA numeric collapses to NaN, bin retained
  expect_true(all(res$sensor_unit == ""))      # all-NA character concatenates to ""
})
