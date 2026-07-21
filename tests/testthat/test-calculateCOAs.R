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
