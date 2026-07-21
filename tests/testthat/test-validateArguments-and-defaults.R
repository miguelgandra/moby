test_that("as_moby validates types and stores metadata", {
  df <- make_detections(ids = c("A", "B"), n_per_id = 3)
  md <- as_moby(df)
  expect_true(is_moby(md))
  expect_s3_class(md, "data.frame")
  expect_equal(mobyMeta(md)$id.col, "ID")

  expect_error(as_moby(df, id.col = "nope"), "not found")
  expect_error(as_moby(df, tagging.dates = "2020-01-01"), "POSIXct")
  expect_error(as_moby(df, epsg.code = "abc"), "EPSG")
})

test_that("empty datasets are rejected with a clear message", {
  empty <- data.frame(ID = factor(character(0)),
                      datetime = as.POSIXct(character(0), tz = "UTC"))
  expect_error(
    suppressWarnings(calculateCOAs(empty, id.col = "ID", timebin.col = "datetime",
                                   lon.col = "lon", lat.col = "lat")),
    "empty|no rows|data frame"
  )
})

test_that(".dropCols never silently drops every column", {
  df <- data.frame(a = 1:2, b = 3:4)
  expect_equal(ncol(moby:::.dropCols(df, "missing")), 2L)
  expect_equal(colnames(moby:::.dropCols(df, "a")), "b")
})

test_that(".dataTZ falls back to UTC and reads the tzone attribute", {
  expect_equal(moby:::.dataTZ(as.POSIXct("2020-01-01", tz = "Australia/Sydney")),
               "Australia/Sydney")
  expect_equal(moby:::.dataTZ(1:3), "UTC")
})
