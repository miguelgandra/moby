test_that("mobyData is a data.frame subclass carrying metadata", {
  df <- make_detections(ids = c("A", "B"), n_per_id = 4)
  md <- as_moby(df, epsg.code = 32629,
                tagging.dates = as.POSIXct("2023-06-01", tz = "UTC"),
                id.groups = list(g1 = c("A", "B")))
  expect_s3_class(md, c("mobyData", "data.frame"))
  expect_true(is.data.frame(md))
  m <- mobyMeta(md)
  expect_equal(m$epsg.code, 32629)
  expect_s3_class(m$tagging.dates, "POSIXct")
  expect_named(m$id.groups, "g1")
})

test_that("re-wrapping a mobyData updates only supplied fields and inherits the rest", {
  df <- make_detections(ids = c("A", "B"), n_per_id = 3)
  md <- as_moby(df, id.col = "ID", epsg.code = 32629)
  md2 <- as_moby(md, tagging.dates = as.POSIXct("2023-06-01", tz = "UTC"))
  expect_equal(mobyMeta(md2)$epsg.code, 32629)            # inherited
  expect_s3_class(mobyMeta(md2)$tagging.dates, "POSIXct") # updated
})

test_that("subsetting preserves the mobyData class and metadata", {
  df <- make_detections(ids = c("A", "B"), n_per_id = 3)
  md <- as_moby(df)
  sub <- md[md$ID == "A", ]
  expect_true(is_moby(sub))
  expect_equal(mobyMeta(sub)$id.col, "ID")
})

test_that("functions resolve column names from mobyData metadata without explicit args", {
  df <- data.frame(
    animal = factor(rep(c("A", "B"), each = 3)),
    tbin = as.POSIXct("2023-01-01", tz = "UTC") + rep(c(0, 3600, 7200), 2),
    x = -8, y = 37, rec = "R1"
  )
  md <- as_moby(df, id.col = "animal", timebin.col = "tbin",
                lon.col = "x", lat.col = "y", station.col = "rec")
  coas <- suppressWarnings(calculateCOAs(md))
  expect_true("animal" %in% colnames(coas))
  expect_true(nrow(coas) > 0)
})

test_that("plain data frames with canonical column names still work", {
  df <- data.frame(
    ID = factor(rep(c("A", "B"), each = 3)),
    timebin = as.POSIXct("2023-01-01", tz = "UTC") + rep(c(0, 3600, 7200), 2),
    lon = -8, lat = 37, station = "R1"
  )
  coas <- suppressWarnings(calculateCOAs(df))
  expect_true(nrow(coas) > 0)
})

test_that("as_moby stores, validates and inherits nominal.delay", {
  d <- as.data.frame(rays)
  # a single value applies to all individuals; a named vector keys per animal (mixed tag families)
  expect_identical(mobyMeta(as_moby(d, nominal.delay = 120))$nominal.delay, 120)
  expect_identical(mobyMeta(as_moby(d, nominal.delay = c(R01 = 60, R02 = 120)))$nominal.delay,
                   c(R01 = 60, R02 = 120))
  # inherited on re-wrap, like the other metadata fields
  expect_identical(mobyMeta(as_moby(as_moby(d, nominal.delay = 120)))$nominal.delay, 120)
  # must be a positive number of seconds
  expect_error(as_moby(d, nominal.delay = -5), "nominal.delay")
  expect_error(as_moby(d, nominal.delay = "fast"), "nominal.delay")
})
