test_that("createWideTable returns a time-bin x individual table with metadata", {
  det <- make_detections(ids = c("A", "B"), n_per_id = 6)
  det$timebin <- getTimeBins(det$datetime, interval = "1 hour")
  w <- createWideTable(det, id.col = "ID", timebin.col = "timebin",
                       value.col = "detections", verbose = FALSE)
  expect_s3_class(w, "data.frame")
  expect_true(all(c("A", "B") %in% colnames(w)))
  expect_equal(attr(w, "ids"), c("A", "B"))
  expect_false(is.null(attr(w, "start.dates")))
})

test_that("createWideTable does not error when all bins are already present", {
  det <- data.frame(
    ID = c("A", "B"),
    timebin = as.POSIXct("2023-01-01 00:00:00", tz = "UTC") + c(0, 3600),
    station = "R1"
  )
  expect_no_error(
    createWideTable(det, id.col = "ID", timebin.col = "timebin",
                    value.col = "detections", verbose = FALSE)
  )
})

test_that("createWideTable preserves the input timezone", {
  det <- make_detections(ids = c("A", "B"), n_per_id = 6,
                         start = "2023-06-01 00:00:00", tz = "Australia/Sydney")
  det$timebin <- getTimeBins(det$datetime, interval = "1 hour")
  w <- createWideTable(det, id.col = "ID", timebin.col = "timebin",
                       value.col = "detections", verbose = FALSE)
  expect_equal(attr(w$timebin, "tzone"), "Australia/Sydney")
})

test_that("createWideTable errors clearly on a single time bin", {
  det <- data.frame(ID = "A",
                    timebin = as.POSIXct("2023-01-01 00:00:00", tz = "UTC"),
                    station = "R1")
  expect_error(
    createWideTable(det, id.col = "ID", timebin.col = "timebin",
                    value.col = "detections", verbose = FALSE),
    "time-bin interval"
  )
})
