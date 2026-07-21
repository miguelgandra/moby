test_that("getTimeBins floors/ceils/rounds and preserves timezone", {
  x <- as.POSIXct("2024-11-15 08:23:45", tz = "UTC")
  expect_equal(getTimeBins(x, "30 mins", "floor"),
               as.POSIXct("2024-11-15 08:00:00", tz = "UTC"))
  expect_equal(getTimeBins(x, "30 mins", "ceiling"),
               as.POSIXct("2024-11-15 08:30:00", tz = "UTC"))
  expect_equal(attr(getTimeBins(x, "1 hour"), "tzone"), "UTC")
})

test_that("getSeason maps months correctly in both hemispheres", {
  jan <- as.POSIXct("2024-01-15", tz = "UTC")
  jul <- as.POSIXct("2024-07-15", tz = "UTC")
  expect_equal(as.character(getSeason(jan, "Northern")), "winter")
  expect_equal(as.character(getSeason(jul, "Northern")), "summer")
  expect_equal(as.character(getSeason(jan, "Southern")), "summer")
  expect_equal(as.character(getSeason(jul, "Southern")), "winter")
})

test_that("getReprodPeriod handles wrap-around periods and year-aware formats", {
  # November-February spawning: January is within the (wrapped) period
  expect_equal(
    as.character(getReprodPeriod(as.POSIXct("2024-01-15", tz = "UTC"), "11", "02", format = "%m")),
    "spawning"
  )
  # July is outside it
  expect_equal(
    as.character(getReprodPeriod(as.POSIXct("2024-07-15", tz = "UTC"), "11", "02", format = "%m")),
    "resting"
  )
  # POSIXct bounds with an explicit year-aware format and a time-of-day component
  expect_equal(
    as.character(getReprodPeriod(as.POSIXct("2024-05-30 06:30:00", tz = "UTC"),
                                 as.POSIXct("2024-04-01", tz = "UTC"),
                                 as.POSIXct("2024-09-30", tz = "UTC"),
                                 format = "%Y-%m-%d")),
    "spawning"
  )
})
