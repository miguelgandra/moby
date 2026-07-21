test_that("summaryTable reports dates in the input timezone (UTC-shift regression)", {
  d <- data.frame(
    ID = "A",
    datetime = as.POSIXct("2020-06-01 09:00:00", tz = "Australia/Sydney") + (0:2) * 3600,
    station = "R1"
  )
  st <- suppressWarnings(suppressMessages(
    summaryTable(d, tagging.dates = as.POSIXct("2020-06-01 00:00:00", tz = "Australia/Sydney"),
                 residency.index = "IR1")
  ))
  expect_true(any(grepl("01/06/2020", st[["Last detection"]])))
})

test_that("summaryTable works with a single individual (single-column apply regression)", {
  d <- make_detections(ids = "A", n_per_id = 5)
  expect_no_error(suppressWarnings(suppressMessages(
    summaryTable(d, tagging.dates = as.POSIXct("2023-06-01 00:00:00", tz = "UTC"),
                 residency.index = "IR1")
  )))
})

test_that("summaryTable preserves a single metadata column (2-column id.metadata regression)", {
  # id.metadata[,-1] used to collapse to a vector when only one non-ID column remained, so aggregate()
  # named the merged column "x" instead of its real name (e.g. "species"), silently losing it.
  d <- make_detections(ids = c("A", "B"), n_per_id = 6)
  md <- data.frame(ID = c("A", "B"), species = c("tope", "tope"), stringsAsFactors = FALSE)
  st <- suppressWarnings(suppressMessages(
    summaryTable(d, tagging.dates = as.POSIXct("2023-06-01 00:00:00", tz = "UTC"),
                 residency.index = "IR1", id.metadata = md)
  ))
  expect_true("species" %in% colnames(st))
  expect_false("x" %in% colnames(st))
  expect_true(all(st[["species"]][st[["ID"]] %in% c("A", "B")] == "tope"))
})
