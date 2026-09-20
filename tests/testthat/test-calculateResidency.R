res_dataset <- function() {
  d <- data.frame(
    ID = factor(c(rep("A", 5), rep("B", 3))),
    datetime = as.POSIXct(c("2023-01-01 10:00", "2023-01-02 10:00", "2023-01-02 12:00",
                            "2023-01-05 10:00", "2023-01-10 10:00",
                            "2023-01-03 10:00", "2023-01-03 11:00", "2023-01-04 10:00"), tz = "UTC"),
    station = c("R1", "R1", "R2", "R1", "R2", "R3", "R3", "R3"))
  as_moby(d, tagging.dates = as.POSIXct(c(A = "2023-01-01", B = "2023-01-01"), tz = "UTC"))
}

test_that("calculateResidency returns a tidy numeric table with correct values", {
  md <- res_dataset()
  res <- calculateResidency(md, tag.durations = c(A = 30, B = 30),
                            residency.index = c("IR1", "IR2", "IWR"))
  expect_true(is.data.frame(res))
  expect_true(is.numeric(res$IR1))
  # A: 4 unique detection days; span Jan1->Jan10 = 10 days; Dt = 30
  expect_equal(res$days_detected[res$ID == "A"], 4)
  expect_equal(res$detection_span[res$ID == "A"], 10)
  expect_equal(res$IR1[res$ID == "A"], 0.4)
  expect_equal(res$IR2[res$ID == "A"], 4 / 30)
  expect_equal(res$IWR[res$ID == "A"], (4 / 30) * (10 / 30))
  expect_true(all(c("tagging_date", "first_detection", "last_detection",
                    "monitoring_end", "monitoring_duration") %in% colnames(res)))
})

test_that("calculateResidency cap argument bounds indices at 1", {
  md <- res_dataset()
  capped <- calculateResidency(md, tag.durations = c(A = 30, B = 30), residency.index = "IR1", cap = TRUE)
  expect_true(all(capped$IR1 <= 1, na.rm = TRUE))
  uncapped <- calculateResidency(md, tag.durations = c(A = 30, B = 30), residency.index = "IR1", cap = FALSE)
  expect_true(is.numeric(uncapped$IR1))
})

test_that("residency days are inclusive calendar dates in the detection time zone", {
  tz <- "America/New_York"
  d <- data.frame(ID = factor("A"), station = "R1",
                  datetime = as.POSIXct(c("2024-03-09 23:50:00", "2024-03-10 00:10:00"), tz = tz))
  md <- as_moby(d, tagging.dates = as.POSIXct("2024-03-09 23:00:00", tz = tz))
  res <- calculateResidency(md, residency.index = "IR1", cap = FALSE, verbose = FALSE)
  expect_equal(res$days_detected, 2)
  expect_equal(res$detection_span, 2)
  expect_equal(res$IR1, 1)

  # The spring DST transition is a 23-hour day, but still one calendar date.
  res_end <- calculateResidency(md, residency.index = c("IR1", "IR2"), cap = FALSE,
                                last.monitoring.date = as.POSIXct("2024-03-11 00:00:00", tz = tz),
                                verbose = FALSE)
  expect_equal(res_end$monitoring_duration, 2)
  expect_equal(res_end$IR2, 1)
})

test_that("start.point retains the release versus first-detection choice", {
  d <- data.frame(ID = factor("A"), station = "R1",
                  datetime = as.POSIXct(c("2023-01-03 11:00:00", "2023-01-04 12:00:00"), tz = "UTC"))
  md <- as_moby(d, tagging.dates = as.POSIXct("2023-01-01 12:00:00", tz = "UTC"))
  release <- calculateResidency(md, residency.index = "IR1", verbose = FALSE)
  first <- calculateResidency(md, residency.index = "IR1", start.point = "first.detection", verbose = FALSE)
  expect_equal(release$detection_span, 4)
  expect_equal(first$detection_span, 2)
  expect_equal(release$IR1, 0.5)
  expect_equal(first$IR1, 1)
})

test_that("unknown monitoring ends remain typed NA and display cleanly", {
  md <- res_dataset()
  res <- calculateResidency(md, residency.index = "IR1", verbose = FALSE)
  expect_s3_class(res$monitoring_end, "POSIXct")
  expect_true(all(is.na(res$monitoring_end)))
  expect_true(all(is.na(res$monitoring_duration)))
  expect_true(all(!is.na(res$IR1)))

  tab <- suppressWarnings(summaryTable(md, residency.index = "IR1", verbose = FALSE))
  expect_true(all(is.na(tab$monitoring_duration_d)))
  expect_true(all(format(tab)$monitoring_duration_d == "-"))
})

test_that("monitoring days count positive calendar-date overlap", {
  d <- data.frame(ID = factor("A"), station = "R1",
                  datetime = as.POSIXct("2023-01-01 13:00:00", tz = "UTC"))
  md <- as_moby(d, tagging.dates = as.POSIXct("2023-01-01 12:00:00", tz = "UTC"))
  midday <- calculateResidency(md, tag.durations = 1, residency.index = "IR2", verbose = FALSE)
  midnight <- calculateResidency(md, tag.durations = 0.5, residency.index = "IR2", verbose = FALSE)
  expect_equal(midday$days_detected, 1)
  expect_equal(midday$detection_span, 1)
  expect_equal(midday$monitoring_duration, 2)
  expect_equal(midnight$monitoring_duration, 1)
  expect_equal(midday$IR2, 0.5)
  expect_equal(midnight$IR2, 1)
})

test_that("the earlier monitoring cutoff is used for each individual", {
  md <- res_dataset()
  ends <- as.POSIXct(c(A = "2023-01-20", B = "2023-02-01"), tz = "UTC")
  res <- calculateResidency(md, tag.durations = c(A = 30, B = 30),
                            last.monitoring.date = ends, residency.index = "IR2", verbose = FALSE)
  expect_equal(res$monitoring_duration[res$ID == "A"], 19)
  expect_equal(res$monitoring_duration[res$ID == "B"], 30)
  expect_equal(as.Date(res$monitoring_end[res$ID == "B"]), as.Date("2023-01-31"))
})

test_that("Dt indices need endpoints and detections must be inside monitoring", {
  md <- res_dataset()
  expect_error(calculateResidency(md, residency.index = "IR2/IR1", verbose = FALSE),
               "monitoring duration")
  expect_error(calculateResidency(md, tag.durations = c(A = NA_real_, B = 30),
                                  residency.index = "IR2", verbose = FALSE),
               "No monitoring endpoint for ID\\(s\\): A")
  expect_error(summaryTable(md, residency.index = "IR2/IR1", verbose = FALSE),
               "monitoring duration|monitoring dates")

  d <- data.frame(ID = factor("A"), station = "R1",
                  datetime = as.POSIXct("2023-01-02 12:00:00", tz = "UTC"))
  valid <- as_moby(d, tagging.dates = as.POSIXct("2023-01-01 12:00:00", tz = "UTC"))
  expect_error(calculateResidency(valid, tag.durations = 1, residency.index = "IR1", verbose = FALSE),
               "outside the tagging-to-monitoring interval")
  expect_error(calculateResidency(valid, last.monitoring.date = as.POSIXct("2023-01-02", tz = "UTC"),
                                  residency.index = "IR1", verbose = FALSE),
               "outside the tagging-to-monitoring interval")

  too_early <- as_moby(d, tagging.dates = as.POSIXct("2023-01-03", tz = "UTC"))
  expect_error(calculateResidency(too_early, residency.index = "IR1", verbose = FALSE),
               "outside the tagging-to-monitoring interval")
})

test_that("calculateResidency computes partial residencies and validates inputs", {
  md <- res_dataset()
  resp <- calculateResidency(md, tag.durations = c(A = 30, B = 30),
                             residency.index = "IR1", residency.by = "station")
  expect_true(any(grepl("^IR1 R", colnames(resp))))
  # IR2/IWR require durations or last monitoring date
  expect_error(calculateResidency(md, residency.index = "IR2"), "tag.durations")
  expect_error(calculateResidency(md, tag.durations = c(A = 30, B = 30),
                                  residency.index = "BOGUS"), "residency.index")
})

test_that("summaryTable residency values match calculateResidency (single source of truth)", {
  md <- res_dataset()
  st <- suppressWarnings(suppressMessages(
    summaryTable(md, tag.durations = c(A = 30, B = 30), residency.index = c("IR1", "IR2"))))
  res <- suppressWarnings(suppressMessages(
    calculateResidency(md, tag.durations = c(A = 30, B = 30), residency.index = c("IR1", "IR2"))))
  ind <- st[st$ID %in% c("A", "B"), ]
  expect_equal(ind[["IR1"]], res$IR1[match(ind$ID, res$ID)])
  expect_equal(as.integer(ind[["n_days_detected"]]), res$days_detected[match(ind$ID, res$ID)])
})
