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
  # craft a case where Dd > Di (detection after the recorded span start) cannot exceed 1 when capped
  md <- res_dataset()
  capped <- calculateResidency(md, tag.durations = c(A = 30, B = 30), residency.index = "IR1", cap = TRUE)
  expect_true(all(capped$IR1 <= 1, na.rm = TRUE))
  uncapped <- calculateResidency(md, tag.durations = c(A = 30, B = 30), residency.index = "IR1", cap = FALSE)
  expect_true(is.numeric(uncapped$IR1))
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
  expect_equal(ind[["IR1"]], sprintf("%.2f", res$IR1[match(ind$ID, res$ID)]))
  expect_equal(as.integer(ind[["N days detected"]]), res$days_detected[match(ind$ID, res$ID)])
})
