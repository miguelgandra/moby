ss_dataset <- function() {
  set.seed(1)
  n <- 600
  data.frame(
    ID = factor(sample(paste0("id", 1:8), n, replace = TRUE)),
    timebin = as.POSIXct("2023-01-01", tz = "UTC") + sample(0:300, n, replace = TRUE) * 3600,
    station = factor(sample(c("Reef-North", "Reef-South", "Channel", "Bay", "Offshore"), n, replace = TRUE)),
    habitat = NA,
    stringsAsFactors = FALSE
  ) -> d
  # one habitat per station (for aggregate.by)
  hab <- c("Reef-North" = "Reef", "Reef-South" = "Reef", "Channel" = "Soft", "Bay" = "Soft", "Offshore" = "Pelagic")
  d$habitat <- factor(unname(hab[as.character(d$station)]))
  d
}

grp <- list(sp1 = paste0("id", 1:4), sp2 = paste0("id", 5:8))

test_that("plotStationStats runs across types, scales, labels and groups", {
  d <- ss_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 7)
  on.exit(dev.off(), add = TRUE)
  ss <- function(...) plotStationStats(d, id.col = "ID", timebin.col = "timebin", station.col = "station", ...)
  expect_no_error(suppressWarnings(suppressMessages(ss())))                                   # default (detections)
  expect_no_error(suppressWarnings(suppressMessages(ss(type = "individuals"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(type = "average detections"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(type = "co-occurrences"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(type = c("detections", "individuals", "co-occurrences")))))
  expect_no_error(suppressWarnings(suppressMessages(ss(scale = "count"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(scale = "proportion"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(station.labels = "rotated"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(station.labels = "numbered"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(annotate = "both", annot.min = 5))))
  expect_no_error(suppressWarnings(suppressMessages(ss(annotate = "none"))))
  expect_no_error(suppressWarnings(suppressMessages(ss(aggregate.by = "habitat"))))
  expect_no_error(suppressWarnings(suppressMessages(
    ss(type = c("detections", "co-occurrences"), id.groups = grp, group.comparisons = "all", main = "By species"))))
  expect_no_error(suppressWarnings(suppressMessages(
    ss(id.groups = grp, group.comparisons = "within"))))
  expect_no_error(suppressWarnings(suppressMessages(
    ss(type = "co-occurrences", id.groups = grp, group.comparisons = "between"))))
})

test_that("plotStationStats handles a single station without error (regression)", {
  d <- ss_dataset(); d$station <- factor("OnlyStation")
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(
    plotStationStats(d, type = "detections", id.col = "ID", timebin.col = "timebin", station.col = "station"))))
  expect_no_error(suppressWarnings(suppressMessages(
    plotStationStats(d, type = c("detections", "co-occurrences"), id.col = "ID", timebin.col = "timebin", station.col = "station"))))
})

test_that("plotStationStats validates inputs and returns a tidy table invisibly", {
  d <- ss_dataset()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_error(
    plotStationStats(d, type = "bogus", id.col = "ID", timebin.col = "timebin", station.col = "station"),
    "Invalid 'type'")
  expect_error(
    plotStationStats(d, id.groups = grp, group.comparisons = "sideways",
                     id.col = "ID", timebin.col = "timebin", station.col = "station"),
    "group.comparisons")
  ret <- suppressWarnings(suppressMessages(
    plotStationStats(d, type = c("detections", "individuals"), id.col = "ID", timebin.col = "timebin", station.col = "station")))
  expect_s3_class(ret, "data.frame")
  expect_true(all(c("series", "type", "location", "count", "proportion") %in% names(ret)))
  expect_equal(nlevels(ret$location), 5L)
})

test_that("plotStationStats detection counts match a direct tabulation", {
  d <- ss_dataset()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  ret <- suppressWarnings(suppressMessages(
    plotStationStats(d, type = "detections", id.col = "ID", timebin.col = "timebin", station.col = "station")))
  det <- ret[ret$type == "detections", ]
  expected <- as.numeric(table(d$station)[as.character(det$location)])
  expect_equal(det$count, expected)
  expect_equal(sum(det$proportion), 1)
})
