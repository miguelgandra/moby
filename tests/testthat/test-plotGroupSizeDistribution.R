gsd_dataset <- function() {
  set.seed(1)
  n <- 4000
  ids <- paste0("id", 1:12)
  tb <- as.POSIXct("2023-01-01", tz = "UTC") + sample(0:(24 * 40), n, replace = TRUE) * 3600
  sp <- setNames(rep(c("sp1", "sp2"), each = 6), ids)
  d <- data.frame(
    ID = factor(sample(ids, n, replace = TRUE), levels = ids),
    timebin = tb,
    station = factor(sample(c("S1", "S2", "S3"), n, replace = TRUE)),
    diel = factor(sample(c("day", "night"), n, replace = TRUE), levels = c("day", "night")),
    stringsAsFactors = FALSE
  )
  d$species <- factor(sp[as.character(d$ID)], levels = c("sp1", "sp2"))
  d
}
gsd_groups <- list(sp1 = paste0("id", 1:6), sp2 = paste0("id", 7:12))

test_that("plotGroupSizeDistribution runs across modes and returns a tidy table", {
  d <- gsd_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 8, height = 7); on.exit(dev.off(), add = TRUE)
  gg <- function(...) plotGroupSizeDistribution(d, id.col = "ID", timebin.col = "timebin",
                                                station.col = "station", ...)
  expect_no_error(suppressWarnings(suppressMessages(gg())))                                # whole study
  expect_no_error(suppressWarnings(suppressMessages(gg(split.by = "diel"))))               # stacked
  expect_no_error(suppressWarnings(suppressMessages(gg(id.groups = gsd_groups))))          # faceted
  expect_no_error(suppressWarnings(suppressMessages(
    gg(id.groups = gsd_groups, group.comparisons = "between", split.by = "diel", main = "x"))))
  expect_no_error(suppressWarnings(suppressMessages(gg(annotate = FALSE, legend = FALSE))))

  ret <- suppressWarnings(suppressMessages(gg(split.by = "diel")))
  expect_s3_class(ret, "data.frame")
  expect_true(all(c("series", "group_size", "split_level", "count", "freq") %in% names(ret)))
  expect_true(all(ret$group_size >= 2))                                                    # a group is 2+ animals
  # per-series frequencies sum to 100%
  agg <- tapply(ret$freq, ret$series, sum)
  expect_true(all(abs(agg[!is.na(agg)] - 100) < 1e-6))
})

test_that(".countCooccurrences: 'between' is per-cluster (regression for the old whole-row bug)", {
  # within-group cluster at S1 (2x sp1) while an unrelated sp2 animal sits alone at S2
  row <- c("S1", "S1", "S2"); grp <- c("sp1", "sp1", "sp2")
  expect_equal(as.numeric(.countCooccurrences(row, "within", grp, "sizes")), 2)            # the S1 pair
  expect_length(.countCooccurrences(row, "between", grp, "sizes"), 0)                       # NOT cross-group
  # a genuine cross-group cluster at S1 (sp1 + 2x sp2)
  row2 <- c("S1", "S1", "S1"); grp2 <- c("sp1", "sp2", "sp2")
  expect_equal(as.numeric(.countCooccurrences(row2, "between", grp2, "sizes")), 3)
  # 'stations' return mode (used by plotStationStats) still works
  expect_equal(.countCooccurrences(row, "within", grp, "stations"), "S1")
})

test_that("plotGroupSizeDistribution validates inputs and writes files without leaking a device", {
  d <- gsd_dataset()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_error(
    plotGroupSizeDistribution(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                              group.comparisons = "sideways"),
    "group.comparisons")
  expect_error(
    plotGroupSizeDistribution(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                              split.by = "nope"),
    "split.by")
  before <- length(dev.list())
  f <- tempfile(fileext = ".pdf")
  expect_no_error(suppressWarnings(suppressMessages(
    plotGroupSizeDistribution(d, id.col = "ID", timebin.col = "timebin", station.col = "station",
                              id.groups = gsd_groups, split.by = "diel", file = f))))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)          # no leaked device
})
