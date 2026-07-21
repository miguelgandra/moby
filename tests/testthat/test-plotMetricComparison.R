mc_data <- function(levs = c("day", "night"), effect = 4, n_ids = 16, seed = 1) {
  set.seed(seed)
  ids <- paste0("id", seq_len(n_ids))
  do.call(rbind, lapply(ids, function(id) {
    base <- rnorm(1, 10, 3)
    data.frame(ID = id, phase = factor(levs, levels = levs),
               home_range = pmax(0.1, base + effect * seq(0, length(levs) - 1) + rnorm(length(levs), 0, 1)),
               activity   = pmax(0, rnorm(length(levs), 0.5, 0.15)),
               stringsAsFactors = FALSE)
  }))
}

test_that("plotMetricComparison runs across levels, plot types and options; returns tidy stats", {
  pdf(tempfile(fileext = ".pdf"), width = 8, height = 6); on.exit(dev.off(), add = TRUE)
  d2 <- mc_data(c("day", "night")); d3 <- mc_data(c("spring", "summer", "autumn"))
  mc <- function(d, ...) plotMetricComparison(d, id.col = "ID", split.by = "phase", ...)
  expect_no_error(suppressWarnings(suppressMessages(mc(d2, metrics = c("home_range", "activity")))))
  expect_no_error(suppressWarnings(suppressMessages(mc(d3, metrics = "home_range"))))            # 3 levels (old crash)
  expect_no_error(suppressWarnings(suppressMessages(mc(d3, metrics = "home_range", plot.type = "violin"))))
  expect_no_error(suppressWarnings(suppressMessages(mc(d3, metrics = "home_range", plot.type = "points", display.n = TRUE))))
  expect_no_error(suppressWarnings(suppressMessages(mc(d2, metrics = "home_range", test = "none"))))
  expect_no_error(suppressWarnings(suppressMessages(mc(d2, metrics = "home_range", paired = FALSE))))
  expect_no_error(suppressWarnings(suppressMessages(mc(d3, metrics = "home_range", flag.significant = TRUE))))

  ret <- suppressWarnings(suppressMessages(mc(d2, metrics = c("home_range", "activity"))))
  expect_s3_class(ret, "data.frame")
  expect_true(all(c("metric", "test", "paired", "n_complete", "statistic", "p_raw", "p_adj", "effect", "effect_type") %in% names(ret)))
  expect_s3_class(attr(ret, "values"), "data.frame")
})

test_that("design-appropriate tests are chosen (repeated measures, not independent groups)", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  r2 <- suppressWarnings(suppressMessages(
    plotMetricComparison(mc_data(c("day", "night")), id.col = "ID", split.by = "phase", metrics = "home_range")))
  expect_equal(r2$test[1], "Wilcoxon signed-rank")          # paired, 2 levels
  r3 <- suppressWarnings(suppressMessages(
    plotMetricComparison(mc_data(c("a", "b", "c")), id.col = "ID", split.by = "phase", metrics = "home_range")))
  expect_equal(r3$test[1], "Friedman")                       # paired, >2 levels (NOT Kruskal-Wallis)
  r2u <- suppressWarnings(suppressMessages(
    plotMetricComparison(mc_data(c("day", "night")), id.col = "ID", split.by = "phase", metrics = "home_range", paired = FALSE)))
  expect_equal(r2u$test[1], "Mann-Whitney")                  # only when paired = FALSE
})

test_that(".rankBiserialPaired has the correct sign and range (footgun pin)", {
  x <- c(1, 2, 3, 4, 5); y <- x + 2                          # y always > x
  expect_equal(.rankBiserialPaired(x, y), 1)                 # +1
  expect_equal(.rankBiserialPaired(y, x), -1)                # symmetric -> -1
  expect_equal(.rankBiserialPaired(x, x), NA_real_)          # all ties -> undefined
  # a genuine paired difference recovers direction
  set.seed(2); a <- rnorm(30); b <- a + rnorm(30, 0.6, 0.3)
  expect_gt(.rankBiserialPaired(a, b), 0)                    # b tends > a
})

test_that("the test is suppressed below min.n complete blocks", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- mc_data(c("day", "night"), n_ids = 4)                # only 4 pairs
  r <- suppressWarnings(suppressMessages(
    plotMetricComparison(d, id.col = "ID", split.by = "phase", metrics = "home_range", min.n = 6)))
  expect_true(is.na(r$p_raw[1]))
  expect_match(r$note[1], "min.n")
})

test_that("plotMetricComparison validates inputs and writes files without leaking a device", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- mc_data()
  expect_error(plotMetricComparison(d, id.col = "ID", split.by = "phase", metrics = "nope"), "not found")
  expect_error(plotMetricComparison(d, id.col = "ID", split.by = "phase", metrics = "home_range",
                                    p.adjust.method = "holm") -> ok, NA)
  # single-level grouping errors clearly
  d1 <- d; d1$phase <- factor("all")
  expect_error(suppressWarnings(plotMetricComparison(d1, id.col = "ID", split.by = "phase", metrics = "home_range")),
               "at least 2 levels")
  before <- length(dev.list()); f <- tempfile(fileext = ".pdf")
  expect_no_error(suppressWarnings(suppressMessages(
    plotMetricComparison(d, id.col = "ID", split.by = "phase", metrics = c("home_range", "activity"), file = f))))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)
})
