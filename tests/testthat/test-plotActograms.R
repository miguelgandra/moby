acto_dataset <- function() {
  set.seed(1)
  ids <- c("A", "B", "C", "D")
  d <- do.call(rbind, lapply(ids, function(id) {
    n <- 300
    dt <- as.POSIXct("2023-03-01", tz = "UTC") + sort(runif(n, 0, 200 * 86400))
    data.frame(ID = id, datetime = dt,
               station = sample(c("R1", "R2", "R3"), n, replace = TRUE),
               stringsAsFactors = FALSE)
  }))
  d$ID <- factor(d$ID, levels = ids)
  d$station <- factor(d$station)
  d
}

coords <- c(-8.9, 37.0)
tags <- rep(as.POSIXct("2023-03-01", tz = "UTC"), 4)

test_that("plotActograms runs across colour, groups, diel and date modes", {
  d <- acto_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 7)
  on.exit(dev.off(), add = TRUE)
  pa <- function(...) plotActograms(d, id.col = "ID", datetime.col = "datetime",
                                    tagging.dates = tags, coords = coords, ...)
  expect_no_error(suppressWarnings(suppressMessages(pa())))                                  # defaults
  expect_no_error(suppressWarnings(suppressMessages(pa(color.by = "station"))))              # colour + legend
  expect_no_error(suppressWarnings(suppressMessages(pa(diel.lines = 4))))                    # 4 diel lines
  expect_no_error(suppressWarnings(suppressMessages(pa(diel.lines = 0))))                    # diel off
  expect_no_error(suppressWarnings(suppressMessages(pa(tag.durations = 120))))               # tag expiry line
  expect_no_error(suppressWarnings(suppressMessages(pa(grid = TRUE, ncol = 1))))             # grid + single col
  expect_no_error(suppressWarnings(suppressMessages(pa(date.interval = 2, date.format = "%b/%y"))))  # manual axis
  expect_no_error(suppressWarnings(suppressMessages(
    pa(id.groups = list(g1 = c("A", "B"), g2 = c("C", "D")), color.by = "station", main = "Diel"))))
})

test_that("plotActograms requires coords only when diel lines are on", {
  d <- acto_dataset()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(
    plotActograms(d, id.col = "ID", datetime.col = "datetime", tagging.dates = tags, diel.lines = 0))))
  expect_error(
    plotActograms(d, id.col = "ID", datetime.col = "datetime", tagging.dates = tags, diel.lines = 4),
    "coords")
})

test_that("plotActograms validates inputs and returns invisibly", {
  d <- acto_dataset()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_error(
    plotActograms(d, id.col = "ID", datetime.col = "datetime", tagging.dates = tags,
                  diel.lines = 0, date.interval = "weekly"),
    "date.interval")
  ret <- suppressWarnings(suppressMessages(
    plotActograms(d, id.col = "ID", datetime.col = "datetime", tagging.dates = tags, diel.lines = 0)))
  expect_null(ret)
})

test_that("plotActograms restores graphical parameters on exit", {
  d <- acto_dataset()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  before <- par("mfrow")
  suppressWarnings(suppressMessages(
    plotActograms(d, id.col = "ID", datetime.col = "datetime", tagging.dates = tags, diel.lines = 0)))
  expect_equal(par("mfrow"), before)
})
