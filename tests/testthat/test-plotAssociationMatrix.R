# Reuse the association-input builder pattern (overlaps + null-model results)
assocmat_rr <- function(ids = c("A", "B", "C", "D", "E"), seed = 1) {
  set.seed(seed)
  n <- 90; tb <- as.POSIXct("2023-01-01", tz = "UTC") + (0:(n - 1)) * 3600
  df <- do.call(rbind, lapply(ids, function(id) data.frame(
    ID = id, timebin = tb, station = sample(c("R1", "R2", "R3"), n, replace = TRUE), stringsAsFactors = FALSE)))
  df$ID <- factor(df$ID)
  md <- as_moby(df, timebin.col = "timebin", station.col = "station", tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  wt <- suppressWarnings(suppressMessages(createWideTable(md, value.col = "station", verbose = FALSE)))
  ov <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  suppressWarnings(suppressMessages(randomizeAssociations(wt, ov, iterations = 80, random.seed = 1)))
}

test_that("plotAssociationMatrix draws significance and returns the id x id matrix", {
  rr <- assocmat_rr()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(plotAssociationMatrix(rr)))
  expect_true(is.matrix(res) || is.data.frame(res))
  # cells hold the recoded significance labels
  expect_true(all(unlist(res) %in% c("+", "-", "ns", NA)))
})

test_that("plotAssociationMatrix type='mean overlap' no longer crashes", {
  rr <- assocmat_rr()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  # previously a guaranteed crash: range() on a whole (mixed-type) data.frame
  expect_no_error(suppressWarnings(suppressMessages(plotAssociationMatrix(rr, type = "mean overlap"))))
  expect_no_error(suppressWarnings(suppressMessages(plotAssociationMatrix(rr, type = "mean overlap", full.scale = TRUE))))
})

test_that("plotAssociationMatrix sorts, restores par when file=NULL, and validates", {
  rr <- assocmat_rr()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 1)); mar0 <- par("mar")
  sizes <- stats::setNames(seq(50, 90, length.out = 5), as.character(attributes(rr)$ids))
  expect_no_error(suppressWarnings(suppressMessages(plotAssociationMatrix(rr, sort.by = sizes, sort.by.title = "Size"))))
  expect_equal(par("mar"), mar0)                              # par-leak fix (was nested in the file block)
  expect_error(plotAssociationMatrix(list(a = 1)), "randomizeAssociations")   # validation + correct fn name
  expect_error(plotAssociationMatrix(rr, color.pal = c("a", "b")), "3 colours")
})

test_that("plotAssociationMatrix accepts the factor and POSIXct sort.by classes it documents", {
  rr <- assocmat_rr()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  ids <- as.character(attributes(rr)$ids)

  # factor: previously errored ('round' not meaningful for factors) via .decimalPlaces()
  sizes <- factor(c("big", "sml", "big", "sml", "big"), levels = c("sml", "big"))
  res <- suppressWarnings(suppressMessages(
    plotAssociationMatrix(rr, sort.by = sizes, sort.by.title = "size")))
  expect_true(is.matrix(res) || is.data.frame(res))

  # dates: previously rendered as raw epoch seconds ("1672531200") rather than formatted dates
  tagged <- as.POSIXct("2023-01-01", tz = "UTC") + (0:4) * 86400
  expect_no_error(suppressWarnings(suppressMessages(
    plotAssociationMatrix(rr, sort.by = tagged, sort.by.title = "tagged"))))
})

test_that("plotAssociationMatrix labels sort.by by class, keeping the underlying sort order", {
  # the header row/column is drawn from these labels, so assert on them directly
  f <- factor(c("big", "sml", "big"), levels = c("sml", "big"))
  expect_equal(.sortByLabels(f, f[1:2]), list(labs1 = c("big", "sml", "big"), labs2 = c("big", "sml")))

  d <- as.POSIXct(c("2023-01-01", "2023-01-03"), tz = "UTC")
  expect_equal(.sortByLabels(d, d)$labs1, format(d))
  expect_false(any(grepl("^[0-9]+$", .sortByLabels(d, d)$labs1)))   # not bare epoch seconds

  # numeric: shared decimal places across both margins
  expect_equal(.sortByLabels(c(1, 2.25), c(3))$labs1, c("1.00", "2.25"))
  expect_equal(.sortByLabels(c(1, 2), c(3))$labs1, c("1", "2"))

  # order() drives the sort, so factors sort by level and dates chronologically
  expect_equal(order(f), c(2L, 1L, 3L))
  expect_equal(order(d), c(1L, 2L))
})

test_that("plotAssociationMatrix writes files without leaking a device", {
  rr <- assocmat_rr()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  before <- length(dev.list()); f <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(plotAssociationMatrix(rr, file = f))))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)
})
