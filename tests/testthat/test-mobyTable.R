# The typed-table contract shared by summaryTable(), movementTable() and transitionsTable():
# the object carries NUMBERS, format() renders them, print() displays them. Before this, the tables
# came back as character with the mean row and the group headings spliced in as data rows.

sum_tbl <- function(...) {
  summaryTable(rays, last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"),
               verbose = FALSE, ...)
}
# `rays` carries id.groups in its metadata, so summaryTable() groups by species unless told
# otherwise; several tests want the plain ungrouped table.
flat <- function(f) f[!grepl("^mean", f[[1]]), , drop = FALSE]
ray_ids <- function() unique(as.character(rays$ID))

test_that("summaryTable returns a typed table you can compute on", {
  tbl <- sum_tbl()
  expect_s3_class(tbl, "mobyTable")
  expect_s3_class(tbl, "data.frame")
  expect_type(tbl$n_detections, "integer")
  expect_true(is.numeric(tbl$IR1))
  expect_s3_class(tbl$tagging_date, "POSIXct")
  expect_s3_class(tbl$last_detection, "POSIXct")
  # the point of all of it
  expect_false(is.na(mean(tbl$n_detections)))
  # no mean row and no group heading rows hiding in the data
  expect_false(any(grepl("^mean", tbl$ID)))
  expect_equal(nrow(tbl), length(unique(tbl$ID)))
})

test_that("columns are stable snake_case names", {
  tbl <- sum_tbl()
  expect_true(all(c("tagging_date", "last_detection", "n_detections", "n_receivers",
                    "monitoring_duration_d", "detection_span_d", "n_days_detected") %in% names(tbl)))
  # the residency indices keep the names calculateResidency() gives them
  expect_true("IR1" %in% names(tbl))
})

test_that("format() returns a rectangular character table fit for write.csv()", {
  tbl <- sum_tbl(id.groups = mobyMeta(rays)$id.groups)
  f <- format(tbl)
  expect_s3_class(f, "data.frame")
  expect_true(all(vapply(f, is.character, logical(1))))
  # no blank separator records: every row carries a label
  expect_false(any(f[[1]] == ""))
  expect_equal(ncol(f), ncol(tbl))
  # missing values render as "-", never NA
  expect_false(anyNA(f))
})

test_that("format() is ASCII by default in every style, and unicode on request", {
  tbl <- sum_tbl()
  for (st in c("internal", "report", "concise")) {
    f <- format(tbl, style = st)
    txt <- c(names(f), unlist(f, use.names = FALSE))
    expect_false(any(grepl("[^\x01-\x7f]", txt)), info = st)
    expect_true(any(grepl("+/-", txt, fixed = TRUE)), info = st)
  }
  fu <- format(tbl, symbols = "unicode")
  expect_true(any(grepl("\u00b1", unlist(fu, use.names = FALSE), fixed = TRUE)))
})

test_that("style changes the headers and nothing else", {
  tbl <- sum_tbl()
  a <- format(tbl, style = "internal"); b <- format(tbl, style = "report")
  c3 <- format(tbl, style = "concise")
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)))
  expect_equal(unname(as.matrix(a)), unname(as.matrix(c3)))
  expect_equal(names(b)[names(a) == "n_detections"], "N Detect")
  expect_equal(names(c3)[names(a) == "n_detections"], "N det.")
})

test_that("display precision is fixed, and `decimals` merges over it", {
  tbl <- sum_tbl()
  f <- format(tbl)
  # counts have no decimals; residency indices have two (mean rows carry the "+/-" and are skipped)
  expect_true(all(grepl("^[0-9]+$", flat(f)$n_detections)))
  expect_true(all(grepl("^[0-9]+\\.[0-9]{2}$", flat(f)$IR1)))
  f2 <- format(tbl, decimals = c(IR1 = 4))
  expect_true(all(grepl("^[0-9]+\\.[0-9]{4}$", flat(f2)$IR1)))
  # naming one column leaves the others alone
  expect_equal(f2$n_detections, f$n_detections)
  # the mean row follows the override too
  expect_true(grepl("\\.[0-9]{4} \\+/- ", utils::tail(f2$IR1, 1)))
})

test_that("`decimals` rejects bad input and points at the right column name", {
  tbl <- sum_tbl()
  expect_error(format(tbl, decimals = c(nope = 1)), "not in this table")
  expect_error(format(tbl, decimals = c(IR1 = -1)), "whole numbers")
  expect_error(format(tbl, decimals = 2), "named numeric")
  expect_error(format(tbl, decimals = c(ID = 1)), "non-numeric")
  # a header read off the formatted table is mapped back to its column
  expect_error(format(tbl, decimals = c(`N Detect` = 0)), "display header")
})

test_that("grouping is presentation: the object is one row per individual either way", {
  tbl <- sum_tbl(id.groups = mobyMeta(rays)$id.groups)
  expect_true("group" %in% names(tbl))
  expect_s3_class(tbl$group, "factor")
  n_ind <- nrow(tbl)

  g <- format(tbl)                       # groups by the `group` column by default
  expect_equal(nrow(g), n_ind + 2L)      # one mean row per group
  expect_equal(sum(grepl("^mean", g$ID)), 2L)
  expect_equal(attr(g, "table.groups")[1], "Raja clavata")

  u <- format(tbl, group.by = FALSE)     # ungrouped: a single mean row
  expect_equal(nrow(u), n_ind + 1L)
  expect_equal(sum(grepl("^mean", u$ID)), 1L)

  expect_error(format(tbl, group.by = "nope"), "must name a column")
})

test_that("a group of one gets no mean row - it would only restate the row", {
  ids <- ray_ids()
  # one group of three, one group of one
  groups <- list(many = ids[1:3], solo = ids[4])
  tbl <- sum_tbl(id.groups = groups)
  f <- format(tbl)
  expect_equal(nrow(f), 4L + 1L)                      # 4 individuals + ONE mean row
  expect_equal(sum(grepl("^mean", f$ID)), 1L)
  # and the mean row belongs to the group that has more than one member
  expect_equal(attr(f, "table.groups")[grepl("^mean", f$ID)], "many")
})

test_that("a single-row table gets no mean row and no '+/- NA'", {
  one <- rays[rays$ID == ray_ids()[1], ]
  tbl <- summaryTable(one, last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"),
                      id.groups = list(solo = ray_ids()[1]), verbose = FALSE)
  f <- format(tbl)
  expect_equal(nrow(f), 1L)
  expect_false(any(grepl("+/-", unlist(f, use.names = FALSE), fixed = TRUE)))
})

test_that("print renders group headings without putting them in the table", {
  tbl <- sum_tbl(id.groups = mobyMeta(rays)$id.groups)
  out <- paste(utils::capture.output(print(tbl)), collapse = "\n")
  expect_match(out, "<mobyTable: summary>")
  expect_match(out, "Raja clavata", fixed = TRUE)
  expect_match(out, "Dasyatis pastinaca", fixed = TRUE)
  expect_match(out, "mean")
  # the grouping column is not repeated on every row under its own heading
  body <- sub("^.*?\n", "", out)                       # drop the banner, which names the group count
  expect_false(grepl("Raja clavata +Raja clavata", body))
})

test_that("transitionsTable splits its composite cells into typed columns", {
  net <- calculateTransitions(rays, spatial.col = "station", verbose = FALSE)
  tt <- transitionsTable(net, verbose = FALSE)
  expect_s3_class(tt, "mobyTable")
  expect_true(all(c("transition", "n_movements", "n_individuals", "pct_individuals",
                    "mean_duration", "error_duration") %in% names(tt)))
  expect_true(is.numeric(tt$n_individuals))
  expect_true(is.numeric(tt$pct_individuals))
  expect_true(is.numeric(tt$mean_duration))
  # a percentage of the group is not something to average across transitions
  f <- format(tt)
  expect_equal(utils::tail(f$pct_individuals, 1), "-")
  expect_true(grepl("+/-", utils::tail(f$mean_duration, 1), fixed = TRUE))
})

test_that("the error statistic is carried on the object and named in the output", {
  tbl <- sum_tbl(error.stat = "se")
  expect_match(utils::tail(format(tbl)$ID, 1), "mean \\+/- se")
  tbl2 <- sum_tbl(error.stat = "sd")
  expect_match(utils::tail(format(tbl2)$ID, 1), "mean \\+/- sd")
})

test_that("an empty table formats and prints without erroring", {
  tbl <- .newMobyTable(data.frame(), kind = "summary")
  expect_equal(nrow(format(tbl)), 0L)
  expect_output(print(tbl), "0 rows")
})
