# Regression tests for the base-R reshape/join helpers that replaced reshape2 and plyr.
# Expectations are frozen (hard-coded), not compared against reshape2/plyr, so the suite has no
# dependency on the removed packages. The .castWide cases pin the fill / drop = FALSE / factor-level
# and class-preservation semantics the association, UD and network matrices depend on.

test_that(".castWide spreads a POSIXct row key x factor col key (sum, fill = 0)", {
  tz <- "Etc/GMT-10"
  d <- data.frame(
    timebin = as.POSIXct(c("2023-06-01 00:00", "2023-06-01 01:00", "2023-06-01 00:00",
                           "2023-06-01 02:00", "2023-06-01 01:00"), tz = tz),
    id      = factor(c("A", "A", "B", "B", "A"), levels = c("A", "B", "C")),
    value   = c(10, 20, 30, 40, 5))
  w <- moby:::.castWide(d, "timebin", "id", "value", fun.aggregate = sum, fill = 0)

  expect_identical(names(w), c("timebin", "A", "B", "C"))
  expect_s3_class(w$timebin, "POSIXct")
  expect_identical(attr(w$timebin, "tzone"), tz)                 # tzone survives the round-trip
  expect_identical(w$timebin,
                   as.POSIXct(c("2023-06-01 00:00", "2023-06-01 01:00", "2023-06-01 02:00"), tz = tz))
  expect_equal(w$A, c(10, 25, 0))                                # 20 + 5 aggregated at 01:00
  expect_equal(w$B, c(30, 0, 40))
  expect_equal(w$C, c(0, 0, 0))                                  # empty factor level kept (drop = FALSE)
})

test_that(".castWide keeps a legitimately-NA present cell rather than filling it", {
  tz <- "Etc/GMT-10"
  assignVal <- function(x) names(which.max(table(x)))
  d <- data.frame(
    timebin = as.POSIXct(c("2023-06-01 00:00", "2023-06-01 00:00", "2023-06-01 01:00"), tz = tz),
    id      = factor(c("A", "A", "B"), levels = c("A", "B")),
    value   = c("x", "x", "y"), stringsAsFactors = FALSE)
  w <- moby:::.castWide(d, "timebin", "id", "value", fun.aggregate = assignVal, fill = NA_character_)

  # A@00:00 present -> "x"; B@00:00 structurally empty -> NA (fill); A@01:00 empty -> NA; B@01:00 -> "y"
  expect_equal(w$A, c("x", NA))
  expect_equal(w$B, c(NA, "y"))
})

test_that(".castWide builds a square padded id1 x id2 matrix (single value per cell)", {
  ids <- c("A", "B", "C")
  pw <- data.frame(id1 = factor(c("A", "A", "B"), levels = ids),
                   id2 = factor(c("B", "C", "C"), levels = ids),
                   association = c(0.5, 0.2, 0.8))
  w <- moby:::.castWide(pw, "id1", "id2", "association", fun.aggregate = function(z) z[1L], fill = NA)

  expect_identical(names(w), c("id1", "A", "B", "C"))
  expect_equal(nrow(w), 3L)                                      # all 3 levels as rows (padded)
  expect_equal(w$B, c(0.5, NA, NA))
  expect_equal(w$C, c(0.2, 0.8, NA))
})

test_that(".castWide fills only structurally-empty cells (max, fill = 0)", {
  pc <- data.frame(day = c("d1", "d1", "d2", "d3"), hour = c("h1", "h2", "h1", "h2"),
                   var = c(3, 7, 2, 9))
  w <- moby:::.castWide(pc, "day", "hour", "var", fun.aggregate = max, fill = 0)
  expect_identical(w$day, c("d1", "d2", "d3"))
  expect_equal(w$h1, c(3, 2, 0))                                 # d3/h1 empty -> 0
  expect_equal(w$h2, c(7, 0, 9))                                 # d2/h2 empty -> 0
})

test_that(".castWide yields NA (not a silent list-column) for an all-NA cell, and errors on multi-value", {
  # most-frequent aggregator returns character(0) when a populated cell's values are all NA; the guard
  # must turn that into a scalar NA so the wide table never becomes a corrupt list-column.
  assignVal <- function(x) names(which.max(table(x)))
  d <- data.frame(tb = c("t1", "t1", "t2"), id = factor(c("A", "A", "B"), levels = c("A", "B")),
                  v = c(NA_character_, NA_character_, "y"), stringsAsFactors = FALSE)
  w <- moby:::.castWide(d, "tb", "id", "v", fun.aggregate = assignVal, fill = NA_character_)
  expect_false(any(vapply(w, is.list, logical(1))))              # no list-columns
  expect_true(is.na(w$A[w$tb == "t1"]))                          # all-NA cell -> NA scalar
  # a genuinely non-scalar aggregator must fail loudly, not silently corrupt
  expect_error(moby:::.castWide(d, "tb", "id", "v", fun.aggregate = function(z) z, fill = NA),
               "single value per cell")
})

test_that(".meltList turns a named list into value / L1 columns", {
  L <- list(grp1 = c("A", "B"), grp2 = "C")
  m <- moby:::.meltList(L)
  expect_identical(names(m), c("value", "L1"))
  expect_identical(m$value, c("A", "B", "C"))
  expect_identical(m$L1, c("grp1", "grp1", "grp2"))
})

test_that(".joinKeep is a left join that preserves x's row order and column layout", {
  x <- data.frame(k = c("c", "a", "b", "a"), v = 1:4, stringsAsFactors = FALSE)
  y <- data.frame(k = c("b", "a"), w = c(99, 88), stringsAsFactors = FALSE)
  z <- moby:::.joinKeep(x, y, by = "k", type = "left")
  expect_identical(names(z), c("k", "v", "w"))                   # x cols then y non-key cols
  expect_identical(z$k, x$k)                                     # x order preserved
  expect_identical(z$v, x$v)
  expect_equal(z$w, c(NA, 88, 99, 88))                           # unmatched "c" -> NA
  expect_error(moby:::.joinKeep(x, y, by = "k", type = "inner"), "only supports")
})

test_that(".joinKeep keeps both columns when x and y share a non-key column name", {
  # shadeSeasons joins a 'start' and an 'end' table that both carry an aggregate column named "x";
  # the result must keep all three columns (key + both values) so the positional colnames() succeeds.
  ss <- data.frame(Group.1 = c("spring", "summer"),
                   x = as.POSIXct(c("2023-03-01", "2023-06-01"), tz = "UTC"))
  se <- data.frame(Group.1 = c("spring", "summer"),
                   x = as.POSIXct(c("2023-05-31", "2023-08-31"), tz = "UTC"))
  j <- moby:::.joinKeep(ss, se, by = "Group.1", type = "left")
  expect_equal(ncol(j), 3L)
  expect_identical(j[[2]], ss$x)                                 # x's value column preserved
  expect_identical(j[[3]], se$x)                                 # y's value column appended, not overwritten
})

test_that(".joinKeep preserves a list-column in x (shadeSeasons' aggregate simplify=FALSE)", {
  # shadeSeasons joins two aggregate(simplify = FALSE) tables whose value column is a LIST-column,
  # then unlist()s it downstream; the join must keep the list-column intact, not expand it.
  ss <- stats::aggregate(list(x = c("a", "b", "c")), by = list(Group.1 = c("g1", "g1", "g2")),
                         FUN = min, simplify = FALSE)
  se <- stats::aggregate(list(x = c("a", "b", "c")), by = list(Group.1 = c("g1", "g1", "g2")),
                         FUN = max, simplify = FALSE)
  expect_true(is.list(ss$x))
  j <- moby:::.joinKeep(ss, se, by = "Group.1", type = "left")
  expect_equal(ncol(j), 3L)
  expect_true(is.list(j[[2]]) && is.list(j[[3]]))                # list-columns intact, not expanded
  expect_identical(unlist(j[[2]]), unlist(ss$x))
  expect_identical(unlist(j[[3]]), unlist(se$x))
})

test_that(".rbindFill unions columns, keeps 0-row frames, drops NULLs", {
  a <- data.frame(x = 1:2, y = c("p", "q"), stringsAsFactors = FALSE)
  b <- data.frame(y = "r", z = 9, stringsAsFactors = FALSE)
  out <- moby:::.rbindFill(a, b)
  expect_identical(names(out), c("x", "y", "z"))
  expect_equal(out$x, c(1, 2, NA))
  expect_identical(out$y, c("p", "q", "r"))
  expect_equal(out$z, c(NA, NA, 9))

  expect_equal(moby:::.rbindFill(list(a, b)), out)               # single-list form
  expect_equal(moby:::.rbindFill(a, NULL, b), out)               # NULLs dropped
  expect_null(moby:::.rbindFill(NULL, NULL))                     # all-NULL -> NULL

  # 0-row frame still contributes its columns to the union
  e <- a[0, ]
  expect_true("x" %in% names(moby:::.rbindFill(e, b)))
})

test_that(".mapValues replaces matched values and leaves the rest unchanged", {
  expect_identical(moby:::.mapValues(c("A", "B", "D", "A"), from = c("A", "B"), to = c("aa", "bb")),
                   c("aa", "bb", "D", "aa"))
  # the plyr call sites pass warn_missing = FALSE; it must be silently absorbed
  expect_identical(moby:::.mapValues("A", "A", "z", warn_missing = FALSE), "z")
})
