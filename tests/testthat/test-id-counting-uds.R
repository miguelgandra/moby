# calculateUDs(): $summary_table must have one row per UD that was actually ESTIMATED, not one row
# per DECLARED factor level. .computeUDs() estimates from .cleanData()'s droplevels()'d frame but
# counts COAs from the raw one (R/calculateUDs.R:700, `table(coords[[id.col]])`), and that count is
# the LEFT side of the left join at :704, so it drives the row set. A stale level therefore buys a
# phantom row carrying "N COAs" = 0 and NA areas - a home range that was never estimated, sitting in
# the table a user tabulates into a paper.
#
# Only the single-group path (multiple = FALSE) is affected; the subset/id.groups path droplevels()
# each group at :333 before estimating. Both are covered so the fix cannot regress the healthy one.
#
# Fixture note: this file does NOT use helper-id-divergence.R's ghost_levels_df(). That fixture puts
# every detection at one coordinate (-8.0/37.0), and adehabitatHR::kernelUD() cannot smooth a cloud
# with no spread ("all z values are equal"). ghost_coas() below is the promised variant: the same
# phantom-level idea, the same last/first/middle placement (which is what proves POSITIONAL
# alignment rather than tolerance of a trailing NA), but three well-separated clouds of COAs.


# Three real animals (A, B, C) in disjoint clouds, plus a level the data never had. `month` gives a
# subset variable under which all three IDs recur (the merge branch of the multiple = TRUE path) and
# `sex` one under which they do not (the rbind branch); both branches must stay immune.
ghost_coas <- function(pos = c("last", "first", "middle"), n = 20, ghost = "GHOST") {
  pos <- match.arg(pos)
  set.seed(42)
  cloud <- function(id, cx, cy, month, sex) {
    data.frame(ID = id,
               timebin = as.POSIXct("2021-01-01", tz = "UTC") + seq_len(n) * 1800,
               x = rnorm(n, cx, 300), y = rnorm(n, cy, 300),
               month = month, sex = sex, stringsAsFactors = FALSE)
  }
  d <- rbind(cloud("A",    0,    0, "Jan", "F"), cloud("A",    0,    0, "Feb", "F"),
             cloud("B", 5000, 5000, "Jan", "F"), cloud("B", 5000, 5000, "Feb", "F"),
             cloud("C",    0, 5000, "Jan", "M"), cloud("C",    0, 5000, "Feb", "M"))
  real <- c("A", "B", "C")
  d$ID <- factor(d$ID, levels = switch(pos,
                                       last   = c(real, ghost),
                                       first  = c(ghost, real),
                                       middle = c("A", ghost, "B", "C")))
  d
}

# A fixed grid keeps the run deterministic and fast (no repeat-loop expansion); values are filled
# because terra::as.points() warns on a valueless SpatRaster.
ud_grid <- function() {
  g <- terra::rast(terra::ext(-3000, 8000, -3000, 8000), res = 250, crs = "EPSG:32629")
  terra::values(g) <- 0
  g
}

run_uds <- function(d, ...) {
  suppressWarnings(calculateUDs(d, method = "kde", bandwidth = 800, spatial.grid = ud_grid(),
                                lon.col = "x", lat.col = "y", epsg.code = 32629,
                                verbose = FALSE, ...))
}

# The rows for the animals that really exist, normalised so they can be compared across runs whose
# ID columns carry different level sets (4 levels before the fix, 3 after).
real_rows <- function(tbl, id.col = "ID") {
  out <- tbl[as.character(tbl[[id.col]]) %in% c("A", "B", "C"), , drop = FALSE]
  out[[id.col]] <- as.character(out[[id.col]])
  rownames(out) <- NULL
  out
}


# ---- single-group path (multiple = FALSE) ------------------------------------------------------

test_that("calculateUDs() summary_table is roster-shaped: one row per declared level", {
  skip_on_cran()
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  d   <- ghost_coas("last")            # levels A, B, C, GHOST; only A, B, C observed
  res <- run_uds(d)
  ref <- run_uds(droplevels(d))        # the same data, truth computed independently

  # REGRESSION: one row per UD ACTUALLY ESTIMATED, not per declared factor level. Before the fix
  # this was 4 (nlevels), one more than the 3 UDs computed.
  expect_equal(nrow(res$summary_table), 3L)
  expect_equal(nrow(res$summary_table), length(res$ud))
  expect_equal(nrow(res$summary_table), nrow(ref$summary_table))
  expect_lt(nrow(res$summary_table), nlevels(d$ID))

  # REGRESSION: no phantom row at all - a home range that was never estimated must not appear in a
  # table a user might publish. Before the fix GHOST was present with "N COAs" = 0 and NA areas.
  expect_false("GHOST" %in% as.character(res$summary_table$ID))
  expect_equal(as.character(res$summary_table$ID), c("A", "B", "C"))
  expect_false(any(is.na(unlist(res$summary_table[, c("UD 50% (Km2)", "UD 95% (Km2)")]))))
  expect_true(all(res$summary_table[["N COAs"]] > 0L))

  # INVARIANT: must be identical before and after the fix. Only three UDs are estimated, and only
  # the real animals reach the contours - summary_table is the single place the roster leaks out.
  expect_equal(length(res$ud), 3L)
  expect_equal(names(res$ud), c("A", "B", "C"))
  expect_equal(sort(as.character(res$K50$id)), c("A", "B", "C"))
  expect_equal(sort(as.character(res$K95$id)), c("A", "B", "C"))
  expect_equal(nrow(res$K95), 3L)

  # INVARIANT: must be identical before and after the fix. Column layout, and every number reported
  # for a real animal, are unchanged by the phantom - proof the fix is surgical.
  expect_equal(colnames(res$summary_table), c("ID", "N COAs", "UD 50% (Km2)", "UD 95% (Km2)"))
  expect_equal(colnames(res$summary_table), colnames(ref$summary_table))
  expect_equal(real_rows(res$summary_table), real_rows(ref$summary_table))
  expect_equal(real_rows(res$summary_table)[["N COAs"]], c(40L, 40L, 40L))   # 2 months x 20 COAs
})


test_that("the phantom row sits at the stale level's position, shifting the real animals down", {
  skip_on_cran()
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  # Placing the phantom first and middle (not merely last) is what distinguishes a positional
  # misalignment from a harmless trailing NA: anything reading summary_table by row index - or
  # pairing it with a level-ordered vector such as tagging.dates - lands on the wrong animal.
  first  <- run_uds(ghost_coas("first"))
  middle <- run_uds(ghost_coas("middle"))

  # REGRESSION: wherever the phantom sat in the level order, it is gone and the real animals keep
  # their positions. Before the fix these read c("GHOST","A","B","C") and c("A","GHOST","B","C"),
  # so anything indexing summary_table by row - or pairing it with a level-ordered vector such as
  # tagging.dates - landed on the wrong animal.
  expect_equal(as.character(first$summary_table$ID),  c("A", "B", "C"))
  expect_equal(as.character(middle$summary_table$ID), c("A", "B", "C"))
  expect_equal(first$summary_table[["N COAs"]],  c(40L, 40L, 40L))
  expect_equal(middle$summary_table[["N COAs"]], c(40L, 40L, 40L))

  # INVARIANT: must be identical before and after the fix. Where the phantom sits in the level
  # order changes nothing about the estimates themselves.
  expect_equal(names(first$ud),  c("A", "B", "C"))
  expect_equal(names(middle$ud), c("A", "B", "C"))
  expect_equal(real_rows(first$summary_table), real_rows(middle$summary_table))
})


test_that("a droplevels()'d input already yields one summary_table row per estimated UD", {
  skip_on_cran()
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  res <- run_uds(droplevels(ghost_coas("middle")))

  # INVARIANT: must be identical before and after the fix. This is the healthy case the fix has to
  # leave completely alone - it is the reference the buggy runs above are compared against.
  expect_equal(nrow(res$summary_table), 3L)
  expect_equal(nrow(res$summary_table), length(res$ud))
  expect_equal(as.character(res$summary_table$ID), c("A", "B", "C"))
  expect_equal(res$summary_table[["N COAs"]], c(40L, 40L, 40L))
  expect_false(anyNA(res$summary_table))
})


# ---- grouped path (multiple = TRUE) - already immune, must stay immune -------------------------

test_that("the subset path droplevels() each group and so is immune to stale levels", {
  skip_on_cran()
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  d <- ghost_coas("middle")

  # subset = "month": A, B and C all recur across groups, so calculateUDs takes the merge branch
  # (R/calculateUDs.R:365) and one row per individual carries per-group columns.
  by_month <- run_uds(d, subset = "month")
  # INVARIANT: must be identical before and after the fix.
  expect_equal(nrow(by_month$summary_table), 3L)
  expect_equal(as.character(by_month$summary_table$ID), c("A", "B", "C"))
  expect_false("GHOST" %in% as.character(by_month$summary_table$ID))
  expect_equal(by_month$summary_table[["N COAs [Jan]"]], c(20L, 20L, 20L))
  expect_equal(by_month$summary_table[["N COAs [Feb]"]], c(20L, 20L, 20L))

  # subset = "sex": the groups hold disjoint IDs (F = A, B; M = C), so the rbind branch
  # (R/calculateUDs.R:375) runs instead and rows are stacked with a `group` column.
  by_sex <- run_uds(d, subset = "sex")
  # INVARIANT: must be identical before and after the fix.
  expect_equal(nrow(by_sex$summary_table), 3L)
  expect_equal(colnames(by_sex$summary_table)[1:2], c("group", "ID"))
  expect_equal(as.character(by_sex$summary_table$group), c("F", "F", "M"))
  expect_equal(as.character(by_sex$summary_table$ID), c("A", "B", "C"))
  expect_equal(by_sex$summary_table[["N COAs"]], c(40L, 40L, 40L))
})


# ---- the mechanism itself ----------------------------------------------------------------------

test_that("table() over a stale-level factor is what drives the phantom rows, not .joinKeep()", {
  # Reproduces R/calculateUDs.R:700-704 in isolation, so whoever applies the fix knows the defect is
  # on the ncoas side. .joinKeep() keeping every row of `x` is contractual (see
  # test-helpers-reshape.R) and must NOT be touched.
  f <- factor(c("A", "A", "B", "B"), levels = c("A", "B", "GHOST"))

  ncoas <- as.data.frame(table(f))
  colnames(ncoas) <- c("ID", "N COAs")
  # table() emits one entry PER LEVEL, including a zero for the level the data never had.
  expect_equal(nrow(ncoas), 3L)
  expect_equal(ncoas[["N COAs"]], c(2L, 2L, 0L))

  # kud_table only ever holds the individuals kernelUD() actually estimated, because .cleanData()
  # droplevels() first.
  kud_table <- data.frame(ID = c("A", "B"), "UD 95% (Km2)" = c("1.00", "2.00"),
                          check.names = FALSE, stringsAsFactors = FALSE)

  # INVARIANT: must be identical before and after the fix. ncoas is the LEFT side, so the join
  # returns nrow(ncoas) rows - the roster - with NA areas for the level that has no UD.
  joined <- .joinKeep(ncoas, kud_table, by = "ID", type = "left")
  expect_equal(nrow(joined), nrow(ncoas))
  expect_equal(as.character(joined$ID), c("A", "B", "GHOST"))
  expect_true(is.na(joined[["UD 95% (Km2)"]][3]))

  # ...and droplevels() on the counted factor is all it takes for the row set to match the UDs.
  fixed <- as.data.frame(table(droplevels(f)))
  colnames(fixed) <- c("ID", "N COAs")
  expect_equal(nrow(.joinKeep(fixed, kud_table, by = "ID", type = "left")), nrow(kud_table))
})
