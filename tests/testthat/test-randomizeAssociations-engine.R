# The matrix-based permutation engine must reproduce the previous implementation's null model and RNG
# stream, and randomized_overlaps is now a [dyads x iterations] matrix. The reference values in
# ref-randomizeAssociations.rds were captured from the pre-optimisation engine.
#
# WHAT IS ASSERTED, AND WHY IT IS NOT UNIFORM:
#
# * p_value / p_adjusted are compared EXACTLY. These carry the reproducibility claim: they are tail
#   counts over the null draws, formatted with sprintf("%.3f", ...) (randomizeAssociations.R:227), so
#   they are integer-derived and bit-identical on every IEEE-754 platform. If the permutation stream
#   drifted at all, they would move.
#
# * mean_null is compared to within ONE unit of its reported precision, because asserting it exactly
#   is not sound on any platform mix. It is rounded to 2 dp (randomizeAssociations.R:327) and comes
#   from rowMeans() (helpers-stats.R:53), whose C backend accumulates in a single pass into LDOUBLE.
#   LDOUBLE's WIDTH is platform-dependent - 53-bit on Apple arm64 and on CRAN's noLD build, 64-bit on
#   x86-64, 113-bit elsewhere - so the same data can give row means differing in the last ulp. Null
#   values here are all k/120*100, so exact .xx5 ties are structural rather than rare, and round(x, 2)
#   turns a ~1e-15 difference into 0.01. Measured: rows 7 and 10 of the `subset` scenario have true
#   means of exactly 32.875 and 33.175, and rounded differently on arm64 vs x86-64.
#
#   Cost of that tolerance, stated plainly: null values are quantised in steps of 100/120 pp, so one
#   altered draw out of 200 moves a row mean by ~0.004 pp. A regression that changed only one or two
#   draws per dyad would slip past this bound - which is precisely what the exact p_value and
#   p_adjusted assertions above are there to catch.

# mean_null equality "to the last reported digit". NOTE: testthat's `tolerance` is RELATIVE, so
# tolerance = 0.011 on values near 33 would silently permit a difference of ~0.36. Compare absolutely.
expect_mean_null <- function(got, want)
  expect_lt(max(abs(got - want), na.rm = TRUE), 0.011)

mk_assoc_fixture <- function() {
  set.seed(11); ids <- c("A", "B", "C", "D", "E"); n <- 240
  tb <- as.POSIXct("2023-01-01", tz = "UTC") + (0:(n - 1)) * 3600
  df <- do.call(rbind, lapply(ids, function(id) data.frame(
    ID = id, timebin = tb, station = sample(c("R1", "R2", "R3"), n, replace = TRUE), stringsAsFactors = FALSE)))
  df$ID <- factor(df$ID)
  md <- as_moby(df, timebin.col = "timebin", station.col = "station",
                tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  wt <- suppressWarnings(suppressMessages(createWideTable(md, value.col = "station", verbose = FALSE)))
  wt$phase <- c("day", "night")[(as.integer(format(wt$timebin, "%H")) %% 2) + 1]
  wt$half  <- ifelse(seq_len(nrow(wt)) <= nrow(wt) / 2, "h1", "h2")
  wt
}
rz <- function(...) suppressWarnings(suppressMessages(
  randomizeAssociations(..., iterations = 200, random.seed = 1, verbose = FALSE)))

test_that("permutation engine reproduces the reference null model and RNG stream", {
  ref <- readRDS(test_path("ref-randomizeAssociations.rds"))
  wt <- mk_assoc_fixture()

  ov1 <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  r <- rz(wt, ov1)
  expect_mean_null(r$pairwise_results$mean_null,  ref$plain$mean_null)
  expect_equal(r$pairwise_results$p_value,    ref$plain$p)
  expect_equal(r$pairwise_results$p_adjusted, ref$plain$padj)

  rc <- rz(wt, ov1, constraint.by = "phase")
  expect_mean_null(rc$pairwise_results$mean_null,  ref$constraint$mean_null)
  expect_equal(rc$pairwise_results$p_value,     ref$constraint$p)
  expect_equal(rc$pairwise_results$p_adjusted,  ref$constraint$padj)

  ov2 <- suppressWarnings(suppressMessages(
    calculateAssociations(wt, id.groups = list(g1 = c("A", "B"), g2 = c("C", "D", "E")))))
  rg <- rz(wt, ov2)
  expect_mean_null(rg$pairwise_results$mean_null,  ref$groups$mean_null)
  expect_equal(rg$pairwise_results$p_adjusted, ref$groups$padj)

  ov3 <- suppressWarnings(suppressMessages(calculateAssociations(wt, subset = "half")))
  rs <- rz(wt, ov3)
  expect_mean_null(rs$pairwise_results$mean_null, ref$subset$mean_null)
  expect_equal(rs$pairwise_results$p_value,   ref$subset$p)
})

test_that("randomized_overlaps is a [dyads x iterations] matrix with dyad row names", {
  wt <- mk_assoc_fixture()
  ov <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  ro <- rz(wt, ov)$randomized_overlaps
  expect_true(is.matrix(ro))
  expect_equal(dim(ro), c(choose(5, 2), 200))
  expect_true(all(rz(wt, ov)$pairwise_results$pair %in% rownames(ro)))

  # subset -> a named list of matrices, one per subset level
  ov3 <- suppressWarnings(suppressMessages(calculateAssociations(wt, subset = "half")))
  ros <- rz(wt, ov3)$randomized_overlaps
  expect_type(ros, "list")
  expect_named(ros, c("h1", "h2"))
  expect_true(all(vapply(ros, is.matrix, logical(1))))
})

test_that("a fixed seed gives reproducible results", {
  wt <- mk_assoc_fixture()
  ov <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  expect_equal(rz(wt, ov)$pairwise_results, rz(wt, ov)$pairwise_results)
})
