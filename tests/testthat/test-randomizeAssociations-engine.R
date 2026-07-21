# The matrix-based permutation engine must reproduce the previous implementation's output bit-for-bit
# (same null model, same RNG stream), and randomized_overlaps is now a [dyads x iterations] matrix.
# The reference values in ref-randomizeAssociations.rds were captured from the pre-optimisation engine.

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

test_that("permutation engine reproduces the reference output bit-for-bit", {
  ref <- readRDS(test_path("ref-randomizeAssociations.rds"))
  wt <- mk_assoc_fixture()

  ov1 <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  r <- rz(wt, ov1)
  expect_equal(r$pairwise_results$mean_null,  ref$plain$mean_null)
  expect_equal(r$pairwise_results$p_value,    ref$plain$p)
  expect_equal(r$pairwise_results$p_adjusted, ref$plain$padj)

  rc <- rz(wt, ov1, constraint.by = "phase")
  expect_equal(rc$pairwise_results$mean_null,  ref$constraint$mean_null)
  expect_equal(rc$pairwise_results$p_value,     ref$constraint$p)
  expect_equal(rc$pairwise_results$p_adjusted,  ref$constraint$padj)

  ov2 <- suppressWarnings(suppressMessages(
    calculateAssociations(wt, id.groups = list(g1 = c("A", "B"), g2 = c("C", "D", "E")))))
  rg <- rz(wt, ov2)
  expect_equal(rg$pairwise_results$mean_null,  ref$groups$mean_null)
  expect_equal(rg$pairwise_results$p_adjusted, ref$groups$padj)

  ov3 <- suppressWarnings(suppressMessages(calculateAssociations(wt, subset = "half")))
  rs <- rz(wt, ov3)
  expect_equal(rs$pairwise_results$mean_null, ref$subset$mean_null)
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
