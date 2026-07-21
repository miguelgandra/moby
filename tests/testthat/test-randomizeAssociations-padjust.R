make_overlap_inputs <- function(seed = 1) {
  set.seed(seed)
  ids <- c("A", "B", "C", "D"); n <- 60
  tb <- as.POSIXct("2023-01-01", tz = "UTC") + (0:(n - 1)) * 3600
  df <- do.call(rbind, lapply(ids, function(id) data.frame(
    ID = id, timebin = tb, station = sample(c("R1", "R2", "R3"), n, replace = TRUE),
    stringsAsFactors = FALSE)))
  df$ID <- factor(df$ID)
  md <- as_moby(df, timebin.col = "timebin", station.col = "station",
                tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  wt <- suppressWarnings(suppressMessages(createWideTable(md, value.col = "station", verbose = FALSE)))
  ov <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  list(table = wt, overlaps = ov)
}

test_that("randomizeAssociations reports both raw and multiple-comparison-adjusted p-values", {
  inp <- make_overlap_inputs()
  r <- suppressWarnings(suppressMessages(
    randomizeAssociations(inp$table, inp$overlaps, iterations = 100, random.seed = 1, p.adjust.method = "fdr")))
  pr <- r$pairwise_results
  expect_true(all(c("p_value", "p_adjusted") %in% colnames(pr)))
  padj <- suppressWarnings(as.numeric(pr$p_adjusted))
  praw <- suppressWarnings(as.numeric(pr$p_value))
  # Benjamini-Hochberg never decreases a p-value
  expect_true(all(padj >= praw - 1e-9, na.rm = TRUE))
})

test_that("p.adjust.method = 'none' leaves p-values unchanged; invalid method errors", {
  inp <- make_overlap_inputs()
  r <- suppressWarnings(suppressMessages(
    randomizeAssociations(inp$table, inp$overlaps, iterations = 100, random.seed = 1, p.adjust.method = "none")))
  pr <- r$pairwise_results
  expect_equal(suppressWarnings(as.numeric(pr$p_adjusted)),
               suppressWarnings(as.numeric(pr$p_value)))
  expect_error(
    suppressWarnings(suppressMessages(randomizeAssociations(inp$table, inp$overlaps, p.adjust.method = "bogus"))),
    "p.adjust.method")
})