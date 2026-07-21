# Build association inputs (overlaps + null-model results) from a small synthetic dataset
assoc_inputs <- function(ids = c("A", "B", "C", "D", "E"), id.groups = NULL, seed = 1) {
  set.seed(seed)
  n <- 90; tb <- as.POSIXct("2023-01-01", tz = "UTC") + (0:(n - 1)) * 3600
  df <- do.call(rbind, lapply(ids, function(id) data.frame(
    ID = id, timebin = tb, station = sample(c("R1", "R2", "R3"), n, replace = TRUE), stringsAsFactors = FALSE)))
  df$ID <- factor(df$ID)
  md <- as_moby(df, timebin.col = "timebin", station.col = "station", tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  wt <- suppressWarnings(suppressMessages(createWideTable(md, value.col = "station", verbose = FALSE)))
  ov <- suppressWarnings(suppressMessages(calculateAssociations(wt, id.groups = id.groups)))
  rr <- suppressWarnings(suppressMessages(randomizeAssociations(wt, ov, iterations = 80, random.seed = 1)))
  list(overlaps = ov, random.results = rr)
}

test_that("plotAssociations draws network / histogram / both and returns tidy per-type stats", {
  skip_if_not_installed("qgraph")
  x <- assoc_inputs()
  pdf(tempfile(fileext = ".pdf"), width = 7, height = 8); on.exit(dev.off(), add = TRUE)
  pa <- function(...) suppressWarnings(suppressMessages(plotAssociations(...)))
  res <- pa(overlaps = x$overlaps)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("type", "n_dyads", "mean_overlap", "mean_binary_degree", "edge_cut") %in% names(res)))
  expect_equal(sum(res$n_dyads), choose(5, 2))               # 10 dyads among 5 individuals
  expect_no_error(pa(random.results = x$random.results))     # histogram only
  expect_no_error(pa(overlaps = x$overlaps, random.results = x$random.results, main = "net"))  # both
  expect_no_error(pa(overlaps = x$overlaps, color.by = "single"))
})

test_that("plotAssociations accepts edge.color without the id.groups-not-found crash", {
  skip_if_not_installed("qgraph")
  x <- assoc_inputs()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  # previously: 'object id.groups not found' whenever edge.color was supplied
  expect_no_error(suppressWarnings(suppressMessages(plotAssociations(overlaps = x$overlaps, edge.color = "grey40"))))
})

test_that("plotAssociations aligns per-panel stats across id.groups comparisons", {
  skip_if_not_installed("qgraph")
  x <- assoc_inputs(ids = LETTERS[1:6], id.groups = list(g1 = c("A", "B", "C"), g2 = c("D", "E", "F")))
  pdf(tempfile(fileext = ".pdf"), width = 8, height = 8); on.exit(dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(plotAssociations(overlaps = x$overlaps, random.results = x$random.results)))
  expect_true(nrow(res) >= 2)                                # one row per comparison type
  expect_true(all(res$n_dyads > 0))                          # each aligned panel has its own dyads
})

# Capture the node colours plotAssociations hands to qgraph, one entry per panel, plus any
# warnings raised along the way. Mocking is the only way in: the colours never reach the return value.
capture_node_colors <- function(...) {
  panels <- list(); warns <- character()
  testthat::local_mocked_bindings(
    qgraph = function(input, ...) {
      panels[[length(panels) + 1]] <<- list(ids = colnames(input), color = list(...)$color)
      invisible(NULL)
    }, .package = "qgraph")
  pdf(tempfile(fileext = ".pdf"), width = 8, height = 8); on.exit(dev.off(), add = TRUE)
  withCallingHandlers(suppressMessages(plotAssociations(...)),
                      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(panels = panels, warnings = warns)
}

test_that("plotAssociations colours nodes by id.group instead of blanking them to NA", {
  skip_if_not_installed("qgraph")
  # regression: node_group is a factor, so mapping ids through its *labels* and coercing back
  # ("A" -> as.numeric -> NA) left every node colour NA and made nodes.color inert
  x <- assoc_inputs(ids = LETTERS[1:6], id.groups = list(g1 = c("A", "B", "C"), g2 = c("D", "E", "F")))
  out <- capture_node_colors(overlaps = x$overlaps)
  pal <- moby:::.okabe_ito_pal(2)

  expect_length(out$panels, 3)                                  # g1<->g1, g1<->g2, g2<->g2
  for (panel in out$panels) {
    expect_false(any(is.na(panel$color)))
    expect_equal(panel$color, pal[ifelse(panel$ids %in% c("A", "B", "C"), 1L, 2L)])
  }
  between <- Filter(function(p) length(p$ids) == 6, out$panels)[[1]]   # the panel holding every id
  expect_equal(between$ids, LETTERS[1:6])
  expect_equal(between$color, pal[c(1, 1, 1, 2, 2, 2)])
  expect_false(any(grepl("NAs introduced by coercion", out$warnings)))
})

test_that("plotAssociations honours an explicit nodes.color across id.groups", {
  skip_if_not_installed("qgraph")
  # regression: nodes.color used to be entirely inert -- the plot was byte-identical whatever was passed
  x <- assoc_inputs(ids = LETTERS[1:6], id.groups = list(g1 = c("A", "B", "C"), g2 = c("D", "E", "F")))
  out <- capture_node_colors(overlaps = x$overlaps, nodes.color = c("red", "blue"))
  between <- Filter(function(p) length(p$ids) == 6, out$panels)[[1]]
  expect_equal(between$color, c("red", "red", "red", "blue", "blue", "blue"))
})

test_that("plotAssociations validates inputs and writes files without leaking a device", {
  skip_if_not_installed("qgraph")
  x <- assoc_inputs()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_error(plotAssociations(), "must be provided")
  before <- length(dev.list()); f <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(plotAssociations(overlaps = x$overlaps, file = f))))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)
})

test_that("plotAssociations histogram panels use valid graphical params (cex.lab/cex.axis regression)", {
  skip_if_not_installed("qgraph")
  x <- assoc_inputs()
  pdf(tempfile(fileext = ".pdf"), width = 7, height = 8); on.exit(dev.off(), add = TRUE)
  # the null-distribution histograms once passed cex_lab=/cex_axis= (underscores) to title()/axis(),
  # which are not graphical parameters and warned on every panel; capture warnings and assert none.
  warns <- character(0)
  withCallingHandlers(
    suppressMessages(plotAssociations(random.results = x$random.results)),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  expect_false(any(grepl("not a graphical parameter", warns)))
})
