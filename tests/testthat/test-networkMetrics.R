test_that("networkMetrics computes directed metrics for movement networks", {
  skip_if_not_installed("igraph")
  d <- data.frame(
    ID = factor(c(rep("A", 6), rep("B", 4), rep("C", 4))),
    datetime = as.POSIXct("2023-01-01", tz = "UTC") +
      c(0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0, 1, 2, 3) * 3600,
    site = c("S1", "S1", "S2", "S2", "S3", "S1", "S2", "S3", "S3", "S2", "S1", "S2", "S3", "S1"))
  md <- as_moby(d, station.col = "site",
                tagging.dates = as.POSIXct(c(A = "2023-01-01", B = "2023-01-01", C = "2023-01-01"), tz = "UTC"))
  net <- calculateTransitions(md, spatial.col = "site")
  m <- networkMetrics(net)
  expect_s3_class(m, "mobyNetworkMetrics")
  expect_true(all(c("in_degree", "out_degree", "in_strength", "out_strength",
                    "betweenness", "community") %in% colnames(m$nodes)))
  expect_true(all(c("density", "modularity", "n_communities", "reciprocity") %in% colnames(m$network)))
  # a through-route node should have positive betweenness (movement corridor)
  expect_true(any(m$nodes$betweenness > 0))
})

test_that("networkMetrics computes undirected metrics for association networks", {
  skip_if_not_installed("igraph")
  set.seed(1)
  ids <- c("A", "B", "C", "D"); n <- 60
  tb <- as.POSIXct("2023-01-01", tz = "UTC") + (0:(n - 1)) * 3600
  df <- do.call(rbind, lapply(ids, function(i)
    data.frame(ID = i, timebin = tb, station = sample(c("R1", "R2", "R3"), n, replace = TRUE))))
  df$ID <- factor(df$ID)
  md <- as_moby(df, timebin.col = "timebin", station.col = "station",
                tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  wt <- suppressWarnings(suppressMessages(createWideTable(md, value.col = "station", verbose = FALSE)))
  net <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  m <- networkMetrics(net)
  expect_true(all(c("degree", "strength", "betweenness", "eigenvector", "clustering", "community") %in% colnames(m$nodes)))
  expect_true("transitivity" %in% colnames(m$network))
  expect_equal(length(unique(m$nodes$node)), 4L)            # all individuals are nodes
})

test_that("networkMetrics validates input and weight column", {
  skip_if_not_installed("igraph")
  expect_error(networkMetrics(data.frame(a = 1)), "mobyNetwork")
  d <- data.frame(ID = factor(c("A", "A", "B")),
                  datetime = as.POSIXct("2023-01-01", tz = "UTC") + c(0, 1, 0) * 3600,
                  site = c("S1", "S2", "S1"))
  net <- calculateTransitions(as_moby(d, station.col = "site",
                                      tagging.dates = as.POSIXct("2023-01-01", tz = "UTC")),
                              spatial.col = "site")
  expect_error(networkMetrics(net, weight = "nope"), "not found")
})
