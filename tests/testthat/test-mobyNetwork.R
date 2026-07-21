make_assoc_inputs <- function(seed = 1) {
  set.seed(seed)
  ids <- c("A", "B", "C", "D"); n <- 50
  tb <- as.POSIXct("2023-01-01", tz = "UTC") + (0:(n - 1)) * 3600
  df <- do.call(rbind, lapply(ids, function(id) data.frame(
    ID = id, timebin = tb, station = sample(c("R1", "R2", "R3"), n, replace = TRUE),
    stringsAsFactors = FALSE)))
  df$ID <- factor(df$ID)
  md <- as_moby(df, timebin.col = "timebin", station.col = "station",
                tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  suppressWarnings(suppressMessages(createWideTable(md, value.col = "station", verbose = FALSE)))
}

test_that("calculateAssociations returns a mobyNetwork that is also a data.frame", {
  wt <- make_assoc_inputs()
  net <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  expect_true(is_mobyNetwork(net))
  expect_s3_class(net, c("mobyNetwork", "data.frame"))
  expect_equal(networkType(net), "association")
  expect_true(is.data.frame(net))                      # back-compat: still a data frame
  expect_true("ids" %in% names(attributes(net)))       # consumers read this attribute
})

test_that("mobyNetwork accessors return edges and nodes", {
  wt <- make_assoc_inputs()
  net <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  edges <- networkEdges(net)
  nodes <- networkNodes(net)
  expect_s3_class(edges, "data.frame")
  expect_false(inherits(edges, "mobyNetwork"))         # plain data frame
  expect_true("association" %in% colnames(edges))      # the association index column
  expect_equal(nrow(nodes), 4L)
  expect_true("ID" %in% colnames(nodes))
})

test_that("randomizeAssociations accepts a mobyNetwork object", {
  wt <- make_assoc_inputs()
  net <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  r <- suppressWarnings(suppressMessages(randomizeAssociations(wt, net, iterations = 50, random.seed = 1)))
  expect_true(is.list(r))
  expect_true("pairwise_results" %in% names(r))
})

test_that("plot.mobyNetwork errors on an unknown network type", {
  wt <- make_assoc_inputs()
  net <- suppressWarnings(suppressMessages(calculateAssociations(wt)))
  fake <- net; attr(fake, "network.type") <- "bogus"
  expect_error(plot(fake), "Unknown network type")
})
