# Build a movement-type mobyNetwork from a small synthetic detection dataset
move_net <- function(id.groups = NULL, seed = 1) {
  set.seed(seed)
  sites <- data.frame(site = paste0("S", 1:5),
                      lon = c(-8, -7.6, -7, -6.6, -6), lat = c(37, 37.5, 38, 38.5, 39))
  ids <- LETTERS[1:6]
  rows <- do.call(rbind, lapply(ids, function(id) {
    k <- sample(6:12, 1); s <- sites[sample(nrow(sites), k, replace = TRUE), ]
    data.frame(ID = id, datetime = as.POSIXct("2023-01-01", tz = "UTC") + seq_len(k) * 3600,
               site = s$site, lon = s$lon, lat = s$lat, stringsAsFactors = FALSE)
  }))
  rows$ID <- factor(rows$ID)
  md <- as_moby(rows, station.col = "site", lon.col = "lon", lat.col = "lat",
                tagging.dates = stats::setNames(as.POSIXct(rep("2023-01-01", 6), tz = "UTC"), ids))
  calculateTransitions(md, spatial.col = "site", id.groups = id.groups)
}

test_that("plotMovements draws a projected map and returns a tidy node table", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  net <- move_net()
  res <- suppressWarnings(suppressMessages(plotMovements(net, epsg.code = 32629)))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("site", "group", "x", "y", "n_individuals") %in% names(res)))
  expect_s3_class(attr(res, "edges"), "data.frame")
  # coordinates were projected to metres (UTM), not left as lon/lat degrees
  expect_gt(max(abs(res$x)), 1000)
})

test_that("plotMovements handles groups, edge.metric and the no-coordinate fallback", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  netg <- move_net(id.groups = list(g1 = c("A", "B", "C"), g2 = c("D", "E", "F")))
  pm <- function(n, ...) suppressWarnings(suppressMessages(plotMovements(n, ...)))
  expect_no_error(pm(netg, epsg.code = 32629, color.nodes.by = "group", ncol = 2))
  expect_no_error(pm(netg, epsg.code = 32629, edge.metric = "individuals", cex = 0.8))
  # a network without coordinates falls back to a force-directed layout (with a warning)
  net_nc <- move_net(); attr(net_nc, "has.coords") <- FALSE
  expect_warning(plotMovements(net_nc), "force-directed")
})

test_that("plotMovements validates inputs and writes files without leaking a device", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  net <- move_net()
  expect_error(plotMovements(data.frame(a = 1)), "movement-type")
  expect_error(plotMovements(net, nodes.size = 0.05), "min and max")
  before <- length(dev.list()); f <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(plotMovements(net, epsg.code = 32629, file = f))))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)                     # no device leaked
})

test_that("plotMovements draws a default coastline when land.shape is NULL and coastline = TRUE", {
  skip_if_not_installed("maps")
  net <- move_net()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(
    plotMovements(net, epsg.code = 32629, coastline = TRUE))))       # projected map + auto coastline
})
