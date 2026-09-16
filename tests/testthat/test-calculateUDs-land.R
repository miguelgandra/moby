land_test_grid <- function() {
  r <- terra::rast(terra::ext(-5000, 5000, -5000, 5000),
                   resolution = 100, crs = "EPSG:32629")
  terra::values(r) <- 0
  r
}


land_test_ud <- function() {
  set.seed(7301)
  xy <- cbind(stats::rnorm(40, 500, 250), stats::rnorm(40, 0, 250))
  pts <- sp::SpatialPointsDataFrame(
    xy, data.frame(ID = factor(rep("fish", nrow(xy)))),
    proj4string = sp::CRS(SRS_string = sf::st_crs(32629)$wkt)
  )
  adehabitatHR::kernelUD(pts, h = 500,
                         grid = .spatRasterToSpatialPixels(land_test_grid()))
}


legacy_subtract_land <- function(uds, land, crs) {
  xy <- as.data.frame(uds[[1]]@coords)
  colnames(xy) <- c("x", "y")
  points <- sf::st_as_sf(xy, coords = c("x", "y"), crs = crs)
  mask <- as.logical(sf::st_intersects(points, sf::st_union(land), sparse = FALSE))

  for (i in seq_along(uds)) {
    density <- uds[[i]]$ud
    total <- sum(density, na.rm = TRUE)
    density[mask] <- 0
    corrected <- sum(density, na.rm = TRUE)
    if (corrected > 0) density <- density * (total / corrected)
    uds[[i]]$ud <- density
  }
  uds
}


test_that("SpatRaster conversion preserves the regular grid without materialising points", {
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  r <- terra::rast(terra::ext(455000, 455800, 673500, 674100),
                   ncol = 8, nrow = 3, crs = "EPSG:3812")
  direct <- .spatRasterToSpatialPixels(r)

  # Reproduce the former conversion on a filled raster: the new path must preserve
  # every numeric part of the grid, differing at most in coordinate dimnames.
  terra::values(r) <- 0
  legacy <- methods::as(
    sf::as_Spatial(sf::st_as_sf(terra::as.points(r[[1]]))),
    "SpatialPixels"
  )

  expect_identical(unname(direct@coords), unname(legacy@coords))
  expect_identical(unname(direct@bbox), unname(legacy@bbox))
  expect_identical(unname(direct@grid@cellcentre.offset),
                   unname(legacy@grid@cellcentre.offset))
  expect_identical(unname(direct@grid@cellsize), unname(legacy@grid@cellsize))
  expect_identical(unname(direct@grid@cells.dim), unname(legacy@grid@cells.dim))
  expect_true(sp::identicalCRS(direct, legacy))
})


test_that("direct and legacy grids give identical KDE values and contours", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7302)
  xy <- cbind(stats::rnorm(30, 0, 250), stats::rnorm(30, 0, 250))
  pts <- sp::SpatialPointsDataFrame(
    xy, data.frame(ID = factor(rep("fish", nrow(xy)))),
    proj4string = sp::CRS(SRS_string = sf::st_crs(32629)$wkt)
  )
  r <- land_test_grid()
  direct_grid <- .spatRasterToSpatialPixels(r)
  legacy_grid <- methods::as(
    sf::as_Spatial(sf::st_as_sf(terra::as.points(r[[1]]))),
    "SpatialPixels"
  )

  direct <- adehabitatHR::kernelUD(pts, h = 500, grid = direct_grid)
  legacy <- adehabitatHR::kernelUD(pts, h = 500, grid = legacy_grid)

  expect_identical(direct[[1]]$ud, legacy[[1]]$ud)
  for (percent in c(50, 95)) {
    expect_identical(
      adehabitatHR::getverticeshr(direct, percent = percent, unin = "m", unout = "km2"),
      adehabitatHR::getverticeshr(legacy, percent = percent, unin = "m", unout = "km2")
    )
  }
})


test_that("sharing invariant estUD grid slots preserves values and copy-on-modify", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7303)
  xy <- rbind(cbind(stats::rnorm(25, -1000, 200), stats::rnorm(25, 0, 200)),
              cbind(stats::rnorm(25,  1000, 200), stats::rnorm(25, 0, 200)))
  pts <- sp::SpatialPointsDataFrame(
    xy, data.frame(ID = factor(rep(c("one", "two"), each = 25))),
    proj4string = sp::CRS(SRS_string = sf::st_crs(32629)$wkt)
  )
  original <- adehabitatHR::kernelUD(
    pts, h = 500, grid = .spatRasterToSpatialPixels(land_test_grid())
  )
  shared <- .shareUDGrid(original)

  expect_identical(shared, original)
  expect_identical(
    adehabitatHR::getverticeshr(shared, percent = 95, unin = "m", unout = "km2"),
    adehabitatHR::getverticeshr(original, percent = 95, unin = "m", unout = "km2")
  )
  expect_identical(unserialize(serialize(shared, NULL)), shared)

  first_before <- shared[[1]]@coords[1, 1]
  shared[[2]]@coords[1, 1] <- shared[[2]]@coords[1, 1] + 1
  expect_identical(shared[[1]]@coords[1, 1], first_before)
})


test_that("multiple contours reuse volume UDs without changing polygon objects", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7304)
  xy <- rbind(cbind(stats::rnorm(25, -1000, 200), stats::rnorm(25, 0, 200)),
              cbind(stats::rnorm(25,  1000, 200), stats::rnorm(25, 0, 200)))
  pts <- sp::SpatialPointsDataFrame(
    xy, data.frame(ID = factor(rep(c("one", "two"), each = 25))),
    proj4string = sp::CRS(SRS_string = sf::st_crs(32629)$wkt)
  )
  ud <- adehabitatHR::kernelUD(
    pts, h = 500, grid = .spatRasterToSpatialPixels(land_test_grid())
  )
  percents <- c(50, 95)
  reference <- lapply(percents, function(percent) {
    adehabitatHR::getverticeshr(ud, percent = percent, unin = "m", unout = "km2")
  })
  names(reference) <- paste0("K", percents)

  expect_identical(.extractKernelContours(ud, percents, grid.supplied = TRUE),
                   reference)
})


test_that("land masking retains st_intersects boundary and hole semantics", {
  skip_if_not_installed("sf")

  outer <- rbind(c(2.5, 2.5), c(8.5, 2.5), c(8.5, 8.5),
                 c(2.5, 8.5), c(2.5, 2.5))
  hole <- rbind(c(4.5, 4.5), c(4.5, 6.5), c(6.5, 6.5),
                c(6.5, 4.5), c(4.5, 4.5))
  island <- rbind(c(0.5, 0.5), c(1.5, 0.5), c(1.5, 1.5),
                  c(0.5, 1.5), c(0.5, 0.5))
  land <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(outer, hole)), sf::st_polygon(list(island)), crs = 3812
  ))
  grid <- terra::rast(terra::ext(0, 10, 0, 10), ncol = 10, nrow = 10,
                      crs = "EPSG:3812")
  xy <- terra::xyFromCell(grid, seq_len(terra::ncell(grid)))
  colnames(xy) <- c("x", "y")
  points <- sf::st_as_sf(as.data.frame(xy), coords = c("x", "y"), crs = 3812)
  reference <- as.logical(sf::st_intersects(
    points, sf::st_union(land), sparse = FALSE
  ))

  expect_identical(
    .landOverlapMask(xy, land, sf::st_crs(3812), chunk.size = 7L),
    reference
  )
  expect_identical(
    .landOverlapMaskSf(xy, sf::st_union(land), sf::st_crs(3812), chunk.size = 7L),
    reference
  )
})


test_that("chunked land clipping is numerically identical to the former implementation", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  ud <- land_test_ud()
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-5000, -5000), c(0, -5000), c(0, 5000),
    c(-5000, 5000), c(-5000, -5000)
  ))), crs = 32629))

  reference <- legacy_subtract_land(ud, land, sf::st_crs(32629))
  result <- .subtractLand(ud, land, sf::st_crs(32629), verbose = FALSE)

  expect_identical(result[[1]]$ud, reference[[1]]$ud)
  expect_identical(sum(result[[1]]$ud), sum(reference[[1]]$ud))
  for (percent in c(50, 95)) {
    expect_identical(
      adehabitatHR::getverticeshr(result, percent = percent, unin = "m", unout = "km2"),
      adehabitatHR::getverticeshr(reference, percent = percent, unin = "m", unout = "km2")
    )
  }
})


test_that("a fully land-masked UD fails before contour extraction with its ID", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  ud <- land_test_ud()
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-6000, -6000), c(6000, -6000), c(6000, 6000),
    c(-6000, 6000), c(-6000, -6000)
  ))), crs = 32629))

  err <- tryCatch(
    .subtractLand(ud, land, sf::st_crs(32629), verbose = FALSE),
    error = identity
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "Land clipping removed all KDE density")
  expect_match(conditionMessage(err), "fish")
  expect_match(conditionMessage(err), "enlarging 'spatial.grid' will not resolve", fixed = TRUE)
  expect_false(grepl("too small", conditionMessage(err), fixed = TRUE))
})


test_that("calculateUDs exposes fully masked UDs as a non-retryable input error", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7305)
  d <- data.frame(
    ID = factor(rep("fish", 12)),
    timebin = as.POSIXct("2021-01-01", tz = "UTC") + seq_len(12) * 3600,
    x = stats::rnorm(12, 0, 100), y = stats::rnorm(12, 0, 100)
  )
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-6000, -6000), c(6000, -6000), c(6000, 6000),
    c(-6000, 6000), c(-6000, -6000)
  ))), crs = 32629))

  err <- tryCatch(
    calculateUDs(
      d, id.col = "ID", timebin.col = "timebin", lon.col = "x", lat.col = "y",
      method = "kde", bandwidth = 500, contour.percent = c(50, 95),
      epsg.code = 32629, spatial.grid = land_test_grid(), land.shape = land,
      verbose = FALSE
    ),
    error = identity
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "Land clipping removed all KDE density")
  expect_false(grepl("too small", conditionMessage(err), fixed = TRUE))
})
