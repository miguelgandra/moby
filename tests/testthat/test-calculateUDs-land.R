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
  expect_length(result$empty.ids, 0L)
  result <- result$uds

  expect_identical(result[[1]]$ud, reference[[1]]$ud)
  expect_identical(sum(result[[1]]$ud), sum(reference[[1]]$ud))
  for (percent in c(50, 95)) {
    expect_identical(
      adehabitatHR::getverticeshr(result, percent = percent, unin = "m", unout = "km2"),
      adehabitatHR::getverticeshr(reference, percent = percent, unin = "m", unout = "km2")
    )
  }
})


test_that("land clipping identifies a fully masked UD without passing it to contours", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  ud <- land_test_ud()
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-6000, -6000), c(6000, -6000), c(6000, 6000),
    c(-6000, 6000), c(-6000, -6000)
  ))), crs = 32629))

  result <- .subtractLand(ud, land, sf::st_crs(32629), verbose = FALSE)
  expect_identical(result$empty.ids, "fish")
  expect_true(all(result$uds[["fish"]]$ud == 0))
})


test_that("calculateUDs returns zero areas and empty contours for fully masked UDs", {
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

  warnings <- character()
  result <- withCallingHandlers(
    calculateUDs(
      d, id.col = "ID", timebin.col = "timebin", lon.col = "x", lat.col = "y",
      method = "kde", bandwidth = 500, contour.percent = c(50, 95),
      epsg.code = 32629, spatial.grid = land_test_grid(), land.shape = land,
      verbose = FALSE
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("Pre-flight land check", warnings, fixed = TRUE)))
  expect_true(any(grepl("Land clipping removed all KDE density", warnings, fixed = TRUE)))
  expect_s3_class(result$ud, "estUDm")
  expect_length(result$ud, 0L)
  expect_identical(attr(result, "on.empty"), "warn")
  expect_identical(attr(result, "empty.ids"), "fish")
  expect_identical(result$summary_table[["UD 50% (Km2)"]], "0.00")
  expect_identical(result$summary_table[["UD 95% (Km2)"]], "0.00")
  expect_true(all(sf::st_is_empty(result$K50)))
  expect_true(all(sf::st_is_empty(result$K95)))
  expect_identical(result$K50$area, 0)
  expect_identical(result$K95$area, 0)

  overlap <- calculateUDOverlap(result, verbose = FALSE)
  expect_equal(nrow(overlap), 0L)

  map_file <- tempfile(fileext = ".png")
  on.exit(unlink(map_file), add = TRUE)
  expect_no_error(plotMaps(
    d, uds = result, id.col = "ID", lon.col = "x", lat.col = "y",
    epsg.code = 32629, coastline = FALSE, verbose = FALSE, file = map_file
  ))
  expect_true(file.exists(map_file))
})


test_that("a fully masked individual does not discard valid KDE results", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7307)
  d <- rbind(
    data.frame(ID = "masked", x = stats::rnorm(20, -2500, 80),
               y = stats::rnorm(20, 0, 80)),
    data.frame(ID = "water", x = stats::rnorm(20, 2500, 80),
               y = stats::rnorm(20, 0, 80))
  )
  d$ID <- factor(d$ID)
  d$timebin <- as.POSIXct("2021-01-01", tz = "UTC") + seq_len(nrow(d)) * 3600
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-5000, -5000), c(0, -5000), c(0, 5000),
    c(-5000, 5000), c(-5000, -5000)
  ))), crs = 32629))

  result <- suppressWarnings(calculateUDs(
    d, id.col = "ID", timebin.col = "timebin", lon.col = "x", lat.col = "y",
    method = "kde", bandwidth = 300, contour.percent = c(50, 95),
    epsg.code = 32629, spatial.grid = land_test_grid(), land.shape = land,
    verbose = FALSE
  ))
  reference <- calculateUDs(
    d, id.col = "ID", timebin.col = "timebin", lon.col = "x", lat.col = "y",
    method = "kde", bandwidth = 300, contour.percent = c(50, 95),
    epsg.code = 32629, spatial.grid = land_test_grid(), verbose = FALSE
  )

  expect_identical(names(result$ud), "water")
  expect_identical(result$ud[["water"]]$ud, reference$ud[["water"]]$ud)
  expect_identical(attr(result, "empty.ids"), "masked")
  masked <- result$summary_table[result$summary_table$ID == "masked", ]
  water <- result$summary_table[result$summary_table$ID == "water", ]
  reference_water <- reference$summary_table[reference$summary_table$ID == "water", ]
  expect_identical(masked[["UD 50% (Km2)"]], "0.00")
  expect_identical(masked[["UD 95% (Km2)"]], "0.00")
  expect_identical(water[["UD 50% (Km2)"]], reference_water[["UD 50% (Km2)"]])
  expect_identical(water[["UD 95% (Km2)"]], reference_water[["UD 95% (Km2)"]])
  expect_true(sf::st_is_empty(result$K95[result$K95$id == "masked", ]))
  expect_false(sf::st_is_empty(result$K95[result$K95$id == "water", ]))
})


test_that("grouped KDE retains an all-empty group alongside valid groups", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7308)
  d <- rbind(
    data.frame(ID = "masked", group = "land", x = stats::rnorm(20, -2500, 80),
               y = stats::rnorm(20, 0, 80)),
    data.frame(ID = "water", group = "water", x = stats::rnorm(20, 2500, 80),
               y = stats::rnorm(20, 0, 80))
  )
  d$ID <- factor(d$ID)
  d$group <- factor(d$group)
  d$timebin <- as.POSIXct("2021-01-01", tz = "UTC") + seq_len(nrow(d)) * 3600
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-5000, -5000), c(0, -5000), c(0, 5000),
    c(-5000, 5000), c(-5000, -5000)
  ))), crs = 32629))

  result <- suppressWarnings(calculateUDs(
    d, id.col = "ID", timebin.col = "timebin", lon.col = "x", lat.col = "y",
    method = "kde", bandwidth = 300, contour.percent = c(50, 95),
    epsg.code = 32629, spatial.grid = land_test_grid(), land.shape = land,
    subset = "group", verbose = FALSE
  ))

  expect_length(result$ud[["land"]], 0L)
  expect_identical(names(result$ud[["water"]]), "water")
  expect_identical(attr(result, "empty.ids"), "masked [land]")
  masked <- result$summary_table[result$summary_table$ID == "masked", ]
  expect_identical(masked[["UD 50% (Km2)"]], "0.00")
  expect_true(sf::st_is_empty(result$K95[result$K95$id == "masked", ]))
})


test_that("on.empty = 'error' retains strict non-retryable behaviour", {
  skip_if_not_installed("adehabitatHR")
  skip_if_not_installed("sp")
  skip_if_not_installed("terra")

  set.seed(7306)
  d <- data.frame(
    ID = factor(rep("fish", 12)),
    timebin = as.POSIXct("2021-01-01", tz = "UTC") + seq_len(12) * 3600,
    x = stats::rnorm(12, 0, 100), y = stats::rnorm(12, 0, 100)
  )
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-6000, -6000), c(6000, -6000), c(6000, 6000),
    c(-6000, 6000), c(-6000, -6000)
  ))), crs = 32629))

  err <- suppressWarnings(tryCatch(
    calculateUDs(
      d, id.col = "ID", timebin.col = "timebin", lon.col = "x", lat.col = "y",
      method = "kde", bandwidth = 500, contour.percent = c(50, 95),
      epsg.code = 32629, spatial.grid = land_test_grid(), land.shape = land,
      on.empty = "error", verbose = FALSE
    ),
    error = identity
  ))

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "Land clipping removed all KDE density")
  expect_match(conditionMessage(err), "fish")
  expect_match(conditionMessage(err), "enlarging 'spatial.grid' will not resolve", fixed = TRUE)
  expect_false(grepl("too small", conditionMessage(err), fixed = TRUE))
})


test_that("pre-flight land checks distinguish expected water from actionable mismatches", {
  coords <- sf::st_as_sf(
    data.frame(ID = factor(c("water", "land")), x = c(5, 15), y = c(5, 5)),
    coords = c("x", "y"), crs = 32629
  )
  land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(10, 0), c(20, 0), c(20, 10), c(10, 10), c(10, 0)
  ))), crs = 32629))

  expect_warning(.checkPositionsAgainstLand(coords, "ID", land),
                 "1 position of 2.*1 individual", perl = TRUE)

  water <- coords[1, ]
  outer <- rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0))
  hole <- rbind(c(2, 2), c(2, 8), c(8, 8), c(8, 2), c(2, 2))
  enclosing_land <- sf::st_sf(
    geometry = sf::st_sfc(sf::st_polygon(list(outer, hole)), crs = 32629)
  )
  expect_no_warning(.checkPositionsAgainstLand(water, "ID", enclosing_land))

  distant_land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(100, 100), c(110, 100), c(110, 110), c(100, 110), c(100, 100)
  ))), crs = 32629))
  expect_warning(.checkPositionsAgainstLand(coords, "ID", distant_land),
                 "position extent does not overlap")
})
