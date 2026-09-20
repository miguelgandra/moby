# Spatial regression harness (Phase 0 of the raster->terra / gdistance migration).
#
# Freezes deterministic fixtures (a strait: west mainland + central island + east mainland,
# with detections in the channels) and GOLDEN outputs captured from the CURRENT raster/gdistance
# implementation, then re-runs the spatial functions and asserts the results stay within
# TOLERANCE BANDS (not byte-identity). The bands exist precisely because a backend swap
# (raster->terra, gdistance->terra+igraph) legitimately perturbs numbers at the grid-cell level:
# focal/modal tie-breaking, projection alignment, nearest-cell vs nearest-boundary relocation,
# and least-cost path discretisation. A regression is a change LARGER than these bands.
#
# Regenerate fixtures/goldens with `Rscript data-raw/build-spatial-fixtures.R` ONLY when the
# intended behaviour changes; never to paper over an unexplained drift.

# --- tolerance bands (with rationale) --------------------------------------------------------
TOL_STEP_ABS  <- 100    # one grid cell (grid.resolution = 100 m): least-cost path discretisation
TOL_STEP_REL  <- 0.05   # 5%: geoCorrection / router differences on longer detours
TOL_COORD_M   <- 150    # relocation target may shift ~1-2 raster cells between backends
TOL_LAND_ABS  <- 100    # one raster cell (grid.resolution = 100 m)
TOL_LAND_REL  <- 0.05

fx     <- readRDS(test_path("_spatial", "fixtures.rds"))
golden <- readRDS(test_path("_spatial", "golden.rds"))
EPSG   <- fx$epsg
gres   <- fx$params$grid.resolution
mdir   <- fx$params$mov.directions

# helper: element-wise "within tolerance band" for two numeric vectors, NA-structure-aware
expect_within_band <- function(current, golden, abs_tol, rel_tol, label) {
  expect_equal(is.na(current), is.na(golden),
               info = paste0(label, ": NA structure changed"))
  ok <- !is.na(golden) & !is.na(current)
  if (any(ok)) {
    band <- pmax(abs_tol, rel_tol * abs(golden[ok]))
    diff <- abs(current[ok] - golden[ok])
    expect_true(all(diff <= band),
                info = paste0(label, ": max diff ", round(max(diff), 2),
                              " exceeds band ", round(max(band), 2),
                              " (current=", paste(round(current[ok], 1), collapse = ","),
                              " | golden=", paste(round(golden[ok], 1), collapse = ","), ")"))
  }
}


test_that("calculateStepDistances (gdistance least-cost) matches frozen golden within band", {
  sd <- calculateStepDistances(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                               grid.resolution = gres, mov.directions = mdir,
                               id.col = "ID", lon.col = "lon", lat.col = "lat",
                               verbose = FALSE)
  expect_within_band(sd$dist_m, golden$step_dist_m,
                     TOL_STEP_ABS, TOL_STEP_REL, "step dist_m")
})


test_that("least-cost path invariants hold (>= straight line; land crossings detour)", {
  sd <- calculateStepDistances(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                               grid.resolution = gres, mov.directions = mdir,
                               id.col = "ID", lon.col = "lon", lat.col = "lat",
                               verbose = FALSE)
  # great-circle straight-line distance between consecutive detections
  pts <- as.matrix(fx$track[, c("lon", "lat")])
  straight <- geosphere::distGeo(pts[-nrow(pts), ], pts[-1, ])   # length n-1
  d <- sd$dist_m[-length(sd$dist_m)]                             # drop trailing NA
  # a least-cost path is never (meaningfully) shorter than the straight line
  expect_true(all(d >= straight - TOL_STEP_ABS),
              info = "some in-water distance is shorter than the great-circle line")
  # segments whose straight line crosses the central island must detour (strictly longer)
  crossing <- c(1L, 3L, 5L)  # W<->E transitions in the fixture track (see build script)
  expect_true(all(d[crossing] > straight[crossing]),
              info = "a land-crossing segment did not detour (least-cost path never engaged)")
})


test_that("correctPositions RASTER path matches frozen golden (relocated set + coords)", {
  land_rast <- terra::rast(test_path("_spatial", "land_raster.tif"))
  cp <- suppressWarnings(suppressMessages(
    correctPositions(fx$onland, spatial.layer = land_rast, raster.type = fx$params$raster.type,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat")))
  # which points are on land is structural — must match exactly
  expect_identical(attr(cp, "points.relocated"), golden$correct_raster$relocated)
  # relocated coordinates within a cell or two (geodesic distance golden<->current)
  gd <- geosphere::distGeo(cbind(golden$correct_raster$lon, golden$correct_raster$lat),
                           cbind(cp$data$lon, cp$data$lat))
  expect_true(all(gd <= TOL_COORD_M),
              info = paste0("raster relocation moved > ", TOL_COORD_M, " m from golden (max ",
                            round(max(gd), 1), " m)"))
})


test_that("correctPositions SF path matches frozen golden (relocated set + coords)", {
  cp <- suppressWarnings(suppressMessages(
    correctPositions(fx$onland, spatial.layer = fx$land_sf,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat")))
  expect_identical(attr(cp, "points.relocated"), golden$correct_sf$relocated)
  gd <- geosphere::distGeo(cbind(golden$correct_sf$lon, golden$correct_sf$lat),
                           cbind(cp$data$lon, cp$data$lat))
  expect_true(all(gd <= TOL_COORD_M),
              info = paste0("sf relocation moved > ", TOL_COORD_M, " m from golden (max ",
                            round(max(gd), 1), " m)"))
})


test_that("polygon relocation clears the calculateUDs land pre-flight check", {
  cp <- suppressWarnings(suppressMessages(
    correctPositions(fx$onland, spatial.layer = fx$land_sf,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                     verbose = FALSE)))
  land <- sf::st_transform(fx$land_sf, EPSG)
  corrected <- sf::st_transform(
    sf::st_as_sf(cp$data, coords = c("lon", "lat"), crs = 4326), EPSG)

  expect_identical(attr(cp, "points.relocated"), 3L)
  expect_true(all(lengths(sf::st_intersects(corrected, land)) == 0L))
  expect_no_warning(.checkPositionsAgainstLand(corrected, "ID", land))

  projected <- fx$onland
  xy <- sf::st_coordinates(sf::st_transform(
    sf::st_as_sf(projected, coords = c("lon", "lat"), crs = 4326), EPSG))
  projected$lon <- xy[, "X"]
  projected$lat <- xy[, "Y"]
  land_projected <- sf::st_transform(fx$land_sf, EPSG)
  cp_projected <- suppressWarnings(suppressMessages(
    correctPositions(projected, spatial.layer = land_projected,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                     verbose = FALSE)))
  corrected_projected <- sf::st_as_sf(cp_projected$data,
                                      coords = c("lon", "lat"), crs = EPSG)
  expect_true(all(lengths(sf::st_intersects(corrected_projected, land_projected)) == 0L))
})


test_that("water-side offsets handle exact boundaries, holes and adjacent land", {
  rect <- function(x0, x1, y0, y1) sf::st_polygon(list(rbind(
    c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
  land <- sf::st_sf(geometry = sf::st_sfc(
    rect(500000, 500001, 4256000, 4256001),
    rect(500001.002, 500002, 4256000, 4256001), crs = EPSG))
  origin <- sf::st_sfc(sf::st_point(c(500000.999, 4256000.5)), crs = EPSG)
  boundary <- sf::st_sfc(sf::st_point(c(500001, 4256000.5)), crs = EPSG)
  moved <- moby:::.nudgeToWater(origin, boundary, land, 1, TRUE)
  expect_false(is.null(moved))
  expect_equal(lengths(sf::st_intersects(moved, land)), 0L)
  expect_equal(lengths(sf::st_intersects(
    sf::st_transform(sf::st_transform(moved, 4326), EPSG), land)), 0L)

  corner <- sf::st_sfc(sf::st_point(c(500001, 4256001)), crs = EPSG)
  moved_corner <- moby:::.nudgeToWater(corner, corner, land, 1, TRUE)
  expect_false(is.null(moved_corner))
  expect_equal(lengths(sf::st_intersects(moved_corner, land)), 0L)

  outer <- rbind(c(500010,4256010), c(500020,4256010),
                 c(500020,4256020), c(500010,4256020), c(500010,4256010))
  hole <- rbind(c(500013,4256013), c(500013,4256017),
                c(500017,4256017), c(500017,4256013), c(500013,4256013))
  island_with_hole <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(outer, hole)), crs = EPSG))
  inside_land <- sf::st_sfc(sf::st_point(c(500012, 4256015)), crs = EPSG)
  hole_edge <- sf::st_sfc(sf::st_point(c(500013, 4256015)), crs = EPSG)
  moved_into_hole <- moby:::.nudgeToWater(inside_land, hole_edge, island_with_hole, 1, TRUE)
  expect_false(is.null(moved_into_hole))
  expect_equal(lengths(sf::st_intersects(moved_into_hole, island_with_hole)), 0L)

  # A radius smaller than the minimum offset must not return a boundary point as a success.
  expect_null(moby:::.nudgeToWater(corner, corner, land, 0.0000001, TRUE))
})


test_that("parallel polygon relocation matches the single-core result", {
  skip_on_cran()
  skip_if_not_installed("foreach")
  skip_if_not_installed("doSNOW")

  serial <- suppressWarnings(suppressMessages(
    correctPositions(fx$onland, spatial.layer = fx$land_sf,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                     cores = 1, verbose = FALSE)))
  parallel <- suppressWarnings(suppressMessages(
    correctPositions(fx$onland, spatial.layer = fx$land_sf,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                     cores = 2, verbose = FALSE)))
  expect_equal(parallel$data[, c("lon", "lat")], serial$data[, c("lon", "lat")])
})


test_that("correctPositions searches the requested radius below 10 km", {
  land_rast <- terra::rast(test_path("_spatial", "land_raster.tif"))
  layers <- list(sf = fx$land_sf, raster = land_rast)

  for (layer in layers) {
    reference <- suppressWarnings(suppressMessages(
      correctPositions(fx$onland, spatial.layer = layer, raster.type = fx$params$raster.type,
                       epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                       max.distance.km = 50, verbose = FALSE)))
    cp <- suppressWarnings(suppressMessages(
      correctPositions(fx$onland, spatial.layer = layer, raster.type = fx$params$raster.type,
                       epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                       max.distance.km = 5, verbose = FALSE)))

    expect_identical(attr(cp, "points.relocated"), 3L)
    expect_identical(attr(cp, "points.skipped"), 0L)
    expect_true(all(cp$summary$distance_m < 5000))
    expect_equal(cp$data[, c("lon", "lat")], reference$data[, c("lon", "lat")])
  }
})


test_that("correctPositions converts a single relocated point back to geographic coordinates", {
  one_on_land <- fx$onland[c(1, 2), , drop = FALSE]
  # Place the first point fractionally inside the island boundary. It is classified as on land,
  # but its sub-millimetre move to the boundary is reported as 0 m after summary rounding.
  one_on_land$lon[1] <- -9.003999999
  one_on_land$lat[1] <- 38.45

  cp <- suppressWarnings(suppressMessages(
    correctPositions(one_on_land, spatial.layer = fx$land_sf,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                     max.distance.km = 5, verbose = FALSE)))

  expect_identical(attr(cp, "points.relocated"), 1L)
  expect_identical(attr(cp, "points.skipped"), 0L)
  expect_equal(nrow(cp$summary), 1L)
  expect_equal(cp$summary$distance_m, 0)
  expect_true(all(abs(cp$data$lon) <= 180))
  expect_true(all(abs(cp$data$lat) <= 90))
  expect_equal(cp$data$lon[2], one_on_land$lon[2])
  expect_equal(cp$data$lat[2], one_on_land$lat[2])
  corrected <- sf::st_transform(
    sf::st_as_sf(cp$data, coords = c("lon", "lat"), crs = 4326), EPSG)
  expect_true(all(lengths(sf::st_intersects(corrected, sf::st_transform(fx$land_sf, EPSG))) == 0L))
})


test_that("correctPositions includes the exact maximum and final partial increment", {
  point <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(0, 0)), crs = EPSG))
  target <- sf::st_sf(geometry = sf::st_sfc(sf::st_point(c(12000, 0)), crs = EPSG))

  expect_null(moby:::.setSearchRegion(point, target, max.distance.km = 11.999))
  expect_s3_class(moby:::.setSearchRegion(point, target, max.distance.km = 12), "sf")
  expect_s3_class(moby:::.setSearchRegion(point, target, max.distance.km = 15), "sf")
})


test_that("correctPositions enforces and validates the maximum radius", {
  cp <- suppressWarnings(suppressMessages(
    correctPositions(fx$onland, spatial.layer = fx$land_sf,
                     epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                     max.distance.km = 0.25, verbose = FALSE)))

  expect_identical(attr(cp, "points.relocated"), 0L)
  expect_identical(attr(cp, "points.skipped"), 3L)
  expect_true(all(is.na(cp$summary$distance_m)))

  for (bad_radius in list(0, -1, NA_real_, Inf, c(1, 2), "5")) {
    expect_error(
      correctPositions(fx$onland, spatial.layer = fx$land_sf,
                       epsg.code = EPSG, lon.col = "lon", lat.col = "lat",
                       max.distance.km = bad_radius, verbose = FALSE),
      "positive, finite",
      fixed = TRUE
    )
  }
})


test_that("calculateLandDists (terra) matches frozen golden within band", {
  ld <- calculateLandDists(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                           id.col = "ID", lon.col = "lon", lat.col = "lat", verbose = FALSE)
  expect_within_band(ld$land_dist, golden$land_dist,
                     TOL_LAND_ABS, TOL_LAND_REL, "land_dist")
})


# --- Phase 2 (terra + igraph productionisation) regression guards --------------------------

test_that("calculateStepDistances() returns only the data + trajectories (no cost graph)", {
  sd <- calculateStepDistances(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                               id.col = "ID", lon.col = "lon", lat.col = "lat", verbose = FALSE)
  expect_s3_class(sd, "data.frame")            # the data itself, never a list
  expect_true("dist_m" %in% names(sd))
  expect_false(is.null(getTrajectories(sd)))   # trajectories are still attached
  # the least-cost graph is an internal implementation detail and is never handed to users
  expect_null(attr(sd, "transition_layer"))
  expect_null(attr(sd, "cost.graph"))
})

test_that("the cost surface is a serialisable mobyCostGraph (no live SpatRaster)", {
  # the graph is internal-only: reach it through the engine, not the public return value
  cg <- moby:::.stepDistances(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                              id.col = "ID", lon.col = "lon", lat.col = "lat",
                              verbose = FALSE)$cost.graph
  expect_s3_class(cg, "mobyCostGraph")
  expect_true(inherits(cg$graph, "igraph"))
  # the grid descriptor is plain R (round-trips through serialize()); no SpatRaster pointer
  expect_type(cg$grid, "list")
  expect_true(all(c("xmin", "ymax", "resx", "resy", "nr", "nc") %in% names(cg$grid)))
  expect_identical(unserialize(serialize(cg$grid, NULL)), cg$grid)
})

test_that("internal cost-graph reuse reproduces a fresh build exactly", {
  # filterDetections() builds one full-extent graph and reuses it across many calls; reuse must be exact
  s1 <- moby:::.stepDistances(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                              id.col = "ID", lon.col = "lon", lat.col = "lat", verbose = FALSE)
  s2 <- moby:::.stepDistances(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                              id.col = "ID", lon.col = "lon", lat.col = "lat",
                              cost.graph = s1$cost.graph, verbose = FALSE)
  expect_equal(s1$data$dist_m, s2$data$dist_m)
})

test_that("ALL land polygons are impassable (multi-feature bug fixed)", {
  # two separate barriers; a west->east segment must detour around the SECOND polygon too,
  # so the routed distance exceeds the straight-line distance (the old code left it passable).
  rect <- function(x0, x1, y0, y1)
    sf::st_polygon(list(rbind(c(x0,y0), c(x1,y0), c(x1,y1), c(x0,y1), c(x0,y0))))
  land2 <- sf::st_sf(id = 1:2, geometry = sf::st_sfc(
    rect(-9.006, -9.004, 38.43, 38.47), rect(-8.996, -8.994, 38.43, 38.47), crs = 4326))
  trk <- data.frame(ID = "A01", datetime = as.POSIXct("2023-05-01", tz = "UTC") + (0:1) * 3600,
                    lon = c(-9.02, -8.98), lat = c(38.449, 38.451))
  d <- suppressWarnings(calculateStepDistances(trk, land.shape = land2, epsg.code = EPSG,
                        id.col = "ID", lon.col = "lon", lat.col = "lat", verbose = FALSE)$dist_m[1])
  straight <- geosphere::distVincentyEllipsoid(c(-9.02, 38.449), c(-8.98, 38.451))
  expect_gt(d, straight + 1)   # detoured around both barriers
})

test_that("a projected epsg.code is still required (geographic CRS rejected)", {
  expect_error(
    calculateStepDistances(fx$track, land.shape = fx$land_sf, epsg.code = 4326,
                           id.col = "ID", lon.col = "lon", lat.col = "lat", verbose = FALSE))
})

test_that("an on-land endpoint routes to a finite distance without error", {
  onland_trk <- data.frame(ID = "A01", datetime = as.POSIXct("2023-05-01", tz = "UTC") + (0:1) * 3600,
                           lon = c(-9.000, -8.988), lat = c(38.450, 38.452))  # first point on the island
  d <- suppressWarnings(calculateStepDistances(onland_trk, land.shape = fx$land_sf, epsg.code = EPSG,
                        id.col = "ID", lon.col = "lon", lat.col = "lat", verbose = FALSE)$dist_m)
  expect_true(is.finite(d[1]))
  expect_true(is.na(d[2]))     # last position always NA
})

test_that("parallel routing (cores = 2) matches single-core exactly", {
  skip_on_cran()
  skip_if_not_installed("doSNOW"); skip_if_not_installed("foreach"); skip_if_not_installed("parallel")
  tr2 <- rbind(transform(fx$track, ID = "A01"),
               transform(fx$track, ID = "B02", lat = fx$track$lat + 0.002))
  d1 <- calculateStepDistances(tr2, land.shape = fx$land_sf, epsg.code = EPSG,
                               id.col = "ID", lon.col = "lon", lat.col = "lat", cores = 1, verbose = FALSE)$dist_m
  d2 <- calculateStepDistances(tr2, land.shape = fx$land_sf, epsg.code = EPSG,
                               id.col = "ID", lon.col = "lon", lat.col = "lat", cores = 2, verbose = FALSE)$dist_m
  expect_equal(d1, d2)
})
