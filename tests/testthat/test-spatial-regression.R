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


test_that("calculateLandDists (terra) matches frozen golden within band", {
  ld <- calculateLandDists(fx$track, land.shape = fx$land_sf, epsg.code = EPSG,
                           id.col = "ID", lon.col = "lon", lat.col = "lat")
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
