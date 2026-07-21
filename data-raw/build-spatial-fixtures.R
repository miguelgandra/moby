#!/usr/bin/env Rscript
## Build deterministic spatial fixtures + golden outputs for the moby spatial regression harness.
## Captures golden outputs from the CURRENT raster/gdistance implementation so a future
## terra/igraph migration can be validated against tolerance bands (not byte-identity).

suppressWarnings(suppressMessages(pkgload::load_all(".", quiet = TRUE)))
set.seed(1)  # only for any incidental use; geometry below is fully explicit

EPSG <- 32629  # UTM 29N (metric), matches the toy-data region off Portugal

## ---- 1. Land polygons (sf, WGS84): a strait ----------------------------------------------
## West mainland + central island "wall" + east mainland. The land extent encloses all
## detections (so land-distance raster covers them), while straight west<->east segments
## still cross the central island, forcing least-cost detours.
rect <- function(xmin, xmax, ymin = 38.420, ymax = 38.480) {
  sf::st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin),
                            c(xmax, ymax), c(xmin, ymax), c(xmin, ymin))))
}
land_sf <- sf::st_sf(
  id = 1:3,
  part = c("west_mainland", "island", "east_mainland"),
  geometry = sf::st_sfc(rect(-9.030, -9.020),   # west mainland
                        rect(-9.004, -8.996),    # central island (the barrier)
                        rect(-8.984, -8.975),    # east mainland
                        crs = 4326))

## ---- 2a. Detections for step-distance / land-distance (all in water) ----------------------
## Single individual weaving west<->east; several consecutive pairs cross the wall.
track <- data.frame(
  ID  = "A01",
  datetime = as.POSIXct("2023-05-01 00:00", tz = "UTC") + (0:7) * 3600,
  lon = c(-9.010, -8.990, -8.988, -9.012, -9.011, -8.991, -8.993, -8.989),
  lat = c( 38.448, 38.452, 38.460, 38.458, 38.440, 38.438, 38.445, 38.455),
  stringsAsFactors = FALSE
)

## ---- 2b. Detections for correctPositions (mix of on-land + water) -------------------------
## On-land points must be relocated; water points must stay put.
onland <- data.frame(
  ID = "A01",
  datetime = as.POSIXct("2023-05-01 00:00", tz = "UTC") + (0:5) * 3600,
  lon = c(-9.000, -9.010, -8.999, -8.990, -9.001, -9.011),  # #1,#3,#5 on the wall
  lat = c( 38.450, 38.448, 38.455, 38.452, 38.445, 38.458),
  stringsAsFactors = FALSE
)

## ---- 3. Land raster (RasterLayer, PROJECTED EPSG:32629): land = 1, water = NA -------------
## Built directly in the analysis CRS so correctPositions receives it already projected and does
## NOT reproject it. (The raster::projectRaster reprojection path is both a migration hotspot and
## currently fragile: on a clean session, a geographic raster + projected epsg.code errors inside
## raster:::.computeRes -> terra::xmax. A projected raster is also what users typically supply.)
## ~60 m cells over a padded bbox; exercises the trim/focal/modal/extract/rasterToPoints path.
bb        <- sf::st_bbox(land_sf)                    # geographic bbox, for the summary print only
land_proj <- sf::st_transform(land_sf, EPSG)
bbp       <- sf::st_bbox(land_proj); padm <- 2000    # metres
tmpl <- raster::raster(xmn = bbp["xmin"] - padm, xmx = bbp["xmax"] + padm,
                       ymn = bbp["ymin"] - padm, ymx = bbp["ymax"] + padm,
                       resolution = 60, crs = "EPSG:32629")
land_rast <- raster::rasterize(sf::as_Spatial(land_proj), tmpl, field = 1, background = NA)

cat("fixtures built:\n")
cat("  land_sf   :", nrow(land_sf), "polygon | bbox",
    paste(round(as.numeric(bb), 4), collapse = ", "), "\n")
cat("  track     :", nrow(track), "detections (all water)\n")
cat("  onland    :", nrow(onland), "detections (on-land + water)\n")
cat("  land_rast :", raster::ncell(land_rast), "cells |",
    sum(!is.na(raster::values(land_rast))), "land cells\n\n")

## ---- 4. Capture GOLDEN outputs from the current implementation ----------------------------
golden <- list()

## (a) calculateStepDistances — the gdistance least-cost path (dist_m per detection)
sd <- calculateStepDistances(track, land.shape = land_sf, epsg.code = EPSG,
                             grid.resolution = 100, mov.directions = 16,
                             id.col = "ID", lon.col = "lon", lat.col = "lat",
                             verbose = FALSE)
golden$step_dist_m <- sd$dist_m
cat("(a) calculateStepDistances dist_m:\n    ",
    paste(ifelse(is.na(sd$dist_m), "NA", round(sd$dist_m, 1)), collapse = ", "), "\n")

## (b) correctPositions — RASTER path (the heavy raster:: code)
cp_r <- correctPositions(onland, spatial.layer = land_rast, raster.type = "land",
                         epsg.code = EPSG, lon.col = "lon", lat.col = "lat")
golden$correct_raster <- list(
  relocated   = attr(cp_r, "points.relocated"),
  lon         = cp_r$data$lon,
  lat         = cp_r$data$lat
)
cat("(b) correctPositions [raster]: relocated =", attr(cp_r, "points.relocated"), "of", nrow(onland), "\n")

## (c) correctPositions — SF path (pure-sf branch, for cross-check)
cp_s <- correctPositions(onland, spatial.layer = land_sf,
                         epsg.code = EPSG, lon.col = "lon", lat.col = "lat")
golden$correct_sf <- list(
  relocated = attr(cp_s, "points.relocated"),
  lon       = cp_s$data$lon,
  lat       = cp_s$data$lat
)
cat("(c) correctPositions [sf]    : relocated =", attr(cp_s, "points.relocated"), "of", nrow(onland), "\n")

## (d) calculateLandDists — already terra (reference idiom); distance to nearest land
ld <- calculateLandDists(track, land.shape = land_sf, epsg.code = EPSG,
                         id.col = "ID", lon.col = "lon", lat.col = "lat")
ld_col <- setdiff(names(ld), names(track))
golden$land_dist <- ld[[tail(ld_col, 1)]]
cat("(d) calculateLandDists (", tail(ld_col,1), "):\n    ",
    paste(round(golden$land_dist, 1), collapse = ", "), "\n")

## ---- 5. Save fixtures + goldens ----------------------------------------------------------
## Fixtures are the FROZEN inputs; goldens are the FROZEN current-implementation outputs.
## The land raster is stored as a GeoTIFF so the correctPositions INPUT is identical whether
## later read by raster::raster() or terra::rast() — isolating the function's behaviour from
## the rasterize step. sf/data.frame fixtures go in an .rds. Nothing here is regenerated by
## the test; the test loads these and compares within tolerance bands.
outdir <- "tests/testthat/_spatial"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
fixtures <- list(epsg = EPSG, land_sf = land_sf, track = track, onland = onland,
                 params = list(grid.resolution = 100, mov.directions = 16, raster.type = "land"))
saveRDS(fixtures, file.path(outdir, "fixtures.rds"))
raster::writeRaster(land_rast, file.path(outdir, "land_raster.tif"), overwrite = TRUE)
saveRDS(golden, file.path(outdir, "golden.rds"))
cat("\nsaved fixtures + goldens to", outdir, "\n")
cat("  files:", paste(list.files(outdir), collapse = ", "), "\n")
