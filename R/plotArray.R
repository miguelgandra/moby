#######################################################################################################
## Plot the receiver array (deployment map) ##########################################################
#######################################################################################################

#' Map the receiver array
#'
#' @description Draws the spatial configuration of a receiver array from a deployment log: one point
#' per station, optionally over a coastline and/or background raster (e.g. bathymetry), with a metric
#' scale bar. It is a pre-analysis quality-control tool - use it to inspect the array layout, spot
#' implausible receiver coordinates, and judge spacing and overall coverage before running any
#' downstream analysis.
#'
#' Semi-transparent circles of a user-defined radius (`detection.range`) can be drawn around each
#' receiver to give a quick visual sense of nominal coverage and potential detection gaps, and the
#' array's convex hull can be outlined as a footprint. Stations can be coloured by any metadata column
#' (`color.by`) or by their deployment status at a chosen date (`status.at`), and - if a detection
#' dataset is supplied - stations that never logged a detection are highlighted.
#'
#' This is the spatial companion to \code{\link{checkDeployments}} (numeric QC) and
#' \code{\link{plotDeployments}} (temporal coverage), and follows the same mapping conventions as
#' \code{\link{plotMovements}} and \code{\link{plotMaps}}.
#'
#' @param deployments A receiver-deployment log (a data frame), as produced by
#' \code{\link{importDeployments}} - one row per deployment, with at least the station and coordinate
#' columns. Multiple rows per station (repeated servicing) are reduced to one point per station.
#' @template deploymentSpatialArgs
#' @template deploymentDateArgs
#' @param detections Optional. A detection dataset (data frame or `mobyData`). When supplied, stations
#' with zero detections are highlighted and counted - a quick check for dead receivers or coverage
#' gaps. The detection station column is taken from the dataset's `mobyData` metadata when present,
#' otherwise from the canonical `"station"`.
#' @param color.by Optional. Name of a `deployments` column used to colour the stations (e.g. habitat,
#' receiver model); a colourblind-safe Okabe-Ito palette and a legend are added. Takes precedence over
#' `status.at`.
#' @param status.at Optional. A single date; each station is classified and coloured as `"active"`,
#' `"recovered"` or `"not yet deployed"` at that instant, from the deployment/recovery dates (a station is
#' active while `deploy <= status.at < recover`). Character/`Date` input is interpreted in the data's
#' timezone (not the session's) for reproducibility.
#' @param detection.range Optional. Nominal detection radius in metres, drawn as a semi-transparent
#' circle around each receiver: a single value (applied to all stations), a numeric vector **named by
#' station**, or the name of a `deployments` column holding a per-receiver radius. Stations with a
#' missing (`NA`) or zero radius get no circle. This is a visual aid for coverage assessment, NOT a
#' modelled detection probability.
#' @param range.color Fill colour of the detection-range circles. If NULL, a muted blue is used.
#' @param range.alpha Opacity of the range circles, from 0 (transparent) to 1 (opaque). Defaults to 0.15.
#' @param range.border Border colour of the range circles. Defaults to NA (no border).
#' @param hull Logical. If TRUE, outline the convex hull of the stations (the array footprint).
#' Defaults to FALSE.
#' @param hull.color,hull.lty,hull.lwd Colour, line type and width of the hull outline.
#' @param label Draw station labels: `FALSE` (default) none, `TRUE` the station name, or the name of a
#' `deployments` column to label by.
#' @param label.cex,label.color Size and colour of the station labels.
#' @param label.wrap Logical. If TRUE, wrap multi-word labels onto separate lines. Defaults to FALSE.
#' @param pch Plotting symbol for the stations (a fillable symbol such as 21-25 uses `pt.bg`).
#' Defaults to 21.
#' @param pt.cex,pt.color,pt.bg,pt.lwd Size, border colour, fill colour and border width of the station
#' symbols. `pt.bg` is overridden when `color.by`/`status.at` is used.
#' @param land.shape Optional. An `sf` object with coastlines/landmasses, drawn underneath the array.
#' @param coastline When the map is projected (an `epsg.code`, `land.shape` or `background.layer` is
#' given) and no `land.shape` is supplied, draw a default coastline for the extent. `TRUE` auto-picks a
#' scale; a string forces an `rnaturalearth` scale; `FALSE` draws none. Defaults to
#' \code{getOption("moby.coastline", TRUE)}. Requires the (Suggested) `rnaturalearth` (with
#' `rnaturalearthdata`) or `maps` package.
#' @param land.color Colour of land areas. Defaults to "gray50".
#' @param epsg.code Optional. A projected coordinate reference system (numeric EPSG code or `crs`
#' object, in metre units) used to project the coordinates and map layers. Required for the scale bar;
#' when omitted (and no map layer is supplied) the array is drawn in geographic (longitude/latitude)
#' space with coordinate axes.
#' @param background.layer Optional. A projected raster displayed in the background (e.g. bathymetry).
#' @param background.color Solid colour for the plot background when no `background.layer` is supplied.
#' Defaults to "#F3F7F7".
#' @param background.pal Colour palette for the `background.layer` raster. If NULL, a muted
#' colourblind-safe viridis palette is used.
#' @param scale.km Length of the scale bar, in kilometres. If NULL, chosen automatically. Projected
#' maps only.
#' @param scale.pos Position of the scale bar (a keyword such as "bottomright"). Defaults to "bottomright".
#' @param scale.inset Inset of the scale bar from the plot edges. A single value or a length-2 (x, y)
#' vector. Defaults to 0.05.
#' @param scale.height Thickness of the scale bar. Defaults to 1.5.
#' @param scale.color Colour of the scale-bar labels. Defaults to "black".
#' @param extent.factor Numeric factor to expand the plotting region around the stations. A value of 1
#' keeps the bounding box; values > 1 add a margin. Defaults to 1.15.
#' @param legend Logical. Draw a legend when `color.by`/`status.at`/`detections` is used. Defaults to TRUE.
#' @param main Optional title. If NULL, "Receiver array" is used; FALSE omits it.
#' @param cex Global expansion factor scaling every text element. Defaults to 1.
#' @template deviceArgs
#' @param verbose Logical; print the array summary to the console. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return Invisibly, a data frame of the plotted stations (station, receiver count, coordinates, and -
#' when requested - the `color.by` value, deployment status and detection flag). Called mainly for its
#' side effect (the array map).
#' @seealso \code{\link{checkDeployments}}, \code{\link{plotDeployments}}, \code{\link{importDeployments}},
#' \code{\link{plotMovements}}, \code{\link{plotMaps}}
#' @examples
#' # geographic quick look (no projection needed), with 800 m nominal detection ranges
#' plotArray(rays_deployments, detection.range = 800, label = TRUE)
#'
#' \donttest{
#' # projected map with an automatic coastline and the array footprint
#' plotArray(rays_deployments, epsg.code = 32629, detection.range = 800,
#'           hull = TRUE, coastline = FALSE)
#' }
#' @export

plotArray <- function(deployments,
                      deployment.station.col = "station",
                      deployment.lon.col = "lon",
                      deployment.lat.col = "lat",
                      deployment.deploy.col = "deploy",
                      deployment.recover.col = "recover",
                      detections = NULL,
                      color.by = NULL,
                      status.at = NULL,
                      detection.range = NULL,
                      range.color = NULL,
                      range.alpha = 0.15,
                      range.border = NA,
                      hull = FALSE,
                      hull.color = "grey40",
                      hull.lty = 2,
                      hull.lwd = 1,
                      label = FALSE,
                      label.cex = 0.7,
                      label.color = "grey15",
                      label.wrap = FALSE,
                      pch = 21,
                      pt.cex = 1.4,
                      pt.color = "black",
                      pt.bg = "white",
                      pt.lwd = 1,
                      land.shape = NULL,
                      coastline = getOption("moby.coastline", TRUE),
                      land.color = "gray50",
                      epsg.code = NULL,
                      background.layer = NULL,
                      background.color = "#F3F7F7",
                      background.pal = NULL,
                      scale.km = NULL,
                      scale.pos = "bottomright",
                      scale.inset = 0.05,
                      scale.height = 1.5,
                      scale.color = "black",
                      extent.factor = 1.15,
                      legend = TRUE,
                      main = NULL,
                      cex = 1,
                      file = NULL,
                      width = NULL,
                      height = NULL,
                      res = 300,
                      verbose = getOption("moby.verbose", TRUE)) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  dep <- as.data.frame(deployments)
  required <- c(deployment.station.col, deployment.lon.col, deployment.lat.col)
  miss <- setdiff(required, colnames(dep))
  if (length(miss) > 0)
    .mobyAbort("'deployments' is missing required column(s): ", paste(miss, collapse = ", "),
               ". See importDeployments().")
  if (!is.null(color.by) && !color.by %in% colnames(dep))
    .mobyAbort("'color.by' column '", color.by, "' not found in 'deployments'.")
  if (is.character(label) && length(label) == 1 && !label %in% colnames(dep))
    .mobyAbort("'label' column '", label, "' not found in 'deployments'.")
  if (!is.null(status.at)) {
    if (!deployment.deploy.col %in% colnames(dep))
      .mobyAbort("'status.at' requires the deployment-date column '", deployment.deploy.col, "', which is not in 'deployments'.")
    # coerce character/Date input in the DATA's timezone (not the session's) so the classification is
    # reproducible and locale-independent; an existing POSIXct is an absolute instant, left as-is
    if (!inherits(status.at, "POSIXct"))
      status.at <- tryCatch(as.POSIXct(status.at, tz = .dataTZ(dep[[deployment.deploy.col]])), error = function(e) NA)
    if (length(status.at) != 1 || is.na(status.at))
      .mobyAbort("'status.at' must be a single date (POSIXct or coercible to one).")
  }
  if (!is.null(detection.range)) {
    if (is.character(detection.range)) {
      if (length(detection.range) != 1 || !detection.range %in% colnames(dep))
        .mobyAbort("'detection.range' given as a column name must be a single existing 'deployments' column.")
    } else if (!is.numeric(detection.range)) {
      .mobyAbort("'detection.range' must be a numeric radius in metres, a numeric vector named by station, ",
                 "or the name of a 'deployments' column.")
    } else if (any(detection.range < 0, na.rm = TRUE)) {
      .mobyAbort("'detection.range' must be non-negative (metres).")
    }
  }

  # accept a raster::RasterLayer background for convenience, but work in terra internally
  if (!is.null(background.layer) && inherits(background.layer, "RasterLayer"))
    background.layer <- terra::rast(background.layer)


  #####################################################################################
  # Reduce the deployment log to one point per station ################################
  #####################################################################################

  dep[[deployment.station.col]] <- as.character(dep[[deployment.station.col]])
  by_station <- split(dep, dep[[deployment.station.col]], drop = TRUE)
  first_val <- function(d, col) if (col %in% names(d)) d[[col]][which(!is.na(d[[col]]))[1]] else NA

  stations <- data.frame(
    station     = names(by_station),
    lon         = vapply(by_station, function(d) mean(d[[deployment.lon.col]], na.rm = TRUE), numeric(1)),
    lat         = vapply(by_station, function(d) mean(d[[deployment.lat.col]], na.rm = TRUE), numeric(1)),
    n_receivers = vapply(by_station, function(d)
      if ("receiver" %in% names(d)) length(unique(stats::na.omit(d[["receiver"]]))) else nrow(d), integer(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  if (!is.null(color.by))
    stations$group <- vapply(by_station, function(d) as.character(first_val(d, color.by)), character(1))

  # surface coordinate errors: warn when a station's rows disagree on position (it is plotted at the
  # mean, which would otherwise silently hide exactly the kind of mistake this QC map is meant to catch)
  coord_spread <- vapply(by_station, function(d) {
    lo <- d[[deployment.lon.col]][is.finite(d[[deployment.lon.col]])]; la <- d[[deployment.lat.col]][is.finite(d[[deployment.lat.col]])]
    max(if (length(lo) >= 2) diff(range(lo)) else 0, if (length(la) >= 2) diff(range(la)) else 0)
  }, numeric(1))
  divergent <- names(coord_spread)[coord_spread > 5e-4]                # ~50 m of longitude/latitude
  if (length(divergent) > 0)
    .mobyWarn(length(divergent), " station(s) have inconsistent coordinates across deployment rows and were ",
              "plotted at their mean position (", paste(utils::head(divergent, 5), collapse = ", "),
              if (length(divergent) > 5) ", ..." else "", "); check the deployment log.")

  # deployment status at a chosen instant: active > recovered > not-yet, reduced across a station's rows.
  # A station is 'active' while deployed and still in the water (deploy <= t < recover); at the exact
  # recover instant it is 'recovered'.
  if (!is.null(status.at)) {
    stations$status <- vapply(by_station, function(d) {
      dp <- d[[deployment.deploy.col]]; rc <- if (deployment.recover.col %in% names(d)) d[[deployment.recover.col]] else as.POSIXct(NA)
      active <- any(!is.na(dp) & dp <= status.at & (is.na(rc) | rc > status.at))
      if (active) "active"
      else if (any(!is.na(dp) & dp <= status.at)) "recovered"
      else "not yet deployed"
    }, character(1))
    stations$status <- factor(stations$status, levels = c("active", "recovered", "not yet deployed"))
  }

  # drop stations with missing/implausible coordinates
  bad <- !is.finite(stations$lon) | !is.finite(stations$lat) |
    abs(stations$lon) > 180 | abs(stations$lat) > 90
  n_dropped <- sum(bad)
  if (n_dropped > 0) {
    .mobyInform("- ", n_dropped, " station(s) with missing/implausible coordinates were dropped from the map.",
                verbose = verbose)
    stations <- stations[!bad, , drop = FALSE]
  }
  if (nrow(stations) == 0) .mobyAbort("No stations with valid coordinates to plot.")
  n_stations <- nrow(stations)

  # resolve detection.range to one value per station, ALIGNED TO `stations` BY NAME (never by input
  # position: split() sorts the stations, so a positional vector would land on the wrong receivers)
  if (!is.null(detection.range)) {
    if (is.character(detection.range)) {                       # a per-receiver column
      detection.range <- vapply(by_station[stations$station],
                                function(d) suppressWarnings(as.numeric(first_val(d, detection.range))), numeric(1))
    } else if (length(detection.range) == 1) {                # one radius for all
      detection.range <- rep(detection.range, n_stations)
    } else if (!is.null(names(detection.range))) {            # named vector -> match by station
      unmatched <- !stations$station %in% names(detection.range)
      detection.range <- unname(detection.range[stations$station])
      if (any(unmatched))                                     # a genuine name miss (not an intentional NA value)
        .mobyWarn(sum(unmatched), " station(s) are not named in 'detection.range'; their circles are omitted.")
    } else {
      .mobyAbort("A per-station 'detection.range' must be named by station (or supply a 'deployments' ",
                 "column name); an unnamed vector cannot be reliably matched to stations.")
    }
  }

  # optional detection cross-check: which stations never logged a detection?
  zero_det <- NULL
  if (!is.null(detections)) {
    det <- as.data.frame(detections)
    meta <- attr(detections, "moby")
    det_station <- if (!is.null(meta) && !is.null(meta$station.col)) meta$station.col else "station"
    if (!det_station %in% colnames(det)) {
      .mobyWarn("Could not find the station column ('", det_station, "') in 'detections'; ",
                "skipping the zero-detection check.")
    } else {
      # trim whitespace on both sides so a trailing space does not flag every receiver as dead
      detected <- unique(trimws(as.character(det[[det_station]])))
      matched <- trimws(stations$station) %in% detected
      zero_det <- stations$station[!matched]
      stations$detected <- matched
      if (all(!matched))
        .mobyWarn("None of the ", n_stations, " station(s) matched a detection; check that the station ",
                  "names in 'detections' correspond to those in 'deployments'.")
    }
  }


  #####################################################################################
  # Prepare spatial objects ###########################################################
  #####################################################################################

  projected <- !is.null(land.shape) || !is.null(epsg.code) || !is.null(background.layer)
  geo <- stations[, c("lon", "lat")]                    # geographic coords retained for range circles

  if (projected) {
    coords <- sf::st_as_sf(stations[, c("station", "lon", "lat")], coords = c("lon", "lat"))
    layers <- .prepareMapLayers(coords, land.shape, background.layer, epsg.code)
    coords <- layers$coords; land.shape <- layers$land.shape
    background.layer <- layers$background.layer; epsg.code <- layers$epsg.code
    xy <- sf::st_coordinates(coords)
    bbox <- sf::st_bbox(coords)
    epsg_label <- if (inherits(epsg.code, "crs") && !is.na(epsg.code$epsg)) epsg.code$epsg else NA
  } else {
    xy <- as.matrix(geo)
    bbox <- c(xmin = min(xy[, 1]), ymin = min(xy[, 2]), xmax = max(xy[, 1]), ymax = max(xy[, 2]))
    epsg_label <- NA
  }
  # pad the extent (also guards against a zero-width box for a single station / one lon or lat)
  if (bbox["xmax"] - bbox["xmin"] == 0) { bbox["xmin"] <- bbox["xmin"] - 1; bbox["xmax"] <- bbox["xmax"] + 1 }
  if (bbox["ymax"] - bbox["ymin"] == 0) { bbox["ymin"] <- bbox["ymin"] - 1; bbox["ymax"] <- bbox["ymax"] + 1 }
  bbox <- .expandBbox(bbox, extent.factor)

  # coastline fallback: on a projected map with no user land.shape, fetch a default coastline
  if (projected && is.null(land.shape) && !isFALSE(coastline))
    land.shape <- .defaultCoastline(bbox, epsg.code, coastline, verbose = verbose)
  if (projected && !is.null(land.shape)) land.shape <- sf::st_crop(sf::st_geometry(land.shape), bbox)
  if (!is.null(background.layer)) background.layer <- terra::crop(background.layer, .bboxToExtent(bbox))

  # background raster palette default (resolved once): colourblind-safe viridis, muted for continuous
  if (!is.null(background.layer) && is.null(background.pal)) {
    if (!isTRUE(terra::is.factor(background.layer)[1])) background.pal <- grDevices::adjustcolor(.viridis_pal(100), 0.6)
    else background.pal <- .viridis_pal(nrow(terra::cats(background.layer)[[1]]))
  }

  # detection-range circles, built geodesically in WGS84 then (if projected) reprojected. Stations with
  # a missing (NA) or zero radius get no circle - in either mode, and without erroring.
  range_geom <- NULL
  if (!is.null(detection.range)) {
    draw_r <- which(is.finite(detection.range) & detection.range > 0)
    if (length(draw_r) > 0) {
      ang <- seq(0, 360, length.out = 65)
      rings <- lapply(draw_r, function(i) {
        ring <- .destPoint(cbind(rep(geo$lon[i], length(ang)), rep(geo$lat[i], length(ang))),
                           d = detection.range[i], b = ang)
        rbind(ring, ring[1, ])
      })
      if (projected) {
        range_geom <- sf::st_sfc(lapply(rings, function(r) sf::st_polygon(list(r))), crs = 4326)
        range_geom <- sf::st_transform(range_geom, sf::st_crs(epsg.code))
      } else {
        range_geom <- rings
      }
    }
  }
  if (is.null(range.color)) range.color <- "#3182BD"

  # station fill colours: color.by > status.at > uniform pt.bg
  fill_col <- rep(pt.bg, n_stations); legend_levels <- NULL; legend_cols <- NULL
  if (!is.null(color.by)) {
    lv <- sort(unique(stations$group)); legend_cols <- .okabe_ito_pal(length(lv))
    fill_col <- legend_cols[match(stations$group, lv)]; legend_levels <- lv
  } else if (!is.null(status.at)) {
    all_lv <- levels(stations$status)
    status_cols <- stats::setNames(c("#2CA25F", "#DE9A2A", "grey70"), all_lv)
    fill_col <- status_cols[as.character(stations$status)]
    legend_levels <- all_lv[all_lv %in% as.character(stations$status)]   # observed levels only
    legend_cols <- status_cols[legend_levels]
  }

  # scale bar length + nearest-neighbour spacing (a QC metric), resolved once
  if (projected && is.null(scale.km)) scale.km <- pretty((bbox["xmax"] - bbox["xmin"]) * 0.2 / 1000)[2]
  nn_m <- NA
  if (n_stations >= 2) {
    dm <- matrix(as.numeric(sf::st_distance(sf::st_as_sf(geo, coords = c("lon", "lat"), crs = 4326))),
                 n_stations, n_stations)
    diag(dm) <- NA
    nn_m <- apply(dm, 1, min, na.rm = TRUE)
  }

  # global text sizes derived from a single 'cex'
  cex_title <- 1.1 * cex; cex_sub <- 0.8 * cex; cex_lab <- label.cex * cex
  cex_scale <- 0.6 * cex; cex_leg <- 0.7 * cex


  #####################################################################################
  # Device + canvas ###################################################################
  #####################################################################################

  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  if (!is.null(file)) {
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 7, slope = 0, n = 1, lo = 5, hi = 12),
                     h.rule = list(base = 6.5, slope = 0, n = 1, lo = 4.5, hi = 12))
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  graphics::par(mar = if (!is.null(background.layer)) c(3, 3, 2.5, 6) else c(3, 3, 2.5, 1))

  if (projected) {
    .newMapCanvas(bbox)
  } else {
    # geographic space: use a latitude-correct aspect so the array is not visually distorted
    asp <- 1 / cos(mean(bbox[c("ymin", "ymax")]) * pi / 180)
    graphics::plot(bbox[c("xmin", "xmax")], bbox[c("ymin", "ymax")], type = "n", asp = asp,
                   xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i",
                   xlim = bbox[c("xmin", "xmax")], ylim = bbox[c("ymin", "ymax")])
  }

  # background: solid colour or raster
  if (is.null(background.layer)) {
    graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],
                   graphics::par("usr")[4], col = background.color, border = NA)
  } else {
    if (!isTRUE(terra::is.factor(background.layer)[1])) {
      terra::plot(background.layer, col = background.pal, legend = TRUE, axes = FALSE, add = TRUE,
                  plg = list(cex = cex_leg, shrink = 0.6))
    } else {
      terra::plot(background.layer, col = background.pal, legend = FALSE, axes = FALSE, add = TRUE)
      cats_df <- terra::cats(background.layer)[[1]]
      graphics::legend("right", legend = rev(cats_df[[ncol(cats_df)]]), xpd = TRUE,
                       fill = rev(background.pal), inset = -0.18, bty = "n", cex = cex_leg)
    }
  }

  # land
  .drawLandOverlay(if (projected) land.shape else NULL, land.color)

  # detection-range circles (under the points)
  if (!is.null(range_geom)) {
    fill <- grDevices::adjustcolor(range.color, alpha.f = range.alpha)
    if (projected) graphics::plot(range_geom, col = fill, border = range.border, add = TRUE)
    else for (r in range_geom) graphics::polygon(r[, 1], r[, 2], col = fill, border = range.border)
  }

  # array footprint (convex hull)
  if (isTRUE(hull) && n_stations >= 3) {
    h <- grDevices::chull(xy)
    graphics::polygon(xy[h, 1], xy[h, 2], border = hull.color, col = NA, lty = hull.lty, lwd = hull.lwd)
  } else if (isTRUE(hull)) {
    .mobyInform("- Convex hull needs at least 3 stations; skipped.", verbose = verbose)
  }

  # zero-detection highlight ring (independent of fill colour)
  if (!is.null(zero_det) && length(zero_det) > 0) {
    zi <- which(stations$station %in% zero_det)
    graphics::points(xy[zi, 1], xy[zi, 2], pch = 1, cex = pt.cex * cex * 1.9,
                     col = "#D6604D", lwd = 1.6 * cex)
  }

  # stations
  graphics::points(xy[, 1], xy[, 2], pch = pch, cex = pt.cex * cex,
                   col = pt.color, bg = fill_col, lwd = pt.lwd)

  # labels
  if (!isFALSE(label)) {
    lab <- if (isTRUE(label)) stations$station else as.character(vapply(by_station[stations$station],
                                                                        function(d) as.character(first_val(d, label)), character(1)))
    if (label.wrap) lab <- gsub(" ", "\n", lab, fixed = TRUE)
    graphics::text(xy[, 1], xy[, 2], labels = lab, pos = 3, offset = 0.5, cex = cex_lab, col = label.color)
  }

  # geographic axes / projected scale bar
  if (projected) {
    .drawScaleBar(bbox, scale.km = scale.km, scale.pos = scale.pos, scale.inset = scale.inset,
                  height = scale.height, cex = cex_scale, color = scale.color)
  } else {
    graphics::axis(1, cex.axis = cex_scale); graphics::axis(2, cex.axis = cex_scale, las = 1)
    graphics::box()
    graphics::title(xlab = "Longitude", ylab = "Latitude", line = 2, cex.lab = cex_sub)
  }

  # title
  if (!isFALSE(main)) {
    graphics::title(main = if (is.null(main)) "Receiver array" else main, cex.main = cex_title,
                    line = 1, xpd = TRUE)
  }

  # legend (fill grouping and/or the zero-detection marker)
  if (legend && (!is.null(legend_levels) || (!is.null(zero_det) && length(zero_det) > 0))) {
    leg_txt <- character(0); leg_pch <- numeric(0); leg_pt.bg <- character(0); leg_col <- character(0)
    if (!is.null(legend_levels)) {
      leg_txt <- legend_levels; leg_pch <- rep(pch, length(legend_levels))
      leg_pt.bg <- legend_cols[seq_along(legend_levels)]; leg_col <- rep(pt.color, length(legend_levels))
    }
    if (!is.null(zero_det) && length(zero_det) > 0) {
      leg_txt <- c(leg_txt, "no detections"); leg_pch <- c(leg_pch, 1)
      leg_pt.bg <- c(leg_pt.bg, NA); leg_col <- c(leg_col, "#D6604D")
    }
    graphics::legend("topright", legend = leg_txt, pch = leg_pch, pt.bg = leg_pt.bg, col = leg_col,
                     pt.cex = pt.cex, bty = "n", cex = cex_leg, inset = 0.02, xpd = NA)
  }


  #####################################################################################
  # Console summary ###################################################################
  #####################################################################################

  if (verbose) {
    kv <- .kv
    .summaryOpen("Receiver array")
    kv("Stations", sprintf("%d (%d receiver%s)", n_stations, sum(stations$n_receivers),
                           if (sum(stations$n_receivers) == 1) "" else "s"))
    kv("Projection", if (projected) sprintf("projected%s", if (!is.na(epsg_label)) paste0(" (EPSG:", epsg_label, ")") else "")
       else "geographic (lon/lat)")
    if (n_stations >= 2)
      kv("Spacing (NN)", sprintf("median %s (range %s-%s)", .fmtDist(stats::median(nn_m)),
                                 .fmtDist(min(nn_m)), .fmtDist(max(nn_m))))
    if (!is.null(detection.range)) {
      rr <- detection.range[is.finite(detection.range)]
      if (length(rr)) kv("Detection range", if (length(unique(rr)) == 1) .fmtDist(rr[1])
                         else sprintf("%s-%s", .fmtDist(min(rr)), .fmtDist(max(rr))))
    }
    if (!is.null(status.at)) {
      tb <- table(stations$status); tb <- tb[tb > 0]
      kv("Status", paste(sprintf("%d %s", tb, names(tb)), collapse = ", "))
    }
    if (!is.null(zero_det)) kv("Zero detections", sprintf("%d station%s", length(zero_det), if (length(zero_det) == 1) "" else "s"))
    if (projected && !is.na(scale.km)) kv("Scale bar", sprintf("%g km", scale.km))
    if (n_dropped > 0) kv("Dropped", sprintf("%d (bad coordinates)", n_dropped))
    .summaryClose()
  }

  invisible(stations)
}


#######################################################################################################
## Internal helper ####################################################################################

#' Format a distance in metres or kilometres for the array summary card
#' @keywords internal
#' @noRd
.fmtDist <- function(m) {
  if (!is.finite(m)) return("NA")
  if (m >= 1000) sprintf("%.1f km", m / 1000) else sprintf("%.0f m", m)
}

#######################################################################################################
#######################################################################################################
