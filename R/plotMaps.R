##################################################################################################
## Plot maps with movement trajectories + home ranges ############################################
##################################################################################################

#' Plot home-range maps with movement trajectories
#'
#' @description Draws one map per individual showing its kernel utilisation distribution (UD home
#' range) as graded utilisation isopleths, optionally over a coastline and a background raster
#' (e.g. bathymetry or temperature), together with inferred movement tracks and time-ordered
#' detection points. Panels can be grouped by `id.groups`.
#'
#' @details The home range is drawn as true **utilisation isopleths**: the UD is converted to a
#' volume distribution with \code{\link[adehabitatHR]{getvolumeUD}}, so `ud.contour = 0.95` really
#' masks everything outside the 95% home range, and the colour ramp and legend encode the isopleth
#' level (core = brightest). All layers (coordinates, coastline, background raster) are reconciled to
#' one projected CRS via the shared spatial pipeline before plotting.
#'
#' @inheritParams as_moby
#' @param data A `mobyData` or data frame of animal positions (centres of activity) with lon/lat.
#' @param uds Utilisation distributions from \code{\link{calculateUDs}} - either KDE (`estUD`/`estUDm`,
#' drawn as graded volume isopleths) or AKDE (drawn from the stored isopleth contours, with a
#' confidence envelope; see `ud.uncertainty`). Pass the whole `calculateUDs()` output; the method is
#' detected automatically.
#' @param animal.tracks The distance-enriched output of \code{\link{calculateStepDistances}} (tracks
#' are read from its `"trajectories"` attribute), or a pre-extracted list of per-individual
#' linestring geometries.
#' @param id.groups Optional named list of ID groups (one block of panels each).
#' @param land.shape Optional `sf` coastline/landmass polygon drawn over the background.
#' @param coastline When no `land.shape` is supplied, draw a default coastline for the map extent.
#' `TRUE` auto-picks a scale; a string (`"small"`/`"medium"`/`"large"`) forces an `rnaturalearth`
#' scale; `FALSE` draws no coastline. Defaults to \code{getOption("moby.coastline", TRUE)} (set that
#' option to switch the fallback off package-wide). Requires the (Suggested) `rnaturalearth` (with
#' `rnaturalearthdata`) or `maps` package; if neither is installed, no coastline is drawn.
#' @param epsg.code Projected CRS (numeric EPSG code or `crs` object, metre units) for the map. If
#' NULL, inferred from `land.shape`, then from `uds`.
#' @param land.color Colour of land areas. Defaults to "gray50".
#' @param background.layer Optional projected raster displayed behind the home range.
#' @param background.color Solid background colour when no `background.layer` is supplied.
#' Defaults to "#F3F7F7".
#' @param background.pal Colour palette for the `background.layer` raster. If NULL, a muted
#' bathymetry palette is used.
#' @param background.title Legend title for the background layer. Defaults to "Layer values".
#' @param discard.missing Logical; drop individuals with fewer than 5 detections (UD needs >= 5).
#' Defaults to TRUE.
#' @param color.pal Colour palette for the home-range isopleths. If NULL, viridis is used.
#' @param ud.contour Outer utilisation isopleth to display, as a proportion in (0, 1]. Cells outside
#' this contour are made transparent. e.g. 0.95 shows the 95% home range. Defaults to 0.95.
#' @param ud.uncertainty For AKDE UDs only: how to draw the confidence envelope of the outer isopleth.
#' `"band"` (default) shades the region between the low and high contours; `"outline"` draws the low
#' and high contours as dashed outlines; `"none"` draws only the estimate. Ignored for KDE (which has
#' no CI).
#' @param tracks.color,tracks.lty,tracks.lwd Colour, line type and width of the movement trajectories.
#' @param plot.detections Logical; overlay detection points. Defaults to TRUE.
#' @param pts.color Colour(s) for detection points: one colour, or three for the first / intermediate
#' / last detection. Defaults to c("#009E73", "white", "#D55E00").
#' @param pts.cex Size(s) for detection points: one value, or three for first / intermediate / last.
#' Defaults to c(1, 0.6, 1).
#' @param ud.legend Logical; draw the home-range isopleth legend. Defaults to TRUE.
#' @param title.color Colour of the per-panel ID label. Defaults to "black".
#' @param title.pos Keyword position of the ID label (e.g. "topleft"). Defaults to "topleft".
#' @param title.inset Inset of the ID label from the plot edges (single value or x/y vector).
#' Defaults to c(-0.08, 0).
#' @param scale.km Scale-bar length in kilometres. If NULL (default), ~20% of the map width.
#' @param scale.color Colour of the scale-bar labels. Defaults to "black".
#' @param scale.pos Keyword position of the scale bar. Defaults to "bottomright".
#' @param scale.inset Inset of the scale bar from the plot edges. Defaults to 0.05.
#' @param extent.factor Bounding-box expansion factor around the positions. Defaults to 1.1.
#' @param main Optional overall title above the panel grid.
#' @param cex Global expansion factor for all plot text (ID labels, legends, scale bar). Defaults to 1.
#' @param ncol Number of panel columns. If NULL, set from the number of individuals.
#' @template deviceArgs
#'
#' @return Invisibly, a tidy data frame with one row per mapped individual: the isopleth level, the
#' home-range area (km^2, when coordinates are metric) and the number of detections.
#' @seealso \code{\link{calculateUDs}}, \code{\link{calculateStepDistances}}, \code{\link{plotMovements}}
#' @export
#'
#' @examples
#' \dontrun{
#' pdf("./home-range-maps.pdf", width = 12, height = 20)
#' plotMaps(data = animal_coas, uds = kud_output, animal.tracks = tracks,
#'          land.shape = coastline, id.col = "ID", lon.col = "lon", lat.col = "lat")
#' dev.off()
#' }

plotMaps <- function(data,
                     uds = NULL,
                     animal.tracks = NULL,
                     id.groups = NULL,
                     id.col = NULL,
                     lon.col = NULL,
                     lat.col = NULL,
                     land.shape = NULL,
                     coastline = getOption("moby.coastline", TRUE),
                     epsg.code = NULL,
                     land.color = "gray50",
                     background.layer = NULL,
                     background.color = "#F3F7F7",
                     background.pal = NULL,
                     background.title = NULL,
                     discard.missing = TRUE,
                     color.pal = NULL,
                     ud.contour = 0.95,
                     ud.uncertainty = c("band", "outline", "none"),
                     tracks.color = "black",
                     tracks.lty = 2,
                     tracks.lwd = 0.2,
                     plot.detections = TRUE,
                     pts.color = c("#009E73", "white", "#D55E00"),
                     pts.cex = c(1, 0.6, 1),
                     ud.legend = TRUE,
                     title.color = "black",
                     title.pos = "topleft",
                     title.inset = c(-0.08, 0),
                     scale.km = NULL,
                     scale.color = "black",
                     scale.pos = "bottomright",
                     scale.inset = 0.05,
                     extent.factor = 1.1,
                     main = NULL,
                     cex = 1,
                     ncol = NULL,
                     file = NULL,
                     width = NULL,
                     height = NULL,
                     res = 300) {

  ##############################################################################
  ## Checks + input normalisation ##############################################
  ##############################################################################

  reviewed_params <- .validateArguments()
  data <- as.data.frame(reviewed_params$data)
  land.shape <- reviewed_params$land.shape

  errors <- c()
  ud.uncertainty <- match.arg(ud.uncertainty)
  # AKDE home ranges are drawn from pre-computed sf isopleth contours (no adehabitatHR needed); KDE
  # home ranges go through the adehabitatHR getvolumeUD() raster path.
  akde <- !is.null(uds) && identical(attr(uds, "method"), "akde")
  if (!is.numeric(ud.contour) || ud.contour <= 0 || ud.contour > 1)
    errors <- c(errors, "'ud.contour' must be a proportion in (0, 1] (e.g. 0.95). If you meant a percentage, divide by 100.")
  if (!is.null(uds) && !akde && !requireNamespace("adehabitatHR", quietly = TRUE))
    errors <- c(errors, "Plotting KDE 'uds' requires the 'adehabitatHR' package.")
  if (length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)

  # normalise the UD input:
  #  - AKDE  -> the sf isopleth contours (K50/K95, each carrying low/est/high CI polygons)
  #  - KDE   -> an estUDm raster (unwrapping the whole calculateUDs() output if that was passed)
  akde_contours <- NULL
  if (!is.null(uds)) {
    if (akde) {
      akde_contours <- uds[grepl("^K[0-9]+$", names(uds))]
      if (length(akde_contours) == 0)
        stop("AKDE 'uds' carry no isopleth contours (K50/K95). Re-run calculateUDs(method = 'akde').", call. = FALSE)
      ud_ids <- unique(unlist(lapply(akde_contours, function(s) as.character(s[[id.col]]))))
      uds <- NULL                                          # disable the KDE raster path
    } else {
      if ("ud" %in% names(uds)) uds <- uds$ud
      if (!inherits(uds, "estUDm")) {
        if (inherits(uds[[1]], "estUDm")) {
          ids <- unlist(lapply(uds, names))
          uds <- unlist(uds); names(uds) <- ids
          class(uds) <- "estUDm"
        } else stop("Supplied 'uds' are in the wrong format. Provide the output of calculateUDs().", call. = FALSE)
      }
      ud_ids <- names(uds)
    }
    if (all(!unique(data[, id.col]) %in% ud_ids))
      stop("Names of 'uds' do not match the supplied data IDs.", call. = FALSE)
  }
  has_ud <- akde || !is.null(uds)

  # normalise animal.tracks (read the 'trajectories' attribute if present)
  if (!is.null(animal.tracks)) {
    tracks_geom <- getTrajectories(animal.tracks)
    if (!is.null(tracks_geom)) animal.tracks <- tracks_geom
    else warning("- No 'trajectories' attribute in 'animal.tracks'; assuming it already holds linestring geometries.", call. = FALSE)
  }

  # global text sizes from a single 'cex'
  cex_title <- 1.5 * cex; cex_legend <- 0.8 * cex; cex_scale <- 0.7 * cex
  if (is.null(color.pal)) color.pal <- .viridis_pal(100)

  ##############################################################################
  ## Cleanup + spatial reconciliation ##########################################
  ##############################################################################

  if (is.null(id.groups)) id.groups <- list(levels(factor(data[, id.col])))

  n_before <- nlevels(factor(data[, id.col])); missing_individuals <- character(0)
  if (discard.missing) {
    counts <- table(data[, id.col])
    missing_individuals <- names(counts[counts < 5])
    data <- data[!data[, id.col] %in% missing_individuals, ]
    data[, id.col] <- droplevels(factor(data[, id.col]))
    id.groups <- lapply(id.groups, function(x) x[!x %in% missing_individuals])
  }
  id.groups <- id.groups[lengths(id.groups) > 0]
  if (nlevels(factor(data[, id.col])) == 0) stop("No individuals left to plot after discarding.", call. = FALSE)

  # infer the map CRS from the UDs when neither epsg.code nor land.shape carries one
  if (is.null(epsg.code) && is.null(land.shape) && has_ud) {
    kcrs <- tryCatch(sf::st_crs(if (akde) akde_contours[[1]] else uds[[1]]), error = function(e) NA)
    if (!is.na(kcrs)) epsg.code <- kcrs
  }

  # accept a raster::RasterLayer background for convenience, but work in terra internally
  if (!is.null(background.layer) && inherits(background.layer, "RasterLayer"))
    background.layer <- terra::rast(background.layer)

  # reconcile coordinates + land.shape + background.layer onto one projected CRS
  coords <- sf::st_as_sf(data, coords = c(lon.col, lat.col))
  layers <- .prepareMapLayers(coords, land.shape, background.layer, epsg.code)
  coords <- layers$coords; land.shape <- layers$land.shape
  background.layer <- layers$background.layer; epsg.code <- layers$epsg.code

  ##############################################################################
  ## Plot parameters ###########################################################
  ##############################################################################

  if (!is.null(background.layer)) {
    if (is.null(background.pal)) background.pal <- rev(.bathy_deep_pal(100)[c(30:100)])
    if (is.null(background.title)) background.title <- "Layer values"
  }
  if (length(pts.color) == 1) pts.color <- rep(pts.color, 3)
  if (length(pts.cex) == 1) pts.cex <- rep(pts.cex, 3)

  # bounding box: union of positions and the UD extent, expanded
  bbox <- sf::st_bbox(coords)
  if (akde) {
    ke <- sf::st_bbox(do.call(rbind, lapply(akde_contours, function(s) s["geometry"])))
    bbox <- c(xmin = min(bbox["xmin"], ke["xmin"]), ymin = min(bbox["ymin"], ke["ymin"]),
              xmax = max(bbox["xmax"], ke["xmax"]), ymax = max(bbox["ymax"], ke["ymax"]))
  } else if (!is.null(uds)) {
    ke <- terra::ext(terra::rast(methods::as(uds[[1]], "SpatialGridDataFrame")))
    bbox <- c(xmin = min(bbox["xmin"], terra::xmin(ke)), ymin = min(bbox["ymin"], terra::ymin(ke)),
              xmax = max(bbox["xmax"], terra::xmax(ke)), ymax = max(bbox["ymax"], terra::ymax(ke)))
  }
  bbox <- .expandBbox(bbox, extent.factor)
  # coastline fallback: if the user supplied no land.shape, fetch a default coastline for the extent
  if (is.null(land.shape) && !isFALSE(coastline))
    land.shape <- .defaultCoastline(bbox, epsg.code, coastline)
  if (!is.null(land.shape)) land.shape <- suppressWarnings(sf::st_crop(sf::st_geometry(land.shape), bbox))
  if (!is.null(background.layer)) background.layer <- terra::crop(background.layer, .bboxToExtent(bbox))
  if (is.null(scale.km)) scale.km <- pretty((bbox["xmax"] - bbox["xmin"]) * 0.2 / 1000)[2]

  # layout: a shared legend slot when a UD or background legend is needed
  add_legend_space <- (has_ud && ud.legend) || !is.null(background.layer)
  if (is.null(ncol)) {
    n_ind <- sum(lengths(id.groups)); ncol <- if (n_ind == 1) 1 else if (n_ind == 2) 2 else 3
  }

  ##############################################################################
  ## Device + layout ###########################################################
  ##############################################################################

  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  if (!is.null(file)) {
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 2, slope = 3.2, n = ncol, lo = 4, hi = 30),
                     h.rule = list(base = 1, slope = 2.8, n = ceiling(sum(lengths(id.groups)) / ncol), lo = 4, hi = 30),
                     crowd.unit = "panels")
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  layout_params <- .setLayout(ncol, id.groups, plots.height = 6, dividers.height = 0.6,
                              legend = add_legend_space, min.legend.plots = 1, expand.legend = FALSE)
  nplots <- max(layout_params$matrix, na.rm = TRUE)
  layout_params$matrix[is.na(layout_params$matrix)] <- nplots + 1
  layout(mat = layout_params$matrix, heights = layout_params$heights)
  par(mar = c(0.4, 0.4, 0.4, 0.4), oma = if (length(id.groups) == 1) c(0, 0, 0, 0) else c(0, 3, 0, 1))
  if (!is.null(main)) par(oma = par("oma") + c(0, 0, 2, 0))

  ##############################################################################
  ## Compute (home-range areas) + summary ######################################
  ##############################################################################

  ids <- levels(factor(data[, id.col]))
  data_individual <- split(data, f = factor(data[, id.col]), drop = FALSE)
  results <- data.frame()
  .printMapsSummary(n_ids = length(ids), n_total = n_before, n_discarded = length(missing_individuals),
                    has_kud = has_ud, has_tracks = !is.null(animal.tracks),
                    epsg = if (inherits(epsg.code, "crs") && !is.na(epsg.code$epsg)) epsg.code$epsg else NA,
                    ud.contour = ud.contour, scale.km = scale.km)

  ##############################################################################
  ## Draw panels ###############################################################
  ##############################################################################

  n_legend_slot <- if (add_legend_space) nplots else NA
  for (i in seq_len(nplots)) {

    if (i <= length(ids)) {
      id <- ids[i]
      .newMapCanvas(bbox)

      # background (raster or solid); drawn inline (legend is hand-drawn via .colorlegend)
      if (!is.null(background.layer)) terra::image(background.layer, col = background.pal, axes = FALSE, add = TRUE)
      else rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = background.color, border = NA)

      # home range: KDE utilisation isopleths (getvolumeUD raster, masked at ud.contour), or AKDE
      # isopleth contours (sf) drawn with their confidence envelope
      plot_kernel <- FALSE
      if (!is.null(uds) && id %in% names(uds) && nrow(data_individual[[id]]) > 0) {
        vol <- adehabitatHR::getvolumeUD(uds[[id]])
        r <- terra::rast(methods::as(vol, "SpatialGridDataFrame"))
        r[r > ud.contour * 100] <- NA
        r <- terra::crop(r, .bboxToExtent(bbox))
        suppressWarnings(terra::image(r, zlim = c(0, ud.contour * 100), col = rev(color.pal), axes = FALSE, add = TRUE))
        cell_km2 <- prod(terra::res(r)) / 1e6
        results <- rbind(results, data.frame(id = id, isopleth = ud.contour,
                                             area_km2 = sum(!is.na(terra::values(r))) * cell_km2,
                                             n_detections = nrow(data_individual[[id]]), stringsAsFactors = FALSE))
        plot_kernel <- TRUE
      } else if (akde && nrow(data_individual[[id]]) > 0) {
        area_km2 <- .drawAkdeContours(akde_contours, id, id.col, ud.contour, color.pal, ud.uncertainty)
        if (!is.na(area_km2)) {
          results <- rbind(results, data.frame(id = id, isopleth = ud.contour, area_km2 = area_km2,
                                               n_detections = nrow(data_individual[[id]]), stringsAsFactors = FALSE))
          plot_kernel <- TRUE
        }
      }

      .drawLandOverlay(land.shape, land.color)

      # movement tracks
      if (!is.null(animal.tracks) && !is.null(animal.tracks[[id]]))
        lines(sf::st_coordinates(animal.tracks[[id]]), col = tracks.color, lwd = tracks.lwd, lty = tracks.lty)

      # time-ordered detection points (first / intermediate / last)
      if (plot.detections && nrow(data_individual[[id]]) > 0) {
        id_coords <- coords[coords[[id.col]] == id, ]
        xy <- sf::st_coordinates(id_coords)
        points(xy, pch = 21, bg = pts.color[2], cex = pts.cex[2], col = "gray15", lwd = 0.25)
        points(xy[1, , drop = FALSE], pch = 21, bg = pts.color[1], cex = pts.cex[1], lwd = 0.25)
        points(xy[nrow(xy), , drop = FALSE], pch = 21, bg = pts.color[3], cex = pts.cex[3], lwd = 0.25)
      }

      legend(title.pos, legend = id, inset = title.inset, cex = cex_title, text.col = title.color, bty = "n", text.font = 2)
      .drawScaleBar(bbox, scale.km = scale.km, scale.pos = scale.pos, scale.inset = scale.inset,
                    cex = cex_scale, color = scale.color)
      box(lty = "solid", lwd = 0.75)

    } else if (!is.na(n_legend_slot) && i == n_legend_slot) {
      # shared legend slot: home-range isopleths + background layer
      plot.new()
      first_legend <- FALSE
      if (has_ud && ud.legend) {
        if (akde) {
          .drawAkdeLegend(akde_contours, ud.contour, color.pal, ud.uncertainty, cex_legend)
        } else {
          iso <- c(25, 50, 75, 95); iso <- iso[iso <= ud.contour * 100]
          .colorlegend(col = rev(color.pal), zlim = c(0, ud.contour * 100), zval = iso, zlab = paste0(iso, "%"),
                       posx = c(0.15, 0.85), posy = c(0.75, 0.8), main = "Home-range isopleth",
                       main.cex = cex_legend + 0.2, cex = cex_legend, horizontal = TRUE)
        }
        first_legend <- TRUE
      }
      if (!is.null(background.layer)) {
        vals <- terra::values(background.layer)
        layer_range <- range(vals, na.rm = TRUE)
        layer_labs <- pretty(vals, min.n = 4); layer_labs <- layer_labs[layer_labs >= layer_range[1] & layer_labs <= layer_range[2]]
        posy_vals <- if (first_legend) c(0.4, 0.45) else c(0.75, 0.8)
        .colorlegend(col = background.pal, zlim = layer_range, zval = layer_labs, main = background.title,
                     posx = c(0.15, 0.85), posy = posy_vals, main.cex = cex_legend + 0.2,
                     digit = max(.decimalPlaces(layer_labs)), cex = cex_legend, horizontal = TRUE)
      }
    }
  }

  # id.group labels (rotated, left outer margin)
  if (length(id.groups) > 1) {
    label_pos <- grconvertY(1 - (layout_params$group_positions / sum(layout_params$heights)), "ndc", "user")
    text(x = grconvertX(0.01, "ndc", "user"), y = label_pos, labels = names(id.groups),
         srt = 90, cex = cex_title, font = 2, xpd = NA, adj = c(0.5, 0.5))
  }
  if (!is.null(main)) mtext(main, side = 3, outer = TRUE, font = 2, cex = cex_title * 0.8, line = 0.2)

  invisible(results)
}


##################################################################################################
## Internal helper ###############################################################################

#' @keywords internal
#' @noRd
.printMapsSummary <- function(n_ids, n_total, n_discarded, has_kud, has_tracks, epsg, ud.contour, scale.km) {
  kv <- .kv
  .summaryOpen("Home-range maps")
  kv("Individuals", sprintf("%d of %d%s", n_ids, n_total, if (n_discarded > 0) sprintf(" (%d discarded, < 5 detections)", n_discarded) else ""))
  kv("Layers", sprintf("UD: %s; tracks: %s", if (has_kud) "yes" else "no", if (has_tracks) "yes" else "no"))
  kv("Projection", if (!is.na(epsg)) sprintf("EPSG:%s", epsg) else "projected")
  if (has_kud) kv("Home range", sprintf("%g%% isopleth", ud.contour * 100))
  kv("Scale bar", sprintf("%g km", scale.km))
  .summaryClose()
}


#' Draw the AKDE isopleth contours for one individual, with its confidence envelope.
#'
#' @description Draws the estimate isopleth(s) (`ci == "est"`) up to `ud.contour`, and — on the
#' outermost level — the AKDE confidence envelope: either the shaded band between the low and high
#' contours (`uncertainty = "band"`) or the low/high contours as dashed outlines
#' (`uncertainty = "outline"`). `uncertainty = "none"` draws only the estimate. Returns the estimate
#' area (km^2) of the outermost drawn isopleth, for the summary table.
#' @keywords internal
#' @noRd
.drawAkdeContours <- function(contours, id, id.col, ud.contour, color.pal, uncertainty) {
  pct  <- as.integer(sub("^K", "", names(contours)))
  keep <- which(pct <= ud.contour * 100 + 1e-9)
  if (!length(keep)) keep <- which.min(abs(pct - ud.contour * 100))
  ord  <- keep[order(pct[keep], decreasing = TRUE)]        # outermost isopleth first
  outer_area <- NA_real_
  for (j in seq_along(ord)) {
    k   <- ord[j]
    sfk <- contours[[k]]; sfk <- sfk[as.character(sfk[[id.col]]) == id, , drop = FALSE]
    if (nrow(sfk) == 0) next
    geo <- function(lab) { g <- sf::st_geometry(sfk[sfk$ci == lab, ]); if (length(g)) g[1] else NULL }
    est <- geo("est"); if (is.null(est)) est <- sf::st_geometry(sfk)[1]
    base_col <- rev(color.pal)[min(pct[k], 100)]
    est_fill <- grDevices::adjustcolor(base_col, 0.6)
    is_outer <- (j == 1)                                   # outermost level carries the CI envelope + area
    low <- high <- NULL
    if (is_outer) {
      outer_area <- tryCatch(as.numeric(sf::st_area(est)) / 1e6, error = function(e) NA_real_)
      if (uncertainty != "none") { low <- geo("low"); high <- geo("high") }
      # a band is drawn UNDER the estimate (which then becomes an outline, so the band shows the fill)
      if (uncertainty == "band" && !is.null(low) && !is.null(high)) {
        band <- tryCatch(sf::st_difference(high, low), error = function(e) high)
        suppressWarnings(plot(band, col = grDevices::adjustcolor(base_col, 0.22), border = NA, add = TRUE))
        est_fill <- NA
      }
    }
    suppressWarnings(plot(est, col = est_fill, border = base_col, lwd = 1.6, add = TRUE))
    # dashed low/high bounds are drawn OVER the filled estimate, so both are visible
    if (is_outer && uncertainty == "outline" && !is.null(low) && !is.null(high)) {
      br <- grDevices::adjustcolor(base_col, 0.95)
      suppressWarnings(plot(high, col = NA, border = br, lty = 2, lwd = 1, add = TRUE))
      suppressWarnings(plot(low,  col = NA, border = br, lty = 2, lwd = 1, add = TRUE))
    }
  }
  outer_area
}

#' Discrete legend for AKDE home ranges (isopleth level(s) + the confidence representation)
#' @keywords internal
#' @noRd
.drawAkdeLegend <- function(contours, ud.contour, color.pal, uncertainty, cex) {
  pct <- as.integer(sub("^K", "", names(contours)))
  pct <- sort(pct[pct <= ud.contour * 100 + 1e-9], decreasing = TRUE)
  if (!length(pct)) pct <- max(as.integer(sub("^K", "", names(contours))))
  base_cols <- rev(color.pal)[pmin(pct, 100)]
  txt  <- paste0(pct, "% home range")
  fill <- grDevices::adjustcolor(base_cols, 0.6); bord <- base_cols
  lty  <- rep(NA_integer_, length(pct))
  if (uncertainty == "band") {
    txt  <- c(txt, sprintf("%g%% CI band", max(pct)))
    fill <- c(fill, grDevices::adjustcolor(rev(color.pal)[min(max(pct), 100)], 0.22))
    bord <- c(bord, NA); lty <- c(lty, NA)
  } else if (uncertainty == "outline") {
    txt  <- c(txt, sprintf("%g%% CI (low / high)", max(pct)))
    fill <- c(fill, NA); bord <- c(bord, NA); lty <- c(lty, 2L)
  }
  legend("top", legend = txt, fill = fill, border = bord, lty = lty, lwd = 1, seg.len = 1.4,
         cex = cex, bty = "n", title = "Home range (AKDE)", title.font = 2)
}

##################################################################################################
##################################################################################################
##################################################################################################
