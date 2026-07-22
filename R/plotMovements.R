#######################################################################################################
## Plot movement (transition) network ################################################################
#######################################################################################################

#' Plot a movement (transition) network on a map
#'
#' @description Plots a spatial movement network from a \code{\link{mobyNetwork}} object of type
#' "movement" (as returned by \code{\link{calculateTransitions}}). Nodes are drawn at the geographic
#' location of each spatial unit (receiver, station, habitat, region) and directed edges (arrows)
#' represent transitions between them. Node size reflects the number of individuals detected at each
#' location, while edge width and labels reflect either the number of movements or the number of
#' transiting individuals.
#'
#' An optional coastline (`land.shape`) and a background raster (`background.layer`, e.g. bathymetry
#' or temperature) can be supplied to render the network over a projected map, complete with a scale
#' bar. When node coordinates are unavailable, a force-directed layout is used instead.
#'
#' This function is also invoked automatically by \code{plot()} on a movement-type
#' \code{\link{mobyNetwork}} (see \code{\link{plot.mobyNetwork}}).
#'
#' @param network A \code{\link{mobyNetwork}} object of type `"movement"`, as returned by
#' \code{\link{calculateTransitions}}. When the network carries ID groups, a separate panel is drawn
#' for each group.
#' @param edge.metric Metric used to scale and label the edges: `"movements"` (number of transitions,
#' the default) or `"individuals"` (number of distinct individuals transiting between sites).
#' @param land.shape Optional. An `sf` object containing coastlines, drawn underneath the network.
#' @param coastline When the map is projected (an `epsg.code`, `land.shape` or `background.layer` is
#' given) and no `land.shape` is supplied, draw a default coastline for the extent. `TRUE` auto-picks
#' a scale; a string forces an `rnaturalearth` scale; `FALSE` draws none. Defaults to
#' \code{getOption("moby.coastline", TRUE)}. Requires the (Suggested) `rnaturalearth` (with
#' `rnaturalearthdata`) or `maps` package.
#' @param land.color Colour of land areas. Defaults to "gray50".
#' @param epsg.code Optional. A projected coordinate reference system (numeric EPSG code or `crs`
#' object, in metre units) used to project the node coordinates and the map layers. Required for the
#' scale bar when coordinates are geographic.
#' @param background.layer Optional. A projected raster displayed in the background (e.g. bathymetry,
#' temperature).
#' @param background.color Solid colour used for the plot background when no `background.layer` is
#' supplied. Defaults to "#F3F7F7".
#' @param background.pal Colour palette for the `background.layer` raster. If NULL (default), a
#' muted colourblind-safe viridis palette is used (continuous layers are drawn semi-transparent).
#' @param color.nodes.by How to colour the nodes: `"detection"` (default) uses a single colour for
#' all visited sites, or `"group"` assigns a different colour per ID group (one per panel).
#' @param nodes.color Colour(s) for the nodes. If NULL (default), a single blue for `"detection"`
#' mode, or a colourblind-friendly Okabe-Ito palette (one colour per group) for `"group"` mode.
#' @param nodes.alpha Numeric transparency for the nodes (0 = transparent, 1 = opaque). Defaults to 0.8.
#' @param nodes.size A length-2 numeric vector giving the minimum and maximum node sizes relative to
#' the x-axis (as proportions). Defaults to c(0.04, 0.08).
#' @param nodes.label Logical. If TRUE (default), node labels (site name and number of individuals)
#' are drawn inside the nodes.
#' @param nodes.label.wrap Logical. If TRUE, splits node labels into multiple lines. Defaults to FALSE.
#' @param nodes.label.color Font colour for node labels. Defaults to "white".
#' @param repel.nodes Logical. If TRUE (default), nodes are repositioned with a repulsion algorithm
#' to reduce overlap (requires the `packcircles` package).
#' @param repel.buffer Controls the spacing between nodes when `repel.nodes` is TRUE. Defaults to 1.1.
#' @param edge.color Colour(s) for the edges (recycled across groups). Defaults to "darkblue".
#' @param edge.curved Curvature of the edges. A numeric value sets the curvature (0 = straight);
#' TRUE means 0.5, FALSE means 0. Defaults to 0.5.
#' @param edge.width A length-2 numeric vector giving the minimum and maximum edge widths.
#' Defaults to c(0.4, 3.5).
#' @param edge.arrow.size Size of the arrowheads. Defaults to 0.5.
#' @param edge.arrow.width Width of the arrowheads. Defaults to 1.5.
#' @param edge.label Logical. If TRUE (default), each edge is labelled with its `edge.metric` value.
#' @param edge.label.color Colour of the edge labels. Defaults to "black".
#' @param edge.label.font Font for the edge labels (1 = plain, 2 = bold, 3 = italic, 4 = bold italic).
#' Defaults to 1.
#' @param scale.km Length of the scale bar, in kilometres. If NULL, it is defined automatically.
#' Only drawn when coordinates are projected.
#' @param scale.pos Position of the scale bar, specified by keyword (e.g. "bottom"). Defaults to "bottom".
#' @param scale.inset Controls how far the scale bar is inset from the plot edges. A single value or a
#' length-2 vector (x, y). Defaults to c(0, 0.05).
#' @param scale.height Thickness of the scale bar. Defaults to 1.5.
#' @param scale.color Colour of the scale-bar labels. Defaults to "black".
#' @param extent.factor Numeric factor to expand the plotting region around the node positions. A
#' value of 1 keeps the bounding box; values > 1 increase the extent. Defaults to 1.1.
#' @param main Optional overall title drawn above the panel grid.
#' @param cex Global expansion factor scaling every text element (titles, node/edge labels, scale
#' bar, legend). Defaults to 1.
#' @param ncol Number of columns in the multi-panel layout. Defaults to 1.
#' @template deviceArgs
#' @param ... Further arguments passed to \code{\link[igraph]{plot.igraph}}.
#'
#' @return Invisibly, a tidy data frame of the plotted nodes (site, group, coordinates, number of
#' individuals and detections). The corresponding edge table is attached as attribute `"edges"`.
#' Called mainly for its side effect (the network plot).
#' @seealso \code{\link{calculateTransitions}}, \code{\link{plot.mobyNetwork}},
#' \code{\link{networkMetrics}}, \code{\link{plotAssociations}}
#' @examples
#' # Build a station-to-station transition network, then draw it on a projected map
#' net <- calculateTransitions(rays, spatial.col = "station")
#' plotMovements(net, epsg.code = 32629, coastline = FALSE)   # coastline = TRUE draws a coastline
#' @export

plotMovements <- function(network,
                          land.shape = NULL,
                          coastline = getOption("moby.coastline", TRUE),
                          epsg.code = NULL,
                          edge.metric = c("movements", "individuals"),
                          land.color = "gray50",
                          background.layer = NULL,
                          background.color = "#F3F7F7",
                          background.pal = NULL,
                          color.nodes.by = c("detection", "group"),
                          nodes.color = NULL,
                          nodes.alpha = 0.8,
                          nodes.size = c(0.04, 0.08),
                          nodes.label = TRUE,
                          nodes.label.wrap = FALSE,
                          nodes.label.color = "white",
                          repel.nodes = TRUE,
                          repel.buffer = 1.1,
                          edge.color = "darkblue",
                          edge.curved = 0.5,
                          edge.width = c(0.4, 3.5),
                          edge.arrow.size = 0.5,
                          edge.arrow.width = 1.5,
                          edge.label = TRUE,
                          edge.label.color = "black",
                          edge.label.font = 1,
                          scale.km = NULL,
                          scale.pos = "bottom",
                          scale.inset = c(0, 0.05),
                          scale.height = 1.5,
                          scale.color = "black",
                          extent.factor = 1.1,
                          main = NULL,
                          ncol = 1,
                          cex = 1,
                          file = NULL,
                          width = NULL,
                          height = NULL,
                          res = 300,
                          ...) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  if (!inherits(network, "mobyNetwork") || !identical(attr(network, "network.type"), "movement")) {
    stop("'network' must be a movement-type 'mobyNetwork' object (see calculateTransitions()).", call. = FALSE)
  }
  edge.metric <- match.arg(edge.metric)
  color.nodes.by <- match.arg(color.nodes.by)

  errors <- c()
  if (!requireNamespace("igraph", quietly = TRUE))
    errors <- c(errors, "The 'igraph' package is required. Install it with install.packages('igraph').")
  if (length(nodes.size) != 2) errors <- c(errors, "'nodes.size' must supply a min and max value.")
  if (length(edge.width) != 2) errors <- c(errors, "'edge.width' must supply a min and max value.")
  if (length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)

  # node repulsion relies on the optional 'packcircles' package; degrade gracefully if absent
  if (repel.nodes && !requireNamespace("packcircles", quietly = TRUE)) {
    warning("The 'packcircles' package is not installed; skipping node repulsion. Install it to enable it.", call. = FALSE)
    repel.nodes <- FALSE
  }

  # global text sizes derived from a single 'cex' (no cex.* sprawl)
  cex_title <- 1.1 * cex; cex_sub <- 0.8 * cex; cex_nodelab <- 0.5 * cex
  cex_edgelab <- 0.6 * cex; cex_scale <- 0.6 * cex; cex_leg <- 0.7 * cex

  # extract edges, nodes and grouping
  edges <- networkEdges(network)
  nodes <- networkNodes(network)
  group_levels <- unique(nodes$group)
  n_groups <- length(group_levels)
  has_coords <- isTRUE(attr(network, "has.coords")) && all(c("lon", "lat") %in% colnames(nodes))
  wcol <- if (edge.metric == "movements") "n_movements" else "n_individuals"
  if (length(edge.color) == 1) edge.color <- rep(edge.color, n_groups)

  # node colours: single blue (detection mode) or an Okabe-Ito palette across groups
  if (is.null(nodes.color))
    nodes.color <- if (color.nodes.by == "group") .okabe_ito_pal(n_groups) else "darkblue"

  # save/restore graphical parameters unconditionally (device hygiene, P3)
  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  if (!is.null(file)) {
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 4.5, slope = 3.2, n = ncol, lo = 4.5, hi = 30),
                     h.rule = list(base = 3, slope = 2.8, n = ceiling(n_groups / ncol), lo = 4.5, hi = 30),
                     crowd.unit = "panels")
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }


  #####################################################################################
  # Fallback: no coordinates -> force-directed layout #################################
  #####################################################################################

  if (!has_coords) {
    warning("This movement network has no node coordinates; using a force-directed layout (no map).", call. = FALSE)
    .printMovementsSummary(n_nodes = nrow(nodes), n_edges = sum(edges$from != edges$to), n_groups = n_groups,
                           projected = FALSE, epsg = NA, edge.metric = edge.metric, scale.km = NA, n_dropped = 0)
    graphics::par(mfrow = grDevices::n2mfrow(n_groups), mar = c(1, 1, 2, 1))
    if (!is.null(main)) graphics::par(oma = c(0, 0, 2, 0))
    for (g in seq_along(group_levels)) {
      gname <- group_levels[g]
      n_g <- nodes[nodes$group == gname, , drop = FALSE]
      e_g <- edges[edges$group == gname & edges$from != edges$to & !is.na(edges[[wcol]]) & edges[[wcol]] > 0, , drop = FALSE]
      if (nrow(n_g) == 0) next
      gr <- igraph::graph_from_data_frame(
        d = if (nrow(e_g) > 0) e_g[, c("from", "to", wcol)] else data.frame(from = character(0), to = character(0)),
        vertices = data.frame(name = n_g$site), directed = TRUE)
      ewidth <- if (nrow(e_g) > 0) .rescale(e_g[[wcol]], to = edge.width) else NULL
      elab <- if (edge.label && nrow(e_g) > 0) paste0("(", e_g[[wcol]], ")") else NA
      igraph::plot.igraph(gr, layout = igraph::layout_with_fr(gr),
                          vertex.size = .rescale(n_g$n_individuals, to = c(8, 26)),
                          vertex.color = grDevices::adjustcolor(nodes.color[((g - 1) %% length(nodes.color)) + 1], alpha.f = nodes.alpha),
                          vertex.frame.color = NA,
                          vertex.label = if (nodes.label) n_g$site else NA,
                          vertex.label.color = nodes.label.color, vertex.label.cex = cex_nodelab,
                          edge.width = ewidth, edge.color = edge.color[g], edge.arrow.size = edge.arrow.size,
                          edge.curved = edge.curved, edge.label = elab, edge.label.cex = cex_edgelab,
                          edge.label.color = edge.label.color, edge.label.font = edge.label.font,
                          main = if (n_groups > 1) gname else "Movement Network", ...)
    }
    if (!is.null(main)) graphics::mtext(main, side = 3, outer = TRUE, font = 2, cex = cex_title, line = 0.2)
    return(invisible(nodes))
  }


  #####################################################################################
  # Prepare spatial objects (projected map) ###########################################
  #####################################################################################

  projected <- !is.null(land.shape) || !is.null(epsg.code) || !is.null(background.layer)

  # accept a raster::RasterLayer background for convenience, but work in terra internally
  if (!is.null(background.layer) && inherits(background.layer, "RasterLayer"))
    background.layer <- terra::rast(background.layer)

  # global site coordinates (mean across groups); drop sites with missing coordinates
  site_coords <- stats::aggregate(nodes[, c("lon", "lat")], by = list(site = nodes$site),
                                   FUN = function(z) mean(z, na.rm = TRUE))
  bad <- !stats::complete.cases(site_coords[, c("lon", "lat")]) | !is.finite(site_coords$lon) | !is.finite(site_coords$lat)
  n_dropped <- sum(bad)
  if (n_dropped > 0) {
    message(sprintf("- %d site(s) with missing coordinates were dropped from the map.", n_dropped))
    site_coords <- site_coords[!bad, , drop = FALSE]
  }
  if (nrow(site_coords) == 0) stop("No sites with valid coordinates to plot.", call. = FALSE)

  if (projected) {
    coords <- sf::st_as_sf(site_coords, coords = c("lon", "lat"))
    layers <- .prepareMapLayers(coords, land.shape, background.layer, epsg.code)
    coords <- layers$coords; land.shape <- layers$land.shape
    background.layer <- layers$background.layer; epsg.code <- layers$epsg.code
    xy <- sf::st_coordinates(coords)
    bbox <- sf::st_bbox(coords)
    epsg_label <- if (inherits(epsg.code, "crs") && !is.na(epsg.code$epsg)) epsg.code$epsg else NA
  } else {
    xy <- as.matrix(site_coords[, c("lon", "lat")])
    bbox <- c(xmin = min(xy[, 1]), ymin = min(xy[, 2]), xmax = max(xy[, 1]), ymax = max(xy[, 2]))
    epsg_label <- NA
  }
  bbox <- .expandBbox(bbox, extent.factor)

  # coastline fallback: on a projected map with no user land.shape, fetch a default coastline
  if (projected && is.null(land.shape) && !isFALSE(coastline))
    land.shape <- .defaultCoastline(bbox, epsg.code, coastline)

  # crop map layers to the (expanded) bounding box; raster uses the axis-correct extent
  if (projected && !is.null(land.shape)) land.shape <- sf::st_crop(sf::st_geometry(land.shape), bbox)
  if (!is.null(background.layer)) background.layer <- terra::crop(background.layer, .bboxToExtent(bbox))

  # raster palette default (resolved ONCE, before the panel loop): colourblind-safe viridis
  # (perceptually uniform), muted for a background layer, replacing the non-CVD-safe terrain.colors
  if (!is.null(background.layer) && is.null(background.pal)) {
    if (!isTRUE(terra::is.factor(background.layer)[1])) background.pal <- grDevices::adjustcolor(.viridis_pal(100), 0.6)
    else background.pal <- .viridis_pal(nrow(terra::cats(background.layer)[[1]]))
  }

  # common edge-width scale + scale-bar length, resolved ONCE
  e_all <- edges[edges$from != edges$to & !is.na(edges[[wcol]]) & edges[[wcol]] > 0, , drop = FALSE]
  val_range <- if (nrow(e_all) > 0) range(e_all[[wcol]]) else c(0, 1)
  if (projected && is.null(scale.km)) scale.km <- pretty((bbox["xmax"] - bbox["xmin"]) * 0.15 / 1000)[2]

  .printMovementsSummary(n_nodes = nrow(site_coords), n_edges = nrow(e_all), n_groups = n_groups,
                         projected = projected, epsg = epsg_label, edge.metric = edge.metric,
                         scale.km = if (projected) scale.km else NA, n_dropped = n_dropped)


  #####################################################################################
  # Configure layout ##################################################################
  #####################################################################################

  rows <- ceiling(n_groups / ncol)
  graphics::par(mfrow = c(rows, ncol), oma = if (is.null(main)) c(0, 0, 0, 0) else c(0, 0, 2, 0))
  if (!is.null(background.layer)) graphics::par(mar = c(1, 2, 2, 6)) else graphics::par(mar = c(1, 2, 2, 1))


  #####################################################################################
  # Plot network(s) ###################################################################
  #####################################################################################

  layout_out <- list()
  for (g in seq_along(group_levels)) {

    gname <- group_levels[g]
    n_g <- nodes[nodes$group == gname & nodes$site %in% site_coords$site, , drop = FALSE]
    e_g <- edges[edges$group == gname & edges$from != edges$to &
                   !is.na(edges[[wcol]]) & edges[[wcol]] > 0, , drop = FALSE]
    if (nrow(n_g) == 0) next

    g_xy <- xy[match(n_g$site, site_coords$site), , drop = FALSE]

    # empty canvas + background + land (background drawing kept inline: raster::plot legend differs)
    .newMapCanvas(bbox)
    if (is.null(background.layer)) {
      graphics::rect(graphics::par("usr")[1], graphics::par("usr")[3], graphics::par("usr")[2],
                     graphics::par("usr")[4], col = background.color, border = "black")
    } else {
      if (!isTRUE(terra::is.factor(background.layer)[1])) {
        # continuous layer: terra draws the colour bar via plg=list()
        terra::plot(background.layer, col = background.pal, legend = TRUE, axes = FALSE, add = TRUE,
                    plg = list(cex = cex_leg, shrink = 0.6))
      } else {
        # categorical layer: hand-drawn legend from the category labels (last cats() column)
        terra::plot(background.layer, col = background.pal, legend = FALSE, axes = FALSE, add = TRUE)
        cats_df <- terra::cats(background.layer)[[1]]
        graphics::legend("right", legend = rev(cats_df[[ncol(cats_df)]]), xpd = TRUE,
                         fill = rev(background.pal), inset = -0.18, bty = "n", cex = cex_leg)
      }
      graphics::box()
    }
    .drawLandOverlay(if (projected) land.shape else NULL, land.color)

    # titles
    graphics::title(main = if (n_groups > 1) gname else "Movement Network", cex.main = cex_title, line = 1.2, xpd = TRUE)
    metric_lab <- if (edge.metric == "movements") "no. of movements" else "no. of transiting individuals"
    graphics::title(main = tools::toTitleCase(metric_lab), cex.main = cex_sub, line = 0.4, font.main = 1, xpd = TRUE)

    # scale bar (projected only)
    if (projected) .drawScaleBar(bbox, scale.km = scale.km, scale.pos = scale.pos, scale.inset = scale.inset,
                                 height = scale.height, cex = cex_scale, color = scale.color)

    # node sizes, labels, colours
    node_size <- .rescale_vertex_igraph(n_g$n_individuals, minmax.relative.size = c(nodes.size[1], nodes.size[2]))
    if (nodes.label) {
      node_titles <- paste0(n_g$site, "\n", "n=", n_g$n_individuals)
      if (nodes.label.wrap) node_titles <- gsub(" ", "\n", node_titles, fixed = TRUE)
    } else node_titles <- rep(NA, nrow(n_g))
    vertex_color <- if (color.nodes.by == "group") rep(nodes.color[((g - 1) %% length(nodes.color)) + 1], nrow(n_g))
                    else rep(nodes.color[1], nrow(n_g))
    vertex_color <- grDevices::adjustcolor(vertex_color, alpha.f = nodes.alpha)

    # node repulsion (avoid overlapping nodes)
    if (repel.nodes) {
      radii <- (graphics::par("usr")[2] - graphics::par("usr")[1]) * c(0.035, 0.07)
      radii <- .rescale(node_size, to = radii)
      repel_data <- cbind(g_xy[, 1:2, drop = FALSE], radius = radii * repel.buffer)
      new_coords <- packcircles::circleRepelLayout(repel_data, xlim = bbox[c("xmin", "xmax")],
                                                   ylim = bbox[c("ymin", "ymax")], sizetype = "radius",
                                                   wrap = FALSE, weights = 1)$layout
      g_xy <- as.matrix(new_coords[, c("x", "y")])
    }

    # build and plot the igraph network on top of the map
    gr <- igraph::graph_from_data_frame(
      d = if (nrow(e_g) > 0) e_g[, c("from", "to", wcol)] else data.frame(from = character(0), to = character(0)),
      vertices = data.frame(name = n_g$site), directed = TRUE)
    ewidth <- if (nrow(e_g) > 0) .rescale(e_g[[wcol]], from = val_range, to = edge.width) else NULL
    elab <- if (edge.label && nrow(e_g) > 0) paste0("(", e_g[[wcol]], ")") else NA
    igraph::plot.igraph(gr, layout = g_xy, add = TRUE, rescale = FALSE,
                        edge.color = edge.color[g], edge.curved = edge.curved, edge.width = ewidth,
                        vertex.size = node_size, vertex.color = vertex_color, vertex.frame.color = NA,
                        vertex.label = node_titles, vertex.label.family = "Helvetica",
                        vertex.label.color = nodes.label.color, vertex.label.cex = cex_nodelab,
                        edge.arrow.size = edge.arrow.size, edge.arrow.width = edge.arrow.width,
                        edge.label = elab, edge.label.family = "Helvetica", edge.label.cex = cex_edgelab,
                        edge.label.color = edge.label.color, edge.label.font = edge.label.font, ...)

    layout_out[[gname]] <- data.frame(site = n_g$site, group = gname, x = g_xy[, 1], y = g_xy[, 2],
                                      n_individuals = n_g$n_individuals,
                                      n_detections = if ("n_detections" %in% names(n_g)) n_g$n_detections else NA,
                                      stringsAsFactors = FALSE)
  }

  if (!is.null(main)) graphics::mtext(main, side = 3, outer = TRUE, font = 2, cex = cex_title, line = 0.2)

  out <- if (length(layout_out)) do.call(rbind, layout_out) else site_coords
  attr(out, "edges") <- e_all
  invisible(out)
}


#######################################################################################################
## Internal helper ####################################################################################

#' @keywords internal
#' @noRd
.printMovementsSummary <- function(n_nodes, n_edges, n_groups, projected, epsg, edge.metric, scale.km, n_dropped) {
  kv <- .kv
  .summaryOpen("Movement network")
  kv("Nodes / edges", sprintf("%d sites, %d edges", n_nodes, n_edges))
  kv("Groups", n_groups)
  kv("Edge metric", edge.metric)
  kv("Projection", if (projected) sprintf("projected%s", if (!is.na(epsg)) paste0(" (EPSG:", epsg, ")") else "") else "force-directed (no map)")
  if (projected && !is.na(scale.km)) kv("Scale bar", sprintf("%g km", scale.km))
  if (n_dropped > 0) kv("Dropped sites", sprintf("%d (missing coordinates)", n_dropped))
  .summaryClose()
}

#######################################################################################################
#######################################################################################################
