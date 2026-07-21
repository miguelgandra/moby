#######################################################################################################
## Calculate movement transitions (movement network) #################################################
#######################################################################################################

#' Build a movement (transition) network between locations
#'
#' @description Builds a **movement network**, in which nodes are spatial units (receivers,
#' stations, habitats or any user-defined spatial class) and directed edges represent transitions
#' - i.e. individuals moving from one location to another. Edge weights summarise how many
#' movements occurred and how many distinct individuals performed them. This is the spatial
#' counterpart to \code{\link{calculateAssociations}}, and the foundation for movement-network
#' metrics, visualisation and (in a forthcoming release) randomisation.
#'
#' Transitions are defined between *consecutive distinct* residence events (visits) in each
#' individual's time-ordered sequence, where visits are segmented by `max.gap` (see
#' \code{\link{calculateVisits}}): consecutive detections at the same location are collapsed into one
#' visit unless separated by an absence longer than `max.gap`, so a transition is genuine movement
#' between two different nodes and a long absence is not mistaken for one continuous stay.
#'
#' @inheritParams as_moby
#' @param data A data frame (or \code{\link{mobyData}}) of detections.
#' @param spatial.col Name of the column defining the network nodes (e.g. receiver, station,
#' habitat, region).
#' @param id.groups Optional named list of ID groups; when supplied, an independent network is
#' built for each group (carried in the `group` column of the node/edge tables).
#' @param max.gap Maximum tolerated gap between successive detections within a single visit (passed
#' to \code{\link{calculateVisits}}); a longer absence ends a stay, so a later return counts as a new
#' visit. Defaults to 48 (hours). Use `Inf` to segment on location changes only.
#' @param max.gap.units Units of `max.gap`: one of `"hours"` (default), `"days"`, `"mins"`, `"secs"`.
#'
#' @return A \code{\link{mobyNetwork}} object of type `"movement"`. Its rows are the directed
#' edges with columns: `group`, `from`, `to`, `n_movements` (number of transitions),
#' `n_individuals` (distinct individuals performing them), and `mean_duration_h` (mean transit time
#' in hours, computed over continuously-observed movements only — transits that spanned a `max.gap`
#' absence are excluded as their timing is unobserved). The node table (\code{\link{networkNodes}})
#' has one row per location with `group`, `site`, `n_detections`, `n_individuals`, `n_residence`
#' (number of residence events / visits), `mean_residence_h` (mean visit duration in hours), and
#' `lon`/`lat` (mean coordinates, when available). The full per-transition records (departure /
#' arrival times, hour, month, `duration_h`, and `crossed_gap`) are stored in the
#' `"transition_records"` attribute.
#'
#' @seealso \code{\link{calculateVisits}}, \code{\link{calculateAssociations}},
#' \code{\link{mobyNetwork}}, \code{\link{transitionsTable}}
#'
#' @examples
#' data(rays)
#' # build a movement network with the receiver stations as nodes
#' trans <- calculateTransitions(rays, spatial.col = "station")
#' trans
#' head(networkEdges(trans))
#'
#' @export

calculateTransitions <- function(data,
                                 spatial.col,
                                 id.col = NULL,
                                 datetime.col = NULL,
                                 id.groups = NULL,
                                 max.gap = 48,
                                 max.gap.units = c("hours", "days", "mins", "secs")) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # capture mobyData metadata before the data is coerced to a plain data.frame
  meta <- attr(data, "moby")

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # residence events (and hence node residency + transit gaps) are segmented on max.gap: an absence
  # longer than this ends a stay, so a later return is a new visit rather than one continuous stay.
  max.gap.units <- match.arg(max.gap.units)
  max.gap.secs  <- .gapToSecs(max.gap, max.gap.units)
  if (is.finite(max.gap.secs))
    message("- Segmenting residence events with max.gap = ", max.gap, " ", max.gap.units,
            "; tune/justify per system (Inf = split on location change only).")

  # resolve (optional) node coordinates from the mobyData metadata / canonical defaults
  lon.col <- if (!is.null(meta) && !is.null(meta$lon.col)) meta$lon.col else "lon"
  lat.col <- if (!is.null(meta) && !is.null(meta$lat.col)) meta$lat.col else "lat"
  has_coords <- lon.col %in% colnames(data) && lat.col %in% colnames(data) &&
    is.numeric(data[[lon.col]]) && is.numeric(data[[lat.col]])

  # ensure the spatial column is a factor (defines node ordering)
  if (!inherits(data[, spatial.col], "factor")) {
    data[, spatial.col] <- as.factor(data[, spatial.col])
  }
  ordered_sites <- levels(data[, spatial.col])

  # group label per row (single "all" group when id.groups is not supplied) and the
  # number of (tagged) individuals per group, used as the denominator for transit percentages
  if (is.null(id.groups)) {
    data$.group <- "all"
    group_levels <- "all"
    group_sizes <- stats::setNames(nlevels(data[, id.col]), "all")
  } else {
    grp_map <- stats::setNames(rep(names(id.groups), lengths(id.groups)), as.character(unlist(id.groups)))
    data$.group <- unname(grp_map[as.character(data[, id.col])])
    group_levels <- names(id.groups)
    group_sizes <- lengths(id.groups)
  }

  tz <- .dataTZ(data[, datetime.col])

  ##############################################################################
  ## Per-group computation #####################################################
  ##############################################################################

  edges_list <- list()
  nodes_list <- list()
  records_list <- list()

  for (g in group_levels) {
    dg <- data[data$.group == g, , drop = FALSE]
    if (nrow(dg) == 0) next

    # ---- residence events (shared segmentation) + transition records ---------
    # segment each individual into visits (breaks on location change OR gap > max.gap); the same
    # events drive node residency below, so counts/durations stay consistent with calculateVisits().
    ev <- .residenceRuns(dg, spatial.col, id.col, datetime.col, max.gap.secs)

    recs <- list()
    for (ind in unique(ev$id)) {
      e <- ev[ev$id == ind, , drop = FALSE]          # already ordered by arrival
      if (nrow(e) < 2) next
      from <- e$site[-nrow(e)]; to <- e$site[-1]
      keep <- from != to                             # gap-separated same-site re-visits are not edges
      if (!any(keep)) next
      recs[[ind]] <- data.frame(group = g, from = from[keep], to = to[keep], id = ind,
                                departure = e$departure[-nrow(e)][keep],
                                arrival   = e$arrival[-1][keep],
                                stringsAsFactors = FALSE)
    }
    recs <- if (length(recs) > 0) do.call(rbind, recs) else
      data.frame(group = character(0), from = character(0), to = character(0),
                 id = character(0), departure = as.POSIXct(character(0), tz = tz),
                 arrival = as.POSIXct(character(0), tz = tz), stringsAsFactors = FALSE)

    if (nrow(recs) > 0) {
      recs$duration_h <- as.numeric(difftime(recs$arrival, recs$departure, units = "hours"))
      # a transit spanning a > max.gap absence was not continuously observed: keep the edge for
      # connectivity, but flag it and exclude its (meaningless) transit time from mean_duration_h.
      recs$crossed_gap <- as.numeric(difftime(recs$arrival, recs$departure, units = "secs")) > max.gap.secs
      recs$departure_hour  <- strftime(recs$departure, "%H", tz = tz)
      recs$arrival_hour    <- strftime(recs$arrival, "%H", tz = tz)
      recs$departure_month <- strftime(recs$departure, "%m", tz = tz)
      recs$arrival_month   <- strftime(recs$arrival, "%m", tz = tz)
    }
    records_list[[g]] <- recs

    # ---- edges (aggregate transitions) ---------------------------------------
    if (nrow(recs) > 0) {
      key <- paste(recs$from, recs$to, sep = "\r")
      open_dur <- ifelse(recs$crossed_gap, NA_real_, recs$duration_h)  # transit time of observed moves
      edge_df <- data.frame(
        group = g,
        from = tapply(recs$from, key, `[`, 1),
        to   = tapply(recs$to, key, `[`, 1),
        n_movements = as.integer(tapply(recs$id, key, length)),
        n_individuals = as.integer(tapply(recs$id, key, function(x) length(unique(x)))),
        mean_duration_h = as.numeric(tapply(open_dur, key, function(x) {
          x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_
        })),
        row.names = NULL, stringsAsFactors = FALSE)
      edge_df <- edge_df[order(factor(edge_df$from, levels = ordered_sites),
                               factor(edge_df$to, levels = ordered_sites)), , drop = FALSE]
      edges_list[[g]] <- edge_df
    }

    # ---- nodes ---------------------------------------------------------------
    n_det <- tapply(seq_len(nrow(dg)), dg[, spatial.col], length)
    n_ind <- tapply(as.character(dg[, id.col]), dg[, spatial.col], function(x) length(unique(x)))
    # residency per node from the residence events: number of visits + mean visit duration
    ev_site <- factor(ev$site, levels = ordered_sites)
    n_res <- tapply(ev$id, ev_site, length)
    res_h <- tapply(ev$residence_h, ev_site, function(x) mean(x, na.rm = TRUE))
    node_df <- data.frame(group = g, site = ordered_sites,
                          n_detections = as.integer(n_det[ordered_sites]),
                          n_individuals = as.integer(n_ind[ordered_sites]),
                          n_residence = as.integer(n_res[ordered_sites]),
                          mean_residence_h = as.numeric(res_h[ordered_sites]),
                          row.names = NULL, stringsAsFactors = FALSE)
    node_df$n_detections[is.na(node_df$n_detections)] <- 0L
    node_df$n_individuals[is.na(node_df$n_individuals)] <- 0L
    node_df$n_residence[is.na(node_df$n_residence)] <- 0L
    if (has_coords) {
      node_df$lon <- as.numeric(tapply(data[[lon.col]][data$.group == g], dg[, spatial.col], mean, na.rm = TRUE)[ordered_sites])
      node_df$lat <- as.numeric(tapply(data[[lat.col]][data$.group == g], dg[, spatial.col], mean, na.rm = TRUE)[ordered_sites])
    }
    # keep only nodes that were visited in this group
    node_df <- node_df[node_df$n_detections > 0, , drop = FALSE]
    nodes_list[[g]] <- node_df
  }

  edges <- if (length(edges_list) > 0) do.call(rbind, edges_list) else
    data.frame(group = character(0), from = character(0), to = character(0),
               n_movements = integer(0), n_individuals = integer(0),
               mean_duration_h = numeric(0), stringsAsFactors = FALSE)
  nodes <- do.call(rbind, nodes_list)
  rownames(edges) <- NULL; rownames(nodes) <- NULL

  ##############################################################################
  ## Assemble mobyNetwork ######################################################
  ##############################################################################

  .mobyNetwork(edges, nodes = nodes, type = "movement",
               spatial.col = spatial.col, id.col = id.col, id.groups = id.groups,
               group.sizes = group_sizes, ordered.sites = ordered_sites,
               has.coords = has_coords, transition_records = records_list,
               processing.date = Sys.time())
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
