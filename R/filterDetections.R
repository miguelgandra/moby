#######################################################################################################
## Filter detections data #############################################################################
#######################################################################################################

#' Filter and quality-control acoustic detections
#'
#' @description Cleans a detection dataset by applying a sequence of quality-control filters, each of
#' which is opt-in (nothing but the temporal bounds is applied unless you ask for it). Removed
#' detections are never discarded silently: they are returned in `data_discarded` with the exact
#' `reason` for removal, and the filtered detections carry a `qc_flag` column.
#'
#' @details Filters are applied in this order (a filter is skipped when its controlling argument is
#' left at its "off" value):
#'
#' \enumerate{
#'   \item \strong{Duplicate removal} (`remove.duplicates`): exact-duplicate records (same animal,
#'     timestamp and station) are dropped up front so they cannot slip past the other filters.
#'   \item \strong{Pre-tagging}: detections before an animal's `tagging.date`.
#'   \item \strong{Cut-off}: detections after an optional `cutoff.date`.
#'   \item \strong{False-detection (min_lag)}: the standard short-interval false-detection filter
#'     (Pincock 2012; Simpfendorfer et al. 2015). For each detection it computes the shortest interval
#'     to the nearest other detection of the \emph{same tag on the same receiver}; a detection whose
#'     min_lag exceeds `min.lag.factor * nominal.delay` (default 30x the transmitter nominal delay), or
#'     that has no same-receiver companion at all, is treated as a spurious (code-collision) decode.
#'     Requires `nominal.delay`; off otherwise.
#'   \item \strong{Isolation} (`isolation.window`): an optional residency / sporadic-visitor tool that
#'     removes a detection when the time gap to BOTH temporal neighbours (across the whole array)
#'     exceeds the window. This is NOT a false-detection filter and is off by default.
#'   \item \strong{Speed}: flags steps exceeding `max.speed` using a minimum inter-detection distance
#'     (`straight/least-cost distance - 2 * acoustic.range`, a conservative lower bound). For each
#'     flagged step the local +/-3-detection window is hierarchically clustered; an isolated spatial
#'     outlier (a singleton cluster) is removed, whereas a corroborated group of `min.corroboration`
#'     or more detections is FLAGGED (`qc_flag = "overspeed_review"`) and retained, so a genuine
#'     relocation or a systematic error (bad coordinates, clock drift, under-estimated range) is
#'     surfaced for review rather than cascaded away. Requires `max.speed`; off otherwise.
#'   \item \strong{Minimum detections / minimum days}: drop individuals with too few surviving
#'     detections or too few days with detections.
#' }
#'
#' When a `land.shape` is supplied, distances (for the speed filter) are shortest in-water (least-cost)
#' distances rather than straight lines. `nominal.delay` and `max.speed` may be a single value applied
#' to all individuals or a named numeric vector keyed by `id.col` (for mixed tag families or
#' multi-species datasets); `nominal.delay` is additionally read from the `mobyData` metadata when not
#' supplied.
#'
#' @inheritParams as_moby
#' @param data A data frame (or `mobyData`) of raw animal detections.
#' @param cutoff.dates Optional. Cut-off date(s) beyond which detections are discarded (tag expiry,
#' last download, etc.). A single POSIXct applied to all IDs, or a named POSIXct vector keyed by `id.col`.
#' @param remove.duplicates Logical; drop exact-duplicate records (same animal ID, timestamp and
#' station) before any other filter. Defaults to TRUE.
#' @param nominal.delay Transmitter nominal (mean) delay, in seconds, used to scale the min_lag
#' false-detection window. A single value applied to all individuals, or a named numeric vector keyed
#' by `id.col` (for mixed tag families). Read from the `mobyData` metadata (`nominal.delay`) when not
#' supplied. `NULL` (and no metadata) disables the false-detection filter.
#' @param min.lag.factor Multiplier defining the min_lag threshold as `min.lag.factor * nominal.delay`
#' (in seconds). Defaults to 30 (Pincock 2012 rule of thumb). Use [plotMinLag()] to check whether this
#' default suits a given dataset.
#' @param min.lag.threshold Optional. Set the min_lag threshold directly (in seconds), overriding
#' `min.lag.factor * nominal.delay`.
#' @param isolation.window Optional residency / sporadic-visitor filter (in hours): a detection is
#' removed when the gap to BOTH temporal neighbours exceeds this window. `NULL` (default) or `FALSE`
#' disables it. This is a coverage tool, NOT a false-detection filter (use `nominal.delay` for that).
#' @param max.speed Maximum plausible swim speed. A single value applied to all individuals, or a named
#' numeric vector keyed by `id.col` (for multi-species / size-specific limits). `NULL` disables the
#' speed filter.
#' @param speed.unit Units of `max.speed`: either "m/s" or "km/h".
#' @param acoustic.range Assumed detection range of the receivers (metres), used in the speed filter's
#' minimum-distance correction. Note this makes the speed filter blind to consecutive detections closer
#' than `2 * acoustic.range`. Defaults to 600.
#' @param min.corroboration Integer; in the speed filter, a flagged spatial outlier is auto-removed only
#' if its cluster has fewer than this many detections. A corroborated group of `min.corroboration` or
#' more is flagged for review instead of deleted (prevents cascading deletion of real relocations).
#' Defaults to 2 (remove isolated singletons only).
#' @param max.iterations Maximum passes of the speed-filter convergence loop. Defaults to 20; use `Inf`
#' for unlimited.
#' @param min.detections Discard individuals with fewer than this many surviving detections. 0 = off.
#' @param min.days Discard individuals detected on fewer than this many distinct days. 0 = off.
#' @param land.shape Optional. Coastline/landmass polygons (`sf`, or `SpatialPolygons*`), enabling
#' shortest in-water distances in the speed filter. Passed to \code{\link{calculateStepDistances}}.
#' @param epsg.code Optional integer EPSG code of a projected (metre-based) CRS used to project
#' positions for the speed filter.
#' @param verbose Logical; print progress and the filtering summary. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#' @param ... Further arguments passed to \code{\link{calculateStepDistances}} (e.g. `grid.resolution`,
#' `mov.directions`, `cores`).
#'
#' @return An object of class `mobyFilter` (a list, with a `print` method) containing:
#' \describe{
#'   \item{`data`}{the filtered detections as a `mobyData`, with an added `qc_flag` column
#'     (`"valid"`, or `"overspeed_review"` for speed-flagged detections retained for review);}
#'   \item{`data_discarded`}{the removed detections, each with a `reason` column;}
#'   \item{`summary`}{a per-individual data frame of raw counts and removals by filter.}
#' }
#' The filtering parameters used are stored as the `"parameters"` attribute.
#'
#' @references
#' Pincock, D.G. (2012) False detections: what they are and how to remove them from detection data.
#' VEMCO Application Note DOC-004691.
#'
#' Simpfendorfer, C.A., Huveneers, C., Steckenreuter, A., Tattersall, K., Hoenner, X., Harcourt, R. &
#' Heupel, M.R. (2015) Ghosts in the data: false detections in VEMCO pulse position modulation acoustic
#' telemetry monitoring equipment. Animal Biotelemetry, 3:55.
#'
#' Kessel, S.T., Cooke, S.J., Heupel, M.R., Hussey, N.E., Simpfendorfer, C.A., Vagle, S. & Fisk, A.T.
#' (2014) A review of detection range testing in aquatic passive acoustic telemetry studies. Reviews in
#' Fish Biology and Fisheries, 24:199-218.
#'
#' @examples
#' data(rays)
#' # default: only the temporal bounds are applied (tagging dates read from the mobyData metadata)
#' filtered <- filterDetections(rays)
#' filtered                     # prints the summary
#' head(filtered$data_discarded[, c("ID", "datetime", "reason")])
#'
#' # enable the short-interval false-detection filter (needs the transmitter nominal delay)
#' filtered2 <- filterDetections(rays, nominal.delay = 120)   # 120 s tags
#'
#' \donttest{
#' # add a movement-speed filter (slower: computes step distances); on a 3-animal subset
#' sub <- rays[rays$ID %in% head(levels(factor(rays$ID)), 3), ]
#' filtered3 <- filterDetections(sub, max.speed = 5, speed.unit = "km/h")
#' }
#' @seealso [plotMinLag()] to check the min_lag threshold empirically.
#' @export

filterDetections <- function(data,
                             tagging.dates = NULL,
                             cutoff.dates = NULL,
                             id.col = NULL,
                             datetime.col = NULL,
                             lon.col = NULL,
                             lat.col = NULL,
                             remove.duplicates = TRUE,
                             nominal.delay = NULL,
                             min.lag.factor = 30,
                             min.lag.threshold = NULL,
                             isolation.window = NULL,
                             max.speed = NULL,
                             speed.unit = "m/s",
                             acoustic.range = 600,
                             min.corroboration = 2L,
                             max.iterations = 20L,
                             min.detections = 0,
                             min.days = 0,
                             land.shape = NULL,
                             epsg.code = NULL,
                             verbose = getOption("moby.verbose", TRUE),
                             ...) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  start.time <- Sys.time()
  .printConsole("Filtering detections", verbose = verbose)

  land_shape_expr <- deparse(substitute(land.shape))
  prev_meta <- attr(data, "moby")

  # read nominal.delay from the mobyData metadata BEFORE .validateArguments (reader only; the upstream
  # population of meta$nominal.delay is handled elsewhere)
  if (is.null(nominal.delay) && !is.null(prev_meta) && !is.null(prev_meta$nominal.delay))
    nominal.delay <- prev_meta$nominal.delay

  # resolve NULL column/date/metadata arguments; coerce id.col to a factor
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  land.shape <- reviewed_params$land.shape
  cutoff.dates <- reviewed_params$cutoff.dates

  # scalar-argument validation (robust off-switches; no fragile `!= FALSE`)
  if (!speed.unit %in% c("m/s", "km/h")) .mobyAbort("Wrong speed unit: choose either 'm/s' or 'km/h'.")
  if (!isFALSE(isolation.window) && !is.null(isolation.window) &&
      (!is.numeric(isolation.window) || length(isolation.window) != 1 || is.na(isolation.window) || isolation.window < 0))
    .mobyAbort("'isolation.window' must be a single non-negative number (hours), NULL, or FALSE.")
  if (!is.numeric(min.corroboration) || length(min.corroboration) != 1 || is.na(min.corroboration) || min.corroboration < 1)
    .mobyAbort("'min.corroboration' must be a single integer >= 1.")
  min.corroboration <- as.integer(min.corroboration)
  if (is.null(max.iterations)) max.iterations <- Inf
  if (!is.numeric(max.iterations) || length(max.iterations) != 1 || is.na(max.iterations) || max.iterations < 1)
    .mobyAbort("'max.iterations' must be a positive number (or Inf).")

  # NA validation on the columns the filters rely on
  if (anyNA(data[, datetime.col]))
    .mobyAbort("'", datetime.col, "' contains ", sum(is.na(data[, datetime.col])),
               " missing value(s); please remove or fix them before filtering.")

  # station column (only used by the min_lag filter); resolved from metadata, defaulting to 'station'
  station_col <- if (!is.null(prev_meta) && !is.null(prev_meta$station.col)) prev_meta$station.col else "station"

  nfish <- nlevels(data[, id.col])
  all_ids <- levels(data[, id.col])

  # per-ID nominal.delay (min_lag): tolerant of unknown tags (NA -> skip that individual)
  do_minlag <- FALSE
  if (!is.null(nominal.delay)) {
    cp <- .checkAnimalParams(nominal.delay, "Nominal delay", expected_class = "numeric", data, id.col, allow.missing = TRUE)
    if (!is.null(cp$errors)) .mobyAbort(paste(cp$errors, collapse = "\n"))
    nominal.delay <- cp$vector
    do_minlag <- TRUE
    if (!station_col %in% names(data)) {
      .mobyWarn("The min_lag false-detection filter needs a station/receiver column ('", station_col,
                "'), which is not present; skipping it.")
      do_minlag <- FALSE
    }
  } else {
    .mobyInform("- No 'nominal.delay' supplied or found in metadata: the min_lag false-detection ",
                "filter is OFF. Supply the transmitter nominal delay (s) to enable it.", verbose = verbose)
  }

  # per-ID max.speed (speed filter): strict (all IDs must be covered by a named vector)
  do_speed <- !is.null(max.speed)
  if (do_speed) {
    cp <- .checkAnimalParams(max.speed, "Maximum speed", expected_class = "numeric", data, id.col)
    if (!is.null(cp$errors)) .mobyAbort(paste(cp$errors, collapse = "\n"))
    max.speed <- cp$vector
    if (anyNA(data[, lon.col]) || anyNA(data[, lat.col]))
      .mobyAbort("The speed filter requires complete coordinates, but '", lon.col, "'/'", lat.col,
                 "' contain missing values.")
  }

  # raw per-individual counts (before any filtering)
  n_individual <- as.integer(table(factor(data[, id.col], levels = all_ids)))
  names(n_individual) <- all_ids
  n_total <- nrow(data)

  # stable per-detection id used to track speed-filter shielding across iterations
  data$.detid <- seq_len(nrow(data))

  # accumulator for every removed row (each chunk carries a 'reason' and an internal '.stage')
  discarded <- list()
  add_discarded <- function(rows, stage) {
    if (is.null(rows) || nrow(rows) == 0) return(invisible(NULL))
    rows$.stage <- stage
    discarded[[length(discarded) + 1L]] <<- rows
  }


  ##############################################################################
  ## Stage 0: duplicate removal ################################################
  ##############################################################################

  if (isTRUE(remove.duplicates)) {
    key_cols <- c(id.col, datetime.col, if (station_col %in% names(data)) station_col)
    dup <- duplicated(data[, key_cols, drop = FALSE])
    if (any(dup)) {
      r <- data[dup, , drop = FALSE]; r$reason <- "duplicate detection"
      add_discarded(r, "Duplicates")
      data <- data[!dup, , drop = FALSE]
    }
  }


  ##############################################################################
  ## Prepare spatial objects (speed filter only) ###############################
  ##############################################################################

  if (do_speed) {
    coords <- sf::st_as_sf(data, coords = c(lon.col, lat.col))
    coords_crs <- .checkProjection(coords)
    if (is.null(land.shape) && coords_crs == "geographic" && is.null(epsg.code)) {
      data$tmp_lon <- data[, lon.col]; data$tmp_lat <- data[, lat.col]
    } else {
      spatial_data <- .processSpatial(coords, land.shape, epsg.code)
      coords <- spatial_data$coords; land.shape <- spatial_data$spatial.layer; epsg.code <- spatial_data$epsg.code
      data$tmp_lon <- sf::st_coordinates(coords)[, 1]; data$tmp_lat <- sf::st_coordinates(coords)[, 2]
    }
  }


  ##############################################################################
  ## Stages 1-5: per-individual temporal / short-interval filters ##############
  ##############################################################################

  .printConsole("Applying detection filters...", verbose = verbose)
  data_individual <- split(data, f = data[, id.col], drop = FALSE)
  pb <- .progressBar(nfish, verbose)

  for (i in seq_len(nfish)) {
    sub <- data_individual[[i]]
    if (nrow(sub) == 0) { .progressSet(pb, i); next }
    sub <- sub[order(sub[, datetime.col]), , drop = FALSE]
    tag_date <- tagging.dates[i]
    off_date <- if (!is.null(cutoff.dates)) cutoff.dates[i] else NULL

    # (1) pre-tagging
    prior <- which(sub[, datetime.col] < tag_date)
    if (length(prior) > 0) {
      r <- sub[prior, ]; r$reason <- "before tagging date"; add_discarded(r, "Before tagging")
      sub <- sub[-prior, , drop = FALSE]
    }

    # (2) cut-off
    if (!is.null(off_date)) {
      after <- which(sub[, datetime.col] > off_date)
      if (length(after) > 0) {
        r <- sub[after, ]; r$reason <- "after cut-off date"; add_discarded(r, "After cut-off")
        sub <- sub[-after, , drop = FALSE]
      }
    }

    # (3) false-detection (min_lag), per receiver, scaled to the nominal delay
    if (do_minlag && !is.na(nominal.delay[i]) && nrow(sub) > 0) {
      thr <- if (!is.null(min.lag.threshold)) min.lag.threshold else min.lag.factor * nominal.delay[i]
      fl <- .minLagFlags(as.numeric(sub[, datetime.col]), as.character(sub[[station_col]]), thr)
      if (any(fl)) {
        r <- sub[fl, ]; r$reason <- paste0("false detection (min_lag > ", round(thr), " s)")
        add_discarded(r, "False detection"); sub <- sub[!fl, , drop = FALSE]
      }
    }

    # (4) isolation (opt-in): gaps recomputed HERE, on the surviving detections only
    if (!isFALSE(isolation.window) && !is.null(isolation.window) && nrow(sub) > 1) {
      hd <- as.numeric(difftime(.lead(sub[, datetime.col]), sub[, datetime.col], units = "hours"))
      iso <- unique(c(which(hd > isolation.window & .lag(hd) > isolation.window),
                      which(is.na(hd) & .lag(hd) > isolation.window),
                      which(hd > isolation.window & is.na(.lag(hd)))))
      if (length(iso) > 0) {
        r <- sub[iso, ]; r$reason <- paste0("isolated (>", isolation.window, "h both sides)")
        add_discarded(r, "Isolated"); sub <- sub[-iso, , drop = FALSE]
      }
    }

    data_individual[[i]] <- sub
    .progressSet(pb, i)
  }
  .progressEnd(pb)

  data_filtered <- do.call("rbind", data_individual)
  data_filtered <- data_filtered[order(data_filtered[, datetime.col], data_filtered[, id.col]), , drop = FALSE]


  ##############################################################################
  ## Stage 6: speed filter (clustering + corroboration shield) #################
  ##############################################################################

  shielded_ids <- integer(0)
  if (do_speed && nrow(data_filtered) > 0) {
    .printConsole("Applying speed filter...", verbose = verbose)
    # build the least-cost graph ONCE, over the full set of surviving positions, and reuse it for every
    # call the speed filter makes below (a graph rebuilt from a small subset would span only that
    # subset's bounding box, silently degrading longer paths to straight lines)
    step <- .stepDistances(data_filtered, land.shape, epsg.code, id.col = id.col,
                           lon.col = "tmp_lon", lat.col = "tmp_lat", verbose = FALSE, ...)
    tracks <- step$data
    cost_graph <- step$cost.graph

    sf_out <- .applySpeedFilter(tracks, id.col = id.col, datetime.col = datetime.col, max.speed = max.speed,
                                speed.unit = speed.unit, acoustic.range = acoustic.range,
                                min.corroboration = min.corroboration, max.iterations = max.iterations,
                                land.shape = land.shape, epsg.code = epsg.code, cost.graph = cost_graph,
                                verbose = verbose, ...)
    data_filtered <- sf_out$data
    shielded_ids <- sf_out$shielded
    if (!is.null(sf_out$removed) && nrow(sf_out$removed) > 0) add_discarded(sf_out$removed, "Speed")
    if (sf_out$small_time_diff > 0)
      .mobyWarn(sf_out$small_time_diff, " overspeed removal(s) involved sub-1-minute gaps; these may be valid ",
                "detections flagged by clock drift or an under-estimated acoustic range - review 'data_discarded'.")
    if (sf_out$iterations >= max.iterations)
      .mobyWarn("The speed filter reached max.iterations (", max.iterations, ") before converging.")
  }


  ##############################################################################
  ## Stage 7: minimum detections / minimum days ################################
  ##############################################################################

  if (min.detections > 0 && nrow(data_filtered) > 0) {
    di <- split(data_filtered, f = data_filtered[, id.col], drop = FALSE)
    for (i in seq_along(di)) {
      if (nrow(di[[i]]) > 0 && nrow(di[[i]]) < min.detections) {
        r <- di[[i]]; r$reason <- paste0("individual with fewer than ", min.detections, " detections")
        add_discarded(r, "Min detections"); di[[i]] <- di[[i]][0, , drop = FALSE]
      }
    }
    data_filtered <- do.call("rbind", di)
  }

  if (min.days > 0 && nrow(data_filtered) > 0) {
    di <- split(data_filtered, f = data_filtered[, id.col], drop = FALSE)
    for (i in seq_along(di)) {
      if (nrow(di[[i]]) == 0) next
      days <- strftime(di[[i]][, datetime.col], "%Y-%m-%d", tz = .dataTZ(di[[i]][, datetime.col]))
      if (length(unique(days)) < min.days) {
        r <- di[[i]]; r$reason <- paste0("individual detected on fewer than ", min.days, " days")
        add_discarded(r, "Min days"); di[[i]] <- di[[i]][0, , drop = FALSE]
      }
    }
    data_filtered <- do.call("rbind", di)
  }


  ##############################################################################
  ## Assemble results ##########################################################
  ##############################################################################

  internal_cols <- c(".detid", "tmp_lon", "tmp_lat", "hour_diff", "dist_m", "dist_min", "speed",
                     "row_id", "row_index", "iteration")

  # filtered detections + qc_flag
  qc <- rep("valid", nrow(data_filtered))
  if (length(shielded_ids) > 0) qc[data_filtered$.detid %in% shielded_ids] <- "overspeed_review"
  data_filtered$qc_flag <- qc
  data_filtered <- .dropCols(data_filtered, internal_cols)
  rownames(data_filtered) <- NULL

  # discarded detections + reason (stable schema even when empty)
  data_discarded <- do.call(.rbindFill, discarded)
  if (is.null(data_discarded)) {
    data_discarded <- data_filtered[0, setdiff(names(data_filtered), "qc_flag"), drop = FALSE]
    data_discarded$reason <- character(0)
  } else {
    stages_of_row <- data_discarded$.stage
    data_discarded <- .dropCols(data_discarded, c(internal_cols, ".stage"))
    ord_cols <- c(setdiff(names(data_discarded), "reason"), "reason")
    data_discarded <- data_discarded[, ord_cols, drop = FALSE]
    data_discarded <- data_discarded[order(data_discarded[, id.col], data_discarded[, datetime.col]), , drop = FALSE]
    rownames(data_discarded) <- NULL
  }

  # per-individual summary: pivot the accumulated removals by stage
  stage_order <- c("Duplicates", "Before tagging", "After cut-off", "False detection",
                   "Isolated", "Speed", "Min detections", "Min days")
  all_disc <- do.call(.rbindFill, discarded)
  removed_mat <- matrix(0L, nrow = nfish, ncol = length(stage_order),
                        dimnames = list(all_ids, stage_order))
  if (!is.null(all_disc)) {
    tb <- table(factor(as.character(all_disc[, id.col]), levels = all_ids),
                factor(all_disc$.stage, levels = stage_order))
    removed_mat[, ] <- as.integer(tb)
  }
  active_stages <- stage_order[colSums(removed_mat) > 0 |
                                 stage_order %in% c("Before tagging")]                 # always show pre-tagging
  n_removed_ind <- rowSums(removed_mat)
  filter_summary <- data.frame(ID = all_ids, "Raw detections" = n_individual, check.names = FALSE,
                               row.names = NULL, stringsAsFactors = FALSE)
  for (st in active_stages) filter_summary[[st]] <- removed_mat[, st]
  pct <- sprintf("%.0f", ifelse(n_individual > 0, n_removed_ind / n_individual * 100, 0))
  filter_summary[["Total removed"]] <- paste0(n_removed_ind, " (", pct, "%)")
  filter_summary[["Total removed"]][n_removed_ind == 0 | n_individual == 0] <- "-"

  # console summary
  if (verbose) {
    n_removed_total <- if (is.null(all_disc)) 0L else nrow(all_disc)
    pct_total <- if (n_total > 0) sprintf("%.0f", n_removed_total / n_total * 100) else "0"
    .mobyInform("Detections removed = ", n_removed_total, " (", pct_total, "%) from a total of ", n_total,
                verbose = verbose)
    fmt_pct <- function(n) { p <- if (n_total > 0) n / n_total * 100 else 0
      if (n > 0 && p < 0.5) " (<1%)" else sprintf(" (%.0f%%)", p) }
    stage_run <- c(Duplicates = isTRUE(remove.duplicates), "Before tagging" = TRUE,
                   "After cut-off" = !is.null(cutoff.dates), "False detection" = do_minlag,
                   Isolated = (!isFALSE(isolation.window) && !is.null(isolation.window)),
                   Speed = do_speed, "Min detections" = min.detections > 0, "Min days" = min.days > 0)
    for (st in stage_order[stage_run[stage_order]]) {
      n <- sum(removed_mat[, st]); ids_st <- sum(removed_mat[, st] > 0)
      .mobyInform("  \u2022 ", st, ": ", n, fmt_pct(n), " from ", ids_st, " individual(s)", verbose = verbose)
    }
    n_flagged <- sum(data_filtered$qc_flag == "overspeed_review")
    if (n_flagged > 0)
      .mobyInform("  \u2022 flagged for review (overspeed, retained): ", n_flagged, verbose = verbose)
    ids_discarded <- sum(n_individual > 0 & table(factor(data_filtered[, id.col], levels = all_ids)) == 0)
    .mobyInform("Individuals fully discarded = ", ids_discarded, " from a total of ",
                sum(n_individual > 0), verbose = verbose)
  }
  .reportRuntime(start.time, verbose)

  # re-attach mobyData metadata so the filtered detections stay chainable
  if (!is.null(prev_meta)) {
    attr(data_filtered, "moby") <- prev_meta
    class(data_filtered) <- unique(c("mobyData", "data.frame"))
  }

  results <- list(data = data_filtered, data_discarded = data_discarded, summary = filter_summary)
  attr(results, "parameters") <- list(
    tagging.dates = tagging.dates, cutoff.dates = cutoff.dates, remove.duplicates = remove.duplicates,
    nominal.delay = nominal.delay, min.lag.factor = min.lag.factor, min.lag.threshold = min.lag.threshold,
    isolation.window = isolation.window, max.speed = max.speed, speed.unit = speed.unit,
    acoustic.range = acoustic.range, min.corroboration = min.corroboration, max.iterations = max.iterations,
    min.detections = min.detections, min.days = min.days,
    land.shape = if (!is.null(land.shape)) land_shape_expr else NULL, epsg.code = epsg.code,
    processing.date = Sys.time())
  class(results) <- c("mobyFilter", "list")
  results
}


#######################################################################################################
## Internal helpers ###################################################################################
#######################################################################################################

#' Per-receiver nearest-neighbour gap (min_lag), in seconds
#'
#' @description For time-ordered detections, returns each detection's min_lag: the time (seconds) to the
#'   nearest other detection of the same tag on the same receiver. `NA` where a detection has no
#'   same-receiver companion (a lone decode). Shared by the `filterDetections()` min_lag filter and the
#'   `plotMinLag()` diagnostic, so both compute min_lag from identical code.
#' @param times Numeric datetimes (seconds), in ascending order.
#' @param stations Character receiver/station id, aligned to `times`.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.minLagValues <- function(times, stations) {
  minlag <- rep(NA_real_, length(times))
  for (s in unique(stations)) {
    idx <- which(stations == s)                       # already in time order
    if (length(idx) >= 2) {
      t <- times[idx]
      minlag[idx] <- pmin(c(NA, diff(t)), c(diff(t), NA), na.rm = TRUE)
    }
  }
  minlag
}


#' Per-receiver short-interval (min_lag) false-detection flags
#'
#' @description For time-ordered detections, returns TRUE where a detection is a likely false
#'   (code-collision) decode: it has no same-receiver companion within `threshold` seconds (its min_lag
#'   exceeds the threshold, or there is no other detection of that tag on that receiver at all).
#' @param times Numeric datetimes (seconds), in ascending order.
#' @param stations Character receiver/station id, aligned to `times`.
#' @param threshold The short-interval threshold in seconds.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.minLagFlags <- function(times, stations, threshold) {
  minlag <- .minLagValues(times, stations)
  is.na(minlag) | minlag > threshold
}


#' Apply the maximum-speed filter with the corroboration-shielded clustering heuristic
#'
#' @description Iteratively flags steps exceeding `max.speed` (a per-individual vector). For each
#'   flagged step, the local +/-3-detection window is hierarchically clustered into two groups: an
#'   isolated spatial outlier (a cluster smaller than `min.corroboration`) is removed, while a
#'   corroborated group (>= `min.corroboration`) is SHIELDED (its detection ids are returned so the
#'   caller can flag them) and never deleted, preventing the cascading deletion of real relocations.
#'   Distances are recomputed only for individuals that changed. Returns
#'   `list(data, removed, shielded, small_time_diff, iterations)`.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.applySpeedFilter <- function(data, id.col, datetime.col, max.speed, speed.unit, acoustic.range,
                              min.corroboration, max.iterations, land.shape, epsg.code, cost.graph,
                              verbose, ...) {

  d_split <- split(data, f = data[, id.col], drop = FALSE)
  removed_list <- vector("list", length(d_split))
  shielded <- integer(0)                                  # .detid of shielded (corroborated) detections
  dirty <- rep(FALSE, length(d_split))                    # initial dist_m is valid; recompute after removals
  small_time_diff <- 0L
  iteration <- 1L

  repeat {
    if (verbose && interactive()) { cat(paste0("\rSpeed filter iteration: ", iteration)); utils::flush.console() }
    removed_this <- 0L; new_shield <- 0L

    for (i in seq_along(d_split)) {
      sub <- d_split[[i]]
      if (nrow(sub) < 2 || is.na(max.speed[i])) next
      if (dirty[i]) {
        sub <- suppressWarnings(.stepDistances(sub, land.shape, epsg.code, id.col = id.col,
                 lon.col = "tmp_lon", lat.col = "tmp_lat", cost.graph = cost.graph, verbose = FALSE, ...)$data)
        dirty[i] <- FALSE
      }
      dist_min <- sub$dist_m - (acoustic.range * 2); dist_min[dist_min < 0] <- 0
      hdt <- as.numeric(difftime(.lead(sub[, datetime.col]), sub[, datetime.col], units = "hours"))
      speed <- if (speed.unit == "m/s") dist_min / (hdt * 3600) else (dist_min / 1000) / hdt
      speed[!is.finite(speed)] <- 0                       # zero-gap or NA steps are never "overspeed"
      overspeed <- which(speed > max.speed[i])
      if (length(overspeed) == 0) { d_split[[i]] <- sub; next }

      to_remove <- integer(0)
      for (k in overspeed) {
        win <- (k - 3):(k + 3); win <- win[win > 0 & win <= nrow(sub)]
        if (length(win) < 2) next
        dm <- .calculateDistMatrix(sub[win, ], land.shape, epsg.code, id.col = id.col,
                                   lon.col = "tmp_lon", lat.col = "tmp_lat", cost.graph = cost.graph)
        if (anyNA(dm) || all(as.numeric(dm) == 0)) next   # cannot resolve a spatial outlier
        cl <- stats::cutree(stats::hclust(dm), k = 2)
        tab <- table(cl)
        minority <- win[cl == names(tab)[which.min(tab)]]
        min_ids <- sub$.detid[minority]
        if (all(min_ids %in% shielded)) next              # already handled -> guarantees termination
        if (length(minority) < min.corroboration) {
          to_remove <- c(to_remove, minority[1])          # isolated singleton outlier -> remove
        } else {
          new_shield <- new_shield + length(setdiff(min_ids, shielded))
          shielded <- union(shielded, min_ids)            # corroborated group -> shield (flag, retain)
        }
      }
      # a detection shielded by ANY window's clustering is never deleted, even if an earlier window
      # (which did not see its corroborating twin) queued it as a singleton -> shield wins over remove
      to_remove <- unique(to_remove)
      to_remove <- to_remove[!(sub$.detid[to_remove] %in% shielded)]
      if (length(to_remove) > 0) {
        # count sub-1-minute gaps among the ACTUALLY removed detections (clock-drift / range heuristic)
        small_time_diff <- small_time_diff + sum(hdt[to_remove] <= (1 / 60), na.rm = TRUE)
        rem <- sub[to_remove, , drop = FALSE]
        rem$reason <- paste0("above max speed (", signif(max.speed[i], 3), " ", speed.unit, ")")
        rem$iteration <- iteration
        removed_list[[i]] <- c(removed_list[[i]], list(rem))
        sub <- sub[-to_remove, , drop = FALSE]
        removed_this <- removed_this + length(to_remove)
        dirty[i] <- TRUE
      }
      d_split[[i]] <- sub
    }

    if ((removed_this == 0L && new_shield == 0L) || iteration >= max.iterations) break
    iteration <- iteration + 1L
  }
  if (verbose && interactive()) cat("\n")

  out_data <- do.call("rbind", d_split)
  removed <- do.call("rbind", lapply(removed_list, function(x) if (length(x)) do.call("rbind", x) else NULL))
  list(data = out_data, removed = removed, shielded = shielded, small_time_diff = small_time_diff,
       iterations = iteration)
}


#' Pairwise shortest-distance dissimilarity matrix for a small window of detections
#'
#' @description Builds an `n x n` dissimilarity (`dist`) object from the shortest (in-water, if a
#'   `land.shape` is supplied) distances between detections. Only the lower triangle is computed - the
#'   half that `as.dist()` reads - halving the number of `calculateStepDistances` calls.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.calculateDistMatrix <- function(data, land.shape, epsg.code, id.col, lon.col, lat.col, cost.graph, ...) {
  n <- nrow(data)
  m <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1L)) {                          # lower triangle only (as.dist reads j < i)
      d <- suppressWarnings(.stepDistances(data[c(i, j), ], land.shape, epsg.code, id.col = id.col,
             lon.col = lon.col, lat.col = lat.col, cost.graph = cost.graph, verbose = FALSE, ...)$data)
      m[i, j] <- d$dist_m[1]
    }
  }
  stats::as.dist(m)
}


#' Print method for filterDetections() results
#'
#' @description Prints a compact overview of a `mobyFilter` object: kept/removed totals and the
#'   per-individual summary table. Inspect `$data`, `$data_discarded` and `$summary` for details.
#' @param x A `mobyFilter` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @keywords internal
#' @export
print.mobyFilter <- function(x, ...) {
  n_kept <- nrow(x$data); n_disc <- nrow(x$data_discarded)
  n_flag <- if ("qc_flag" %in% names(x$data)) sum(x$data$qc_flag == "overspeed_review") else 0L
  cat("<mobyFilter> filtered acoustic detections\n")
  cat("  ", n_kept, " retained | ", n_disc, " removed", sep = "")
  if (n_flag > 0) cat(" | ", n_flag, " flagged for review", sep = "")
  cat("\n")
  if (n_disc > 0) {
    tb <- sort(table(x$data_discarded$reason), decreasing = TRUE)
    for (nm in names(tb)) cat(sprintf("    - %-40s %d\n", nm, tb[[nm]]))
  }
  cat("  Per-individual summary:\n")
  print(x$summary, row.names = FALSE)
  cat("  Inspect $data, $data_discarded and attr(, \"parameters\").\n")
  invisible(x)
}

#######################################################################################################
#######################################################################################################
