#######################################################################################################
## Residence events (visits) ##########################################################################
#######################################################################################################

#' Extract residence events (visits)
#'
#' @description Segments each individual's time-ordered detections into discrete **residence events**
#' (also called *visits*): episodes of continuous presence at a location, bounded by arrival and
#' departure. A visit ends when the animal moves to a different location \emph{or} when the gap since
#' its previous detection exceeds `max.gap` — i.e. a long absence is treated as the animal having left
#' and (possibly) returned, rather than as one uninterrupted stay.
#'
#' This is the event-resolution counterpart to \code{\link{calculateResidency}} (which measures
#' day-scale residency as a presence/absence ratio) and the residence layer underlying
#' \code{\link{calculateTransitions}} (which summarises the same visits as movement-network nodes).
#' Both functions derive their residence segmentation from the same internal engine, so counts and
#' durations are consistent across them.
#'
#' @details The `max.gap` threshold is the one biological choice here: it encodes how long an absence
#' you treat as "left and came back" versus "still present but temporarily undetected". There is no
#' universal value — it depends on the species, the detection range and the array layout — so it is
#' exposed as a first-class, reportable parameter. A principled way to pick it is to inspect the
#' distribution of intervals between successive detections, which is often bimodal, with a natural
#' trough separating within-visit gaps (brief range dropouts) from between-visit gaps (true
#' excursions). Set `max.gap = Inf` to disable gap-splitting entirely (a visit then ends only on a
#' change of location).
#'
#' @inheritParams as_moby
#' @param data A data frame (or \code{\link{mobyData}}) of detections.
#' @param spatial.col Name of the column defining the locations visited (e.g. receiver, station,
#' habitat, region). If `NULL` (default), it is taken from the `mobyData` station column (or the
#' canonical `"station"`); set it to track visits to any other spatial unit.
#' @param id.groups Optional named list of ID groups; when supplied, visits are computed within each
#' group and the group is carried in the `group` column.
#' @param max.gap Maximum tolerated gap between successive detections within a single visit. Absences
#' longer than this split the sequence into separate residence events. Defaults to 48 (hours). Use
#' `Inf` to segment on location changes only.
#' @param max.gap.unit Units of `max.gap`: one of `"hours"` (default), `"days"`, `"mins"`, `"secs"`.
#'
#' @return A data frame with one row per residence event (visit), ordered by group, individual and
#' arrival time:
#' \item{group}{the ID group (`"all"` when `id.groups` is not supplied).}
#' \item{id}{the individual.}
#' \item{site}{the location visited (`spatial.col` value).}
#' \item{arrival, departure}{first and last detection times of the visit (POSIXct).}
#' \item{residence_h}{visit duration in hours (`departure - arrival`; 0 for a single-detection visit).}
#' \item{n_detections}{number of detections recorded during the visit.}
#'
#' @references
#' Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., & Abecasis, D. (2023).
#' Residency and space use estimation methods based on passive acoustic telemetry data.
#' Movement Ecology, 11(1), 12.
#'
#' See also the residence-event extractors \code{glatos::residence_events} and
#' \code{VTrack::RunResidenceExtraction} for the same concept in neighbouring toolkits.
#'
#' @seealso \code{\link{calculateResidency}}, \code{\link{calculateTransitions}}
#'
#' @examples
#' data(rays)
#' # discrete visits to each receiver station (48 h gap threshold)
#' visits <- calculateVisits(rays, spatial.col = "station")
#' head(visits)
#'
#' @export

calculateVisits <- function(data,
                            spatial.col = NULL,
                            id.col = NULL,
                            datetime.col = NULL,
                            id.groups = NULL,
                            max.gap = 48,
                            max.gap.unit = c("hours", "days", "mins", "secs")) {

  ##############################################################################
  ## Validate ##################################################################
  ##############################################################################

  # the location column defaults to the dataset's station column (from the mobyData metadata, else
  # the canonical "station"); pass spatial.col explicitly to track visits to any other unit
  meta <- attr(data, "moby")
  if (is.null(spatial.col))
    spatial.col <- if (!is.null(meta) && !is.null(meta$station.col)) meta$station.col else "station"

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  max.gap.unit <- match.arg(max.gap.unit)
  max.gap.secs  <- .gapToSecs(max.gap, max.gap.unit)
  if (is.finite(max.gap.secs))
    message("- Defining residence events with max.gap = ", max.gap, " ", max.gap.unit,
            " (a longer absence starts a new visit); tune/justify per system.")

  # group label per row (single "all" group unless id.groups is supplied)
  if (is.null(id.groups)) {
    data$.group <- "all"
    group_levels <- "all"
  } else {
    grp_map <- stats::setNames(rep(names(id.groups), lengths(id.groups)), as.character(unlist(id.groups)))
    data$.group <- unname(grp_map[as.character(data[, id.col])])
    group_levels <- names(id.groups)
  }

  ##############################################################################
  ## Per-group segmentation ####################################################
  ##############################################################################

  out <- list()
  for (g in group_levels) {
    dg <- data[data$.group == g, , drop = FALSE]
    if (nrow(dg) == 0) next
    ev <- .residenceRuns(dg, spatial.col, id.col, datetime.col, max.gap.secs)
    if (nrow(ev) == 0) next
    out[[g]] <- cbind(group = g, ev, stringsAsFactors = FALSE)
  }

  tz <- .dataTZ(data[, datetime.col])
  result <- if (length(out)) do.call(rbind, out) else
    data.frame(group = character(0), id = character(0), site = character(0),
               arrival = as.POSIXct(character(0), tz = tz),
               departure = as.POSIXct(character(0), tz = tz),
               residence_h = numeric(0), n_detections = integer(0),
               stringsAsFactors = FALSE)
  result <- result[order(result$group, result$id, result$arrival), , drop = FALSE]
  rownames(result) <- NULL
  result
}


################################################################################
# Internal helpers #############################################################
################################################################################

#' Convert a max.gap value + units into seconds
#' @keywords internal
#' @noRd
.gapToSecs <- function(x, units) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0)
    stop("'max.gap' must be a single positive number (or Inf).", call. = FALSE)
  x * switch(units, secs = 1, mins = 60, hours = 3600, days = 86400)
}

#' Segment detections into residence events (visits).
#'
#' Shared segmentation engine used by both \code{calculateVisits()} and \code{calculateTransitions()}
#' so the two agree on what a visit is. A new run starts whenever the location changes OR the gap to
#' the previous detection exceeds `max.gap.secs` (Inf = location changes only). Operates on one group
#' of detections; callers handle id.groups.
#'
#' @return A data frame (one row per run) with columns id, site, arrival, departure, n_detections,
#'   residence_h, ordered by (id, arrival).
#' @keywords internal
#' @noRd
.residenceRuns <- function(d, spatial.col, id.col, datetime.col, max.gap.secs) {
  tz  <- .dataTZ(d[, datetime.col])
  ids <- as.character(d[, id.col])
  uids <- unique(ids)
  out <- vector("list", length(uids))
  names(out) <- uids

  for (id in uids) {
    di <- d[ids == id, , drop = FALSE]
    di <- di[order(di[, datetime.col]), , drop = FALSE]
    s  <- as.character(di[, spatial.col])
    t  <- di[, datetime.col]
    n  <- length(t)

    gap     <- c(Inf, as.numeric(difftime(t[-1], t[-n], units = "secs")))
    new_run <- c(TRUE, s[-1] != s[-n]) | (gap > max.gap.secs)
    run_id  <- cumsum(new_run)

    first_i <- !duplicated(run_id)                       # first detection of each run (run order)
    last_i  <- !duplicated(run_id, fromLast = TRUE)      # last detection of each run
    arrival   <- t[first_i]
    departure <- t[last_i]

    out[[id]] <- data.frame(
      id           = id,
      site         = s[first_i],
      arrival      = arrival,
      departure    = departure,
      n_detections = as.integer(tabulate(run_id)),
      residence_h  = as.numeric(difftime(departure, arrival, units = "hours")),
      stringsAsFactors = FALSE)
  }

  res <- do.call(rbind, out)
  if (is.null(res))
    return(data.frame(id = character(0), site = character(0),
                      arrival = as.POSIXct(character(0), tz = tz),
                      departure = as.POSIXct(character(0), tz = tz),
                      n_detections = integer(0), residence_h = numeric(0),
                      stringsAsFactors = FALSE))
  res <- res[order(res$id, res$arrival), , drop = FALSE]
  rownames(res) <- NULL
  res
}
