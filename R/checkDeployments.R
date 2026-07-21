#######################################################################################################
## Receiver-deployment / metadata quality control ####################################################
#######################################################################################################

# distance (in metres) between two coordinate pairs, using great-circle distance when the values
# look geographic (degrees) and Euclidean distance otherwise (already projected, metres for the
# metric CRSs moby works with)
.coordDist <- function(lon1, lat1, lon2, lat2) {
  geographic <- all(abs(c(lon1, lon2)) <= 180, na.rm = TRUE) && all(abs(c(lat1, lat2)) <= 90, na.rm = TRUE)
  if (geographic) {
    geosphere::distHaversine(cbind(lon1, lat1), cbind(lon2, lat2))
  } else {
    sqrt((lon1 - lon2)^2 + (lat1 - lat2)^2)
  }
}

# render a distance in metres for a human-readable QC message: metres below 1 km, kilometres above.
# .coordDist already returns metres (great-circle for lon/lat, or a metric projection), so the QC
# messages report an interpretable unit rather than a bare "units" count.
.formatDistance <- function(d) {
  if (length(d) != 1 || is.na(d)) return("an undetermined distance")
  if (d >= 1000) sprintf("%.1f km", d / 1000) else sprintf("%.0f m", d)
}

# fill a cleaned recover date: for a receiver's non-final deployments missing 'recover',
# use the next deployment's start; for the final one, use the current time (still active)
.cleanRecover <- function(dep) {
  dep$recover_clean <- dep$recover
  for (rec in unique(dep$receiver)) {
    idx <- which(dep$receiver == rec)
    n <- length(idx)
    if (n > 1) {
      for (i in seq_len(n - 1)) {
        if (is.na(dep$recover[idx[i]])) dep$recover_clean[idx[i]] <- dep$deploy[idx[i + 1]]
      }
    }
    if (is.na(dep$recover[idx[n]])) dep$recover_clean[idx[n]] <- Sys.time()
  }
  dep
}

.qcRow <- function(type, receiver = NA, station = NA, first = as.POSIXct(NA),
                   last = as.POSIXct(NA), n_detections = NA_integer_,
                   n_individuals = NA_integer_, details = NA_character_) {
  data.frame(type = type, receiver = receiver, station = station,
             first = first, last = last, n_detections = n_detections,
             n_individuals = n_individuals, details = details,
             stringsAsFactors = FALSE)
}


#' Check receiver deployment metadata (quality control)
#'
#' @description Quality-controls a receiver-deployment / station log and (optionally) cross-checks it
#' against a detection dataset, flagging the kinds of issues that commonly corrupt acoustic
#' telemetry analyses. This complements \code{\link{filterDetections}} (which cleans the
#' detections themselves) by auditing the *metadata*. It generalises the receiver-log curation
#' step that is otherwise done by hand. The function only *reports* problems (it never
#' edits data or aborts on data issues), so it is a safe, explicit preprocessing step; the companion
#' \code{\link{matchDeployments}} applies the corrections.
#'
#' Metadata-internal checks: missing deploy/recover dates; invalid date ranges
#' (recover before deploy); overlapping deployments on the same receiver; duplicate
#' deployment records; missing/implausible coordinates; inconsistent station naming
#' (one station name with divergent coordinates, or near-identical coordinates under
#' different names); and gaps in a receiver's temporal coverage.
#'
#' Detection-vs-metadata checks (when `detections` is supplied): receivers present in the
#' detections but absent from the metadata; detections occurring before a receiver's first
#' deployment, after its last recovery, or within a gap between deployments; and station-name
#' mismatches (a known receiver associated with a station it was never deployed at).
#'
#' @param deployments A receiver-deployment data frame, e.g. from \code{\link{importDeployments}},
#' with `receiver` and `station` columns plus deployment / recovery date-times and (optionally)
#' longitude / latitude. The date and coordinate columns may carry non-canonical names, resolved via
#' `deploy.col`/`recover.col`/`lon.col`/`lat.col`.
#' @param detections Optional. A detection dataset (`mobyData` or data frame) with `receiver`,
#' `station` and a date-time column, used for the detection-vs-metadata checks.
#' @param id.col,datetime.col,station.col Column names in `detections` (resolved from the
#' `mobyData` metadata or canonical defaults when `NULL`).
#' @param deploy.col,recover.col,lon.col,lat.col Names of the deployment / recovery date-time and
#' the longitude / latitude columns in `deployments`. Default to the canonical
#' `"deploy"`/`"recover"`/`"lon"`/`"lat"` (as produced by \code{\link{importDeployments}}); set them
#' when a log uses other names (e.g. `deploy.col = "deploy_date"`, `lon.col = "Longitude"`), so a
#' non-canonical log need not be renamed first. The `receiver` and `station` columns are the
#' canonical identifiers this function keys on and are always taken as-is.
#' @param checks Character vector selecting which check groups to run (any of, or `"all"`, the
#' default): `"dates"` (deploy/recover date integrity), `"overlaps"` (duplicate and overlapping
#' deployment records), `"gaps"` (coverage gaps between consecutive deployments - often benign
#' servicing, so the easiest to drop), `"coordinates"` (implausible coordinates and
#' station-vs-coordinate naming consistency; skipped if `lon`/`lat` are absent) and `"detections"`
#' (cross-check detections against deployment windows; requires the `detections` argument).
#' @param scope Character; how much of the deployment log the *metadata-internal* checks cover.
#' `"all"` (default) audits every receiver in `deployments`. `"detected"` restricts those checks to
#' receivers that appear in `detections` (i.e. that recorded at least one detection), which removes
#' warnings about receivers that are irrelevant to the data being analysed - typically the bulk of
#' the coverage-gap and duplicate-station noise. Each retained receiver keeps its full deployment
#' timeline, so gap and overlap logic stays correct. `scope` only affects the metadata-internal
#' checks; the `"detections"` cross-checks are inherently detection-scoped and always run against the
#' complete metadata (so a receiver present in the detections but missing from the log is still
#' flagged). `"detected"` needs `detections`; without it, the function warns and audits everything.
#' @param coord.tolerance Numeric. Distance (in the coordinate units of the data; metres for
#' geographic coordinates) beyond which coordinates sharing a station name are flagged as
#' inconsistent, and below which differently-named stations are flagged as possible duplicates.
#' Defaults to 500.
#' @param gap.tolerance Numeric. Minimum gap (in days) between consecutive deployments of a
#' receiver to report as a coverage gap. Defaults to 1.
#' @param verbose Logical; print a summary to the console. Defaults to TRUE.
#'
#' @return An object of class `mobyQC`: a list with
#' \item{report}{A tidy data frame of flagged issues (`type`, `receiver`, `station`, `first`,
#' `last`, `n_detections`, `n_individuals`, `details`). Columns that do not apply to a given issue
#' are `NA` (kept typed, so the frame stays usable programmatically); for a more readable manual
#' export, write with an explicit placeholder, e.g. `write.csv(x$report, "report.csv",
#' row.names = FALSE, na = "-")`.}
#' \item{deployments}{The input deployments with a cleaned `recover_clean` column.}
#' \item{counts}{Named integer vector of issue counts by type.}
#'
#' @seealso \code{\link{importDeployments}}, \code{\link{importDetections}}, \code{\link{filterDetections}}
#' @examples
#' # quality-control a receiver-deployment log
#' checkDeployments(rays_deployments)
#'
#' # also cross-check the detections against the deployment windows
#' data(rays)
#' checkDeployments(rays_deployments, detections = rays)
#'
#' # restrict the metadata checks to receivers that actually recorded detections
#' checkDeployments(rays_deployments, detections = rays, scope = "detected")
#'
#' @export

checkDeployments <- function(deployments,
                          detections = NULL,
                          id.col = NULL,
                          datetime.col = NULL,
                          station.col = NULL,
                          deploy.col = "deploy",
                          recover.col = "recover",
                          lon.col = "lon",
                          lat.col = "lat",
                          checks = "all",
                          scope = c("all", "detected"),
                          coord.tolerance = 500,
                          gap.tolerance = 1,
                          verbose = TRUE) {

  scope <- match.arg(scope)
  dep <- as.data.frame(deployments)
  required <- c("receiver", "station", deploy.col)
  miss <- setdiff(required, colnames(dep))
  if (length(miss) > 0) {
    stop(paste0("'deployments' is missing required column(s): ", paste(miss, collapse = ", "),
                ". See importDeployments()."), call. = FALSE)
  }
  # standardise the deploy/recover/coordinate columns to the canonical names used by the checks below,
  # so a log whose columns are named differently need not be renamed by the user first
  dep$deploy <- dep[[deploy.col]]
  dep$recover <- if (recover.col %in% colnames(dep)) dep[[recover.col]] else as.POSIXct(NA)
  has_coords <- all(c(lon.col, lat.col) %in% colnames(dep))
  # an explicitly named coordinate column that is absent is almost always a typo, not "no coordinates"
  if (!has_coords && (!identical(lon.col, "lon") || !identical(lat.col, "lat")) &&
      any(c(lon.col, lat.col) %in% colnames(dep)))
    message("- coordinate column(s) '", paste(setdiff(c(lon.col, lat.col), colnames(dep)), collapse = "', '"),
            "' not found in 'deployments'; skipping the coordinate checks.")
  if (has_coords) { dep$lon <- dep[[lon.col]]; dep$lat <- dep[[lat.col]] }

  # resolve which check groups to run (default "all"); "detections" needs a detection dataset
  checks <- match.arg(checks, c("all", "dates", "overlaps", "gaps", "coordinates", "detections"), several.ok = TRUE)
  if ("all" %in% checks) checks <- c("dates", "overlaps", "gaps", "coordinates", "detections")
  do_dates <- "dates" %in% checks; do_overlaps <- "overlaps" %in% checks
  do_gaps <- "gaps" %in% checks; do_coords <- "coordinates" %in% checks
  do_det <- "detections" %in% checks
  if (do_det && is.null(detections))
    message("- 'detections' checks requested but no 'detections' supplied; skipping the detection-vs-metadata checks.")

  dep$receiver <- as.character(dep$receiver)
  dep$station <- as.character(dep$station)
  dep <- dep[order(dep$receiver, dep$deploy), , drop = FALSE]
  dep <- .cleanRecover(dep)

  report <- list()
  add <- function(r) report[[length(report) + 1]] <<- r

  # Scope of the metadata-internal checks (dates / overlaps / coordinates / gaps). With
  # scope = "detected" they run only on receivers that appear in the detections, i.e. that recorded
  # at least one detection, so issues affecting receivers absent from the data (a common source of
  # benign coverage-gap and duplicate-station noise) are not reported. The whole receiver timeline is
  # kept for each retained receiver, so gap/overlap logic stays correct. The detection-vs-metadata
  # checks further below are inherently detection-scoped and always run on the FULL metadata.
  detected_receivers <- NULL
  if (!is.null(detections)) {
    det_rec <- as.data.frame(detections)[["receiver"]]
    if (!is.null(det_rec)) detected_receivers <- unique(as.character(det_rec))
  }
  if (scope == "detected" && is.null(detected_receivers))
    message("- scope = \"detected\" needs a 'detections' dataset with a 'receiver' column; ",
            "auditing all deployments instead.")
  dep_scope <- if (scope == "detected" && !is.null(detected_receivers))
    dep[dep$receiver %in% detected_receivers, , drop = FALSE] else dep

  ##############################################################################
  ## Metadata-internal checks ##################################################
  ##############################################################################

  if (do_dates) {
  # missing deploy dates
  for (i in which(is.na(dep_scope$deploy))) {
    add(.qcRow("Missing deploy date", dep_scope$receiver[i], dep_scope$station[i],
               details = "Deployment date is missing."))
  }

  # missing recover dates (non-final deployments of a receiver)
  for (rec in unique(dep_scope$receiver)) {
    idx <- which(dep_scope$receiver == rec)
    if (length(idx) > 1) {
      for (i in idx[-length(idx)]) {
        if (is.na(dep_scope$recover[i])) {
          add(.qcRow("Missing recover date", rec, dep_scope$station[i], first = dep_scope$deploy[i],
                     details = "Non-final deployment is missing a recovery date."))
        }
      }
    }
  }

  # invalid date ranges
  bad_range <- which(!is.na(dep_scope$recover) & !is.na(dep_scope$deploy) & dep_scope$recover < dep_scope$deploy)
  for (i in bad_range) {
    add(.qcRow("Invalid date range", dep_scope$receiver[i], dep_scope$station[i],
               first = dep_scope$deploy[i], last = dep_scope$recover[i],
               details = "Recovery date precedes deployment date."))
  }
  }  # end 'dates' checks

  if (do_overlaps) {
  # duplicate deployment records
  dup_keys <- paste(dep_scope$receiver, dep_scope$station, dep_scope$deploy, sep = "_|_")
  for (k in unique(dup_keys[duplicated(dup_keys)])) {
    i <- which(dup_keys == k)[1]
    add(.qcRow("Duplicate deployment", dep_scope$receiver[i], dep_scope$station[i], first = dep_scope$deploy[i],
               details = paste0("Deployment record appears ", sum(dup_keys == k), " times.")))
  }

  # overlapping deployments on the same receiver
  for (rec in unique(dep_scope$receiver)) {
    idx <- which(dep_scope$receiver == rec)
    if (length(idx) > 1) {
      for (i in seq_len(length(idx) - 1)) {
        curr <- idx[i]; nxt <- idx[i + 1]
        if (!is.na(dep_scope$deploy[nxt]) && dep_scope$deploy[nxt] < dep_scope$recover_clean[curr] &&
            !is.na(dep_scope$recover[curr])) {
          add(.qcRow("Overlapping deployments", rec,
                     paste(dep_scope$station[curr], "/", dep_scope$station[nxt]),
                     first = dep_scope$deploy[nxt], last = dep_scope$recover[curr],
                     details = "Next deployment starts before the previous one ends."))
        }
      }
    }
  }

  }  # end 'overlaps' checks

  # implausible / missing coordinates + station<->coordinate naming consistency
  if (do_coords && has_coords) {
    lon <- suppressWarnings(as.numeric(dep_scope$lon)); lat <- suppressWarnings(as.numeric(dep_scope$lat))
    bad_coord <- which(is.na(lon) | is.na(lat) | lon == 0 | lat == 0 |
                       abs(lon) > 180 | abs(lat) > 90)
    for (i in bad_coord) {
      add(.qcRow("Implausible coordinates", dep_scope$receiver[i], dep_scope$station[i], first = dep_scope$deploy[i],
                 details = "Coordinate is missing, zero, or outside valid bounds."))
    }
    # inconsistent station naming: same station name -> divergent coordinates
    for (st in unique(dep_scope$station)) {
      idx <- which(dep_scope$station == st & !is.na(lon) & !is.na(lat))
      if (length(idx) > 1) {
        d <- .coordDist(lon[idx], lat[idx], rep(lon[idx[1]], length(idx)), rep(lat[idx[1]], length(idx)))
        if (any(d > coord.tolerance, na.rm = TRUE)) {
          add(.qcRow("Inconsistent station coordinates", paste(unique(dep_scope$receiver[idx]), collapse = " | "),
                     st, details = sprintf("Coordinates for this station span %s (> %s tolerance).",
                                           .formatDistance(max(d, na.rm = TRUE)), .formatDistance(coord.tolerance))))
        }
      }
    }
    # possible duplicate station names: near-identical coordinates under different names
    uniq <- dep_scope[!is.na(lon) & !is.na(lat), c("station", "lon", "lat")]
    uniq <- uniq[!duplicated(uniq$station), , drop = FALSE]
    if (nrow(uniq) > 1) {
      for (i in seq_len(nrow(uniq) - 1)) for (j in (i + 1):nrow(uniq)) {
        dd <- .coordDist(uniq$lon[i], uniq$lat[i], uniq$lon[j], uniq$lat[j])
        if (!is.na(dd) && dd < coord.tolerance) {
          add(.qcRow("Possible duplicate stations", NA,
                     paste(uniq$station[i], "/", uniq$station[j]),
                     details = sprintf("Different station names %s apart (< %s tolerance).",
                                       .formatDistance(dd), .formatDistance(coord.tolerance))))
        }
      }
    }
  }

  # coverage gaps in a receiver's activity
  if (do_gaps) {
  for (rec in unique(dep_scope$receiver)) {
    idx <- which(dep_scope$receiver == rec)
    if (length(idx) > 1) {
      for (i in seq_len(length(idx) - 1)) {
        gap <- as.numeric(difftime(dep_scope$deploy[idx[i + 1]], dep_scope$recover_clean[idx[i]], units = "days"))
        if (!is.na(gap) && gap > gap.tolerance) {
          add(.qcRow("Coverage gap", rec, dep_scope$station[idx[i]],
                     first = dep_scope$recover_clean[idx[i]], last = dep_scope$deploy[idx[i + 1]],
                     details = sprintf("%.0f-day gap before the next deployment.", gap)))
        }
      }
    }
  }
  }  # end 'gaps' checks

  ##############################################################################
  ## Detection-vs-metadata checks ##############################################
  ##############################################################################

  if (do_det && !is.null(detections)) {
    det <- as.data.frame(detections)
    args <- .resolveArgs(detections, list(id.col = id.col, datetime.col = datetime.col, station.col = station.col))
    id.col <- args$id.col; datetime.col <- args$datetime.col; station.col <- args$station.col
    if (!"receiver" %in% colnames(det)) stop("'detections' must contain a 'receiver' column.", call. = FALSE)
    if (!datetime.col %in% colnames(det)) stop(paste0("Datetime column ('", datetime.col, "') not found in 'detections'."), call. = FALSE)
    det$receiver <- as.character(det$receiver)
    det_station <- if (station.col %in% colnames(det)) as.character(det[[station.col]]) else rep(NA_character_, nrow(det))
    det_time <- det[[datetime.col]]
    det_id <- if (id.col %in% colnames(det)) as.character(det[[id.col]]) else as.character(det$receiver)

    # match each detection to a valid receiver+station+time window
    matched <- rep(FALSE, nrow(det))
    for (i in seq_len(nrow(dep))) {
      in_win <- det$receiver == dep$receiver[i] &
                (is.na(det_station) | det_station == dep$station[i]) &
                det_time >= dep$deploy[i] & det_time <= dep$recover_clean[i]
      in_win[is.na(in_win)] <- FALSE
      matched[in_win] <- TRUE
    }

    receivers_meta <- unique(dep$receiver)
    meta_pairs <- unique(paste(dep$receiver, dep$station, sep = "_|_"))
    unmatched <- which(!matched)

    if (length(unmatched) > 0) {
      # categorise each unmatched detection
      cat_type <- character(length(unmatched))
      for (k in seq_along(unmatched)) {
        i <- unmatched[k]
        rec <- det$receiver[i]; stn <- det_station[i]; tm <- det_time[i]
        if (!rec %in% receivers_meta) {
          cat_type[k] <- "Receiver missing from metadata"
        } else if (!is.na(stn) && !(paste(rec, stn, sep = "_|_") %in% meta_pairs)) {
          cat_type[k] <- "Station name mismatch"
        } else {
          sub <- dep[dep$receiver == rec & (is.na(stn) | dep$station == stn), ]
          if (nrow(sub) == 0) sub <- dep[dep$receiver == rec, ]
          if (!is.na(tm) && tm < min(sub$deploy, na.rm = TRUE)) cat_type[k] <- "Detection before deployment"
          else if (!is.na(tm) && tm > max(sub$recover_clean, na.rm = TRUE)) cat_type[k] <- "Detection after recovery"
          else cat_type[k] <- "Detection in deployment gap"
        }
      }
      # aggregate unmatched detections by type + receiver + station
      agg <- data.frame(type = cat_type, receiver = det$receiver[unmatched],
                        station = det_station[unmatched], time = det_time[unmatched],
                        id = det_id[unmatched], stringsAsFactors = FALSE)
      grp <- interaction(agg$type, agg$receiver, agg$station, drop = TRUE, sep = "_|_")
      for (g in levels(grp)) {
        rows <- agg[grp == g, ]
        add(.qcRow(rows$type[1], rows$receiver[1], rows$station[1],
                   first = min(rows$time, na.rm = TRUE), last = max(rows$time, na.rm = TRUE),
                   n_detections = nrow(rows), n_individuals = length(unique(rows$id)),
                   details = "Detections not matched to any valid deployment window."))
      }
    }
  }

  ##############################################################################
  ## Assemble result ###########################################################
  ##############################################################################

  report_df <- if (length(report) > 0) do.call(rbind, report) else .qcRow(NA_character_)[0, , drop = FALSE]
  if (nrow(report_df) > 0) report_df <- report_df[order(report_df$type, report_df$receiver), , drop = FALSE]
  rownames(report_df) <- NULL
  counts <- if (nrow(report_df) > 0) table(report_df$type) else integer(0)

  result <- list(report = report_df, deployments = dep, counts = counts)
  class(result) <- "mobyQC"

  if (verbose) print(result)
  invisible(result)
}


#' Match detections to receiver deployments and back-fill metadata
#'
#' @description Joins a receiver-deployment log (see \code{\link{importDeployments}}) to a
#' detection dataset, matching each detection to the deployment window
#' (`receiver` + time within `[deploy, recover]`) it falls in, and back-filling missing
#' coordinates and station names from the receiver log. When a receiver was moved between
#' stations and a detection matches more than one window, the window whose station name agrees
#' with the detection is preferred. Coordinate disagreements between the detection export and the
#' metadata (a common VUE artefact) are flagged. This applies the corrections for the issues
#' reported by \code{\link{checkDeployments}}, and is a natural follow-up to it.
#'
#' @param detections A detection dataset (`mobyData` or data frame) with a `receiver` column.
#' @param deployments A receiver-deployment data frame from \code{\link{importDeployments}}
#' (`receiver`, `station`, `lon`, `lat`, `deploy`, and optionally `recover`).
#' @param datetime.col,station.col,lon.col,lat.col Column names in `detections` (resolved from
#' the `mobyData` metadata or canonical defaults when `NULL`).
#' @param deploy.col,recover.col Names of the deployment and recovery date-time columns in
#' `deployments`. Default to the canonical `"deploy"`/`"recover"`; set them when a log uses other
#' names (e.g. `"deploy_date"`).
#' @param coord.tolerance Numeric. Distance (metres for geographic coordinates) above which a
#' detection's own coordinates are flagged as disagreeing with the deployment metadata. The
#' metadata coordinates are treated as authoritative for back-filling. Defaults to 500.
#' @param fill.coords,fill.station Logical; back-fill missing detection coordinates / station
#' names from the matched deployment. Default TRUE.
#' @param drop.unmatched Logical; drop detections that fall outside every deployment window.
#' Defaults to FALSE (they are retained and flagged).
#' @param verbose Logical; print a short summary. Defaults to TRUE.
#'
#' @return A \code{\link{mobyData}} object (the detections) with back-filled `lon`/`lat`/station
#' where missing, plus two logical flag columns: `deployment_matched` (detection fell within a
#' valid window) and `coord_mismatch` (own vs metadata coordinates differ beyond
#' `coord.tolerance`).
#'
#' @seealso \code{\link{checkDeployments}}, \code{\link{importDeployments}}, \code{\link{importDetections}}
#' @examples
#' # match detections to deployment windows, back-filling station/coordinates
#' matched <- matchDeployments(rays_detections, rays_deployments, station.col = "station")
#' table(matched$deployment_matched)
#'
#' @export

matchDeployments <- function(detections,
                             deployments,
                             datetime.col = NULL,
                             station.col = NULL,
                             lon.col = NULL,
                             lat.col = NULL,
                             deploy.col = "deploy",
                             recover.col = "recover",
                             coord.tolerance = 500,
                             fill.coords = TRUE,
                             fill.station = TRUE,
                             drop.unmatched = FALSE,
                             verbose = TRUE) {

  args <- .resolveArgs(detections, list(datetime.col = datetime.col, station.col = station.col,
                                        lon.col = lon.col, lat.col = lat.col))
  datetime.col <- args$datetime.col; station.col <- args$station.col
  lon.col <- args$lon.col; lat.col <- args$lat.col

  prev_meta <- attr(detections, "moby")
  det <- as.data.frame(detections)
  dep <- as.data.frame(deployments)

  if (!"receiver" %in% colnames(det)) stop("'detections' must contain a 'receiver' column.", call. = FALSE)
  if (!datetime.col %in% colnames(det)) stop(paste0("Datetime column ('", datetime.col, "') not found in 'detections'."), call. = FALSE)
  for (req in c("receiver", "station", deploy.col)) {
    if (!req %in% colnames(dep)) stop(paste0("'deployments' is missing required column '", req, "'. See importDeployments()."), call. = FALSE)
  }
  # standardise the deploy/recover date columns to the canonical names used below
  dep$deploy <- dep[[deploy.col]]
  dep$recover <- if (recover.col %in% colnames(dep)) dep[[recover.col]] else as.POSIXct(NA)
  dep_has_coords <- all(c("lon", "lat") %in% colnames(dep))

  det$receiver <- as.character(det$receiver)
  dep$receiver <- as.character(dep$receiver); dep$station <- as.character(dep$station)
  dep <- dep[order(dep$receiver, dep$deploy), , drop = FALSE]
  dep <- .cleanRecover(dep)

  det_time <- det[[datetime.col]]
  det_station <- if (station.col %in% colnames(det)) as.character(det[[station.col]]) else rep(NA_character_, nrow(det))

  # choose the deployment row for each detection, preferring a station-name match
  dep_idx <- rep(NA_integer_, nrow(det))
  station_matched <- rep(FALSE, nrow(det))
  for (i in seq_len(nrow(dep))) {
    in_win <- det$receiver == dep$receiver[i] & det_time >= dep$deploy[i] & det_time <= dep$recover_clean[i]
    in_win[is.na(in_win)] <- FALSE
    st_match <- in_win & !is.na(det_station) & det_station == dep$station[i]
    assign_now <- in_win & (is.na(dep_idx) | (st_match & !station_matched))
    dep_idx[assign_now] <- i
    station_matched[assign_now] <- st_match[assign_now]
  }

  matched <- !is.na(dep_idx)
  det$deployment_matched <- matched

  # coordinate mismatch flag + back-filling from the (authoritative) metadata
  det$coord_mismatch <- FALSE
  n_filled <- 0L
  if (dep_has_coords) {
    meta_lon <- dep$lon[dep_idx]; meta_lat <- dep$lat[dep_idx]
    if (!lon.col %in% colnames(det)) det[[lon.col]] <- NA_real_
    if (!lat.col %in% colnames(det)) det[[lat.col]] <- NA_real_
    own_lon <- suppressWarnings(as.numeric(det[[lon.col]]))
    own_lat <- suppressWarnings(as.numeric(det[[lat.col]]))
    both <- matched & !is.na(own_lon) & !is.na(own_lat) & !is.na(meta_lon) & !is.na(meta_lat)
    if (any(both)) {
      d <- .coordDist(own_lon[both], own_lat[both], meta_lon[both], meta_lat[both])
      det$coord_mismatch[which(both)] <- !is.na(d) & d > coord.tolerance
    }
    if (fill.coords) {
      need <- matched & (is.na(own_lon) | is.na(own_lat)) & !is.na(meta_lon) & !is.na(meta_lat)
      det[[lon.col]][need] <- meta_lon[need]
      det[[lat.col]][need] <- meta_lat[need]
      n_filled <- sum(need)
    }
  }

  # back-fill missing station names from the matched deployment
  if (fill.station) {
    meta_station <- dep$station[dep_idx]
    if (!station.col %in% colnames(det)) det[[station.col]] <- NA_character_
    need_st <- matched & is.na(det[[station.col]]) & !is.na(meta_station)
    det[[station.col]][need_st] <- meta_station[need_st]
  }

  if (verbose) {
    cat(sprintf("matchDeployments: %d/%d detections matched a deployment window (%d unmatched).\n",
                sum(matched), nrow(det), sum(!matched)))
    if (dep_has_coords) cat(sprintf("   %d coordinate(s) back-filled from metadata; %d coordinate mismatch(es) > %g.\n",
                                    n_filled, sum(det$coord_mismatch, na.rm = TRUE), coord.tolerance))
  }

  if (drop.unmatched) det <- det[matched, , drop = FALSE]
  rownames(det) <- NULL

  base_meta <- if (!is.null(prev_meta)) prev_meta else list()
  do.call(as_moby, c(list(det), base_meta))
}


#' @export
print.mobyQC <- function(x, ...) {
  n <- nrow(x$report)
  cat("<mobyQC> deployment metadata quality-control report\n")
  cat("  ", nrow(x$deployments), " deployment records | ",
      length(unique(x$deployments$receiver)), " receivers | ",
      length(unique(x$deployments$station)), " stations\n", sep = "")
  if (n == 0) {
    cat("  No issues flagged.\n")
  } else {
    cat("  ", n, " issue(s) flagged:\n", sep = "")
    ct <- sort(x$counts, decreasing = TRUE)
    for (nm in names(ct)) cat(sprintf("    - %-32s %d\n", nm, ct[[nm]]))
    cat("  Inspect $report for details.\n")
  }
  invisible(x)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
