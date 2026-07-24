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

# For each (lon, lat), the distance in metres it lies INSIDE a land polygon (0 in water). Uses a
# binary point-in-polygon test, then measures each on-land point's distance to the nearest coastline
# so the caller can apply a tolerance (a genuine near-shore receiver a coarse coastline overlaps sits
# a few metres inside; a metadata error - swapped lon/lat, a sign flip - lands kilometres inside).
# Returns NULL if the test cannot be run safely (no sf, or no usable CRS), so the audit can skip it
# without aborting. st_distance returns metres for a geographic CRS (via s2) and map units for a
# projected one, matching the metric convention used elsewhere in the package.
.distanceIntoLand <- function(lon, lat, land.shape, epsg.code = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) return(NULL)
  land <- tryCatch(sf::st_as_sf(land.shape), error = function(e) NULL)
  if (is.null(land)) return(NULL)
  land_crs <- sf::st_crs(land)
  if (is.na(land_crs)) return(NULL)                       # cannot measure metres without a land CRS
  geographic <- all(abs(lon) <= 180, na.rm = TRUE) && all(abs(lat) <= 90, na.rm = TRUE)
  pt_crs <- if (!is.null(epsg.code)) sf::st_crs(epsg.code) else if (geographic) sf::st_crs(4326) else land_crs
  out <- tryCatch({
    pts <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon", "lat"), crs = pt_crs)
    pts <- sf::st_transform(pts, land_crs)
    on_land <- lengths(sf::st_intersects(pts, land)) > 0
    d <- rep(0, length(lon))
    if (any(on_land)) {
      coast <- sf::st_union(sf::st_boundary(sf::st_geometry(land)))
      d[on_land] <- as.numeric(suppressWarnings(sf::st_distance(pts[on_land, ], coast)))
    }
    d
  }, error = function(e) NULL)
  out
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
#' @description Performs quality control on receiver deployment metadata by identifying
#' inconsistencies that can compromise acoustic telemetry analyses. The function checks the internal
#' consistency of deployment records and, optionally, cross-validates deployment metadata against a
#' detection dataset.
#'
#' The function is intentionally non-destructive: it only reports potential issues and never modifies
#' input data. Reported issues can then be reviewed and, where appropriate, resolved using
#' [matchDeployments()].
#'
#' @details Internal metadata checks include missing or invalid deployment dates, overlapping
#' deployments, duplicate records, coordinate inconsistencies, station naming issues, and gaps in
#' receiver coverage. When a land layer is supplied (`land.shape`), the function additionally flags
#' receiver positions that fall on land - a common consequence of a coordinate-entry error (swapped
#' longitude/latitude, a wrong sign, a decimal typo).
#'
#' When detections are provided, the function additionally checks whether detections are consistent
#' with the deployment history, including unknown receivers, detections outside deployment periods,
#' and station mismatches.
#'
#' Setting `min.active.days` additionally flags stations monitored for less than that total operational
#' duration. This is offered as a **flag only**: excluding low-effort receivers is a common way to
#' reduce monitoring-effort heterogeneity, but it is a deliberate, analysis-specific decision, not a
#' QC fix. Deleting short-duration stations can trade one bias for others - short-lived stations are
#' rarely random in space or time (a peripheral or late-added array, a station lost to a storm in a
#' high-use habitat), so their removal can bias space-use extent, residency and network structure, and
#' erase seasonal signal. It is usually preferable to *account for* effort (e.g. the effort-normalised
#' residency indices from \code{\link{calculateResidency}}) than to delete data. The flag surfaces the
#' candidates and their cost (duration, and detections/individuals held when `detections` is supplied)
#' so the choice can be made explicitly; `checkDeployments()` never removes them.
#'
#' @param deployments A receiver-deployment data frame, e.g. from \code{\link{importDeployments}},
#' with a `receiver` column plus station, deployment / recovery date-times and (optionally)
#' longitude / latitude. Those columns may carry non-canonical names, resolved via the
#' `deployment.*` arguments below.
#' @param detections Optional. A detection dataset (`mobyData` or data frame) with `receiver`,
#' `station` and a date-time column, used for the detection-vs-metadata checks.
#' @param id.col,datetime.col,station.col Column names in the **detection** dataset (`detections`),
#' resolved from its `mobyData` metadata or canonical defaults when `NULL`. Bare `*.col` arguments
#' always refer to the detections; the deployment log's columns use the `deployment.*` arguments.
#' @template deploymentSpatialArgs
#' @template deploymentDateArgs
#' @param checks Character vector selecting which check groups to run (any of, or `"all"`, the
#' default): `"dates"` (deploy/recover date integrity), `"overlaps"` (duplicate and overlapping
#' deployment records), `"gaps"` (coverage gaps between consecutive deployments - often benign
#' servicing, so the easiest to drop), `"coordinates"` (implausible coordinates,
#' station-vs-coordinate naming consistency, and - when `land.shape` is supplied - receiver positions
#' on land; skipped if `lon`/`lat` are absent) and `"detections"` (cross-check detections against
#' deployment windows; requires the `detections` argument).
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
#' @param land.shape Optional `sf` (or `SpatialPolygons*`) polygon layer of landmasses. When
#' supplied, the `"coordinates"` group additionally flags receiver positions that fall on land.
#' Off by default; if `NULL`, it is taken from the `detections` `mobyData` metadata when present (i.e.
#' from `as_moby(detections, land.shape = ...)`). The layer must carry a coordinate reference system;
#' if it lacks one, or cannot be reconciled with the deployment coordinates, the on-land check is
#' skipped with a message (the audit never aborts).
#' @param epsg.code Optional EPSG code for the deployment coordinates, used only by the on-land check.
#' Needed only when the coordinates are projected (not longitude/latitude); geographic coordinates are
#' assumed to be WGS84.
#' @param land.tolerance Numeric. How far inside the coastline (in metres) a receiver position must
#' lie to be reported as on land. This spares genuine near-shore receivers that a coarse coastline
#' overlaps by a small margin; gross coordinate-entry errors (swapped lon/lat, a sign flip) fall
#' kilometres inland and are still flagged. Defaults to 500.
#' @param min.active.days Optional numeric. When set, stations whose **total operational time** - the
#' sum of their deployment windows in days, so servicing gaps do not inflate it - falls below this
#' value are flagged as `"Short monitoring duration"`. Off by default (`NULL`); there is no universal
#' threshold (it is study-specific). The flag is **report-only**: it never removes anything. See the
#' Details for why deletion is a deliberate, bias-prone choice best left to the analyst.
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
#' # additionally flag (never remove) stations monitored for under ~6 months
#' checkDeployments(rays_deployments, detections = rays, min.active.days = 30 * 6)
#'
#' @export

checkDeployments <- function(deployments,
                             detections = NULL,
                             id.col = NULL,
                             datetime.col = NULL,
                             station.col = NULL,
                             deployment.station.col = "station",
                             deployment.lon.col = "lon",
                             deployment.lat.col = "lat",
                             deployment.deploy.col = "deploy",
                             deployment.recover.col = "recover",
                             land.shape = NULL,
                             epsg.code = NULL,
                             checks = "all",
                             scope = c("all", "detected"),
                             coord.tolerance = 500,
                             gap.tolerance = 1,
                             land.tolerance = 500,
                             min.active.days = NULL,
                             verbose = TRUE) {

  scope <- match.arg(scope)
  if (!is.null(min.active.days) && (!is.numeric(min.active.days) || length(min.active.days) != 1 ||
                                    is.na(min.active.days) || min.active.days <= 0))
    stop("'min.active.days' must be a single positive number (or NULL).", call. = FALSE)
  # the on-land check is opt-in: use an explicit land.shape, otherwise borrow one the user already
  # declared on the detections mobyData (as_moby(det, land.shape = ...)). Never required.
  if (is.null(land.shape) && !is.null(detections)) land.shape <- attr(detections, "moby")$land.shape
  dep <- as.data.frame(deployments)
  required <- c("receiver", deployment.station.col, deployment.deploy.col)
  miss <- setdiff(required, colnames(dep))
  if (length(miss) > 0) {
    stop(paste0("'deployments' is missing required column(s): ", paste(miss, collapse = ", "),
                ". See importDeployments()."), call. = FALSE)
  }
  # standardise the station/deploy/recover/coordinate columns to the canonical names used by the checks
  # below, so a log whose columns are named differently need not be renamed by the user first
  dep$station <- dep[[deployment.station.col]]
  dep$deploy <- dep[[deployment.deploy.col]]
  dep$recover <- if (deployment.recover.col %in% colnames(dep)) dep[[deployment.recover.col]] else as.POSIXct(NA)
  has_coords <- all(c(deployment.lon.col, deployment.lat.col) %in% colnames(dep))
  # an explicitly named coordinate column that is absent is almost always a typo, not "no coordinates"
  if (!has_coords && (!identical(deployment.lon.col, "lon") || !identical(deployment.lat.col, "lat")) &&
      any(c(deployment.lon.col, deployment.lat.col) %in% colnames(dep)))
    message("- coordinate column(s) '",
            paste(setdiff(c(deployment.lon.col, deployment.lat.col), colnames(dep)), collapse = "', '"),
            "' not found in 'deployments'; skipping the coordinate checks.")
  if (has_coords) { dep$lon <- dep[[deployment.lon.col]]; dep$lat <- dep[[deployment.lat.col]] }

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
  station_impact <- NULL   # per-station detection tally, filled by the detection block for the short-duration flag

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

    # receiver positions that fall on land (opt-in: only when a land layer is available). Dedup by
    # (station, coordinate) so a receiver moved between an on-land and an off-land position under one
    # name is still caught. A distance-into-land tolerance spares genuine near-shore receivers.
    if (!is.null(land.shape)) {
      ok <- !is.na(lon) & !is.na(lat)
      pos <- dep_scope[ok, c("receiver", "station"), drop = FALSE]
      pos$lon <- lon[ok]; pos$lat <- lat[ok]
      pos <- pos[!duplicated(pos[c("station", "lon", "lat")]), , drop = FALSE]
      din <- if (nrow(pos) > 0) .distanceIntoLand(pos$lon, pos$lat, land.shape, epsg.code) else numeric(0)
      if (is.null(din)) {
        message("- could not run the on-land check: 'land.shape' lacks a usable CRS or is ",
                "incompatible with the coordinates; skipping it.")
      } else {
        for (i in which(din > land.tolerance)) {
          likely_error <- din[i] > 5000
          add(.qcRow("Coordinates on land", pos$receiver[i], pos$station[i],
                     details = paste0("Station coordinates fall ", .formatDistance(din[i]), " inside land",
                       if (likely_error)
                         "; likely a metadata error (e.g. swapped lon/lat or a sign error). Verify against the deployment records."
                       else
                         paste0(" (beyond the ", .formatDistance(land.tolerance),
                                " tolerance). Verify against the deployment records, or use correctPositions() for a genuine near-shore point."))))
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

    # per-station detection tally (so the short-duration flag can show what removing a station costs)
    if (any(!is.na(det_station)))
      station_impact <- list(n_det = tapply(seq_along(det_station), trimws(det_station), length),
                             n_ind = tapply(det_id, trimws(det_station), function(x) length(unique(x))))

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
  ## Short monitoring duration (opt-in) ########################################
  ##############################################################################

  # Flag stations whose total OPERATIONAL time - the sum of their deployment windows, so servicing
  # gaps do not inflate it - falls below min.active.days. Report-only by design: short-effort stations
  # are removed in some workflows to reduce effort heterogeneity, but deletion trades that bias for
  # others (lost spatial/seasonal coverage, selection effects, network artefacts), so the choice must
  # be the user's, made deliberately - not a side effect of an audit. The row reports the actual
  # duration and, when detections are supplied, how many detections/individuals the station holds.
  if (!is.null(min.active.days)) {
    active_days <- as.numeric(difftime(dep_scope$recover_clean, dep_scope$deploy, units = "days"))
    dur_by_station <- tapply(active_days, dep_scope$station, sum, na.rm = TRUE)
    short <- names(dur_by_station)[!is.na(dur_by_station) & dur_by_station < min.active.days]
    for (st in short) {
      rows <- dep_scope[dep_scope$station == st, , drop = FALSE]
      nd <- if (!is.null(station_impact) && st %in% names(station_impact$n_det)) as.integer(station_impact$n_det[[st]]) else NA_integer_
      ni <- if (!is.null(station_impact) && st %in% names(station_impact$n_ind)) as.integer(station_impact$n_ind[[st]]) else NA_integer_
      add(.qcRow("Short monitoring duration", paste(unique(rows$receiver), collapse = " | "), st,
                 first = min(rows$deploy, na.rm = TRUE), last = max(rows$recover_clean, na.rm = TRUE),
                 n_detections = nd, n_individuals = ni,
                 details = sprintf(paste0("Active for %.0f day(s) across %d deployment(s), below the %s-day ",
                                          "threshold. Consider whether removing it is appropriate for your ",
                                          "analysis (it can bias space use, residency and networks; effort-aware ",
                                          "analysis is often preferable to deletion)."),
                                   dur_by_station[[st]], nrow(rows), format(min.active.days))))
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
#' (`receiver` plus station, coordinates and deploy/recover date-times; column names resolved via the
#' `deployment.*` arguments).
#' @param datetime.col,station.col,lon.col,lat.col Column names in the **detection** dataset
#' (`detections`), resolved from its `mobyData` metadata or canonical defaults when `NULL`. Bare
#' `*.col` arguments always refer to the detections; the deployment log's columns use the
#' `deployment.*` arguments.
#' @template deploymentSpatialArgs
#' @template deploymentDateArgs
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
                             deployment.station.col = "station",
                             deployment.lon.col = "lon",
                             deployment.lat.col = "lat",
                             deployment.deploy.col = "deploy",
                             deployment.recover.col = "recover",
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
  for (req in c("receiver", deployment.station.col, deployment.deploy.col)) {
    if (!req %in% colnames(dep)) stop(paste0("'deployments' is missing required column '", req, "'. See importDeployments()."), call. = FALSE)
  }
  # standardise the station/deploy/recover/coordinate columns to the canonical names used below
  dep$station <- dep[[deployment.station.col]]
  dep$deploy <- dep[[deployment.deploy.col]]
  dep$recover <- if (deployment.recover.col %in% colnames(dep)) dep[[deployment.recover.col]] else as.POSIXct(NA)
  dep_has_coords <- all(c(deployment.lon.col, deployment.lat.col) %in% colnames(dep))
  if (dep_has_coords) { dep$lon <- dep[[deployment.lon.col]]; dep$lat <- dep[[deployment.lat.col]] }

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
