#######################################################################################################
## Calculate residency indices ########################################################################
#######################################################################################################

#' Calculate residency indices
#'
#' @description Computes individual residency indices and the temporal building blocks they are
#' derived from, returning a typed table (one row per animal) suitable for
#' plotting and downstream statistical analysis. This is the numeric core used internally by
#' \code{\link{summaryTable}} (which formats these values for publication); use
#' `calculateResidency()` directly when you need the raw values rather than a formatted table.
#'
#' The following indices, used in acoustic telemetry studies (Kraft et al. 2023; Appert et al. 2023),
#' are available:
#' \itemize{
#'   \item \strong{IR1} = Dd / Di, the proportion of days detected over the detection span
#'   (days at liberty); a maximum residency, considering only the period the animal was known
#'   to be alive and the tag operational.
#'   \item \strong{IR2} = Dd / Dt, the proportion of days detected over the full study interval
#'   (release to last data download or tag expiration); a minimum residency.
#'   \item \strong{IWR} = (Dd / Dt) * (Di / Dt), the IR2 index weighted by the ratio of the
#'   detection interval to the study interval.
#'   \item \strong{IR2/IR1}, the ratio of IR2 to IR1 (gap between last detection and end of monitoring).
#' }
#' where Dd = number of days with detections, Di = detection span (days at liberty,
#' first/release to last detection, inclusive), and Dt = study interval (release to monitoring end).
#' Dd and Di use calendar dates in the detection timestamps' time zone. Dt counts calendar
#' dates with positive monitoring time between tagging and the earlier of tag expiration and
#' the end of receiver monitoring. A cutoff exactly at midnight does not count that new day.
#' Detections before tagging or at/after a supplied monitoring cutoff cause an error; review
#' the dates or filter those records before calling this function.
#'
#' @inheritParams as_moby
#' @param data A data frame (or \code{\link{mobyData}}) of detections.
#' @param tagging.dates A POSIXct vector of tagging/release dates (single value or named by ID).
#' Inherited from the `mobyData` metadata when available.
#' @param tag.durations Optional numeric vector of tag battery durations (in days), used (with
#' `last.monitoring.date`) to define the study interval Dt. Required (together with, or instead of,
#' `last.monitoring.date`) when `IR2`, `IWR`, or `IR2/IR1` are requested.
#' @param last.monitoring.date Optional POSIXct value or named vector giving the end of
#' receiver monitoring. The value is a cutoff timestamp, not an automatically counted final
#' day (to include all of December 31, use January 1 at midnight). When both this and
#' `tag.durations` are supplied, the earlier cutoff is used per individual.
#' @param residency.index Character vector of indices to compute; any of `"IR1"`, `"IR2"`,
#' `"IWR"`, `"IR2/IR1"`. Defaults to `c("IR1", "IR2", "IWR")`.
#' @param start.point Starting point for the detection span: `"release"` (default) or
#' `"first.detection"`.
#' @param residency.by Optional column name used to additionally compute partial (spatially
#' structured) residencies, one column per level (e.g. per habitat or array).
#' @param cap Logical; cap index values at 1 (their theoretical maximum). Defaults to TRUE.
#' Set to FALSE to retain raw values (useful for diagnosing edge effects).
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A data frame with one row per individual containing: the ID column, `tagging_date`,
#' `first_detection`, `last_detection`, `monitoring_end` (POSIXct); `days_detected` (Dd),
#' `detection_span` (Di) and `monitoring_duration` (Dt) in days; one numeric column per requested
#' index; and, if `residency.by` is set, additional `"<index> <level>"` partial-residency columns.
#' If neither monitoring cutoff is provided, `monitoring_end` and `monitoring_duration` are `NA`;
#' IR1 remains calculable, while indices using Dt require a known cutoff.
#'
#' @references
#' Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., & Abecasis, D. (2023).
#' Residency and space use estimation methods based on passive acoustic telemetry data.
#' Movement Ecology, 11(1), 12. https://doi.org/10.1186/s40462-022-00364-z
#'
#' Appert, C., Udyawer, V., Simpfendorfer, C. A., et al. (2023). Use, misuse, and ambiguity of
#' indices of residence in acoustic telemetry studies. Marine Ecology Progress Series, 714, 27-44.
#'
#' @seealso \code{\link{summaryTable}}
#' @examples
#' data(rays)
#' # residency indices, using a fixed last-monitoring date to define the study interval (Dt)
#' res <- calculateResidency(rays,
#'          last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"))
#' head(res)
#' @export

calculateResidency <- function(data,
                               id.col = NULL,
                               datetime.col = NULL,
                               tagging.dates = NULL,
                               tag.durations = NULL,
                               last.monitoring.date = NULL,
                               residency.index = c("IR1", "IR2", "IWR"),
                               start.point = "release",
                               residency.by = NULL,
                               cap = TRUE,
                               verbose = getOption("moby.verbose", TRUE)) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  tag.durations <- reviewed_params$tag.durations
  last.monitoring.date <- reviewed_params$last.monitoring.date

  errors <- c()
  if (length(residency.index) == 0 || !all(residency.index %in% c("IR1", "IR2", "IWR", "IR2/IR1"))) {
    errors <- c(errors, "Invalid 'residency.index'. Choose from: 'IR1', 'IR2', 'IWR', 'IR2/IR1'.")
  }
  if (!is.character(start.point) || length(start.point) != 1 || !start.point %in% c("release", "first.detection")) {
    errors <- c(errors, "Invalid 'start.point'. Must be 'release' or 'first.detection'.")
  }
  if (!is.null(residency.by) && !residency.by %in% colnames(data)) {
    errors <- c(errors, "Variable used to calculate partial residencies ('residency.by') not found in the data.")
  }
  if (!is.logical(cap) || length(cap) != 1 || is.na(cap)) errors <- c(errors, "'cap' must be a single logical value.")
  if (any(residency.index %in% c("IR2", "IWR", "IR2/IR1")) && is.null(tag.durations) && is.null(last.monitoring.date)) {
    errors <- c(errors, "Indices using the monitoring duration require 'tag.durations' or 'last.monitoring.date'.")
  }
  if (length(errors) > 0) {
    stop(paste0("\n", paste0("- ", errors, collapse = "\n")), call. = FALSE)
  }

  # ---- header -----------------------------------------------------------------------------------
  # Criteria are the choices that change what the numbers MEAN: which indices are computed, where the
  # detection span starts (release vs first detection moves Di and indices using it), whether
  # values are capped at their theoretical maximum, and the variable partial residencies are broken
  # down by. Everything else this function produces is read straight off the returned table.
  crit <- c(indices = paste(residency.index, collapse = paste0(" ", .mobyGlyph("mid"), " ")))
  crit["span start"] <- if (start.point == "first.detection") "first detection" else "release date"
  crit["cap"] <- if (isTRUE(cap)) "values capped at 1" else "off (raw values retained)"
  if (!is.null(residency.by)) crit["partial residency"] <- paste0("by '", residency.by, "'")

  .mobyHeader("calculateResidency()", "Computing residency indices per individual",
              input = paste0(.fmtCount(.nDetections(data), "detection"), " ", .mobyGlyph("mid"), " ",
                             .fmtCount(.nObserved(data[, id.col]), "individual")),
              criteria = crit, verbose = verbose)

  ##############################################################################
  ## Temporal building blocks ##################################################
  ##############################################################################

  tz <- .dataTZ(data[, datetime.col])
  ids <- levels(data[, id.col])

  # last detection per individual
  last_detections <- tapply(X = data[, datetime.col], INDEX = data[, id.col], FUN = max, na.rm = TRUE)
  last_detections <- as.POSIXct(last_detections[ids], origin = "1970-01-01", tz = tz)

  # first detection per individual
  first_detections <- tapply(X = data[, datetime.col], INDEX = data[, id.col], FUN = min, na.rm = TRUE)
  first_detections <- as.POSIXct(first_detections[ids], origin = "1970-01-01", tz = tz)

  # start of the detection span
  start_dates <- if (start.point == "first.detection") first_detections else tagging.dates

  # Monitoring end: shortest of tag expiration / last monitoring date, if known.
  if (!is.null(tag.durations)) {
    end_dates <- as.POSIXct(rep(NA, length(tagging.dates)), tz = tz)
    tag_expiration_dates <- tagging.dates + tag.durations * 60 * 60 * 24
    for (e in seq_along(tagging.dates)) {
      if (!is.null(last.monitoring.date)) {
        bounds <- c(tag_expiration_dates[e], last.monitoring.date[e])
        if (!all(is.na(bounds))) end_dates[e] <- min(bounds, na.rm = TRUE)
      }
      else end_dates[e] <- tag_expiration_dates[e]
    }
  } else if (!is.null(last.monitoring.date)) {
    end_dates <- last.monitoring.date
  } else {
    end_dates <- as.POSIXct(rep(NA, length(tagging.dates)), tz = tz)
  }

  # The daily index counts calendar dates in the detection data's time zone. An animal's
  # detection dates must fall within the interval when its tag and the array could be active.
  row_id <- match(as.character(data[, id.col]), ids)
  before_tagging <- !is.na(data[, datetime.col]) & !is.na(tagging.dates[row_id]) &
    data[, datetime.col] < tagging.dates[row_id]
  after_monitoring <- !is.na(data[, datetime.col]) & !is.na(end_dates[row_id]) &
    data[, datetime.col] >= end_dates[row_id]
  invalid_ids <- unique(as.character(data[before_tagging | after_monitoring, id.col]))
  if (length(invalid_ids)) {
    stop("Detections outside the tagging-to-monitoring interval for ID(s): ",
         paste(invalid_ids, collapse = ", "), ". Check tagging dates, tag durations, ",
         "and the last monitoring date.", call. = FALSE)
  }
  if (any(residency.index %in% c("IR2", "IWR", "IR2/IR1")) && anyNA(end_dates)) {
    stop("No monitoring endpoint for ID(s): ", paste(ids[is.na(end_dates)], collapse = ", "),
         ". Supply 'tag.durations' or 'last.monitoring.date'.", call. = FALSE)
  }

  # Days with detections (Dd), including each date at most once per individual.
  data$.date <- as.Date(data[, datetime.col], tz = tz)
  Dd <- stats::aggregate(data$.date, by = list(data[, id.col]),
                         function(x) length(unique(x[!is.na(x)])), drop = FALSE)$x

  # Detection span (Di): inclusive calendar dates, matching Dd's unit.
  Di <- as.integer(as.Date(last_detections, tz = tz) - as.Date(start_dates, tz = tz)) + 1L
  Di[!is.na(Di) & Di < 1] <- NA_integer_

  # Monitoring days (Dt): dates with positive monitoring time in [tagging, end).
  # If the tag/array stops exactly at midnight, that new date has no monitoring time.
  end_day <- as.Date(end_dates, tz = tz)
  end_clock <- as.POSIXlt(end_dates, tz = tz)
  Dt <- as.integer(end_day - as.Date(tagging.dates, tz = tz)) +
    as.integer(end_clock$hour > 0 | end_clock$min > 0 | end_clock$sec > 0)
  Dt[is.na(end_dates) | is.na(tagging.dates) |
       (!is.na(end_dates) & !is.na(tagging.dates) & end_dates <= tagging.dates)] <- NA_integer_

  # index calculator
  residencyValue <- function(Dd, Di, Dt, metric) {
    val <- switch(metric,
                  "IR1" = Dd / Di,
                  "IR2" = Dd / Dt,
                  "IWR" = (Dd / Dt) * (Di / Dt),
                  "IR2/IR1" = (Dd / Dt) / (Dd / Di))
    if (cap) val <- pmin(val, 1)
    val
  }

  ##############################################################################
  ## Assemble output ###########################################################
  ##############################################################################

  out <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
                    ID = ids,
                    tagging_date = tagging.dates,
                    first_detection = first_detections,
                    last_detection = last_detections,
                    monitoring_end = end_dates,
                    days_detected = Dd,
                    detection_span = Di,
                    monitoring_duration = Dt)
  colnames(out)[1] <- id.col

  for (index in residency.index) out[[index]] <- residencyValue(Dd, Di, Dt, index)

  # partial residencies (per level of residency.by)
  if (!is.null(residency.by)) {
    data_grouped <- split(data, f = data[, residency.by])
    for (index in residency.index) {
      for (i in seq_along(data_grouped)) {
        d_sub <- data_grouped[[i]]
        Dd_partial <- stats::aggregate(d_sub$.date, by = list(factor(d_sub[, id.col], levels = ids)),
                                       function(x) length(unique(x)), drop = FALSE)$x
        out[[paste(index, names(data_grouped)[i])]] <- residencyValue(Dd_partial, Di, Dt, index)
      }
    }
  }

  rownames(out) <- NULL

  # ---- the one thing the returned table cannot show ---------------------------------------------
  # An individual can end up with no computable residency at all (never detected, or an interval that
  # could not be resolved). The NAs are in the table, but only a row-by-row scan reveals them, so the
  # count is named once here. No completion line: the table itself is the summary.
  n_missing <- sum(rowSums(!is.na(out[, residency.index, drop = FALSE])) == 0)
  if (n_missing > 0) {
    .mobyBlank(verbose)
    .mobyNote(.fmtCount(n_missing, "individual"), " with no computable residency ",
              "(no detections, or an unresolvable interval)", verbose = verbose)
  }

  out
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
