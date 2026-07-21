#######################################################################################################
## Calculate residency indices ########################################################################
#######################################################################################################

#' Calculate residency indices
#'
#' @description Computes individual residency indices and the temporal building blocks they are
#' derived from, returning a tidy, fully numeric table (one row per animal) suitable for
#' plotting and downstream statistical analysis. This is the numeric core used internally by
#' \code{\link{summaryTable}} (which formats these values for publication); use
#' `calculateResidency()` directly when you need the raw values rather than a formatted table.
#'
#' Three indices, widely used in acoustic telemetry studies (Kraft et al. 2023; Appert et al. 2023),
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
#'
#' @inheritParams as_moby
#' @param data A data frame (or \code{\link{mobyData}}) of detections.
#' @param tagging.dates A POSIXct vector of tagging/release dates (single value or named by ID).
#' Inherited from the `mobyData` metadata when available.
#' @param tag.durations Optional numeric vector of tag battery durations (in days), used (with
#' `last.monitoring.date`) to define the study interval Dt. Required (together with, or instead of,
#' `last.monitoring.date`) when `IR2` or `IWR` are requested.
#' @param last.monitoring.date Optional POSIXct value or named vector giving the last date data
#' could be retrieved (last download / receivers operational). When both this and `tag.durations`
#' are supplied, the shorter of the two defines the monitoring end per individual.
#' @param residency.index Character vector of indices to compute; any of `"IR1"`, `"IR2"`,
#' `"IWR"`, `"IR2/IR1"`. Defaults to `c("IR1", "IR2", "IWR")`.
#' @param start.point Starting point for the detection span: `"release"` (default) or
#' `"first.detection"`.
#' @param residency.by Optional column name used to additionally compute partial (spatially
#' structured) residencies, one column per level (e.g. per habitat or array).
#' @param cap Logical; cap index values at 1 (their theoretical maximum). Defaults to TRUE.
#' Set to FALSE to retain raw values (useful for diagnosing edge effects).
#'
#' @return A data frame with one row per individual containing: the ID column, `tagging_date`,
#' `first_detection`, `last_detection`, `monitoring_end` (POSIXct); `days_detected` (Dd),
#' `detection_span` (Di) and `monitoring_duration` (Dt) in days; one numeric column per requested
#' index; and, if `residency.by` is set, additional `"<index> <level>"` partial-residency columns.
#'
#' @references
#' Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., & Abecasis, D. (2023).
#' Residency and space use estimation methods based on passive acoustic telemetry data.
#' Movement Ecology, 11(1), 12. https://doi.org/10.1186/s40462-023-00349-y
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
                               tagging.dates = NULL,
                               tag.durations = NULL,
                               id.col = NULL,
                               datetime.col = NULL,
                               last.monitoring.date = NULL,
                               residency.index = c("IR1", "IR2", "IWR"),
                               start.point = "release",
                               residency.by = NULL,
                               cap = TRUE) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  tag.durations <- reviewed_params$tag.durations
  last.monitoring.date <- reviewed_params$last.monitoring.date

  errors <- c()
  if (!all(residency.index %in% c("IR1", "IR2", "IWR", "IR2/IR1"))) {
    errors <- c(errors, "Invalid 'residency.index'. Choose from: 'IR1', 'IR2', 'IWR', 'IR2/IR1'.")
  }
  if (!is.character(start.point) || length(start.point) != 1 || !start.point %in% c("release", "first.detection")) {
    errors <- c(errors, "Invalid 'start.point'. Must be 'release' or 'first.detection'.")
  }
  if (!is.null(residency.by) && !residency.by %in% colnames(data)) {
    errors <- c(errors, "Variable used to calculate partial residencies ('residency.by') not found in the data.")
  }
  if (!is.logical(cap) || length(cap) != 1) errors <- c(errors, "'cap' must be a single logical value.")
  if (any(residency.index %in% c("IR2", "IWR")) && is.null(tag.durations) && is.null(last.monitoring.date)) {
    errors <- c(errors, "Indices 'IR2'/'IWR' require 'tag.durations' or 'last.monitoring.date' to define the study interval.")
  }
  if (length(errors) > 0) {
    stop(paste0("\n", paste0("- ", errors, collapse = "\n")), call. = FALSE)
  }

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

  # monitoring end (Dt): shortest of tag expiration / last monitoring date, else last detection
  if (!is.null(tag.durations)) {
    end_dates <- as.POSIXct(rep(NA, length(tagging.dates)), tz = tz)
    tag_expiration_dates <- tagging.dates + tag.durations * 60 * 60 * 24
    for (e in seq_along(tagging.dates)) {
      if (!is.null(last.monitoring.date)) end_dates[e] <- min(tag_expiration_dates[e], last.monitoring.date[e], na.rm = TRUE)
      else end_dates[e] <- tag_expiration_dates[e]
    }
  } else if (!is.null(last.monitoring.date)) {
    end_dates <- last.monitoring.date
  } else {
    end_dates <- rep(max(data[, datetime.col], na.rm = TRUE), length(tagging.dates))
    warning(paste(strwrap(paste("- Neither `tag.durations` nor `last.monitoring.date` were provided.",
                                "Using the last detection date in the dataset as the monitoring end for all individuals."),
                          width = getOption("width")), collapse = "\n"), call. = FALSE)
  }

  # days with detections (Dd)
  data$.date <- strftime(data[, datetime.col], format = "%d-%m-%Y", tz = tz)
  Dd <- stats::aggregate(data$.date, by = list(data[, id.col]), function(x) length(unique(x)), drop = FALSE)$x

  # detection span (Di) - days at liberty, inclusive of both endpoints
  Di <- as.integer(difftime(last_detections, start_dates, units = "days")) + 1
  Di[Di < 1] <- NA

  # study interval (Dt) - the monitoring / tag-active duration, in days. Non-positive intervals
  # (e.g. a monitoring end or tag expiration before tagging) become NA, so a bad input yields NA
  # rather than a silently negative residency index that cap = TRUE would not floor.
  Dt <- as.integer(difftime(end_dates, tagging.dates, units = "days"))
  Dt[Dt <= 0] <- NA

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
  out
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
