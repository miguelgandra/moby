#######################################################################################################
## Calculate rates of movement #######################################################################
#######################################################################################################

#' Calculate rates of movement (and total distance travelled)
#'
#' @description Summarises each animal's stepwise distances (as returned by
#' \code{\link{calculateStepDistances}}) into a tidy, fully numeric table — one row per individual —
#' suitable for plotting and downstream statistical analysis. This is the numeric core used
#' internally by \code{\link{movementTable}} (which formats these values, joins home-range areas
#' and adds summary rows for publication); use `calculateROM()` directly when you need
#' the raw values rather than a formatted table.
#'
#' Two related quantities are derived from the per-step distance column:
#' \itemize{
#'   \item \strong{Total distance} — the sum of consecutive in-water (or linear) step distances,
#'   i.e. the total path length travelled by each individual.
#'   \item \strong{Rate of movement (ROM)} — the per-step distance rescaled to an hourly rate using
#'   the (modal) time-bin interval. Both the mean and maximum ROM are reported, in metres per hour.
#' }
#' When the detection series contains irregular time-bin intervals (more than two distinct bin
#' widths), distances are first reconciled to a common interval with
#' \code{\link{interpolateDistances}} so that rates are comparable across individuals.
#'
#' @inheritParams as_moby
#' @param data A data frame (or \code{\link{mobyData}}) of binned detections with stepwise
#' distances, as returned by \code{\link{calculateStepDistances}}.
#' @param dist.col Name of the column containing the (stepwise) distance values, in metres.
#' Defaults to `"dist_m"`.
#'
#' @return A data frame with one row per individual containing: the ID column, `n_steps`
#' (number of movement steps, i.e. non-missing distances), `total_distance_m` (total path length
#' in metres), and `mean_rom` / `max_rom` (mean and maximum rate of movement in metres per hour).
#' The (modal) time-bin interval used to convert distances to rates, in minutes, is stored in the
#' `"interval"` attribute.
#'
#' @seealso \code{\link{calculateStepDistances}}, \code{\link{calculateLinearityIndex}},
#' \code{\link{movementTable}}, \code{\link{interpolateDistances}}
#' @examples
#' data(rays)
#'
#' # build per-time-bin tracks with stepwise distances
#' coas <- calculateCOAs(rays)
#' tracks <- calculateStepDistances(coas, verbose = FALSE)
#'
#' # summarise total distance travelled and rate of movement per individual
#' calculateROM(tracks)
#'
#' @export

calculateROM <- function(data,
                                     id.col = NULL,
                                     timebin.col = NULL,
                                     dist.col = "dist_m") {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  if (!dist.col %in% colnames(data)) {
    stop(paste0("Distance column ('", dist.col, "') not found in the data. ",
                "Did you run calculateStepDistances() first?"), call. = FALSE)
  }

  # reject duplicated id/time-bin combinations (ambiguous step distances)
  if (any(duplicated(data[, c(id.col, timebin.col)]))) {
    stop("Duplicated detection timestamps found. Please ensure all detections are correctly time-binned before proceeding.", call. = FALSE)
  }

  ##############################################################################
  ## Determine the time-bin interval (and interpolate if irregular) ############
  ##############################################################################

  .binInterval <- function(d) {
    by_ind <- split(d, f = d[, id.col])
    diffs <- lapply(by_ind, function(x) difftime(x[, timebin.col], .lag(x[, timebin.col]), units = "min"))
    diffs <- unlist(lapply(diffs, unique))
    unique(diffs[!is.na(diffs)])
  }

  interval <- .binInterval(data)
  # irregular series (more than one distinct bin width): regularise by interpolation so that every
  # step spans one nominal interval before distances are converted to rates of movement
  if (length(interval) > 1) {
    data <- interpolateDistances(data, id.col = id.col, timebin.col = timebin.col,
                                 dist.col = dist.col, keep.intermediate = TRUE)
  }
  # reduce to a single scalar = the MODAL (most common) step interval, in minutes, computed across
  # all consecutive steps. This is the true nominal bin width (after regularisation, the dominant
  # width); using the minimum would inflate rates for any gap-spanning step.
  step_mins <- unlist(lapply(split(data, data[, id.col]), function(x)
    as.numeric(difftime(x[, timebin.col], .lag(x[, timebin.col]), units = "min"))))
  step_mins <- step_mins[!is.na(step_mins) & step_mins > 0]
  interval <- as.numeric(names(sort(table(step_mins), decreasing = TRUE))[1])

  ##############################################################################
  ## Per-individual aggregation ################################################
  ##############################################################################

  ids <- levels(data[, id.col])
  d_by_id <- data[, dist.col]
  f_id <- data[, id.col]

  total_dist <- tapply(d_by_id, f_id, function(x) sum(x, na.rm = TRUE))[ids]
  mean_step  <- tapply(d_by_id, f_id, function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))[ids]
  max_step   <- tapply(d_by_id, f_id, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))[ids]
  n_steps    <- tapply(d_by_id, f_id, function(x) sum(!is.na(x)))[ids]

  # convert per-step distances to hourly rates of movement (m/h)
  mean_rom <- as.numeric(mean_step) * 60 / interval
  max_rom  <- as.numeric(max_step) * 60 / interval

  ##############################################################################
  ## Assemble output ###########################################################
  ##############################################################################

  out <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
                    ID = ids,
                    n_steps = as.integer(n_steps),
                    total_distance_m = as.numeric(total_dist),
                    mean_rom = mean_rom,
                    max_rom = max_rom)
  colnames(out)[1] <- id.col
  rownames(out) <- NULL
  attr(out, "interval") <- interval
  out
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
