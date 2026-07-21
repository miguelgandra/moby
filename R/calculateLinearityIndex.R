#######################################################################################################
## Calculate movement linearity index ################################################################
#######################################################################################################

#' Calculate the linearity index of individual movements
#'
#' @description Computes a movement **linearity (directness) index** for each animal, returning a
#' tidy, fully numeric table (one row per individual). The index is the ratio between the net
#' displacement (the shortest in-water path between an individual's first and last positions) and
#' the total path length actually travelled (the sum of consecutive step distances):
#' \deqn{LI = net\ displacement / total\ distance.}
#' Values approach 1 for highly directional (near straight-line) movements and approach 0 for
#' convoluted, back-and-forth or resident movements.
#'
#' Net displacements are obtained with \code{\link{calculateStepDistances}}, so when a `land.shape` is
#' supplied they follow the shortest in-water route (consistent with the total distance); otherwise
#' great-circle distances are used. This is the numeric core behind the `LI` column of
#' \code{\link{movementTable}}.
#'
#' @inheritParams as_moby
#' @param data A data frame (or \code{\link{mobyData}}) of binned detections with stepwise
#' distances, as returned by \code{\link{calculateStepDistances}}.
#' @param land.shape Optional. A projected coastline/landmass layer (`sf` or convertible). When
#' supplied, net displacements follow the shortest in-water path; otherwise linear (great-circle)
#' distances are used.
#' @param epsg.code Optional. Coordinate reference system used to project positions. If not
#' supplied, the CRS is taken from `land.shape` (when available).
#' @param dist.col Name of the column containing the (stepwise) distance values, in metres.
#' Defaults to `"dist_m"`.
#' @param ... Additional arguments passed to \code{\link{calculateStepDistances}} (e.g. `grid.resolution`,
#' `mov.directions`, `cores`), used when computing net displacements.
#'
#' @return A data frame with one row per individual containing: the ID column, `net_distance_m`
#' (shortest distance between the first and last positions, in metres), `total_distance_m` (total
#' path length, in metres) and `linearity_index` (their ratio, between 0 and 1).
#'
#' @seealso \code{\link{calculateStepDistances}}, \code{\link{calculateROM}},
#' \code{\link{movementTable}}
#' @examples
#' data(rays)
#'
#' # build per-time-bin tracks with stepwise distances
#' coas <- calculateCOAs(rays)
#' tracks <- calculateStepDistances(coas, verbose = FALSE)
#'
#' # movement directness per individual (net displacement / total path length)
#' calculateLinearityIndex(tracks)
#'
#' @export

calculateLinearityIndex <- function(data,
                                    land.shape = NULL,
                                    epsg.code = NULL,
                                    id.col = NULL,
                                    timebin.col = NULL,
                                    lon.col = NULL,
                                    lat.col = NULL,
                                    dist.col = "dist_m",
                                    ...) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  land.shape <- reviewed_params$land.shape

  if (!dist.col %in% colnames(data)) {
    stop(paste0("Distance column ('", dist.col, "') not found in the data. ",
                "Did you run calculateStepDistances() first?"), call. = FALSE)
  }

  ##############################################################################
  ## Total path length (Di denominator) ########################################
  ##############################################################################

  ids <- levels(data[, id.col])
  total_dist <- tapply(data[, dist.col], data[, id.col], function(x) sum(x, na.rm = TRUE))[ids]
  total_named <- stats::setNames(as.numeric(total_dist), ids)

  ##############################################################################
  ## Net displacement (first -> last position) #################################
  ##############################################################################

  # first and last centre-of-activity per individual (in time order)
  first_coas <- by(data, data[, id.col], function(x) x[which.min(x[, timebin.col]), c(lon.col, lat.col)], simplify = FALSE)
  first_coas <- do.call("rbind", first_coas)
  first_coas[[id.col]] <- rownames(first_coas)
  first_coas$type <- "start"
  last_coas <- by(data, data[, id.col], function(x) x[which.max(x[, timebin.col]), c(lon.col, lat.col)], simplify = FALSE)
  last_coas <- do.call("rbind", last_coas)
  last_coas[[id.col]] <- rownames(last_coas)
  last_coas$type <- "end"

  start_end_coas <- rbind(first_coas, last_coas)
  rownames(start_end_coas) <- NULL
  start_end_coas$type <- factor(start_end_coas$type, levels = c("start", "end"))
  start_end_coas <- start_end_coas[order(start_end_coas[[id.col]], start_end_coas$type), ]

  # shortest distance between the start and end positions of each individual
  start_end_dists <- suppressWarnings(calculateStepDistances(start_end_coas, land.shape = land.shape, id.col = id.col,
                                                      lon.col = lon.col, lat.col = lat.col,
                                                      epsg.code = epsg.code, verbose = FALSE, ...))
  net_rows <- start_end_dists[!is.na(start_end_dists$dist_m), , drop = FALSE]
  net_named <- stats::setNames(as.numeric(net_rows$dist_m), as.character(net_rows[[id.col]]))

  ##############################################################################
  ## Assemble output ###########################################################
  ##############################################################################

  net_dist <- net_named[ids]
  li <- as.numeric(net_dist) / total_named[ids]
  li[!is.finite(li)] <- NA_real_

  out <- data.frame(check.names = FALSE, stringsAsFactors = FALSE,
                    ID = ids,
                    net_distance_m = as.numeric(net_dist),
                    total_distance_m = as.numeric(total_named[ids]),
                    linearity_index = as.numeric(li))
  colnames(out)[1] <- id.col
  rownames(out) <- NULL
  out
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
