##################################################################################################
## Function to interpolate (average) distance across all timebins ################################
##################################################################################################

#' Distance interpolation

#' @description Function to interpolate (average) distances across all timebins.
#' If an animal goes undetected for large periods of time, large distance spikes
#' may occur after running the \code{\link{calculateStepDistances}} function.
#' Hence, to avoid large biases when averaging distances per hour or day, this
#' function 'dilutes' traveled distances by the total number of time-bins
#' passed between each two consecutive positions.
#'
#' @inheritParams as_moby
#' @param data Output \code{\link{calculateStepDistances}} function: a data frame
#' with animal IDS, positions and distances.
#' @param dist.col Name of the column containing distances. Defaults to 'dist_m',
#' the output of \code{\link{calculateStepDistances}}.
#' @param keep.intermediate Boolean indicating if intermediate distances (assigned
#' to time-bins without detections) should be kept or discarded. Defaults to false.
#' @return Original data frame plus missing time-bins (with diluted distances)
#' @examples
#' data(rays)
#'
#' # build per-time-bin tracks with stepwise distances
#' coas <- calculateCOAs(rays)
#' tracks <- calculateStepDistances(coas, verbose = FALSE)
#'
#' # dilute distances across time-bins with no detections
#' interp <- interpolateDistances(tracks)
#' head(interp)
#'
#' @export


interpolateDistances <- function(data,
                                 id.col = NULL,
                                 timebin.col = NULL,
                                 dist.col = "dist_m",
                                 keep.intermediate = FALSE){

  # capture mobyData metadata so a detection-like output can be re-wrapped (chainable pipeline)
  prev_meta <- attr(data, "moby")

  # perform argument checks and return reviewed parameters (resolves timebin.col from metadata)
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # print message
  .printConsole("Interpolating distances")

  # split data by individual
  data$origin_temp <- 1
  data_fish <- split(data, f=data[,id.col])

  # set progress bar
  pb <- txtProgressBar(min=0, max=length(data_fish), initial=0, style=3)

  # loop through each individual
  data_list <- list()
  for(i in seq_along(data_fish)) {

    setTxtProgressBar(pb,i)
    data_out <- data_fish[[i]]
    if(nrow(data_out)<=1) {data_list[[i]] <- data_out; next}

    # add missing (in-between) timebins
    interval <- difftime(data_out[[timebin.col]], .lag(data_out[[timebin.col]]), units="min")
    interval <- as.numeric(min(interval[interval>0], na.rm=TRUE))
    complete_seq <- seq.POSIXt(from=min(data_out[[timebin.col]]), to=max(data_out[[timebin.col]]), by=interval*60)
    complete_seq <- data.frame(complete_seq, unique(data_out[,id.col]), stringsAsFactors=FALSE)
    colnames(complete_seq) <- c(timebin.col, id.col)
    data_out <- merge(data_out, complete_seq, by=c(timebin.col, id.col), all=TRUE)
    data_out <- data_out[order(data_out[[timebin.col]]),]

    # interpolate distances
    missing <- which(is.na(data_out[,dist.col][-length(data_out[,dist.col])]))
    to_replace <- split(missing, cumsum(c(1, diff(missing) != 1)))
    to_replace <- lapply(to_replace, function(x) return(c(x[1]-1,x)))
    replace_vals <- lapply(to_replace, function(x) data_out[,dist.col][x])
    replace_vals <- lapply(replace_vals, function(x) rep(x[1]/length(x), length(x)))
    data_out[,dist.col][unlist(to_replace)] <- unlist(replace_vals)
    data_out[,dist.col][nrow(data_out)] <- NA
    if(!keep.intermediate) {data_out[,dist.col][is.na(data_out$origin_temp)]<-NA}
    data_list[[i]] <- data_out
  }
  close(pb)

  # return results
  result <- do.call("rbind", data_list)
  result <- .dropCols(result, "origin_temp")
  if("detections" %in% colnames(data)) {result$detections[is.na(result$detections)] <- 0}
  if("stations" %in% colnames(data)) {result$stations[is.na(result$stations)] <- 0}
  if("transmitter" %in% colnames(data)) {result$transmitter <- .naLocf(result$transmitter)}

  # re-attach mobyData metadata so the (detection-like) output stays chainable
  if(!is.null(prev_meta)){
    attr(result, "moby") <- prev_meta
    class(result) <- unique(c("mobyData", "data.frame"))
  }

  return(result)
}
##################################################################################################
##################################################################################################
##################################################################################################
