##################################################################################################
## Function to interpolate (average) distance across all timebins ################################
##################################################################################################

#' Distance interpolation

#' @description Function to interpolate (average) distances across all timebins.
#' If an animal goes undetected for large periods of time, large distance spikes
#' may occur after running the \code{\link{calculateDistances}} function.
#' Hence, to avoid large biases when averaging distances per hour or day, this
#' function 'dilutes' traveled distances by the total number of time-bins
#' passed between each two consecutive positions.
#'
#' @param data Output \code{\link{calculateDistances}} function: a data frame
#' with animal IDS, positions and distances.
#' @param id.col Name of the column containing animal IDs Defaults to 'ID'.
#' @param dist.col Name of the column containing distances. Defaults to 'dist_m',
#' the output of \code{\link{calculateDistances}}.
#' @param keep.intermediate Boolean indicating if intermediate distances (assigned
#' to time-bins without detections) should be kept or discarded. Defaults to false.
#' @return Original data frame plus missing time-bins (with diluted distances)
#' @export


interpolateDistances <- function(data, id.col="ID", dist.col="dist_m", keep.intermediate=F){

  # print message
  cat("Interpolating distances\n")

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # split data by individual
  data$origin_temp <- 1
  data_fish <- split(data, f=data[,id.col])

  # set progress bar
  pb <- txtProgressBar(min=0, max=length(data_fish), initial=0, style=3)

  # loop through each individual
  data_list <- list()
  for(i in 1:length(data_fish)) {

    setTxtProgressBar(pb,i)
    data_out <- data_fish[[i]]
    if(nrow(data_out)<=1) {data_list[[i]] <- data_out; next}

    # add missing (in-between) timebins
    interval <- difftime(data_out$timebin, dplyr::lag(data_out$timebin), units="min")
    interval <- as.numeric(min(interval[interval>0], na.rm=T))
    complete_seq <- seq.POSIXt(from=min(data_out$timebin), to=max(data_out$timebin), by=interval*60)
    complete_seq <- data.frame("timebin"=complete_seq, "ID"=unique(data_out[,id.col]))
    colnames(complete_seq)[2] <- id.col
    data_out <- merge(data_out, complete_seq, by=c("timebin",id.col), all=T)
    data_out <- data_out[order(data_out$timebin),]

    # interpolate distances
    missing <- which(is.na(data_out[,dist.col][-length(data_out[,dist.col])]))
    to_replace <- split(missing, cumsum(c(1, diff(missing) != 1)))
    to_replace <- lapply(to_replace, function(x) return(c(x[1]-1,x)))
    replace_vals <- lapply(to_replace, function(x) data_out[,dist.col][x])
    replace_vals <- lapply(replace_vals, function(x) rep(x[1]/length(x), length(x)))
    data_out[,dist.col][unlist(to_replace)] <- unlist(replace_vals)
    data_out[,dist.col][nrow(data_out)] <- NA
    if(keep.intermediate==F) {data_out[,dist.col][is.na(data_out$origin_temp)]<-NA}
    data_list[[i]] <- data_out
  }
  close(pb)

  # return results
  result <- do.call("rbind", data_list)
  result <- result[,-which(colnames(result)=="origin_temp")]
  if("detections" %in% colnames(data)) {result$detections[is.na(result$detections)] <- 0}
  if("stations" %in% colnames(data)) {result$stations[is.na(result$stations)] <- 0}
  if("transmitter" %in% colnames(data)) {result$transmitter <- zoo::na.locf(result$transmitter)}


  return(result)
}
##################################################################################################
##################################################################################################
##################################################################################################
