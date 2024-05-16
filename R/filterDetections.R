#######################################################################################################
## Filter detections data #############################################################################
#######################################################################################################

#' Filter data
#'
#' @description Function to remove spurious detections, based on tagging date and/or
#' isolated occurrences.
#'
#' @param data A data frame containing raw animal detections.
#' @param tagging.dates A POSIXct vector containing the tag/release date of each animal.
#' @param id.col Name of the column containing animal identifications. Defaults to "ID".
#' @param datetime.col Name of the column containing datetimes (in POSIXct format). Defaults to "datetime".
#' @param min.detections Optional. Discard individuals with less than this number of detections.
#' @param min.days Optional. Minimum number of days an individual must be detected for it to be included;
#' individuals with fewer days of detections will be discarded.
#' @param hours.threshold Discard single detections occurring alone within a predefined time interval.
#' @param max.speed Maximum allowed speed between detections. If null, no speed-based checks are performed.
#' @param speed.unit Units of supplied maximum speed. Either 'm/s' or 'km/h'.
#' @param acoustic.range Maximum assumed detection range of the acoustic receivers (in metres).
#' Used in the speed filter to account for uncertainties in animal positions within the radius of the receivers.
#' @param ... Further arguments passed to the \code{\link{calculateDistances}} function
#' (e.g. land.shape, epsg.code, grid.resolution and mov.directions).
#' @return A list containing both the filtered data and the rejected detections,
#' as well as a detailed summary of the filtered out data by individual.
#' @importFrom dplyr %>%
#' @export


filterDetections <- function(data, tagging.dates, id.col="ID", datetime.col="datetime",
                             min.detections=0, min.days=0, hours.threshold=24,
                             max.speed=NULL, speed.unit="m/s", acoustic.range=600,
                              land.shape=NULL, epsg.code=NULL,... ) {

  ##############################################################################
  ## Prepare data ##############################################################

  cat("Filtering out spurious data\n")

  if(!id.col %in% colnames(data)) {
    stop("ID column not found. Please specify the right name using the 'id.col' argument")
  }

  nfish <- nlevels(data[,id.col])
  n_total <- nrow(data)

  if(length(tagging.dates) != nfish) {
    stop("Number of animal IDs and tagging dates don't match")
    }

  if(!speed.unit %in% c("m/s", "km/h")) {
    stop("Wrong speed unit. Please select either 'm/s' or 'km/h'")
  }

  if(!datetime.col %in% colnames(data)) {
    stop("Datetime column not found. Please specify the right name using the 'datetime.col' argument")
  }

  # initialize variables
  rejected_start <- vector(mode="list", length=nfish)
  rejected_isolated <- vector(mode="list", length=nfish)
  rejected_speed <- vector(mode="list", length=nfish)
  rejected_min <- vector(mode="list", length=nfish)
  rejected_days <- vector(mode="list", length=nfish)

  # split data by individual
  data_individual <- split(data, f=data[,id.col])
  n_individual <- unlist(lapply(data_individual, nrow))


  ##############################################################################
  ## Apply temporal filters ####################################################

  # print to console and set progress bar
  cat("Applying temporal filters...\n")
  pb <- txtProgressBar(min=1, max=nfish, initial=0, style=3)

  # run loop and apply filters
  for (i in 1:nfish) {
    data_subset <- data_individual[[i]]
    data_subset <- data_subset[order(data_subset[,datetime.col]),]
    tag_date <- tagging.dates[i]

    # calculate time differences between consecutive detections (in hours)
    data_subset$hour_diff <- difftime(data.table::shift(data_subset[,datetime.col], type="lead"), data_subset[,datetime.col], units="hours")
    data_subset$hour_diff <- as.numeric(data_subset$hour_diff)

    # filter out detections occurring before tagging date
    detections_prior <- which(data_subset[,datetime.col]<tag_date)
    if(length(detections_prior)>1){
      rejected_start[[i]] <- data_subset[detections_prior,]
      rejected_start[[i]]$reason <- "before tagging date"
      data_subset <- data_subset[-detections_prior,]
    }

    # filter out detections occurring alone in periods > hours.threshold
    if(hours.threshold!=F){
      detections_isolated <- which(data_subset$hour_diff>hours.threshold & data.table::shift(data_subset$hour_diff, type="lag")>hours.threshold)
      detections_isolated <- c(detections_isolated, which(is.na(data_subset$hour_diff) & data.table::shift(data_subset$hour_diff, type="lag")>hours.threshold))
      detections_isolated <- c(detections_isolated, which(data_subset$hour_diff>hours.threshold & is.na(data.table::shift(data_subset$hour_diff, type="lag"))))
      detections_isolated <- unique(detections_isolated)
      if(length(detections_isolated)>0) {
        rejected_isolated[[i]] <- data_subset[detections_isolated,]
        rejected_isolated[[i]]$reason <- paste0("isolated within a ", hours.threshold, "-h period")
        data_subset <- data_subset[-detections_isolated,]
      }
    }

    # return data and show progress in console
    data_individual[[i]] <- data_subset
    setTxtProgressBar(pb,i)

  }

  # close progress bar
  close(pb)

  # reassemble data
  data_filtered <- do.call("rbind", data_individual)
  data_filtered <- data_filtered[order(data_filtered$datetime, data_filtered[,id.col]),]


  ##############################################################################
  ## Apply speed filter ########################################################

  # if max.speed is defined, calculate distance between each two consecutive detections
  if(!is.null(max.speed)){
    cat("Applying speed filter...\n")
    data_distances <- moby::calculateDistances(data_filtered, land.shape, epsg.code, ...)$data
    data_distances <- base::split(data_distances, f=data_distances[,id.col])
    # filter out detections above the maximum defined speed
    for(i in 1:length(data_distances)){
      data_subset <- data_distances[[i]]
      # account for detection range in estimated travelled distances
      data_subset$dist_min <- data_subset$dist_m - (acoustic.range*2)
      data_subset$dist_min[data_subset$dist_min<0] <- 0
      # re-calculate time differences between consecutive detections (in hours)
      data_subset$hour_diff <- difftime(data.table::shift(data_subset[,datetime.col], type="lead"), data_subset[,datetime.col],  units="hours")
      data_subset$hour_diff <- as.numeric(data_subset$hour_diff)
      # calculate speed
      if(speed.unit=="m/s"){data_subset$speed <- data_subset$dist_min / (data_subset$hour_diff*3600)}
      if(speed.unit=="km/h"){data_subset$speed <- (data_subset$dist_min/1000) / data_subset$hour_diff}
      detections_overspeed <- which(data_subset$speed>max.speed)
      if(length(detections_overspeed)>0){
        rejected_speed[[i]] <- data_subset[detections_overspeed,]
        rejected_speed[[i]]$reason <- paste("discarded based on a max speed threshold of", max.speed, speed.unit)
        data_subset <- data_subset[-detections_overspeed,]
      }
      data_distances[[i]] <- data_subset
    }
    # reassemble data
    data_filtered <- do.call("rbind", data_distances)
    data_filtered <- data_filtered[order(data_filtered[,datetime.col], data_filtered[,id.col]),]
  }

  ##############################################################################
  ## Apply minimum detections filter  ##########################################

  # if number detections < detections threshold remove animal from dataset
  if(min.detections>0){
    data_individual <- split(data_filtered, f=data_filtered[,id.col])
    for (i in 1:nfish) {
      if(nrow(data_individual[[i]])<min.detections & nrow(data_individual[[i]])>0){
        rejected_min[[i]] <- data_individual[[i]]
        rejected_min[[i]]$reason <- paste("individual with less than", min.detections, "detections")
        data_individual[[i]] <- data_filtered[0,]
      }
    }
    data_filtered <- do.call("rbind", data_individual)
  }


  ##############################################################################
  ## Apply minimum days filter  ################################################

  # if number logged dats < min days threshold remove animal from dataset
  if(min.days>0){
    data_individual <- split(data_filtered, f=data_filtered[,id.col])
    for (i in 1:nfish) {
      data_individual[[i]]$day <- strftime(data_individual[[i]][,datetime.col], "%d-%m-%Y", tz="UTC")
      ndays <- length(unique(data_individual[[i]]$day))
      #ndays <- difftime(max(data_individual[[i]][,datetime.col], na.rm=T), min(data_individual[[i]][,datetime.col], na.rm=T),  units="days")
      if(ndays<min.days & nrow(data_individual[[i]])>0){
        rejected_days[[i]] <- data_individual[[i]]
        rejected_days[[i]]$reason <- paste("individual with detections on fewer than", min.days)
        data_individual[[i]] <- data_filtered[0,]
      }
    }
    data_filtered <- do.call("rbind", data_individual)
  }


  ##############################################################################
  ## Format and return results #################################################

  data_filtered <- data_filtered[,-which(colnames(data_filtered) %in% c("hour_diff", "dist_m", "dist_min", "speed"))]
  rownames(data_filtered) <- NULL
  n_removed_start <- unlist(lapply(rejected_start, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_isolated <- unlist(lapply(rejected_isolated, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_speed <- unlist(lapply(rejected_speed, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_min <- unlist(lapply(rejected_min, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_days <- unlist(lapply(rejected_days, function(x) ifelse(!is.null(x), nrow(x), 0)))
  rejected_start <- do.call("rbind", rejected_start)
  rejected_isolated <- do.call("rbind", rejected_isolated)
  rejected_speed <- do.call("rbind", rejected_speed)
  rejected_min <- do.call("rbind", rejected_min)
  rejected_days <- do.call("rbind", rejected_days)
  data_rejected <- do.call(plyr::rbind.fill, list(rejected_start, rejected_isolated, rejected_speed, rejected_min, rejected_days))
  if(!is.null(data_rejected)) {
    data_rejected <- data_rejected %>% dplyr::select(-reason, reason)
  }else{
    data_rejected <- data_filtered[0,]
  }

  data_rejected <- data_rejected[order(data_rejected[,id.col], data_rejected[,datetime.col]),]
  rownames(data_rejected) <- NULL

  # calculate stats
  n_removed_total <- nrow(data_rejected)
  percent_total <- sprintf("%.0f", n_removed_total/n_total*100)
  percent_start <- sprintf("%.0f", nrow(rejected_start)/n_total*100)
  percent_start <- ifelse(length(percent_start)>0, paste0(" (", percent_start, "%)"), "")
  percent_isolated <- sprintf("%.0f", nrow(rejected_isolated)/n_total*100)
  percent_isolated <- ifelse(length(percent_isolated)>0, paste0(" (", percent_isolated, "%)"), "")
  percent_speed <- sprintf("%.0f", nrow(rejected_speed)/n_total*100)
  percent_speed <- ifelse(length(percent_speed)>0, paste0(" (", percent_speed, "%)"), "")
  percent_min <- sprintf("%.0f",  nrow(rejected_min)/n_total*100)
  percent_min <- ifelse(length(percent_min)>0, paste0(" (", percent_min, "%)"), "")
  percent_days <- sprintf("%.0f",  nrow(rejected_days)/n_total*100)
  percent_days <- ifelse(length(percent_days)>0, paste0(" (", percent_days, "%)"), "")

  ids_total <- length(which(unlist(lapply(data_individual, nrow))<n_individual))
  ids_start <- length(which(n_removed_start>0))
  ids_isolated <- length(which(n_removed_isolated>0))
  ids_speed <- length(which(n_removed_speed>0))
  ids_min <- length(which(n_removed_min>0))
  ids_days <- length(which(n_removed_days>0))
  ids_discarded <- nlevels(data[,id.col]) - length(which(table(data_filtered[,id.col])>0))
  percent_discarded <- sprintf("%.0f", ids_discarded/nfish*100)
  percent_discarded <- ifelse(length(percent_discarded)>0, paste0(" (", percent_discarded, "%)"), "")

  # print stats to console
  cat("\n")
  cat(paste0("Detections removed = ", n_removed_total, " (", percent_total, "%) from a total of ", n_total, "\n"))
  cat(paste0("  • before tagging: ", sum(n_removed_start), percent_start, " from ", ids_start, " individuals\n"))
  if(hours.threshold!=F){
    cat(paste0("  • isolated ", hours.threshold, "h radius: ", sum(n_removed_isolated),  percent_isolated, " from ", ids_isolated, " individuals\n"))}
  if(!is.null(max.speed)){
    cat(paste0("  • above max speed (", max.speed, " ", speed.unit, "): ", sum(n_removed_speed), percent_speed, " from ", ids_speed, " individuals\n"))}
  if(min.detections>0){
    cat(paste0("  • < ", min.detections, " detections: ", sum(n_removed_min), percent_min, " from ", ids_min, " individuals\n"))}
  if(min.days>0){
    cat(paste0("  • < ", min.days, " days with detection(s): ", sum(n_removed_days), percent_days, " from ", ids_days, " individuals\n"))}
  cat(paste0("Individuals discarded = ", ids_discarded, percent_discarded, " from a total of ", nfish, "\n"))

  # create detailed summary table
  n_removed_invidual <- n_removed_start + n_removed_isolated + n_removed_speed + n_removed_min
  filter_summary <- data.frame("ID"=levels(data[,id.col]), "Raw Detections"=n_individual,
                        "Before tagging"=n_removed_start, "Isolated"=n_removed_isolated,
                        "Speed Threshold"=n_removed_speed, "Min Detections Threshold"=n_removed_min,
                        "Min Days Threshold"=n_removed_days, "Total removed"=n_removed_invidual, check.names=F, row.names=NULL)
  percent_removed <- sprintf("%.0f", n_removed_invidual/n_individual*100)
  filter_summary$`Total removed` <- paste0(filter_summary$`Total removed`, " (",  percent_removed, "%)")
  filter_summary$`Total removed`[percent_removed=="0"] <- "-"

  # return results
  results <- list("data"=data_filtered, "data_rejected"=data_rejected, "summary"=filter_summary)
  return(results)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
