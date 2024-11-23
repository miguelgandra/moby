#######################################################################################################
## Filter detections data #############################################################################
#######################################################################################################

#' Filter Detections
#'
#' @description This function filters out spurious animal detections based on
#' multiple criteria to ensure data accuracy. It helps users identify and remove
#' false positives and discard individuals with few or no valid detections.
#' It implements six distinct criteria:
#'
#' 1) detections occurring before the animal's tagging date
#' 2) detections occurring after an optional cut-off date
#' 3) isolated detections within a specified time interval
#' 4) detections implying unrealistic movement speeds
#' 5) the total number of detections
#' 6) the total number of logged days
#'
#' For speed filtering, the function iteratively removes detections that exceed the
#' maximum speed threshold, recalculating distances until no violations remain or the
#' maximum iteration limit is reached. When an unrealistic speed is detected, a hierarchical
#' clustering approach is applied to the spatial coordinates of neighboring detections.
#' This helps identify the outlier detection, likely responsible for triggering the overspeed
#' threshold, by considering its context within surrounding detections.
#'
#' Additionally, `filterDetections` accounts for land topographies and complex
#' coastlines when calculating animal speeds if a `land.shape` (shapefile containing
#' coastlines or landmasses) is provided, leading to more accurate distance estimates.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing raw animal detections.
#' @param cutoff.dates Optional. A POSIXct vector containing the estimated expiration dates of the tags
#' or any other cut-off date beyond which detections should be discarded. The length of this vector
#' should match the number of unique animal IDs. Alternatively, if a single value is provided, it will be applied to all IDs.
#' @param min.detections Optional. Discard individuals with fewer than this number of detections.
#' @param min.days Optional. Minimum number of days an individual must be detected for it to be included;
#' individuals with fewer days of detections will be discarded.
#' @param hours.threshold Discard single detections occurring alone within a predefined time interval.
#' @param max.speed Maximum allowed speed between detections. If null, no speed-based checks are performed.
#' @param speed.unit Units of supplied maximum speed. Either 'm/s' or 'km/h'.
#' @param max.iterations An integer specifying the maximum number of iterations
#' for the speed filtering process. The filtering loop will stop early if this number
#' of iterations is reached, even if detections are still being removed. The default is
#' `Inf` (unlimited iterations).
#' @param acoustic.range Maximum assumed detection range of the acoustic receivers (in metres).
#' Used in the speed filter to account for uncertainties in animal positions within the radius of the receivers,
#' when estimating minimum consecutive distances.
#' @param land.shape Optional. A shapefile containing coastlines or landmasses. It can be supplied as
#' an 'sf' object or as an object of class 'SpatialPolygonsDataFrame' or 'SpatialPolygons'.
#' If the provided object is not of class 'sf', the function will attempt to
#' convert it to an 'sf' object for compatibility with subsequent spatial operations.
#' If not supplied, the function will calculate linear (great-circle) distances instead of shortest in-water
#' distances. Linear distances are much faster to compute but may be less realistic or accurate in cases
#' where positions are separated by landmasses or other spatial obstacles. This parameter
#' is passed to the \code{\link{calculateTracks}} function and is only used when
#' a \code{max.speed} is defined.
#' @param epsg.code Coordinate reference system used to project positions (class 'CRS').
#' If not supplied, CRS is assumed to be the same as in land.shape.
#' @param ... Further arguments passed to the \code{\link{calculateTracks}} function
#' (e.g., `grid.resolution`, `mov.directions` and `cores`).
#' @return A list containing both the filtered data and the rejected detections,
#' as well as a detailed summary of the filtered out data by individual.
#' @importFrom dplyr %>%
#' @export


filterDetections <- function(data,
                             tagging.dates = getDefaults("tagging.dates"),
                             cutoff.dates = NULL,
                             id.col = getDefaults("ID"),
                             datetime.col = getDefaults("datetime"),
                             lon.col = getDefaults("lon"),
                             lat.col = getDefaults("lat"),
                             min.detections = 0,
                             min.days = 0,
                             hours.threshold = 24,
                             max.speed = NULL,
                             speed.unit = "m/s",
                             max.iterations = Inf,
                             acoustic.range = 600,
                             land.shape = NULL,
                             epsg.code = getDefaults("epsg"),
                             ... ) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # print message
  .printConsole("Filtering detections")
  cat("\n")

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  land.shape <- reviewed_params$land.shape

  # check if speed.unit is valid
  if(!speed.unit %in% c("m/s", "km/h")) stop("Wrong speed unit. Please select either 'm/s' or 'km/h'", call.=FALSE)

  # check and replicate cutoff.dates if it is a single value
  if(!is.null(cutoff.dates) && length(cutoff.dates)==1){
    cutoff.dates <- rep(cutoff.dates, nlevels(data[,id.col]))
  }

  # capture the name of the 'land.shape' object
  if(!is.null(land.shape)){
    land_shape_name <- deparse(substitute(land.shape))
    grid_resolution_val <- formals(calculateTracks)$grid.resolution
    mov_directions_val <- formals(calculateTracks)$mov.directions
  # if no land shape is provided, set related parameters to NULL
  }else{
    land_shape_name <- NULL
    grid_resolution_val <- NULL
    mov_directions_val <- NULL
  }

  # set the maximum number of iterations (if not defined)
  if(is.null(max.iterations)) max.iterations <- Inf


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # manage spatial objects
  coords <- sf::st_as_sf(data, coords=c(lon.col, lat.col))
  coords_crs <- .checkProjection(coords)
  spatial_data <- .processSpatial(coords, land.shape, epsg.code)
  coords <- spatial_data$coords
  land.shape <- spatial_data$spatial.layer
  epsg.code <- spatial_data$epsg.code

  # extract processed longitude and latitude into temporary columns
  data$tmp_lon <- sf::st_coordinates(coords)[,1]
  data$tmp_lat <- sf::st_coordinates(coords)[,2]

  # add a unique identifier for each row
  data$row_index <- 1:nrow(data)

  # save nº of individuals and rows
  nfish <- nlevels(data[,id.col])
  n_total <- nrow(data)

  # initialize variables
  rejected_start <- vector(mode="list", length=nfish)
  rejected_cutoff <- vector(mode="list", length=nfish)
  rejected_isolated <- vector(mode="list", length=nfish)
  rejected_speed <- vector(mode="list", length=nfish)
  rejected_min <- vector(mode="list", length=nfish)
  rejected_days <- vector(mode="list", length=nfish)

  # split data by individual
  data_individual <- split(data, f=data[,id.col])
  n_individual <- unlist(lapply(data_individual, nrow))


  ##############################################################################
  ## Apply temporal filters ####################################################
  ##############################################################################

  # print to console and set progress bar
  .printConsole("Applying temporal filters...")
  pb <- txtProgressBar(min=1, max=nfish, initial=0, style=3)

  # run loop and apply filters
  for (i in 1:nfish) {
    data_subset <- data_individual[[i]]
    data_subset <- data_subset[order(data_subset[,datetime.col]),]
    tag_date <- tagging.dates[i]
    off_date <- cutoff.dates[i]

    # calculate time differences between consecutive detections (in hours)
    data_subset$hour_diff <- difftime(dplyr::lead(data_subset[,datetime.col]), data_subset[,datetime.col], units="hours")
    data_subset$hour_diff <- as.numeric(data_subset$hour_diff)

    # filter out detections occurring before tagging date
    detections_prior <- which(data_subset[,datetime.col]<tag_date)
    if(length(detections_prior)>1){
      rejected_start[[i]] <- data_subset[detections_prior,]
      rejected_start[[i]]$reason <- "before tagging date"
      data_subset <- data_subset[-detections_prior,]
    }

    # filter out detections occurring after cut-off date
    detections_after <- which(data_subset[,datetime.col]>off_date)
    if(length(detections_after)>1){
      rejected_cutoff[[i]] <- data_subset[detections_after,]
      rejected_cutoff[[i]]$reason <- "after cut-off date"
      data_subset <- data_subset[-detections_after,]
    }

    # filter out detections occurring alone in periods > hours.threshold
    if(hours.threshold!=FALSE){
      detections_isolated <- which(data_subset$hour_diff>hours.threshold & dplyr::lag(data_subset$hour_diff)>hours.threshold)
      detections_isolated <- c(detections_isolated, which(is.na(data_subset$hour_diff) & dplyr::lag(data_subset$hour_diff)>hours.threshold))
      detections_isolated <- c(detections_isolated, which(data_subset$hour_diff>hours.threshold & is.na(dplyr::lag(data_subset$hour_diff))))
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
  cat("\n")

  # reassemble data
  data_filtered <- do.call("rbind", data_individual)
  data_filtered <- data_filtered[order(data_filtered[,datetime.col], data_filtered[,id.col]),]


  ##############################################################################
  ## Apply speed filter ########################################################
  ##############################################################################

  # if max.speed is defined, calculate distance between each two consecutive detections
  if (!is.null(max.speed)) {

    # initialize run counter
    iteration <- 1

    # print message to console
    .printConsole("Applying speed filter...")

    # counter for removed detections and discarded detections with small time differences
    removed_detections <- 0
    small_time_diff <- 0

    # calculate initial distances
    tracks <- moby::calculateTracks(data_filtered, land.shape, epsg.code, id.col=id.col, lon.col="tmp_lon", lat.col="tmp_lat", ...)

    # extract the 'data' and 'transition_layer' components
    data_distances <- tracks$data
    if(!is.null(land.shape)){
      transition_layer <- tracks$transition_layer
    }else{
      transition_layer <- NULL
    }

    # split tracks by individual
    data_distances <- base::split(data_distances, f=data_distances[,id.col])

    # loop until no detections above the max speed are found
    repeat {

      # print message to console with iteration number
      #cat(paste0("Iteration ",iteration, "\n"))
      cat(paste0("\rIteration: ", iteration))
      flush.console()

      # counter for removed detections
      removed_detections <- 0

      # iterate through each track
      for(i in 1:length(data_distances)) {
        # get current subset
        data_subset <- data_distances[[i]]
        if(nrow(data_subset)==0) next
        # add a unique identifier for each row
        data_subset$row_id <- 1:nrow(data_subset)
        # account for detection range in estimated traveled distances
        data_subset$dist_min <- data_subset$dist_m - (acoustic.range * 2)
        data_subset$dist_min[data_subset$dist_min < 0] <- 0
        # calculate time differences
        data_subset$hour_diff <- difftime(dplyr::lead(data_subset[, datetime.col]), data_subset[, datetime.col], units = "hours")
        data_subset$hour_diff <- as.numeric(data_subset$hour_diff)
        # calculate speed
        if(speed.unit == "m/s") data_subset$speed <- data_subset$dist_min / (data_subset$hour_diff * 3600)
        if(speed.unit == "km/h") data_subset$speed <- (data_subset$dist_min / 1000) / data_subset$hour_diff
        data_subset$speed[data_subset$hour_diff==0] <- 0
        # identify overspeed detections
        detections_overspeed <- which(data_subset$speed > max.speed)
        # if there any overspeed detection is identified
        if(length(detections_overspeed) > 0) {
          # count overspeed detections resulting from a time difference <= 1 min
          small_time_flagged <- which(data_subset$hour_diff[detections_overspeed] <= (1/60))
          if(length(small_time_flagged) > 0) small_time_diff <- small_time_diff + length(small_time_flagged)
          # iterate over each overspeed indice and check problematic detections
          for (flagged_detec in detections_overspeed) {
            # define the range of indices around the current overspeed detection for context
            overspeed_context <- (flagged_detec - 3):(flagged_detec + 3)
            overspeed_context <- overspeed_context[overspeed_context > 0 & overspeed_context <= nrow(data_subset)]
            # extract the context coordinates for the current overspeed detection
            context_coords <- sf::st_as_sf(data_subset[overspeed_context, ], coords=c("tmp_lon", "tmp_lat"), crs=epsg.code)
            # create distance matrix
            dist_matrix <- .calculateDistMatrix(data_subset[overspeed_context,], land.shape, epsg.code, id.col=id.col,
                                                lon.col="tmp_lon", lat.col="tmp_lat", transition.layer=transition_layer)
            # perform hierarchical clustering
            hclust_result <- stats::hclust(dist_matrix)
            # specify the number of clusters (n=2)
            clusters <- stats::cutree(hclust_result, k=2)
            # find the outlier (least frequent cluster)
            outlier_indice <- which(clusters %in% names(which.min(table(clusters))))[1]
            # update detections_overspeed with the outlier detection
            detections_overspeed[detections_overspeed==flagged_detec] <- overspeed_context[outlier_indice]
          }
          detections_overspeed <- unique(detections_overspeed)
          # append the removed detections to the corresponding individual's list
          removed_data <- data_subset[detections_overspeed, ]
          removed_data$reason <- paste("discarded based on a max speed threshold of", max.speed, speed.unit)
          removed_data$iteration <- iteration
          rejected_speed[[i]] <- append(rejected_speed[[i]], list(removed_data))
          # identify affected indices using row_id
          affected_indices <- lapply(detections_overspeed, function(x) c(x-1, x+1))
          affected_indices <- lapply(affected_indices, function(x) x[x>0 & x<=nrow(data_subset)])
          recalc_indices <- unlist(lapply(affected_indices, function(x) x[1]))
          affected_indices <- unlist(affected_indices)
          # remove overspeed detections
          data_subset <- data_subset[-detections_overspeed, ]
          removed_detections <- removed_detections + length(detections_overspeed)
          # update distances for affected indices
          recalculated <- suppressWarnings(moby::calculateTracks(data_subset[data_subset$row_id %in% affected_indices, ], land.shape, epsg.code,
                                                                 id.col=id.col, lon.col="tmp_lon", lat.col="tmp_lat",
                                                                 transition.layer=transition_layer, verbose=FALSE, ...)$data)
          data_subset[data_subset$row_id %in% recalc_indices, "dist_m"] <- recalculated[recalculated$row_id %in% recalc_indices, "dist_m"]
        }

        # update the data subset in the list
        data_distances[[i]] <- data_subset
      }

      # stop if no detections were removed or if the maximum number of iterations is reached
      if (removed_detections == 0 || iteration >= max.iterations) break

      # increment iteration counter
      iteration <- iteration + 1
    }

    # reassemble data
    data_filtered <- do.call("rbind", data_distances)
    data_filtered <- data_filtered[order(data_filtered[, datetime.col], data_filtered[, id.col]), ]

    # reassemble flagged data
    rejected_speed <- lapply(rejected_speed, function(x) {if(!is.null(x)) do.call("rbind", x) else NULL})

    # issue a warning about overspeed detections with small time differences
    if(small_time_diff>0) {
      warning_str <- paste0("A total of ", small_time_diff, " discarded overspeed detections had small datetime differences (< 1 min). ",
                            "These may represent valid detections that were flagged due to an underestimated acoustic range or receiver clock drifts. ",
                            "Please review the flagged detections in the 'data_discarded' data frame included in the results.")
      warning(paste(strwrap(warning_str, width=getOption("width")), collapse="\n"), call. = FALSE)
    }

    # print a empty line
    cat("\n")
  }


  ##############################################################################
  ## Apply minimum detections filter  ##########################################
  ##############################################################################

  # if number detections < detections threshold remove animal from dataset
  if(min.detections>0){

    # print message to console
    .printConsole("Applying min detections filter...")

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
  ##############################################################################

  # if number logged dates < min days threshold remove animal from dataset
  if(min.days>0){

    # print message to console
    .printConsole("Applying min days filter...")

    data_individual <- split(data_filtered, f=data_filtered[,id.col])
    for (i in 1:nfish) {
      data_individual[[i]]$day <- strftime(data_individual[[i]][,datetime.col], "%d-%m-%Y", tz="UTC")
      ndays <- length(unique(data_individual[[i]]$day))
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
  ##############################################################################

  data_filtered <- data_filtered[,-which(colnames(data_filtered) %in% c("tmp_lon", "tmp_lat", "hour_diff", "dist_m", "dist_min", "speed", "row_id"))]
  rownames(data_filtered) <- NULL
  n_removed_start <- unlist(lapply(rejected_start, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_cutoff <- unlist(lapply(rejected_cutoff, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_isolated <- unlist(lapply(rejected_isolated, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_speed <- unlist(lapply(rejected_speed, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_min <- unlist(lapply(rejected_min, function(x) ifelse(!is.null(x), nrow(x), 0)))
  n_removed_days <- unlist(lapply(rejected_days, function(x) ifelse(!is.null(x), nrow(x), 0)))
  rejected_start <- do.call("rbind", rejected_start)
  rejected_cutoff <- do.call("rbind", rejected_cutoff)
  rejected_isolated <- do.call("rbind", rejected_isolated)
  rejected_speed <- do.call("rbind", rejected_speed)
  rejected_min <- do.call("rbind", rejected_min)
  rejected_days <- do.call("rbind", rejected_days)
  data_discarded <- do.call(plyr::rbind.fill, list(rejected_start, rejected_cutoff, rejected_isolated, rejected_speed, rejected_min, rejected_days))
  if(!is.null(data_discarded)) {
    data_discarded <- data_discarded %>% dplyr::select(-reason, reason)
  }else{
    data_discarded <- data_filtered[0,]
  }
  data_discarded <- data_discarded[order(data_discarded[,id.col], data_discarded[,datetime.col]),]
  data_discarded <- data_discarded[,-which(colnames(data_discarded) %in% c("tmp_lon", "tmp_lat", "dist_min",  "speed", "row_id"))]
  data_discarded$hour_diff <- round(data_discarded$hour_diff, 3)
  rownames(data_discarded) <- NULL

  # calculate stats
  n_removed_total <- nrow(data_discarded)
  percent_total <- sprintf("%.0f", n_removed_total/n_total*100)
  percent_start <- sprintf("%.0f", nrow(rejected_start)/n_total*100)
  percent_start <- ifelse(length(percent_start)>0, paste0(" (", percent_start, "%)"), "")
  percent_cutoff <- sprintf("%.0f", nrow(rejected_cutoff)/n_total*100)
  percent_cutoff <- ifelse(length(percent_cutoff)>0, paste0(" (", percent_cutoff, "%)"), "")
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
  ids_cutoff <- length(which(n_removed_cutoff>0))
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
  cat(paste0("  \u2022 before tagging: ", sum(n_removed_start), percent_start, " from ", ids_start, " individuals\n"))
  if(!is.null(cutoff.dates)){
    cat(paste0("  \u2022 after cut-off date: ", sum(n_removed_cutoff), percent_cutoff, " from ", ids_cutoff, " individuals\n"))}
  if(hours.threshold!=FALSE){
    cat(paste0("  \u2022 isolated ", hours.threshold, "h interval: ", sum(n_removed_isolated),  percent_isolated, " from ", ids_isolated, " individuals\n"))}
  if(!is.null(max.speed)){
    cat(paste0("  \u2022 above max speed (", max.speed, " ", speed.unit, "): ", sum(n_removed_speed), percent_speed, " from ", ids_speed, " individuals\n"))}
  if(min.detections>0){
    cat(paste0("  \u2022 < ", min.detections, " detections: ", sum(n_removed_min), percent_min, " from ", ids_min, " individuals\n"))}
  if(min.days>0){
    cat(paste0("  \u2022 < ", min.days, " days with detection(s): ", sum(n_removed_days), percent_days, " from ", ids_days, " individuals\n"))}
  cat(paste0("Individuals discarded = ", ids_discarded, percent_discarded, " from a total of ", nfish, "\n\n"))

  # print time taken
  time.taken <- Sys.time() - start.time
  cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n\n"))


  # create detailed summary table
  n_removed_invidual <- n_removed_start + n_removed_isolated + n_removed_speed + n_removed_min
  filter_summary <- data.frame("ID"=levels(data[,id.col]), "Raw Detections"=n_individual,
                        "Before tagging"=n_removed_start, "After cut-off date"=n_removed_cutoff,
                        "Isolated"=n_removed_isolated, "Speed Threshold"=n_removed_speed,
                        "Min Detections Threshold"=n_removed_min, "Min Days Threshold"=n_removed_days,
                        "Total removed"=n_removed_invidual, check.names=F, row.names=NULL)
  percent_removed <- sprintf("%.0f", n_removed_invidual/n_individual*100)
  filter_summary$`Total removed` <- paste0(filter_summary$`Total removed`, " (",  percent_removed, "%)")
  filter_summary$`Total removed`[percent_removed=="0"] <- "-"
  filter_summary$`Total removed`[n_individual=="0"] <- "-"

  # aggregate results
  results <- list("data"=data_filtered, "data_discarded"=data_discarded, "summary"=filter_summary)

  # create new attributes to save relevant params
  attr(results, 'min.detections') <- min.detections
  attr(results, 'min.days') <- min.days
  attr(results, 'hours.threshold') <- hours.threshold
  attr(results, 'max.speed') <- max.speed
  attr(results, 'speed.unit') <- speed.unit
  attr(results, 'acoustic.range') <- acoustic.range
  attr(results, 'land.shape') <- land_shape_name
  attr(results, 'epsg.code') <- epsg.code
  attr(results, 'grid.resolution') <- grid_resolution_val
  attr(results, 'mov.directions') <- mov_directions_val
  attr(results, 'processing.date') <- Sys.time()

  # return results
  return(results)
}


################################################################################
# Define calculateDistMatrix helper function ###################################
################################################################################

#' @description Helper function to calculate shortest in-water dissimilarity matrix
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd


.calculateDistMatrix <- function(data, land.shape, epsg.code, id.col, lon.col, lat.col,
                                 transition.layer, ...) {

  # number of rows
  n <- nrow(data)

  # initialize an empty matrix for dissimilarity (n x n)
  dissimilarity_matrix <- matrix(NA, nrow=n, ncol=n)

  # loop over all unique row pairs to compute distances
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i != j) {
        # create a small dataset containing just the two rows
        pair <- data[c(i, j), ]

        # calculate shortest in-water distance for this pair
        dist_result <- suppressWarnings(moby::calculateTracks(pair, land.shape, epsg.code, id.col=id.col,
                                                              lon.col=lon.col, lat.col=lat.col, transition.layer=transition.layer,
                                                              verbose=FALSE, ...)$data)

        # extract the calculated distance (dist_m) and store it
        dissimilarity_matrix[i, j] <- dist_result$dist_m[1]
      } else {
        # set the diagonal to 0
        dissimilarity_matrix[i, j] <- 0
      }
    }
  }

  # convert to a dist object for compatibility with clustering functions
  dissimilarity_dist <- as.dist(dissimilarity_matrix)
  return(dissimilarity_dist)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
