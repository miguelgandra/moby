#######################################################################################################
# Function to calculate distance to nearest shore/land feature ########################################
#######################################################################################################
#'
#' Calculate distance to nearest shore/land feature and other related metrics.
#'
#' @description This function calculates the distance to the nearest shore or land feature for each
#' animal position. Additionally, it computes movement-related metrics such as the direction of movement
#' relative to the shore (inshore, offshore, or alongshore) and the proportion of movement
#' directed towards or away from land, if step lengths between consecutive positions are provided.
#'
#' @inheritParams setDefaults
#' @param data A data frame with animal detections/positions, containing longitude and latitude values (unprojected).
#' @param land.shape Optional. A shapefile containing coastlines or landmasses. It can be supplied as
#' an 'sf' object or as an object of class 'SpatialPolygonsDataFrame' or 'SpatialPolygons'.
#' If the provided object is not of class 'sf', the function will attempt to
#' convert it to an 'sf' object for compatibility with subsequent spatial operations.
#' @param mov.threshold Movement threshold (0 - 1). It establishes the proportion of a step length that must be
#' perpendicular to shore (away or approaching) in order for a movement to be classified as inshore or offshore
#' (e.g., with the default mov.threshold of 0.5: if the difference in distance to land between each two
#' consecutive positions corresponds to more than 50% of the total step length, the movement
#' is classified either as inshore or offshore, according to its direction).
#' @param dist.col Name of the column containing step lengths in meters (distance between consecutive positions).
#' This column can be generated using the \code{calculateDistances} function. If left NULL,
#' the function will only compute distances to land without calculating additional movement metrics.
#' @return A data frame containing distances to nearest shore and multiple related metrics, including
#' direction of movement (inshore/offshore/alongshore) and percentage of movement
#' perpendicular to shore (land_ratio).
#'
#' @return A data frame that includes the original input data along with the following additional columns:
#' \itemize{
#'   \item \code{land_dist}: The distance (in meters) to the nearest shore or land feature for each position.
#'   \item \code{land_diff}: The difference in distance to land between consecutive positions.
#'   \item \code{land_ratio}: The ratio of the absolute change in distance to land relative to the total distance traveled (if \code{dist.col} is provided).
#'   \item \code{direction}: The direction of movement relative to land (inshore, offshore, alongshore).
#' }
#'
#' @export


calculateLandDists <- function(data,
                               land.shape,
                               epsg.code = getDefaults("epsg"),
                               mov.threshold = 0.5,
                               id.col = getDefaults("id"),
                               lon.col = getDefaults("lon"),
                               lat.col = getDefaults("lat"),
                               dist.col = NULL){


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # print to console
  .printConsole("Calculating distances to nearest land")

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  land.shape <- reviewed_params$land.shape

  # check if data contains dist.col
  if(!is.null(dist.col)){
    if(!dist.col %in% colnames(data)){
      stop("Distance column not found. Please assign the correct column name with 'dist.col'.", call.=FALSE)
    }
  }else{
    warning("The 'dist.col' parameter was not supplied, additional movement metrics will not be calculated.", call.=FALSE)
  }

  # check crs of spatial objects
  coords <- sf::st_as_sf(data, coords=c(lon.col, lat.col))
  coords_crs <- .checkProjection(coords)
  land_shape_crs <- .checkProjection(land.shape)

  # convert spatial objects if required
  if(coords_crs=="projected"){
    spatial_data <- .processSpatial(coords, land.shape, epsg.code)
    coords <- spatial_data$coords
    coords <- sf::st_transform(coords, crs=4326)
  }

  # convert land.shape to geographic format (EPSG:4326).
  if(land_shape_crs=="projected"){
    land.shape <- sf::st_transform(land.shape, crs=4326)
  }

  ##############################################################################
  ## Calculate distance to nearest land ########################################
  ##############################################################################

  # convert 'sf' objects into coordinate matrices.
  coord_matrix <- sf::st_coordinates(coords)
  land_shape_matrix <- sf::st_coordinates(land.shape)[,c("X","Y")]

  # calculate the distance from each point in 'coord_matrix' to every point in 'land_shape_matrix'.
  land_dists <- terra::distance(coord_matrix, land_shape_matrix, lonlat=TRUE, pairwise=FALSE)

  # for each point, retrieve the minimum distance to the nearest landmass
  data$land_dist <- apply(land_dists, 1, min, na.rm=T)


  ##############################################################################
  ## Calculate additional metrics ##############################################
  ##############################################################################

  # if distances between consecutive points have been calculated (i.e., 'dist.col' exists),
  # we proceed to calculate additional metrics such as land movement ratios and movement direction
  if(!is.null(dist.col)){

    # split data by individual
    data_individual <- split(data, f=data[,id.col], drop=T)

    # initialize the progress bar
    pb <- txtProgressBar(min=1, max=length(data_individual), initial=0, style=3)

    # loop over each individual
    for(i in 1:length(data_individual)){
      # select data
      data_animal <- data_individual[[i]]
      # calculate the difference in distance to land between consecutive records
      data_animal$land_diff <- dplyr::lead(data_animal$land_dist) - data_animal$land_dist
      data_animal$land_diff <- round(data_animal$land_diff, 2)
      # calculate the ratio of the absolute land distance difference to the distance traveled
      data_animal$land_ratio <- round(abs(data_animal$land_diff)/data_animal[,dist.col], 2)
      data_animal$land_ratio[data_animal[,dist.col]==0] <- 0
      data_animal$direction <- NA
      # classify movement direction based on the 'land_ratio' and 'land_diff' values
      data_animal$direction[data_animal$land_ratio>=mov.threshold & data_animal$land_diff>0] <- "offshore"
      data_animal$direction[data_animal$land_ratio>=mov.threshold & data_animal$land_diff<0] <- "inshore"
      data_animal$direction[data_animal$land_ratio<mov.threshold] <- "alongshore"
      # set movement direction to NA where no distance was traveled
      data_animal$direction[data_animal[,dist.col]==0] <- NA
      # store the processed data back into the list
      data_individual[[i]] <- data_animal
      # update the progress bar
      setTxtProgressBar(pb, i)
    }

    # close the progress bar once the loop is complete
    close(pb)

    # combine the list of individual datasets back into a single data frame
    data <- do.call("rbind", data_individual)
    rownames(data) <- NULL

    # convert direction column to factor
    data$direction <- factor(data$direction, levels=c("inshore", "alongshore", "offshore"))
  }


  ##############################################################################
  ## Return results ############################################################
  ##############################################################################

  # print the total time taken for execution
  time.taken <- Sys.time() - start.time
  cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n"))

  # return the final processed data
  return(data)

}


#######################################################################################################
#######################################################################################################
#######################################################################################################
