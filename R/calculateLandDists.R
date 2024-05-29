#######################################################################################################
# Function to calculate distance to nearest shore/land feature ########################################
#######################################################################################################
#'
#' Calculate distance to nearest shore/land feature and other related metrics.
#'
#' @description Function to calculate distance to nearest land and additional related metrics,
#' including direction of movement (inshore/offshore/alongshore).

#' @param data A data frame with animal detections/positions, containing longitude and latitude values (unprojected).
#' @param land.shape A shapefile containing coastlines or land surfaces.
#' @param mov.threshold Movement threshold (0 - 1). It establishes the proportion of a step length that must be
#' perpendicular to shore (away or approaching) in order for a movement to be classified as inshore or offshore
#' (e.g., with the default mov.threshold of 0.5: if the difference in distance to land between each two
#' consecutive positions corresponds to more than 50% of the total step length, the movement
#' is classified either as inshore or offshore, according to its direction).
#' @param id.col Name of the column containing animal IDs Defaults to 'ID'.
#' @param lon.col Name of the column containing longitude values. Defaults to 'lon'.
#' @param lat.col Name of the column containing latitude values. Defaults to 'lat'.
#' @param dist.col Name of the column containing step lengths in meters (distance between each consecutive position).
#' If left null, the function only calculates distance to land and no additional metrics.
#' @return A data frame containing distances to nearest shore and multiple related metrics, including
#' direction of movement (inshore/offshore/alongshore) and percentage of movement
#' perpendicular to shore (land_ratio).
#' @export


calculateLandDists <- function(data, land.shape, mov.threshold=0.5, id.col="ID",
                                   lon.col="lon", lat.col="lat", dist.col=NULL){


  ############################################################################
  ## Initial checks ##########################################################
  ############################################################################

  # measure running time
  start.time <- Sys.time()

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # check if data contains lon.col and lat.col
  if(any(!c(lon.col, lat.col) %in% colnames(data))){
    stop("Longitude/latitude columns not found. Please assign the correct column names with 'lon.col' and 'lat.col'")
  }

  # check if data contains dist.col
  if(!is.null(dist.col)){
    if(!dist.col %in% colnames(data)){
      stop("Distance column not found. Please assign the correct column name with 'dist.col'")
    }
  }else{
    cat("Warning: dist.col not supplied, additional movement metrics will not be calculated\n")
  }

  # check if dataset coordinates are in geographic format (unprojected)
  geographic_coords <-  all(data[,lon.col]>=(-180) & data[,lon.col]<=180 & data[,lat.col]>=(-90) & data[,lat.col]<=90)
  if(is.na(geographic_coords)){
    stop("Longitudes/latitudes values are missing and/or in the wrong format\n")
  }


  #############################################################################
  ## Calculate distance to nearest land #######################################
  #############################################################################

  cat("Calculating distances to nearest land...\n")

  coords <- SpatialPoints(cbind(data[,lon.col], data[,lat.col]))
  projection(coords) <- CRS("+proj=longlat +datum=WGS84")
  coords <- sp::spTransform(coords, land.shape@proj4string)

  data$land_dist <- NA
  data$land_dist <- apply(rgeos::gDistance(coords, land.shape, byid=T), 2, min)

  data_individual <- split(data, f=data[,id.col], drop=T)
  pb <- txtProgressBar(min=1, max=length(data_individual), initial=0, style=3)
  for(i in 1:length(data_individual)){
    data_animal <- data_individual[[i]]
    data_animal$land_diff <- dplyr::lead(data_animal$land_dist) - data_animal$land_dist
    data_animal$land_diff <- round(data_animal$land_diff, 2)
    if(!is.null(dist.col)){
      data_animal$land_ratio <- round(abs(data_animal$land_diff)/data_animal[,dist.col], 2)
      data_animal$land_ratio[data_animal[,dist.col]==0] <- 0
      data_animal$direction <- NA
      data_animal$direction[data_animal$land_ratio>=mov.threshold & data_animal$land_diff>0] <- "offshore"
      data_animal$direction[data_animal$land_ratio>=mov.threshold & data_animal$land_diff<0] <- "inshore"
      data_animal$direction[data_animal$land_ratio<mov.threshold] <- "alongshore"
      data_animal$direction[data_animal[,dist.col]==0] <- NA
    }
    data_individual[[i]] <- data_animal
    setTxtProgressBar(pb, i)
  }
  close(pb)
  data <- do.call("rbind", data_individual)
  rownames(data) <- NULL

  # print time taken
  print(Sys.time() - start.time)

  # return results
  return(data)

}


#######################################################################################################
#######################################################################################################
#######################################################################################################
