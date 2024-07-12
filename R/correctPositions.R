#######################################################################################################
# Function to correct positions on land and relocate them to the nearest marine cell ##################
#######################################################################################################

#' Relocate points on land to the nearest marine cell.
#'
#' @description This function relocates positions on land to the nearest marine cell
#' using either a coastline shapefile or a raster containing land surfaces or bathymetry values.
#'
#' @param positions A data frame containing longitude and latitude values (unprojected lat-long coordinates),
#' corresponding to the inferred animals' positions.
#' @param layer A spatial object representing the geographic layer used to determine land and marine areas.
#' This can be either:\cr
#' - A 'RasterLayer' containing binary values or bathymetric information (negative values for underwater).
#' - A 'SpatialPolygons' object representing the coastline.
#' @param layer.type Either 'land', 'ocean' or 'bathy'. If type 'land' (binary) the raster is assumed
#' to contain NA values in water surfaces, if type 'ocean' (binary) the raster is assumed to
#' contain NA values in land surfaces, else if type 'bathy' the raster is assumed to
#' range between 0 (land) to max depth (positive or negative values).
#' @param depth.margin A numeric value. If set, shortest in-water paths will be calculated only
#' for cells with depths >= this threshold. Only takes effect when raster
#' is of type 'bathy'. Defaults to 0 (water vs land).
#' @param lon.col Name of the column containing longitude values. Defaults to 'lon'.
#' @param lat.col Name of the column containing latitude values. Defaults to 'lat'.
#' @param plot Logical, if TRUE, the function will plot the results. Defaults to FALSE.
#'
#' @return A data frame with corrected positions where points on land are relocated
#' to the nearest marine cell.
#'
#' @export

correctPositions <- function(positions, layer, layer.type="land", depth.margin=0,
                             lon.col=getDefaults("lon"), lat.col=getDefaults("lat"), plot=F) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # initial checks
  errors <- c()
  if (!inherits(positions, "data.frame")) errors <- c(errors, "positions must be a data frame")
  if(!lon.col %in% colnames(positions)) errors <- c(errors, "Longitude column not found.")
  if(!lat.col %in% colnames(positions)) errors <- c(errors, "Latitude column not found.")
  if (!layer.type %in% c("land", "ocean", "bathy")) errors <- c(errors, "layer.type must be either 'land', 'ocean' or 'bathy'")
  if (!inherits(layer, "RasterLayer") && !inherits(layer, "SpatialPolygons")) errors <- c(errors, "layer must be either a RasterLayer or a SpatialPolygons object")
  if (!inherits(depth.margin, "numeric")) errors <- c(errors, "depth.margin must be numeric")
  if (!inherits(plot, "logical"))  errors <- c(errors, "plot must be a logical value")
  if(length(errors)>0){
      stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
      stop(stop_message, call.=FALSE)
  }


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  if (inherits(layer, "RasterLayer")) {
    type <- "raster"
    proj <- raster::crs(layer)
  } else if (inherits(layer, "SpatialPolygons")) {
    type <- "shape"
    proj <- sf::st_crs(layer)
  }


  # convert values based on layer.type
  if (type == "raster") {
    if (layer.type == "land") {
      values(layer)[!is.na(values(layer))] <- 1
      values(layer)[is.na(values(layer))] <- 0
    } else if (layer.type == "ocean") {
      values(layer)[!is.na(values(layer))] <- 0
      values(layer)[is.na(values(layer))] <- 1
    } else if (layer.type == "bathy") {
      bathy_vals <- ifelse(mean(values(layer) > 0, na.rm = TRUE) > 0.5, "positive", "negative")
      if(bathy_vals=="positive") cat("Negative values in raster will be converted to NA\n")
      if(bathy_vals=="negative") cat("Positive values in raster will be converted to NA\n")
      if (bathy_vals == "negative") layer <- layer * -1
      values(layer)[values(layer) < depth.margin] <- 1
      values(layer)[values(layer) >= depth.margin] <- 0
      values(layer)[is.na(values(layer))] <- 1
      if(depth.margin!=0){
        cat(paste0("Positions will be relocated to depths > ", depth.margin, "\n"))}
    }
  }


  # check raster type and convert values if required
  if(type=="raster"){
    if(!layer.type %in% c("land","ocean","bathy")){
      stop("Raster type must be either 'land', 'ocean' or 'bathy'\n")
    } else if (layer.type %in% c("land","ocean")){
      cat(paste("Using raster of type", layer.type, "- depth margin ignored\n"))
    } else if (layer.type=="bathy") {
      bathy_vals <- ifelse(length(which(raster::values(layer)>0))/raster::ncell(layer) > 0.5, "positive", "negative")
      cat("Using raster of type bathy\n")
      if(bathy_vals=="positive"){cat("Negative values in raster will be converted to NA\n")}
      if(bathy_vals=="negative"){cat("Positive values in raster will be converted to NA\n"); layer <- layer*-1}
      if(depth.margin!=0){cat(paste0("Positions will be relocated to depths > ", depth.margin, "\n"))}
    }
  }

  # check if dataset coordinates are in geographic format (unprojected)
  geographic_coords <- all(positions[,lon.col]>=(-180) & positions[,lon.col]<=180 & positions[,lat.col]>=(-90) & positions[,lat.col]<=90)
  if(is.na(geographic_coords)){stop("Missing longitudes/latitudes\n")}
  if(geographic_coords==F){stop("Coordinates should be supplied in lat-long format (unprojected)\n")}


  ############################################################################
  # Create transition matrix #################################################
  ############################################################################

  # format layers
  if(type=="raster" & layer.type=="land"){
    raster::values(layer)[!is.na(raster::values(layer))] <- 1
    raster::values(layer)[is.na(raster::values(layer))] <- 0
  }
  if(type=="raster" & layer.type=="ocean"){
    raster::values(layer)[!is.na(raster::values(layer))] <- 0
    raster::values(layer)[is.na(raster::values(layer))] <- 1
  }
  if(type=="raster" & layer.type=="bathy"){
    raster::values(layer)[raster::values(layer)<depth.margin] <- 1
    raster::values(layer)[raster::values(layer)>depth.margin] <- 0
    raster::values(layer)[is.na(raster::values(layer))] <- 1
  }

  # create spatialPoints object and get land overlap indexes
  points <- sf::st_as_sf(positions, coords=c(lon.col, lat.col), crs=4326)
  points <- sf::st_transform(points, crs=sf::st_crs(proj))

  if(type=="raster"){
    point_values <- raster::extract(layer, points)
    pointsOnLand_indexes <- which(point_values==1)
    ocean_mask <- raster::rasterToPoints(layer, fun=function(x){x==0}, spatial=T)
  }else if(type=="shape"){
    pointsOnLand_indexes <- which(!is.na(sf::st_intersects(points, layer, sparse=FALSE)))
    boundary <- sf::st_as_sfc(sf::st_bbox(layer))
    boundary <- sf::st_set_crs(boundary, sf::st_crs(proj))
    ocean_mask <- sf::st_difference(boundary, layer)
    ocean_mask <- sf::st_cast(ocean_mask, "LINESTRING")
  }

  # if no overlaps found, close
  if(length(pointsOnLand_indexes)==0){
    cat("No land overlaps found\n")
    return(positions)
  }

  # subset points
  points <- methods::as(points,"SpatialPoints")
  pointsOnLand <- points[pointsOnLand_indexes]

  # get nearest location on water
  pointsCorrected <- list()
  distances <- c()
  pb <- txtProgressBar(min=0, max=length(pointsOnLand), initial=0, style=3)
  for (i in 1:length(pointsOnLand)) {
    point <- pointsOnLand[i]
    pointBoundingBox <- raster::as.matrix(raster::extent(point))
    pointBoundingBox[,1] <- pointBoundingBox[,1] - pointBoundingBox[,1]*0.05
    pointBoundingBox[,2] <- pointBoundingBox[,2] + pointBoundingBox[,2]*0.05
    pointBoundingBox <- raster::extent(pointBoundingBox)
    pointLandMask <- raster::crop(ocean_mask, pointBoundingBox)
    if(type=="raster"){
      distancesToLand <- geosphere::distHaversine(point, pointLandMask)
      distances[i] <- min(distancesToLand, na.rm=T)
      pointsCorrected[[i]] <- sp::coordinates(pointLandMask)[which.min(distancesToLand),]
    }
    if(type=="shape"){
      nearestPoint <- geosphere::dist2Line(point, line=pointLandMask, distfun=geosphere::distHaversine)
      pointsCorrected[[i]] <- nearestPoint[,c("lon","lat")]
      distances[i] <- nearestPoint[,"distance"]
    }
      setTxtProgressBar(pb, i)
  }
  pointsCorrected <- do.call("rbind", pointsCorrected)
  close(pb)

  if(plot==T){
    par(mgp=c(2,0.7,0))
    info <- paste0("Points relocated (n=", nrow(pointsCorrected), ")")
    raster::plot(layer, main=info, col="grey90", cex.main=0.8, cex.axis=0.8, las=1)
    points(positions[pointsOnLand_indexes, c(lon.col, lat.col)], pch=16, col=adjustcolor("red",alpha.f=0.5), cex=0.7)
    points(pointsCorrected, pch=16, col=adjustcolor("blue", alpha.f=0.5), cex=0.7)
  }

  # assign new locations and return
  cat(paste0("Points relocated: ", length(pointsOnLand), "\n"))
  cat(paste0("Min relocation distance: ", round(min(unlist(distances), na.rm=T)), " m\n"))
  cat(paste0("Max relocation distance: ", round(max(unlist(distances), na.rm=T)), " m\n"))
  positions[pointsOnLand_indexes, lon.col] <- pointsCorrected[,1]
  positions[pointsOnLand_indexes, lat.col] <- pointsCorrected[,2]
  return(positions)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
