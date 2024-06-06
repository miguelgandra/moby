#######################################################################################################
# Function to correct positions on land and relocate them to the nearest marine cell ##################
#######################################################################################################

#' Relocate points on land to the nearest marine cell.
#'
#' Either a coastline shapefile or a raster containing land surfaces or bathymetry values
#' can be provided.
#'
#' @description Function to relocate positions on land to nearest marine cell
#'
#' @param positions A data frame containing longitude and latitude values
#' (unprojected lat-long coordinates), corresponding to the inferred animals' positions.
#' @param layer A raster containing land portions. Or alternatively a raster containing
#' elevation values (negative below water).
#' @param layer.type Either 'land', 'ocean' or 'bathy'. If type 'land' (binary) the raster is assumed
#' to contain NA values in water surfaces, if type 'ocean' (binary) the raster is assumed to
#' contain NA values in land surfaces, else if type 'bathy' the raster is assumed to
#' range between 0 (land) to max depth (positive or negative values).
#' @param depth.margin If set, shortest in-water paths will be calculated only
#' for cells with depths >= than this threshold. Only takes effect when raster
#' is of type bathy. Defaults to 0 (water vs land).
#' @param dist.margin If set, shortest in-water paths will be calculated only
#' for cells with depths <= than this threshold. Only takes effect when raster
#' is of type bathy. Disabled by default.
#' @param lon.col Name of the column containing longitude values. Defaults to 'lon'.
#' @param lat.col Name of the column containing latitude values. Defaults to 'lat'.
#' @export

correctPositions <- function(positions, layer, layer.type="land", depth.margin=0,
                             lon.col="lon", lat.col="lat", plot=F) {

  ############################################################################
  ## Initial checks ##########################################################
  ############################################################################

  # measure running time
  start.time <- Sys.time()

  # check layer type
  if(class(layer)=="RasterLayer"){
    type <- "raster"
    proj <- layer@crs
  }else if(grepl("SpatialPolygon", class(layer))){
    type <- "shape"
    proj <- layer@proj4string
  }else{
    stop("Layer: Please supply either a Raster or a SpatialPolygon object")
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
  points <- sp::SpatialPoints(positions[,c(lon.col, lat.col)], proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
  points <- sp::spTransform(points, proj)
  if(type=="raster"){
    point_values <- raster::extract(layer, points)
    pointsOnLand_indexes <- which(point_values==1)
    ocean_mask <- raster::rasterToPoints(layer, fun=function(x){x==0}, spatial=T)
  }else if(type=="shape"){
    pointsOnLand_indexes <- which(!is.na(sp::over(points, layer)))
    boundary <- sp::Polygon(raster::extent(layer))
    boundary <- sp::SpatialPolygons(list(sp::Polygons(list(boundary), ID=1)))
    boundary@proj4string <- proj
    ocean_mask <- rgeos::gDifference(boundary, layer)
    ocean_mask <- as(ocean_mask, "SpatialLines")
  }

  # if no overlaps found, close
  if(length(pointsOnLand_indexes)==0){
    cat("No land overlaps found\n")
    return(positions)
  }

  # subset points
  points <- as(points,"SpatialPoints")
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
