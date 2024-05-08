#######################################################################################################
# Function to calculate shortest path in water between consecutive positions ##########################
#######################################################################################################

#' Estimate minimum traveled distances between all positions (pairwise)
#'
#' @description Function to calculate the shortest path in water between consecutive positions.
#' If no land.shape is provided, linear ('great-circle') distances are calculated instead.
#' @param data A data frame with animal positions (COAs), containing longitude and latitude values
#' @param land.shape A  shape file containing coastlines.
#' @param epsg.code Coordinate reference system used to project positions (class 'CRS').
#' If not supplied, CRS is assumed to be the same as in land.shape.
#' @param grid.resolution Grid cell size in meters  over which shortest paths are going to be estimated.
#' @param mov.directions Number of directions allowed for shortest path calculation.
#' Passed to \link[gdistance]{transition}.

#' @export


calculatePairDistances <- function(data, land.shape, epsg.code=NULL, grid.resolution=1000, mov.directions=8,
                                   lon.col="lon", lat.col="lat", verbose=FALSE){


  ############################################################################
  ## Initial checks ##########################################################
  ############################################################################

  # measure running time
  start.time <- Sys.time()


  # check if data contains lon.col and lat.col
  if(any(!c(lon.col, lat.col) %in% colnames(data))){
    stop("Longitude/latitude columns not found. Please assign the correct column names with 'lon.col' and 'lat.col'")
  }

  # check if dataset coordinates are in geographic format (unprojected)
  geographic_coords <-  all(data[,lon.col]>=(-180) & data[,lon.col]<=180 & data[,lat.col]>=(-90) & data[,lat.col]<=90)

  if(is.na(geographic_coords)){stop("Missing longitudes/latitudes\n")}


  #############################################################################
  ## Prepare data for linear-distance calculation #############################
  ############################################################################

  # If no land.shape was provided
  if(is.null(land.shape)){
    if(verbose) {cat("No land shape provided\n")}
    if(geographic_coords==F){
      if(is.null(epsg.code)){
        stop("Please supply longitude and latitude in a geographic CRS /
           unprojected format (WGS84) or provide an epsg.code")
      }
      coords <- sp::SpatialPoints(cbind(data[,lon.col], data[,lat.col]))
      raster::projection(coords) <- epsg_code
      coords <- sp::spTransform(coords, sp::CRS("+proj=longlat +datum=WGS84"))
      data$lon_wgs84_tmp <- coords@coords[,1]
      data$lat_wgs84_tmp <- coords@coords[,2]
    } else {
      epsg.code <- sp::CRS("+proj=longlat +datum=WGS84")
      data$lon_wgs84_tmp <- data[,lon.col]
      data$lat_wgs84_tmp <- data[,lat.col]
    }
  }

  ############################################################################
  ## Prepare data for least-cost distance calculation ########################
  ############################################################################

  # If a land.shape was provided
  if(!is.null(land.shape)){
    # retrieve epsg.code if not provided
    if(is.null(epsg.code)){
      if(!grepl("+units=m", land.shape@proj4string, fixed=T)){
        stop("Please supply a projected land.shape (in metres)")
      }else{
        epsg.code <- land.shape@proj4string
        if(verbose) {warning(paste0("Assuming CRS projection '", epsg.code, "'\n"))}
      }
    } else {
      land.shape <- sp::spTransform(land.shape, epsg.code)
    }

    # convert coordinates
    if(geographic_coords==T){
      data$lon_wgs84_tmp <- data[,lon.col]
      data$lat_wgs84_tmp <- data[,lat.col]
      coords <- sp::SpatialPoints(cbind(data$lon_wgs84_tmp, data$lat_wgs84_tmp))
      raster::projection(coords) <- sp::CRS("+proj=longlat +datum=WGS84")
      coords <- sp::spTransform(coords, epsg.code)
      data$lon_m_tmp <- coords@coords[,1]
      data$lat_m_tmp <- coords@coords[,2]
    } else {
      data$lon_m_tmp <- data[,lon.col]
      data$lat_m_tmp <- data[,lat.col]
      coords <- sp::SpatialPoints(cbind(data$lon_m_tmp, data$lat_m_tmp), epsg.code)
      coords <- sp::spTransform(coords, sp::CRS("+proj=longlat +datum=WGS84"))
      data$lon_wgs84_tmp <- coords@coords[,1]
      data$lat_wgs84_tmp <- coords@coords[,2]
    }

    # check if any coordinate is on land
    pts <- sp::SpatialPoints(as.matrix(data[,c("lon_m_tmp","lat_m_tmp")]), proj4string=epsg.code)

    # create transition layer
    if(verbose) {cat(paste0("Creating transition layer (", grid.resolution, "m grid | ", mov.directions, " directions)\n"))}
    projected_coords <- sp::SpatialPoints(cbind(data$lon_m_tmp, data$lat_m_tmp))
    raster::projection(projected_coords) <- epsg.code
    template_raster <- raster::raster(extent(projected_coords)*1.2, res=grid.resolution, crs=epsg.code)
    template_raster[] <- 0
    land_raster <- raster::rasterize(land.shape, template_raster, update=T)
    raster::values(land_raster)[raster::values(land_raster)==1] <- 10^8
    raster::values(land_raster)[raster::values(land_raster)==0] <- 1
    trCost <- gdistance::transition(1/land_raster, transitionFunction=mean, directions=mov.directions)
    trCost <- gdistance::geoCorrection(trCost, type="c")
  }


  ############################################################################
  ## Calculate trajectories and estimate distances ###########################
  ############################################################################

  # get all pairwise combinations
  pairs <- t(combn(nrow(data), m=2))
  npairs <- nrow(pairs)

  # first calculate all linear distances
  dist_matrix <- distm(data[,c("lon_wgs84_tmp", "lat_wgs84_tmp")], fun=distHaversine)
  dist_matrix[upper.tri(dist_matrix, diag=T)] <- NA

  # then, identify which pair segments intersect land and
  # calculate shortest in-water distances

  # set progress bar
  if(verbose){
    cat("Identifying segments overlapping land\n")
    pb <- txtProgressBar(min=1, max=npairs, initial=0, style=3)
  }
  land_overlaps <- matrix(NA, nrow=nrow(dist_matrix), ncol=ncol(dist_matrix))
  for(i in 1:nrow(pairs)){
    index1 <- pairs[i,1]
    index2 <- pairs[i,2]
    p1 <- as.matrix(data[index1, c("lon_m_tmp", "lat_m_tmp")])
    p2 <- as.matrix(data[index2, c("lon_m_tmp", "lat_m_tmp")])
    segment <- raster::spLines(rbind(p1, p2), crs=epsg.code)
    is_pt <- all(p1 == p2)
    in_land <- rgeos::gIntersects(segment, land.shape)
    linear_dist <- dist_matrix[lower.tri(dist_matrix)][i]
    if(in_land==T & is_pt==F & linear_dist>grid.resolution){
      segment <- gdistance::shortestPath(trCost, p1, p2, output="SpatialLines")
      segment@proj4string <- epsg.code
      segment_wgs84 <- sp::spTransform(segment, sp::CRS("+proj=longlat +datum=WGS84"))
      dist_matrix[lower.tri(dist_matrix)][i] <- sum(geosphere::distVincentyEllipsoid(segment_wgs84@lines[[1]]@Lines[[1]]@coords))
    }
    # update progress bar
    if(verbose) {setTxtProgressBar(pb, i)}
  }
  if(verbose) {close(pb)}


  # return pairwise distance matrix
  dist_matrix[is.infinite(dist_matrix)] <- NA
  return(dist_matrix)

}


#######################################################################################################
#######################################################################################################
#######################################################################################################
