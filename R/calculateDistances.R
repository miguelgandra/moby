#######################################################################################################
# Function to calculate shortest path in water between consecutive positions ##########################
#######################################################################################################

#' Estimate minimum traveled distances
#'
#' @description Function to calculate the shortest path in water between consecutive positions.
#' If no land.shape is provided, linear ('great-circle') distances are calculated instead.
#'
#' @inheritParams setDefaults
#' @param data A data frame with animal positions, containing longitude and latitude values (unprojected).
#' @param land.shape A  shape file containing coastlines.
#' @param epsg.code Coordinate reference system used to project positions (class 'CRS').
#' If not supplied, CRS is assumed to be the same as in land.shape.
#' @param grid.resolution Grid cell size in meters  over which shortest paths are going to be estimated.
#' @param mov.directions Number of directions allowed for shortest path calculation.
#' Passed to \link[gdistance]{transition}.
#' @param verbose Output process info and progress bar to console? Defaults to TRUE.
#' @return A data frame containing (minimum) distances traveled between each consecutive animal position,
#' calculated separately for each individual.
#' @export


calculateDistances <- function(data, land.shape=NULL, epsg.code=NULL, grid.resolution=100, mov.directions=16,
                               id.col=getDefaults("id"), lon.col=getDefaults("lon"), lat.col=getDefaults("lat"),
                               verbose=TRUE){


  ############################################################################
  ## Initial checks ##########################################################
  ############################################################################

  # measure running time
  start.time <- Sys.time()

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # check if dataset coordinates are in geographic format (unprojected)
  geographic_coords <-  all(data[,lon.col]>=(-180) & data[,lon.col]<=180 & data[,lat.col]>=(-90) & data[,lat.col]<=90)
  if(is.na(geographic_coords)){
    stop("Longitudes/latitudes values are missing and/or in the wrong format\n")
  }


  #############################################################################
  ## Prepare data for linear-distance calculation #############################
  ############################################################################

  # If no land.shape was provided
  if(is.null(land.shape)){
    if(verbose) {cat("No land shape provided\n")}
    if(geographic_coords==F){
      if(is.null(epsg.code)){
        stop("Please supply longitude and latitude in a geographic CRS /
           unprojected format (WGS84) or provide an epsg.code", call.=FALSE)
      }
      coords <- sp::SpatialPoints(cbind(data[,lon.col], data[,lat.col]))
      raster::projection(coords) <- epsg.code
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
        stop("Please supply a projected land.shape (in metres)", call.=FALSE)
      }else{
        epsg.code <- land.shape@proj4string
        if(verbose) {warning(paste0("Assuming CRS projection '", epsg.code), call.=FALSE)}
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
    pts <- sf::st_multipoint(cbind(data$lon_m_tmp, data$lat_m_tmp))
    pts <- sf::st_sfc(pts, crs=sf::st_crs(epsg.code))
    if(lengths(sf::st_intersects(pts, sf::st_as_sf(land.shape), sparse=T))>0){
      warning("Some of the coordinates overlap the supplied land shape", call.=FALSE)
    }

    # create transition layer
    if(verbose) {cat(paste0("Creating transition layer (", grid.resolution, "m grid | ", mov.directions, " directions)\n"))}
    projected_coords <- sf::st_multipoint(cbind(data$lon_m_tmp, data$lat_m_tmp))
    projected_coords <- sf::st_sfc(projected_coords, crs=sf::st_crs(epsg.code))
    template_raster <- raster::raster(raster::extent(sf::st_bbox(projected_coords))*1.2, res=grid.resolution, crs=epsg.code)
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

  # split COAs from different individuals
  data_individual <- split(data, f=data[,id.col], drop=F)
  ids <- levels(data[,id.col])

  # set progress bar
  if(verbose){
    pb <- txtProgressBar(min=1, max=length(data_individual), initial=0, style=3)
  }

  # initialize variables
  final_trajectories <- vector("list", length=length(data_individual))
  names(final_trajectories) <- names(data_individual)
  lines_skipped <- 0

  ##########################################################
  # if land shape was not provided calculate linear distances
  if(is.null(land.shape)){
    if(verbose) {cat("Calculating linear path between consecutive positions...\n")}
    for (i in 1:length(data_individual)) {
      # if individual doesn't have any detections, jump to next
      if(nrow(data_individual[[i]])<=0){
        if(verbose) {setTxtProgressBar(pb, i)}
        next
      }
      # else, if individual has only a single detection, return point geometry
      if(nrow(data_individual[[i]])==1){
        final_trajectories[[i]] <- NULL
        data_individual[[i]]$dist_m <- NA
        if(verbose) {setTxtProgressBar(pb, i)}
        next
      }
      # else, calculate great circle distances
      coords <- as.matrix(data_individual[[i]][,c("lon_wgs84_tmp","lat_wgs84_tmp")])
      dists <- geosphere::distVincentyEllipsoid(coords)
      data_individual[[i]]$dist_m <- c(dists, NA)
      final_trajectories[[i]] <- raster::spLines(coords, crs=epsg.code)
      if(verbose) {setTxtProgressBar(pb, i)}
    }
  }

  ##########################################################
  # else calculate shortest in-water paths (whenever land is intersected)
  if(!is.null(land.shape)){

    # iterate through each individual
    if(verbose) {cat("Calculating shortest in-water path between consecutive positions...\n")}
    for (i in 1:length(data_individual)) {

      # if individual doesn't have any detections, jump to next
      if(nrow(data_individual[[i]])==0){
        if(verbose) {setTxtProgressBar(pb, i)}
        next
      }

      # else, if individual has only a single detection, return point geometry
      if(nrow(data_individual[[i]])==1){
        final_trajectories[[i]] <- NULL
        data_individual[[i]]$dist_m <- NA
        if(verbose) {setTxtProgressBar(pb, i)}
        next
      }

      # else, prepare data and calculate linear ("great circle") distances
      coords <- as.matrix(data_individual[[i]][,c("lon_wgs84_tmp","lat_wgs84_tmp")])
      coords_m <- as.matrix(data_individual[[i]][,c("lon_m_tmp","lat_m_tmp")])
      dists <- geosphere::distVincentyEllipsoid(coords)
      trajectories <- list()

      # then, split the track into consecutive segments
      for(r in 1:(nrow(coords)-1)){
        segment <- sp::SpatialLines(list(sp::Lines(sp::Line(rbind(coords_m[r,], coords_m[r+1,])), ID="a")), proj4string=epsg.code)
        segment <- sf::st_as_sf(segment)
        is_pt <- all(coords_m[r,] == coords_m[r+1,])
        in_land <- lengths(sf::st_intersects(segment, sf::st_as_sf(land.shape), sparse=T))>0
        if(in_land==T & is_pt==F & floor(dists[r])<=grid.resolution){
          lines_skipped <- lines_skipped+1
        }
        # if segment overlaps land and is not a point, calculate shortest in-water path
        if(in_land==T & is_pt==F & floor(dists[r])>grid.resolution){
          segment <- gdistance::shortestPath(trCost, coords_m[r,], coords_m[r+1,], output="SpatialLines")
          segment <- sf::st_as_sf(segment)
          sf::st_crs(segment) <- sf::st_crs(epsg.code)
          segment_wgs84 <- sf::st_transform(segment, crs=sf::st_crs("+proj=longlat +datum=WGS84"))
          dists[r] <- sum(geosphere::distVincentyEllipsoid(sf::st_coordinates(segment_wgs84)[,1:2]))
        }
        # save segment
        trajectories[[r]] <- segment
      }

      # merge individual's results
      data_individual[[i]]$dist_m <- c(dists, NA)
      final_trajectories[[i]] <- do.call(rbind, trajectories)

      # update progress bar
      if(verbose) {setTxtProgressBar(pb, i)}
    }
  }


  ############################################################################
  ## Return results ##########################################################

  # close progress bar and print time taken
  if(verbose) {close(pb)}
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if(verbose)   cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n"))

  # print warning
  if(verbose & lines_skipped>0) {
  cat("Warning: ", lines_skipped, " track segment(s) overlapping land but with length <= than the defined grid resolution, original (straight) line(s) preserved\n")
  }

  # reassemble data
  data_individual <- data_individual[unlist(lapply(data_individual, nrow))>0]
  final_data <- do.call("rbind", data_individual)
  rownames(final_data) <- NULL

  # clean out helper columns
  out_cols <- c("hour_diff", "lon_wgs84_tmp", "lat_wgs84_tmp", "lon_m_tmp", "lat_m_tmp")
  final_data <- final_data[,-which(colnames(final_data) %in% out_cols)]

  # return distance results
  if(is.null(land.shape)){
    return(list("data"=final_data, "trajectories"=final_trajectories))}
  else{
    return(list("data"=final_data, "trajectories"=final_trajectories, "transition_layer"=trCost))}

}


#######################################################################################################
#######################################################################################################
#######################################################################################################
