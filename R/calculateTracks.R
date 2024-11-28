#######################################################################################################
# Function to calculate shortest path in water between consecutive positions ##########################
#######################################################################################################

#' Estimate minimum traveled distances
#'
#' @description This function calculates the shortest path in water between
#' consecutive positions for multiple individuals. If no land shapefile is provided,
#' the function defaults to calculating linear (great-circle) distances, which are
#' faster to compute. When a land shapefile is provided, the function can utilize
#' parallel computing to expedite least-cost path estimation, depending on the
#' number of CPU cores specified. However, depending on the chosen grid resolution
#' and the spatial extent of the provided positions, this process may take a long time
#' to run, even with parallel computing enabled.
#'
#' @inheritParams setDefaults
#' @param data A data frame with animal positions, containing longitude and latitude values.
#' @param land.shape Optional. A shapefile containing coastlines or landmasses. It can be supplied as
#' an 'sf' object or as an object of class 'SpatialPolygonsDataFrame' or 'SpatialPolygons'.
#' If the provided object is not of class 'sf', the function will attempt to
#' convert it to an 'sf' object for compatibility with subsequent spatial operations.
#' If not supplied, the function will calculate linear (great-circle) distances instead of shortest in-water
#' distances. Linear distances are much faster to compute but may be less realistic or accurate in cases
#' where positions are separated by landmasses or other spatial obstacles.
#' @param epsg.code The EPSG code (integer) representing the coordinate reference system (CRS) to be used
#' for projecting the positions. If not specified, the function will attempt to use the CRS from the
#' provided land.shape (if available).
#' @param grid.resolution The grid cell size (in meters) used to estimate shortest paths. A higher
#' resolution leads to more precise path calculations but may increase computation time.
#' @param mov.directions Number of movement directions allowed for shortest path calculations, passed
#' to the \link[gdistance]{transition} function from the `gdistance` package.
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations. This parameter
#' is only relevant for estimating least-cost paths, specifically when a \code{land.shape} is provided.
#' Run \code{parallel::detectCores()} to check the number of available cores.
#' @param transition.layer Optional. A pre-computed transition layer object, created using the
#' \link[gdistance]{transition} function from the `gdistance` package. If provided, it can
#' partially reduce computation time for repetitive function calls, as it avoids the need to
#' generate a new transition layer for each iteration (as internally implemented in some `moby` functions).
#' If not supplied, the function will create a transition layer based on the input \code{land.shape} and other relevant parameters.
#' @param verbose Logical. Should the function output process information and display a progress bar?
#' Defaults to TRUE.
#'
#' @note This function relies on the **gdistance** and **geosphere** packages for calculating
#' the shortest in-water paths and great-circle distances, respectively.
#'
#' @return A list with the following components:
#'   \item{data}{A data frame containing the calculated (minimum) distances between consecutive animal positions.}
#'   \item{trajectories}{A list of spatial objects representing the paths between positions.}
#'   \item{transition_layer}{(If land.shape is provided) The transition cost layer used for shortest path calculations.}
#'
#' @export


calculateTracks <- function(data,
                            land.shape = NULL,
                            epsg.code = getDefaults("epsg"),
                            grid.resolution = 100,
                            mov.directions = 16,
                            id.col = getDefaults("ID"),
                            lon.col = getDefaults("lon"),
                            lat.col = getDefaults("lat"),
                            cores = 1,
                            transition.layer = NULL,
                            verbose = TRUE){

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # capture the name of the 'land.shape' object
  if(!is.null(land.shape)){
    land_shape_name <- deparse(substitute(land.shape))
  # if no land shape is provided, set related parameters to NULL
  }else{
    land_shape_name <- NULL
    grid.resolution <- NULL
    mov.directions <- NULL
  }

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  land.shape <- reviewed_params$land.shape

  # manage spatial objects
  coords <- sf::st_as_sf(data, coords=c(lon.col, lat.col))
  coords_crs <- .checkProjection(coords)
  spatial_data <- .processSpatial(coords, land.shape, epsg.code)
  coords <- spatial_data$coords
  land.shape <- spatial_data$spatial.layer
  epsg.code <- spatial_data$epsg.code


  ##############################################################################
  ## Prepare data for linear-distance calculation ##############################
  ##############################################################################

  if(is.null(land.shape)){

    # convert coordinates to geographic format (EPSG:4326)
    if(coords_crs=="projected"){
      # transform to WGS84 (EPSG:4326)
      coords_wgs84 <- sf::st_transform(coords, crs=4326)
      # extract the transformed coordinates
      data$lon_wgs84_tmp <- sf::st_coordinates(coords_wgs84)[,1]
      data$lat_wgs84_tmp <- sf::st_coordinates(coords_wgs84)[,2]
    } else {
      epsg.code <- 4326
      data$lon_wgs84_tmp <- data[,lon.col]
      data$lat_wgs84_tmp <- data[,lat.col]
    }
  }

  ##############################################################################
  ## Prepare data for least-cost distance calculation ##########################
  ##############################################################################

   if(!is.null(land.shape)){

    if(coords_crs=="geographic"){
      # if the coordinates are already in geographic (WGS84)
      data$lon_wgs84_tmp <- data[,lon.col]
      data$lat_wgs84_tmp <- data[,lat.col]
      # extract the transformed coordinates
      data$lon_m_tmp <- sf::st_coordinates(coords)[,1]
      data$lat_m_tmp <- sf::st_coordinates(coords)[,2]
    } else {
      # if the coordinates are already projected (in meters)
      data$lon_m_tmp <- data[,lon.col]
      data$lat_m_tmp <- data[,lat.col]
      # transform back to WGS84 (EPSG:4326)
      coords_wgs84 <- sf::st_transform(coords, crs=4326)
      data$lon_wgs84_tmp <- sf::st_coordinates(coords_wgs84)[,1]
      data$lat_wgs84_tmp <- sf::st_coordinates(coords_wgs84)[,2]
    }

    # check if any coordinates overlap land
    pts <- sf::st_multipoint(cbind(data$lon_m_tmp, data$lat_m_tmp))
    pts <- sf::st_sfc(pts, crs=sf::st_crs(epsg.code))
    overlapping_pts <- lengths(sf::st_intersects(pts, sf::st_as_sf(land.shape), sparse=TRUE))
    if(overlapping_pts > 0){
      num_overlapping <- sum(overlapping_pts > 0)  # count of points overlapping land
      warning_str <- paste0("- Some coordinates (n=", num_overlapping,") overlap with the supplied land shape. ",
      "Consider using the 'correctPositions()' function to relocate these points to the nearest marine cell, ",
      "and then rerun the current function with the updated positions.")
      warning(paste(strwrap(warning_str, width=getOption("width")), collapse="\n"), call. = FALSE)
    }

    # validate transition.layer (if supplied)
    if(!is.null(transition.layer)) {
      if(!inherits(transition.layer, "TransitionLayer")) {
        stop("'transition.layer' must be a valid 'TransitionLayer' object from the gdistance package.", call.=FALSE)
      }
      if(verbose) cat("Using user-supplied transition layer.\n")
      trCost <- transition.layer
      grid.resolution <- raster::res(trCost)[1]
    # generate transition layer
    } else {
      if(verbose) {cat(paste0("Creating transition layer (", grid.resolution, "m grid | ", mov.directions, " directions)\n"))}
      template_raster <- raster::raster(raster::extent(sf::st_bbox(coords))*1.2, res=grid.resolution, crs=epsg.code)
      template_raster[] <- 0
      land.shape.sp <- sf::as_Spatial(land.shape)
      land_raster <- raster::rasterize(land.shape.sp, template_raster, update=TRUE)
      raster::values(land_raster)[raster::values(land_raster)==1] <- 10^8
      raster::values(land_raster)[raster::values(land_raster)==0] <- 1
      trCost <- gdistance::transition(1/land_raster, transitionFunction=mean, directions=mov.directions)
      trCost <- gdistance::geoCorrection(trCost, type="c")
    }

  }


  ##############################################################################
  ## Split data by individual and initialize variables #########################
  ##############################################################################

  # split COAs by individual
  data_individual <- split(data, f=data[,id.col], drop=FALSE)

  # initialize progress bar
  if(verbose){
    pb <- txtProgressBar(min=1, max=length(data_individual), initial=0, style=3)
  }

  # initialize trajectories list variable
  final_trajectories <- vector("list", length=length(data_individual))
  names(final_trajectories) <- names(data_individual)

  # initialize a variable to track the number of skipped line segments
  lines_skipped <- 0


  ##############################################################################
  ## A - Calculate great-circle distances (land.shape not provided) ############
  ##############################################################################

  if(is.null(land.shape)){

    # output to console
    if(verbose) cat("Calculating linear paths between consecutive positions...\n")

    # loop through each individual
    for (i in 1:length(data_individual)) {

      # (1) if individual doesn't have any detections, return empty
      if(nrow(data_individual[[i]])<=0){
        final_trajectories[[i]] <- list(NULL)

      # (2) else if individual has only a single detection, return point geometry
      }else if(nrow(data_individual[[i]])==1){
        data_individual[[i]]$dist_m <- NA
        final_trajectories[[i]] <- list(NULL)

      # (3) else, calculate great circle distances
      }else{
        coords <- as.matrix(data_individual[[i]][,c("lon_wgs84_tmp","lat_wgs84_tmp")])
        dists <- geosphere::distVincentyEllipsoid(coords)
        data_individual[[i]]$dist_m <- c(dists, NA)
        lines <- sf::st_linestring(coords)
        final_trajectories[[i]] <- sf::st_sf(geometry=sf::st_sfc(lines, crs=epsg.code))
      }

      # update progress bar
      if(verbose) setTxtProgressBar(pb, i)
    }
  }

  ##############################################################################
  ## B - Calculate least-cost paths (single core) ##############################
  ##############################################################################

  if(!is.null(land.shape) && cores==1){

    # output to console
    if(verbose) cat("Calculating least-cost paths between consecutive positions...\n")

    # loop through each individual
    for (i in 1:length(data_individual)) {

      # run function
      results <- calculateLeastCost(data_individual[[i]], land.shape, trCost, grid.resolution, epsg.code)

      # assign results
      data_individual[[i]] <- results$data
      final_trajectories[[i]] <- results$trajectories
      lines_skipped <- lines_skipped + results$lines_skipped

      # update progress bar
      if(verbose) setTxtProgressBar(pb, i)
    }
  }


  ##############################################################################
  ## C - Calculate least-cost paths (parallel computation) #####################
  ##############################################################################

  if(!is.null(land.shape) && cores>1){

    # print information to console if verbose mode is enabled
    if (verbose) cat(paste0("Starting parallel computation: ", cores, " cores\n"))
    if (verbose) cat("Calculating least-cost paths between consecutive positions...\n")

    # register parallel backend with the specified number of cores
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)

    # ensure the cluster is properly stopped when the function exits
    on.exit(parallel::stopCluster(cl))

    # define the `%dopar%` operator locally for parallel execution
    `%dopar%` <- foreach::`%dopar%`

    # set progress bar options based on verbose mode
    if(verbose) opts <- list(progress = function(n) setTxtProgressBar(pb, n))
    else opts <- NULL

    # perform parallel computation over each individual's data using foreach
    results_list <- foreach::foreach(i=1:length(data_individual), .options.snow=opts, .packages=c("sf", "gdistance", "geosphere")) %dopar% {
      calculateLeastCost(data_individual[[i]], land.shape, trCost, grid.resolution, epsg.code)
    }

    # aggregate results from parallel processing
    for (i in 1:length(results_list)) {
      data_individual[[i]] <- results_list[[i]]$data
      final_trajectories[[i]] <- results_list[[i]]$trajectories
      lines_skipped <- lines_skipped + results_list[[i]]$lines_skipped
    }

  }


  ##############################################################################
  ## Return results ############################################################
  ##############################################################################

  # close progress bar and print time taken
  if(verbose) close(pb)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if(verbose) cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n"))

  # issue a warning if any track segments overlapped land and were skipped
  if(verbose & lines_skipped>0) {
    warning(paste0("- ", lines_skipped,
                   " track segment(s) were found overlapping land but had a length <= to the defined grid resolution. ",
                   "These segments were retained as straight lines in the final output."), call. = FALSE)
  }

  # combine all individual data frames into one final data frame
  data_individual <- data_individual[unlist(lapply(data_individual, nrow))>0]
  final_data <- do.call("rbind", data_individual)
  rownames(final_data) <- NULL

  # clean out temporary helper columns that were used during processing
  out_cols <- c("lon_wgs84_tmp", "lat_wgs84_tmp", "lon_m_tmp", "lat_m_tmp")
  final_data <- final_data[,-which(colnames(final_data) %in% out_cols)]

  # aggregate results into a list, including data and trajectories
  if(is.null(land.shape)){
    results <- list("data"=final_data, "trajectories"=final_trajectories)
  }else{
    results <- list("data"=final_data, "trajectories"=final_trajectories, "transition_layer"=trCost)
  }

  # add relevant metadata as attributes to the results

  attr(results, 'land.shape') <- land_shape_name
  attr(results, 'epsg.code') <- epsg.code
  attr(results, 'grid.resolution') <- grid.resolution
  attr(results, 'mov.directions') <- mov.directions
  attr(results, 'cores') <- cores
  attr(results, 'processing.date') <- Sys.time()

  # return the final results (data, trajectories, and optional transition layer)
  return(results)
}



################################################################################
# Define least-cost path function ##############################################
################################################################################

#' Calculate the Least-Cost Path for Individual Detections
#'
#' This function computes the least-cost path for individual animal tracks based on
#' their detection coordinates. It handles cases where individuals may have
#' zero, one, or multiple detections. The function uses the geodesic distance
#' calculations and considers land overlaps to determine the shortest in-water path.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

calculateLeastCost <- function(data_individual, land.shape, trCost, grid.resolution, epsg.code){

  # Initialize variable
  lines_skipped <- 0

  # 1 - If individual doesn't have any detections, return NULL
  if(nrow(data_individual)==0){
    final_trajectories  <- list(NULL)

  # 2 - Else if individual has only a single detection, return point geometry
  }else if(nrow(data_individual)==1){
    data_individual$dist_m <- NA
    final_trajectories  <- list(NULL)

  # 3 - Else, calculate in-water distances
  }else{
    # prepare coordinate matrices for calculations
    coords <- as.matrix(data_individual[,c("lon_wgs84_tmp","lat_wgs84_tmp")])
    coords_m <- as.matrix(data_individual[,c("lon_m_tmp","lat_m_tmp")])
    # calculate initial geodesic distances between detection points
    dists <- geosphere::distVincentyEllipsoid(coords)
    # initialize list to store trajectory segments
    trajectories <- list()
    # loop through each consecutive pair of points in the track
    for(r in 1:(nrow(coords)-1)){
      # create a segment (line) between the current point and the next point
      segment <- sf::st_sfc(sf::st_linestring(rbind(coords_m[r, ], coords_m[r + 1, ])), crs=epsg.code)
      # check if the current segment is a point (i.e., the same start and end)
      is_pt <- all(coords_m[r,] == coords_m[r+1,])
      # check if the segment intersects with land
      in_land <- lengths(sf::st_intersects(segment, sf::st_as_sf(land.shape), sparse=TRUE))>0
      # if the segment overlaps land and is not a point, and the distance is less than the grid resolution
      if(in_land & !is_pt & floor(dists[r])<=grid.resolution){
        lines_skipped <- lines_skipped+1
      }
      # if the segment overlaps land and is not a point, calculate shortest in-water path
      if(in_land & !is_pt & floor(dists[r])>grid.resolution){
        # calculate the shortest path around land
        segment <- gdistance::shortestPath(trCost, coords_m[r,], coords_m[r+1,], output="SpatialLines")
        # convert the segment to a spatial object and set its CRS
        segment <- sf::st_as_sf(segment)
        sf::st_crs(segment) <- sf::st_crs(epsg.code)
        # transform segment coordinates to WGS84 for distance calculation
        segment_wgs84 <- sf::st_transform(segment, crs=sf::st_crs("+proj=longlat +datum=WGS84"))
        # update the distance with the calculated distance along the segment
        dists[r] <- sum(geosphere::distVincentyEllipsoid(sf::st_coordinates(segment_wgs84)[,1:2]))
      }
      # save the segment to the trajectories list
      trajectories[[r]] <- segment
    }

    # merge the distances into the individual data frame
    data_individual$dist_m <- c(dists, NA)

    # Convert segments into a simple feature geometry collection
    final_trajectories <- lapply(trajectories, sf::st_geometry)
    final_trajectories <- do.call(c, final_trajectories)
  }

  # return the updated individual data and calculated trajectories
  return(list("data"=data_individual, "trajectories"=final_trajectories, "lines_skipped"=lines_skipped))

}

#######################################################################################################
#######################################################################################################
#######################################################################################################
