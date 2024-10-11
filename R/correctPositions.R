#######################################################################################################
# Function to correct positions on land and relocate them to the nearest marine cell ##################
#######################################################################################################

#' Relocate points on land to the nearest marine cell.
#'
#' @description This function relocates positions on land to the nearest marine cell
#' using either a coastline shapefile or a raster containing land surfaces or bathymetry values.
#'
#' @inheritParams setDefaults
#' @param data A data frame with animal positions, containing longitude and latitude values.
#' @param spatial.layer A spatial object used to determine land and marine areas. This can be one of the following:
#' - A `RasterLayer` with binary values or bathymetric data.
#' - An `sf` or `SpatialPolygons` object representing coastlines or landmasses.
#' @param raster.type A character string indicating the type of raster when `spatial.layer` is a `Raster`.
#' This parameter can be either 'land', 'ocean', or 'bathy'.
#' \itemize{
#'   \item If 'land', the raster is assumed to contain NA values in water surfaces.
#'   \item If 'water', the raster is assumed to contain NA values in land surfaces.
#'   \item If 'bathy', the raster is assumed to range between 0 (land) to maximum depth (positive or negative values).
#' }
#' @param depth.threshold A numeric value. If set, shortest in-water paths will be calculated only
#' for cells with depths >= this threshold. Only takes effect when raster
#' is of type 'bathy'. Defaults to 0 (water vs land).
#' @param max.distance.km A numeric value specifying the maximum distance (in kilometers) to consider when relocating points.
#' This parameter limits the search radius for the nearest marine cell, ensuring that only cells within the specified distance are evaluated.
#' Points that are further than this distance from the nearest water will have their coordinates set to NA.
#' @param plot A logical value indicating whether to generate a plot of the relocation process.
#' If set to TRUE, the function will plot the `spatial.layer` along with the original and updated positions.
#' If set to FALSE, no plot will be generated.
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations.
#' Run \code{parallel::detectCores()} to check the number of available cores.

#' @return A list with two elements:
#' \describe{
#'   \item{data}{The original data frame with updated positions for the points that were relocated from land.}
#'   \item{summary}{A summary of the input data and changes made to the positions.}
#' }
#'
#' @export

correctPositions <- function(data,
                             spatial.layer,
                             raster.type = "land",
                             epsg.code = getDefaults("epsg"),
                             lon.col = getDefaults("lon"),
                             lat.col=getDefaults("lat"),
                             depth.threshold = 0,
                             max.distance.km = 50,
                             plot = FALSE,
                             cores = 1) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # print message
  .printConsole("Relocating positions on land to the nearest marine cell")

  # capture the original spatial.layer name
  spatial_layer_name <- deparse(substitute(spatial.layer))

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # validate additional parameters
  errors <- c()
  if (!inherits(spatial.layer, c("RasterLayer", "sf", "SpatialPolygons"))) errors <- c(errors, "The 'spatial.layer' argument must be a 'RasterLayer', 'sf', or 'SpatialPolygons' object.")
  if (!raster.type %in% c("land", "water", "bathy")) errors <- c(errors, "The 'raster.type' must be either 'land', 'water' or 'bathy'")
  if (!inherits(depth.threshold, "numeric")) errors <- c(errors, "The 'depth.threshold' argument must be numeric")
  if (!inherits(plot, "logical"))  errors <- c(errors, "The 'plot' argument must be a logical value (TRUE or FALSE).")
  if(length(errors)>0){
      stop_message <- sapply(errors, function(x) paste(strwrap(x), collapse="\n"))
      stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
      stop(stop_message, call.=FALSE)
  }

  # convert spatial.layer to sf if required
  if(inherits(spatial.layer, "SpatialPolygons")){
    spatial.layer <- sf::st_as_sf(spatial.layer)

  # convert all NAs in the raster to a temporary dummy value for processing
  }else if(inherits(spatial.layer, "Raster")){
    raster::values(spatial.layer)[is.na(raster::values(spatial.layer))] <- 999999
  }

  # save the original spatial.layer extent
  spatial_layer_bbox <- sf::st_bbox(spatial.layer)

  # manage spatial objects
  coords <- sf::st_as_sf(data, coords=c(lon.col, lat.col))
  coords_crs <- .checkProjection(coords)
  spatial_data <- .processSpatial(coords, spatial.layer, epsg.code)
  coords <- spatial_data$coords
  spatial.layer <- spatial_data$spatial.layer
  epsg.code <- spatial_data$epsg.code


  ##############################################################################
  ## Trim Raster ###############################################################
  ##############################################################################

  # remove potential outer NA values from raster after reprojection
  if(inherits(spatial.layer, "Raster")){
     # trim the raster to remove outer rows and columns that contain NA values
     spatial.layer <- raster::trim(spatial.layer)
     # check if there are any remaining NA values in the trimmed raster
     if(any(is.na(raster::values(spatial.layer)))){
       # initialize a variable to track whether NA values remain
       na_remaining <- TRUE
       # iterate until there are no NA values left in the raster
       while(na_remaining) {
         # fill the NA cells based on nearest neighbors (3x3 window)
         spatial.layer <- raster::focal(spatial.layer, w=matrix(1,3,3), fun=raster::modal, na.rm=TRUE, NAonly=TRUE, pad=99999)
         # check if there are still NA values remaining
         na_remaining <- any(is.na(raster::getValues(spatial.layer)))
       }
     }
     # convert all temporary dummy values back to NA
     raster::values(spatial.layer)[raster::values(spatial.layer)==999999] <- NA
   }

  ##############################################################################
  ## Process Raster Values #####################################################
  ##############################################################################

  # convert values based on raster.type
  if (inherits(spatial.layer, "Raster")) {
    raster_layer <- spatial.layer
    if (raster.type == "land") {
      if(!any(is.na(raster::values(spatial.layer)))) stop("The spatial layer contains no NA values. Ensure the layer includes water surfaces as NA for proper processing.", call.=FALSE)
      cat(paste(strwrap(paste("Raster type set to 'land' - all values will be interpreted as representing land surfaces, while NA values will be treated as water."), width=getOption("width")), collapse="\n"), "\n")
      raster::values(raster_layer)[!is.na(raster::values(raster_layer))] <- 1
      raster::values(raster_layer)[is.na(raster::values(raster_layer))] <- 0
    } else if (raster.type == "water") {
      if(!any(is.na(raster::values(spatial.layer)))) stop("The spatial layer contains no NA values. Ensure the layer includes land surfaces as NA for proper processing.", call.=FALSE)
      cat(paste(strwrap(paste("Raster type set to 'water' - all values will be interpreted as representing water surfaces, while NA values will be treated as land"), width=getOption("width")), collapse="\n"), "\n")
      raster::values(raster_layer)[!is.na(raster::values(raster_layer))] <- 0
      raster::values(raster_layer)[is.na(raster::values(raster_layer))] <- 1
    } else if (raster.type == "bathy") {
      if(length(unique(raster::values(spatial.layer)))<=3) stop("Raster type set to 'bathy', but the supplied values seem binary.", call.=FALSE)
      bathy_vals <- ifelse(mean(raster::values(raster_layer) > 0, na.rm = TRUE) > 0.5, "positive", "negative")
      if(bathy_vals=="positive") {
        cat(paste(strwrap(paste("Depth values seem to range from 0 to", round(max(raster::values(raster_layer), na.rm=TRUE)),
                                "m. Negative values will be converted to NA (land)."), width=getOption("width")), collapse="\n"), "\n")
      }else{
        cat(paste(strwrap(paste("Depth values seem to range from 0 to", round(min(raster::values(raster_layer), na.rm=TRUE)),
                                "m. Positive values will be converted to NA (land)."), width=getOption("width")), collapse="\n"), "\n")
        raster_layer <- raster_layer * -1
      }
      raster::values(raster_layer)[raster::values(raster_layer) <= depth.threshold] <- 1
      raster::values(raster_layer)[raster::values(raster_layer) > depth.threshold] <- 0
      raster::values(raster_layer)[is.na(raster::values(raster_layer))] <- 1
      if(depth.threshold>0) cat(paste0("Positions will be relocated to a min depth of ", depth.threshold, " meters.\n"))
    }
  }


  ############################################################################
  ## Identify Points on Land #################################################
  ############################################################################

  # identify points that are on land based on the spatial layer
  if(inherits(spatial.layer, "Raster")){
    point_values <- raster::extract(raster_layer, coords)
    pointsOnLand_indexes <- which(point_values==1)
    land_mask <- raster::rasterToPoints(raster_layer, fun=function(x){x==0}, spatial=TRUE)
    land_mask <- sf::st_as_sf(land_mask)
    land_mask <- sf::st_transform(land_mask, crs=epsg.code)
  }else if(inherits(spatial.layer, "sf")){
    spatial.layer <- sf::st_union(spatial.layer)
    pointsOnLand_indexes <- which(sf::st_intersects(coords, spatial.layer, sparse=FALSE))
    boundary <- sf::st_as_sfc(sf::st_bbox(spatial.layer))
    boundary <- sf::st_set_crs(boundary, epsg.code)
    land_mask <- sf::st_difference(boundary, spatial.layer)
    land_mask <- sf::st_cast(land_mask, "LINESTRING")
  }

  # if no land overlaps are found, return the original data
  if(length(pointsOnLand_indexes)==0){
    cat("No land overlaps were detected. Returning the original data.\n")

    # return the original dataset
    results <- list("data"=data)

    # add attributes to the results list to save relevant parameters
    attr(results, 'points.relocated') <- length(pointsOnLand_indexes)
    attr(results, 'points.skipped') <- 0
    attr(results, 'spatial.layer.name') <- spatial_layer_name
    attr(results, 'spatial.layer.bbox') <- spatial_layer_bbox
    attr(results, 'depth.threshold') <- depth.threshold
    attr(results, 'max.distance.km') <- max.distance.km
    attr(results, 'epsg.code') <- epsg.code
    attr(results, 'processing.date') <- Sys.time()

    # return results
    return(results)
  }


  ##############################################################################
  ## Correct Positions #########################################################
  ##############################################################################

  # retrieve spurious positions (on land)
  pointsOnLand <- coords[pointsOnLand_indexes, ]

  # initialize list to store results
  results <- vector("list", nrow(pointsOnLand))

  # initialize progress bar
  pb <- txtProgressBar(min=1, max=nrow(pointsOnLand), initial=0, style=3)


  ###########################################################
  # single core ############################################
  if(cores==1){

    # iterate over each spurious point to find the nearest in-water location
    for (i in 1:nrow(pointsOnLand)) {
      results[[i]] <- .relocatePoint(i, pointsOnLand, land_mask, max.distance.km)
      setTxtProgressBar(pb, i)
    }


  ###########################################################
  # parallel computation
  }else if(cores>1){

      # print information to console
      cat(paste0("Starting parallel computation: ", cores, " cores\n"))

      # register parallel backend with the specified number of cores
      cl <- parallel::makeCluster(cores)
      doSNOW::registerDoSNOW(cl)

      # ensure the parallel cluster is stopped when the function exits
      on.exit(parallel::stopCluster(cl))

      # define the `%dopar%` operator locally for parallel execution
      `%dopar%` <- foreach::`%dopar%`

      # set progress bar options based on verbose mode
      opts <- list(progress = function(n) setTxtProgressBar(pb, n))

      # perform parallel computation over each individual's data using foreach
      results <- foreach::foreach(i=1:nrow(pointsOnLand), .options.snow=opts, .packages=c("sf"),
                                  .export=c(".relocatePoint", ".setSearchRegion")) %dopar% {
        .relocatePoint(i, pointsOnLand, land_mask, max.distance.km)
      }
  }


  # close progress bar
  close(pb)

  # separate the results into pointsCorrected and distances
  pointsCorrected <- lapply(results, `[[`, "coords")
  pointsCorrected <- do.call("rbind", pointsCorrected)
  distances <- unlist(lapply(results, `[[`, "dist"))

  # identify relocated and skipped indices
  relocated_indices <- which(!is.na(pointsCorrected[,1]) & !is.na(pointsCorrected[,2]))
  skipped_indices <- which(is.na(pointsCorrected[,1]) | is.na(pointsCorrected[,2]))


  ##############################################################################
  ## Generate plot (if required) ###############################################
  ##############################################################################

  # check if the plot should be generated
  if(plot){

    # set up a layout for two plots side by side: original positions and the legend
    layout(matrix(c(1,2), nrow=1), widths=c(3, 1.5))

    # adjust plotting parameters for margins and axis label positioning
    par(mar=c(2,2,2,0))

    # create a bounding box for all points
    plot_points <- rbind(sf::st_coordinates(coords)[pointsOnLand_indexes, 1:2], pointsCorrected[relocated_indices,])
    plot_points <- sf::st_as_sf(as.data.frame(plot_points), coords=c("X", "Y"))
    plot_bbox <- sf::st_bbox(plot_points)

    # expand bounding box by 5% in all directions
    dx <- plot_bbox["xmax"]-plot_bbox["xmin"]
    dy <- plot_bbox["ymax"]-plot_bbox["ymin"]
    plot_bbox <- plot_bbox + c(dx, dy, dx, dy) * c(-0.05, -0.05, 0.05, 0.05)

    ###########################################################
    # crop and plot the raster layer ##########################
    if(inherits(spatial.layer, "Raster")){
      # crop the raster to the expanded bounding box
      spatial.layer <- raster::crop(spatial.layer, raster::extent(plot_bbox))
      # adjust raster values if using bathymetry (depth) and invert negative values if specified
      if(raster.type=="bathy" && bathy_vals=="negative") spatial.layer <- spatial.layer*(-1)
      raster::values(spatial.layer)[raster::values(spatial.layer)<0] <- 0
      if(raster.type=="land") spatial.layer <- raster::calc(spatial.layer, fun=function(x) ifelse(is.na(x), 1, NA))
      # create a mask for NA values in the raster
      na_mask <- is.na(spatial.layer)
      # convert the NA mask into polygons for contouring
      na_polygons <- raster::rasterToPolygons(na_mask, fun=function(x) x==1, dissolve=TRUE)
      # set up the color palette
      if(raster.type=="bathy"){
        color_pal <- rev(.bathy_deep_pal(100)[c(10:95)])
        color_pal <- colorRampPalette(color_pal)(100)
      }else{
        color_pal <- "grey96"
      }
      # plot the raster image
      raster::image(spatial.layer, col=color_pal, main="Results",
                    axes=FALSE, xlab="", ylab="", cex.main=0.8, cex.axis=0.8, las=1, asp=1)
      # overlay the NA areas with a grey color for better visual distinction
      raster::image(na_mask, col=c(NA, "grey"), add=TRUE, legend=FALSE)
      # add a contour around the NA areas
      raster::plot(na_polygons, border="black", lwd=0.2, add=TRUE)
      # add a contour at the depth threshold, if specified,
      if(depth.threshold>0) raster::contour(spatial.layer, levels=depth.threshold, col="black", lwd=0.6, add=TRUE)

    ##########################################################
    # crop and plot the sf layer #############################
    }else{
      spatial.layer <- sf::st_crop(spatial.layer, plot_bbox)
      plot(spatial.layer, main="Results", col="grey90", cex.main=0.8, cex.axis=0.8, las=1)
    }

    ##########################################################
    # plot points ############################################
    pt_colors <- adjustcolor(c("red3","coral","blue"), alpha.f=0.5)
    if(length(skipped_indices)>0) points(coords[pointsOnLand_indexes[skipped_indices],], pch=16, col=pt_colors[1], cex=0.7)
    points(coords[pointsOnLand_indexes[relocated_indices],], pch=16, col=pt_colors[2], cex=0.7)
    points(pointsCorrected[relocated_indices,], pch=16, col=pt_colors[3], cex=0.7)

    ##########################################################
    # draw box ###############################################
    box()

    ##########################################################
    # create an empty plot area for the legend ###############
    par(mar=c(2,0.4,2,0))
    plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE, main="", xaxs="i", yaxs="i")

    ##########################################################
    # add the legend     #####################################
    legend_string <- paste0("Land positions (", nrow(pointsOnLand),")")
    legend_string <- c(legend_string, paste0("Relocated positions (", length(relocated_indices),")"))
    if(length(skipped_indices)>0) legend_string <- c(legend_string, paste0("Skipped positions (", length(skipped_indices),")"))
    legend("left", legend=legend_string, pt.cex=0.8, y.intersp=1.2,
           pch=16, col=pt_colors[c(2,3,1)], bty="n", cex=0.7, xpd=TRUE)

    # reset layout
    par(mar=c(5,4,4,2)+0.1)
    layout(1)

  }


  ##############################################################################
  ## Return results ############################################################
  ##############################################################################

  # output final summary of the relocation process
  cat(paste0("Points relocated: ", length(relocated_indices), "\n"))
  if(length(skipped_indices)>0) cat(paste0("Points skipped: ", length(skipped_indices), "\n"))
  mean_distance <- suppressWarnings(round(mean(unlist(distances), na.rm=T)))
  min_distance <- suppressWarnings(round(min(unlist(distances), na.rm=T)))
  max_distance <- suppressWarnings(round(max(unlist(distances), na.rm=T)))
  if(!is.infinite(mean_distance) && !is.na(mean_distance)){
    cat(paste0("Mean relocation distance: ", mean_distance, " m (from ", min_distance," m to ", max_distance," m)\n"))
  }

  # print time taken
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n"))

  # display the warning message if any points were skipped during processing
  if(length(skipped_indices)>0){
    num_skipped <- length(skipped_indices)
    point_word <- ifelse(num_skipped == 1, "point", "points")
    this_or_these <- ifelse(num_skipped == 1, "this point", "these points")
    warning_msg <- paste0("Maximum search radius of ", max.distance.km, " km reached without finding water areas for ",
                         num_skipped, " ", point_word, ". Coordinates for ", this_or_these, " have been set to NA.")
    warning(paste(strwrap(warning_msg), collapse="\n"), call.=FALSE)
  }


  # convert coordinates back to geographic (if needed)
  if(coords_crs=="geographic" && length(relocated_indices)>0){
    warning("Coordinates were initially in a geographic CRS. They have been projected for spatial processing and converted back to geographic coordinates.", call.=FALSE)
    pointsProjected <- sf::st_as_sf(as.data.frame(pointsCorrected[relocated_indices,]), coords=c(1,2), crs=epsg.code)
    pointsProjected <- sf::st_transform(pointsProjected, crs = 4326)
    pointsProjected <- sf::st_coordinates(pointsProjected)
    pointsCorrected[relocated_indices,] <- pointsProjected
  }

  # create a summary table with relocated points, their original and new coordinates, and the relocation distances
  summary <- data.frame("index"=pointsOnLand_indexes,
                        "lon"=data[pointsOnLand_indexes, lon.col],
                        "lat"=data[pointsOnLand_indexes, lat.col],
                        "new.Lon"=pointsCorrected[,1],
                        "new.Lat"=pointsCorrected[,2],
                        "distance_m"=round(distances,1))
  # rename columns to reflect original and new coordinates
  colnames(summary)[2:3] <- paste0("original.", c(lon.col, lat.col))
  colnames(summary)[4:5] <- paste0("new.", c(lon.col, lat.col))

  # replace original coordinates with corrected positions
  data[pointsOnLand_indexes, lon.col] <- pointsCorrected[,1]
  data[pointsOnLand_indexes, lat.col] <- pointsCorrected[,2]

  # prepare the final results to return, including the updated dataset and the summary table
  results <- list("data"=data, "summary"=summary)

  # add attributes to the results list to save relevant parameters
  attr(results, 'points.relocated') <- length(relocated_indices)
  attr(results, 'points.skipped') <- length(skipped_indices)
  attr(results, 'spatial.layer.name') <- spatial_layer_name
  attr(results, 'spatial.layer.bbox') <- spatial_layer_bbox
  attr(results, 'depth.threshold') <- depth.threshold
  attr(results, 'max.distance.km') <- max.distance.km
  attr(results, 'epsg.code') <- epsg.code
  attr(results, 'processing.date') <- Sys.time()

  # return the updated dataset and summary
  return(results)
}


################################################################################
# Define search region function ################################################
################################################################################

#' Internal function to define a search region for distance calculations
#'
#' This function helps speed up distance calculations by incrementally expanding the search radius
#' around a given position to find the nearest water body. It stops when a water cell is detected
#' or when the search radius exceeds 50 km. The function supports both `RasterLayer` and `sf` objects as spatial inputs.
#'
#' @param point An `sf` object representing the coordinates (longitude, latitude)
#' of the point of interest.
#' @param land.mask An `sf` object that represents the spatial layer
#' containing information about land and water. A value of 0 indicates water.
#' @return A cropped spatial layer of class `sf`, centered around the identified water region.
#' If no water is found within the maximum search radius of 50 kilometers, an error is raised.
#'
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.setSearchRegion <- function(point, land.mask, max.distance.km){

  # initialize variables
  water_found <- FALSE

  # define initial search radius and step size (in meters)
  search_radius <- 0
  search_step <- 10 * 1000

  # loop until a water cell is found or maximum radius is reached
  while(!water_found){

    # increase the search radius
    search_radius <- search_radius + search_step

    # stop if the search radius exceeds 50 km
    if(search_radius >= max.distance.km*1000){
      return()
    }

    # create a buffer around the point to extract from the sf object
    buffered_point <- sf::st_buffer(point, dist=search_radius)

    # check if any features in the sf object intersect with the buffered point
    water_found <- any(sf::st_intersects(buffered_point, land.mask, sparse=FALSE))
  }

  # crop the sf object to the bounding box
  cropped_layer <- suppressWarnings(sf::st_crop(land.mask, sf::st_bbox(buffered_point)))

  # return the cropped spatial layer (RasterLayer or sf)
  return(cropped_layer)
}


################################################################################
# Define relocate point function ###############################################
################################################################################

#' Internal function to relocate a point to the nearest water feature
#'
#' This internal function takes a point and relocates it to the nearest
#' water feature within a specified search region. If no water feature is
#' found within the defined limits, it returns an empty point.
#' @param i An integer representing the index of the point in the `pointsOnLand` object to be relocated.
#' @param pointsOnLand An `sf` object containing points that may be located on land.
#' @param land_mask An `sf` object containing land geometries used to identify the nearest water feature.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.relocatePoint <- function(i, pointsOnLand, land.mask, max.distance.km){

  # retrieve current point
  point <- pointsOnLand[i,]

  # reduce search region to speed up calculations
  land_mask_i <- .setSearchRegion(point, land.mask, max.distance.km)

  # if no water if found within 50 km, return an empty point
  if(length(land_mask_i)==0) return(list(coords=data.frame("X"=NA, "Y"=NA), dist=NA))

  # get the nearest feature (line or pt) index
  nearest_index <- sf::st_nearest_feature(point, land_mask_i)

  # if a nearest feature is found:
  if(length(nearest_index)>0) {

    # extract the nearest geometry using the index
    nearest_geom <- land_mask_i[nearest_index, ]

    # check if the nearest geometry is a LINESTRING or POINT
    if (grepl("LINESTRING", sf::st_geometry_type(nearest_geom), fixed=TRUE)) {
      # calculate the nearest point on the line to the original point
      nearest_geom <- sf::st_nearest_points(point, nearest_geom)
      # the second point is the nearest point on the line
      nearest_point <- sf::st_point(sf::st_coordinates(nearest_geom)[2,])
      nearest_point <- sf::st_sfc(nearest_point, crs=sf::st_crs(nearest_geom))
    } else {
      # it's already a point
      nearest_point <- nearest_geom
    }

    # return the nearest point coordinates and distance
    list(coords=sf::st_coordinates(nearest_point)[,c("X","Y")],
         dist=as.numeric(sf::st_distance(point, nearest_point)))

    # if no nearest feature found, return NULL
  } else {
    return(list(coords=data.frame("X"=NA, "Y"=NA), dist=NA))
  }
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
