#######################################################################################################
# Function to calculate shortest path in water between consecutive positions ##########################
#######################################################################################################

#' Estimate step distances and reconstruct movement tracks
#'
#' @description This function calculates the distance travelled between consecutive positions
#' (step distances) for multiple individuals, reconstructing the underlying movement tracks
#' (shortest in-water paths). If no land shapefile is provided,
#' the function defaults to calculating linear (great-circle) distances, which are If no land shapefile is provided,
#' the function defaults to calculating linear (great-circle) distances, which are
#' faster to compute. When a land shapefile is provided, the function can utilize
#' parallel computing to expedite least-cost path estimation, depending on the
#' number of CPU cores specified. However, depending on the chosen grid resolution
#' and the spatial extent of the provided positions, this process may take a long time
#' to run, even with parallel computing enabled.
#'
#' @inheritParams as_moby
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
#' @param mov.directions Size of the movement neighbourhood used when building the least-cost graph:
#' `4` (rook), `8` (rook + bishop), or `16` (adds knight moves, the default; smoother least-cost paths).
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations. This parameter
#' is only relevant for estimating least-cost paths, specifically when a \code{land.shape} is provided.
#' Run \code{parallel::detectCores()} to check the number of available cores.
#' @param verbose Logical. Should the function output process information and display a progress bar?
#' Defaults to TRUE.
#'
#' @note Least-cost in-water paths are computed on a **terra**-rasterised cost surface routed with
#' **igraph** (Dijkstra); great-circle distances use **geosphere**. The graph is built once and reused
#' across a call's segments (and across calls when the cache is supplied).
#'
#' @return The input data (a \code{\link{mobyData}} when one was supplied, otherwise a data frame)
#' with an added `dist_m` column giving the distance from each position to the next consecutive
#' position of the same individual (in metres; `NA` for each individual's last position). The
#' movement trajectories (one spatial-line geometry per individual) are attached as the
#' `"trajectories"` attribute and retrieved with \code{\link{getTrajectories}}. The distance-enriched
#' output pipes directly into \code{\link{calculateROM}} and
#' \code{\link{calculateLinearityIndex}}.
#'
#' @seealso \code{\link{getTrajectories}}, \code{\link{calculateROM}},
#' \code{\link{calculateLinearityIndex}}, \code{\link{plotMaps}}
#'
#' @examples
#' data(rays)
#'
#' # great-circle step distances between consecutive positions
#' # (no land shape supplied: fast linear paths)
#' rays_dist <- calculateStepDistances(rays)
#' head(rays_dist$dist_m)
#'
#' \donttest{
#' # shortest in-water paths around a coastline (least-cost routing)
#' land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
#'   c(-8.99, 38.44), c(-8.96, 38.44), c(-8.96, 38.47),
#'   c(-8.99, 38.47), c(-8.99, 38.44)))), crs = 4326))
#' rays_lc <- calculateStepDistances(rays[1:60, ], land.shape = land,
#'                                   grid.resolution = 200)
#' }
#'
#' @export


calculateStepDistances <- function(data,
                            land.shape = NULL,
                            epsg.code = NULL,
                            grid.resolution = 100,
                            mov.directions = 16,
                            id.col = NULL,
                            lon.col = NULL,
                            lat.col = NULL,
                            cores = 1,
                            verbose = TRUE){

  # The least-cost graph is a purely internal implementation detail: it is built on demand, used, and
  # discarded. Only the distance-enriched data (with its trajectories) is returned.
  .stepDistances(data = data, land.shape = land.shape, epsg.code = epsg.code,
                 grid.resolution = grid.resolution, mov.directions = mov.directions,
                 id.col = id.col, lon.col = lon.col, lat.col = lat.col, cores = cores,
                 land.shape.name = if(!is.null(land.shape)) deparse(substitute(land.shape)) else NULL,
                 cost.graph = NULL, verbose = verbose)$data
}


#' Internal step-distance engine
#'
#' @description Shared implementation behind \code{\link{calculateStepDistances}}. Returns a list with
#'   the distance-enriched `data` and the `cost.graph` that was used (`NULL` when no land shape was
#'   supplied), so that internal callers - notably `filterDetections()`, whose iterative speed filter
#'   makes many repeated calls - can build the least-cost graph ONCE over the full set of positions and
#'   reuse it. Reusing one full-extent graph is not merely an optimisation: a graph rebuilt from a
#'   two-point subset spans only that pair's bounding box, so paths that need to route around larger
#'   landmasses would silently fall back to straight-line distances. The cost graph is never returned
#'   to users.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.stepDistances <- function(data,
                           land.shape = NULL,
                           epsg.code = NULL,
                           grid.resolution = 100,
                           mov.directions = 16,
                           id.col = NULL,
                           lon.col = NULL,
                           lat.col = NULL,
                           cores = 1,
                           land.shape.name = NULL,
                           cost.graph = NULL,
                           verbose = TRUE){

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # capture mobyData metadata (re-attached to the enriched output at the end)
  meta <- attr(data, "moby")

  # the cost graph in use (stays NULL unless a land shape is supplied)
  trCost <- NULL

  land_shape_name <- land.shape.name
  # if no land shape is provided, set related parameters to NULL
  if(is.null(land.shape)){
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
  # great-circle distances on geographic coordinates require no projection: skip .processSpatial
  # (which would otherwise demand an 'epsg.code') when no land.shape is involved
  if(is.null(land.shape) && coords_crs=="geographic" && is.null(epsg.code)){
    sf::st_crs(coords) <- 4326
    epsg.code <- 4326
  }else{
    spatial_data <- .processSpatial(coords, land.shape, epsg.code)
    coords <- spatial_data$coords
    land.shape <- spatial_data$spatial.layer
    epsg.code <- spatial_data$epsg.code
  }


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

    # check if any coordinates overlap land (per-point, so the count is accurate)
    pts <- sf::st_as_sf(data.frame(lon=data$lon_m_tmp, lat=data$lat_m_tmp),
                        coords=c("lon","lat"), crs=sf::st_crs(epsg.code))
    overlapping_pts <- lengths(sf::st_intersects(pts, sf::st_as_sf(land.shape), sparse=TRUE))
    if(any(overlapping_pts > 0)){
      num_overlapping <- sum(overlapping_pts > 0)  # count of points overlapping land
      warning_str <- paste0("- Some coordinates (n=", num_overlapping,") overlap with the supplied land shape. ",
      "Consider using the 'correctPositions()' function to relocate these points to the nearest marine cell, ",
      "and then rerun the current function with the updated positions.")
      warning(paste(strwrap(warning_str, width=getOption("width")), collapse="\n"), call. = FALSE)
    }

    # reuse a cost graph handed down by an internal caller, or build one from the land shape
    if(!is.null(cost.graph)) {
      trCost <- cost.graph
      grid.resolution <- trCost$grid$resx
    # build the least-cost graph (terra raster -> igraph); replaces gdistance transition/geoCorrection
    } else {
      if(verbose) {cat(paste0("Building least-cost graph (", grid.resolution, "m grid | ", mov.directions, " directions)\n"))}
      trCost <- .buildCostGraph(land.shape, coords, grid.resolution, mov.directions, epsg.code)
    }

  }


  ##############################################################################
  ## Split data by individual and initialize variables #########################
  ##############################################################################

  # split COAs by individual
  data_individual <- split(data, f=data[,id.col], drop=FALSE)

  # initialize progress bar
  if(verbose){
    pb <- txtProgressBar(min=0, max=length(data_individual), initial=0, style=3)
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
    for (i in seq_along(data_individual)) {

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
    for (i in seq_along(data_individual)) {

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
    # trCost (an igraph graph + a plain-R grid descriptor) and calculateLeastCost (a self-contained
    # closure) both serialise to workers; no non-serialisable SpatRaster crosses the process boundary.
    results_list <- foreach::foreach(i=seq_along(data_individual), .options.snow=opts, .packages=c("sf", "igraph", "geosphere")) %dopar% {
      calculateLeastCost(data_individual[[i]], land.shape, trCost, grid.resolution, epsg.code)
    }

    # aggregate results from parallel processing
    for (i in seq_along(results_list)) {
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
  final_data <- .dropCols(final_data, out_cols)

  # attach the movement trajectories as an attribute, so the primary return value is the
  # distance-enriched data frame: consistent with moby's other column-adding functions and directly
  # pipeable into calculateROM() / calculateLinearityIndex().
  attr(final_data, "trajectories") <- final_trajectories

  # processing metadata
  attr(final_data, "land.shape") <- land_shape_name
  attr(final_data, "epsg.code") <- epsg.code
  attr(final_data, "grid.resolution") <- grid.resolution
  attr(final_data, "mov.directions") <- mov.directions
  attr(final_data, "cores") <- cores
  attr(final_data, "processing.date") <- Sys.time()

  # restore the mobyData class/metadata when the input carried it
  if(!is.null(meta)){
    attr(final_data, "moby") <- meta
    class(final_data) <- unique(c("mobyData", "data.frame"))
  }

  # return the distance-enriched data (trajectories retrieved via getTrajectories()) alongside the cost
  # graph, for internal reuse only - calculateStepDistances() hands users just the data.
  return(list(data = final_data, cost.graph = trCost))
}


################################################################################
# Accessor: movement trajectories ##############################################
################################################################################

#' Extract movement trajectories from a calculateStepDistances() result
#'
#' @description Retrieves the per-individual movement trajectories (spatial-line geometries)
#' attached by \code{\link{calculateStepDistances}} to its distance-enriched output.
#'
#' @param x A distance-enriched data frame returned by \code{\link{calculateStepDistances}}.
#' @return A named list of trajectory geometries (one element per individual), or `NULL` if the
#' object carries no trajectories.
#' @seealso \code{\link{calculateStepDistances}}, \code{\link{plotMaps}}
#'
#' @examples
#' data(rays)
#' rays_dist <- calculateStepDistances(rays)
#' traj <- getTrajectories(rays_dist)
#' length(traj)
#'
#' @export

getTrajectories <- function(x) attr(x, "trajectories")



################################################################################
# Build the least-cost cost-graph (terra raster -> igraph) #####################
################################################################################

#' Build a least-cost cost-graph from a land shape
#'
#' @description Rasterises the land shape with **terra** and builds a weighted **igraph** graph whose
#' edge weights reproduce gdistance's `geoCorrection(type = "c")` conductances exactly
#' (`weight = centre-to-centre distance / mean(cell conductance)`). The returned `mobyCostGraph`
#' bundles the graph with a plain-R grid descriptor (extent origin, per-axis cell size, dimensions,
#' CRS) sufficient to map cell indices to centre coordinates arithmetically. It carries no live
#' `SpatRaster`, so it serialises cleanly to parallel workers. All land polygons are treated as
#' impassable (`field = 1`); this fixes a latent bug in the former gdistance path, where
#' `raster::rasterize()` with no field assigned feature indices and only the first polygon became a
#' hard barrier.
#' @param land.shape An `sf` land layer, already projected to `epsg.code`.
#' @param coords Projected `sf` detection points (used only for the grid extent).
#' @param grid.resolution Grid cell size in metres.
#' @param directions Movement neighbourhood: 4, 8, or 16.
#' @param epsg.code Projected CRS of the grid.
#' @return A `mobyCostGraph`: a list with `graph` (igraph) and `grid` (descriptor).
#' @keywords internal
#' @noRd

.buildCostGraph <- function(land.shape, coords, grid.resolution, directions, epsg.code){

  # grid extent = bbox(coords) * 1.2, snapped so the resolution stays EXACT (anchor at xmin/ymax and
  # expand to whole cells - the raster convention). terra::rast(resolution=) would instead keep the
  # extent and adjust the resolution (anisotropic), which shifts the geoCorrection distances.
  bb <- sf::st_bbox(coords)
  hx <- (bb[["xmax"]] - bb[["xmin"]]) / 2; cx <- (bb[["xmax"]] + bb[["xmin"]]) / 2
  hy <- (bb[["ymax"]] - bb[["ymin"]]) / 2; cy <- (bb[["ymax"]] + bb[["ymin"]]) / 2
  xmin <- cx - 1.2 * hx; ymax <- cy + 1.2 * hy
  nc <- max(1L, round((2 * 1.2 * hx) / grid.resolution))
  nr <- max(1L, round((2 * 1.2 * hy) / grid.resolution))
  ext <- terra::ext(xmin, xmin + nc * grid.resolution, ymax - nr * grid.resolution, ymax)
  r <- terra::rast(ext, resolution = grid.resolution, crs = sf::st_crs(epsg.code)$wkt)
  terra::values(r) <- 0

  # rasterise ALL land as impassable (field = 1); conductance basis: water = 1, land = 1e-8
  land_v <- terra::vect(sf::st_geometry(land.shape))
  r <- terra::rasterize(land_v, r, field = 1, update = TRUE, background = 0)
  rv <- terra::values(r)[, 1]
  v <- ifelse(!is.na(rv) & rv == 1, 1e-8, 1)

  nr <- terra::nrow(r); nc <- terra::ncol(r); N <- nr * nc
  resx <- terra::res(r)[1]; resy <- terra::res(r)[2]
  grid <- list(xmin = terra::xmin(r), ymax = terra::ymax(r), resx = resx, resy = resy,
               nr = nr, nc = nc, epsg = epsg.code)

  # neighbour offsets (canonical half; antipodes complete the undirected neighbourhood):
  # 4 = rook, 8 = rook + bishop, 16 = + knight moves
  off <- switch(as.character(directions),
                "4"  = rbind(c(1, 0), c(0, 1)),
                "8"  = rbind(c(1, 0), c(0, 1), c(1, 1), c(1, -1)),
                "16" = rbind(c(1, 0), c(0, 1), c(1, 1), c(1, -1), c(1, 2), c(1, -2), c(2, 1), c(2, -1)),
                stop("'mov.directions' must be 4, 8, or 16.", call. = FALSE))

  cell <- seq_len(N); row <- ((cell - 1L) %/% nc) + 1L; col <- ((cell - 1L) %% nc) + 1L
  parts <- lapply(seq_len(nrow(off)), function(k){
    dr <- off[k, 1]; dc <- off[k, 2]; nrow_ <- row + dr; ncol_ <- col + dc
    ok <- nrow_ >= 1 & nrow_ <= nr & ncol_ >= 1 & ncol_ <= nc
    from <- cell[ok]; to <- (nrow_[ok] - 1L) * nc + ncol_[ok]
    d <- sqrt((dc * resx)^2 + (dr * resy)^2)
    list(from = from, to = to, w = d / ((v[from] + v[to]) / 2))
  })
  from <- unlist(lapply(parts, `[[`, "from")); to <- unlist(lapply(parts, `[[`, "to"))
  w    <- unlist(lapply(parts, `[[`, "w"))
  g <- igraph::add_edges(igraph::make_empty_graph(n = N, directed = FALSE), rbind(from, to), weight = w)
  structure(list(graph = g, grid = grid), class = "mobyCostGraph")
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
    for(r in seq_len(nrow(coords)-1)){
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
        # least-cost in-water path via the terra+igraph cost graph (replaces gdistance::shortestPath).
        # Cell<->coordinate mapping is done arithmetically from the plain-R grid descriptor so this
        # closure stays self-contained (no live SpatRaster), and therefore serialises to parallel workers.
        g <- trCost$grid
        xy2 <- rbind(coords_m[r, ], coords_m[r + 1, ])
        gcol <- floor((xy2[, 1] - g$xmin) / g$resx) + 1L
        grow <- floor((g$ymax - xy2[, 2]) / g$resy) + 1L
        inside <- gcol >= 1 & gcol <= g$nc & grow >= 1 & grow <= g$nr
        if(all(inside)){
          ci <- (grow - 1L) * g$nc + gcol
          vp <- igraph::shortest_paths(trCost$graph, ci[1], ci[2],
                                       weights = igraph::E(trCost$graph)$weight, output = "vpath")$vpath[[1]]
          if(length(vp) >= 2){
            vr <- ((as.integer(vp) - 1L) %/% g$nc) + 1L; vc <- ((as.integer(vp) - 1L) %% g$nc) + 1L
            xy <- cbind(g$xmin + (vc - 0.5) * g$resx, g$ymax - (vr - 0.5) * g$resy)
            segment <- sf::st_sfc(sf::st_linestring(xy), crs = sf::st_crs(epsg.code))
            segment_wgs84 <- sf::st_transform(segment, crs = sf::st_crs("+proj=longlat +datum=WGS84"))
            # distance along the least-cost polyline (great-circle length of the cell-centre route)
            dists[r] <- sum(geosphere::distVincentyEllipsoid(sf::st_coordinates(segment_wgs84)[, 1:2]))
          }
        }
        # (if an endpoint falls outside the graph extent or the goal is unreachable, keep the
        #  straight-line 'segment' and great-circle 'dists[r]' as a safe fallback)
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
