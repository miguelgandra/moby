#######################################################################################################
# Calculate KUDs - Kernel Utilization Distributions  ##################################################
#######################################################################################################

#' Calculate animals' kernel utilization areas

#' @description This function calculates animal utilization areas using kernel
#' density estimation, based on the \code{\link[adehabitatHR]{kernelUD}} function
#' (`adehabitatHR` package). It allows for estimating utilization areas corresponding to
#' specified utilization contours (isopleths), with default values of 50% and 95%.
#' Additionally, the function subtracts portions of these areas that overlap with
#' landmasses and can independently estimate kernel utilization distributions (KUDs)
#' for different groups or time periods.
#'
#' @inheritParams setDefaults
#' @param data A data frame with animal positions (COAs), containing 'longitude' and 'latitude' columns.
#' @param bandwidth Numeric. The smoothing parameter (h) that controls the width of the kernel
#' used in the UD estimation. The bandwidth controls how much smoothing is applied to the data:
#' smaller values result in a more detailed and localized estimation of space use,  while larger
#' values produce smoother, more broad distributions. This parameter is passed to the `h` argument
#' of the \code{\link[adehabitatHR]{kernelUD}} function. For further information, see
#' the "Details" section below.
#' @param spatial.grid Optional. A `Raster` or `SpatialPixels` object representing the
#' grid over which the animal kernel utilization distributions (KUDs) will be estimated
#' (see the `grid` argument in \code{\link[adehabitatHR]{kernelUD}}). If set to `NULL`,
#' the function will automatically generate an appropriate grid based on the spatial extent
#' of the supplied animal's positions.
#' @param subset Optional. A variable used to subset the data, allowing KUDs to be calculated independently
#' for each level of this variable. This should be the name of a column in the provided dataset.
#' If left `NULL`, KUDs are calculated for the whole monitoring period.
#' @param id.groups Optional. A named list where each element represents a group (e.g., species, sex, or age class),
#' containing a vector of IDs for that group. If supplied, KUDs will be calculated
#' independently for each group. The names of the list correspond to the group labels.
#' @param land.shape Optional. A shapefile containing coastlines or landmasses,
#' provided either as an 'sf' object or as a 'SpatialPolygonsDataFrame'/'SpatialPolygons' object.
#' If the input is not in 'sf' format, the function will automatically convert it to 'sf'
#' to ensure compatibility with subsequent spatial operations. Used to clip and exclude
#' any portions of the estimated areas that overlap with landmasses.
#' @param contour.percent Numeric vector. The percentages for which isopleths (contour areas)
#' are calculated. Defaults to 50% and 95%, representing core and total areas of utilization.
#' @param verbose Logical. If TRUE, the function will print detailed processing
#' information. Defaults to FALSE.
#'
#' @return
#' A list containing:
#' \item{kernel_density}{An object of class "estUDm" with the estimated kernel utilization distributions.
#' If a subset variable was provided, this will contain a list of "estUDm" objects, one for each group.}
#' \item{summary_table}{A data frame summarizing the area estimates for the specified contours (isopleths).}
#' \item{contours}{One 'sf' object for each specified contour. If a subset variable was provided, an additional
#' column is included to distinguish between different groups.}
#'
#' The results list also contains multiple attributes to store relevant metadata,
#' such as function options and processing details. These attributes might be useful
#' for tracking parameters and ensuring reproducibility of the analysis.
#'
#'
#' @details
#'
#' This function estimates kernel utilization distributions (KUDs) based on animal location data,
#' allowing for quick analyses of space-use patterns. For a comprehensive overview of other
#' home-range estimation methods, check out Kraft et al. (2023) (full reference below in the References section).
#' The function also includes options for handling landmasses and for grouping data by
#' subsets or groups for independent analysis.
#'
#' **Land clipping**: Land clipping is applied post-hoc, after kernel density estimation.
#' If you need to account for physical barriers like land during UD estimation, consider alternative
#' methods (e.g. dynamic Brownian Bridge Movement Models as provided in the `RSP` package;
#' Niella et al. 2020).
#'
#' **Bandwidth (h)**: The smoothing factor, or bandwidth (h), is a critical parameter in
#' kernel utilization distribution (KUD) analysis, representing the standard deviation of the kernel.
#' It defines the extent to which a location can influence the home range estimation and
#' the overall density estimate. The choice of bandwidth significantly impacts the results
#' of the analysis. The bandwidth can either be fixed (using a single value for all data points)
#' or variable (adapting based on point density). The current function only implements a fixed bandwidth.
#'
#' - A larger bandwidth increases the influence of more distant data points on the home range estimation,
#' leading to a wider utilization distribution (UD) and a larger overall home range size.
#' This increased smoothing can help mitigate sampling errors but may obscure finer details,
#' retaining only the most prominent features of the spatial data.
#'
#' - Conversely, a smaller bandwidth allows for greater detail at smaller spatial scales
#' but might result in more fragmented home range estimates.
#'
#' ### Bandwidth Selection Methods
#' The optimal bandwidth is not universally defined and depends on the specific context and dataset.
#' Several methods are available for selecting the bandwidth, including:
#' - **Reference Bandwidth (Href):** A constant based on data variance, assuming normally
#' distributed data. It provides consistent smoothing but may overestimate home range size in sparse datasets.
#' - **LSCV (Least-Squares Cross-Validation):** Minimizes error between observed and
#' predicted data, but may undersmooth in sparse datasets.
#' - **Ad-Hoc Choice:** Researchers manually select a bandwidth based on prior
#' knowledge or exploratory data analysis.
#' - **Direct Plug-In:** Estimates bandwidth by solving equations based on the second
#' derivative of the kernel density. Computationally intensive but can yield robust results, especially in complex datasets.
#'
#' In some cases, bandwidth is adjusted based on the detection range in the array, particularly when spatial estimates produce smaller, disconnected isopleths.
#' Adjusting the bandwidth can help form continuous areas that better reflect biological relevance.
#'
#' @references
#' Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., & Abecasis, D. (2023).
#' Residency and space use estimation methods based on passive acoustic telemetry data.
#' Movement Ecology, 11(1), 12.
#' https://doi.org/10.1186/s40462-023-00349-y
#'
#' Niella, Y., Flávio, H., Smoothey, A. F., Aarestrup, K., Taylor, M. D., Peddemors,
#' V. M., & Harcourt, R. (2020). Refined Shortest Paths (RSP): Incorporation of
#' topography in space use estimation from node‐based telemetry data.
#' Methods in Ecology and Evolution, 11(12), 1733-1742.
#' https://doi.org/10.1111/2041-210X.13484
#'
#' Worton, B. J. (1989).
#' Kernel methods for estimating the utilization distribution in home‐range studies.
#' Ecology, 70(1), 164-168.
#' https://doi.org/10.2307/1938423
#'
#' @seealso \code{\link[adehabitatHR]{kernelUD}}
#' @export


################################################################################
## Main function - Estimate Kernel Density #####################################

calculateKUDs <- function(data,
                          bandwidth,
                          spatial.grid = NULL,
                          subset = NULL,
                          id.groups = NULL,
                          land.shape = NULL,
                          id.col = getDefaults("ID"),
                          lon.col = getDefaults("lon"),
                          lat.col = getDefaults("lat"),
                          epsg.code = getDefaults("epsg"),
                          contour.percent = c(50,95),
                          verbose = FALSE) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # capture the original spatial.layer name
  if(!is.null(land.shape)){
    land_shape_name <- deparse(substitute(land.shape))
  }else{
    land_shape_name <- NULL
  }

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  land.shape <- reviewed_params$land.shape

  # validate additional parameters
  errors <- c()
  if(!requireNamespace("adehabitatHR", quietly=TRUE)) errors <- c(errors, "The 'adehabitatHR' package is required for this function but is not installed. Please install 'adehabitatHR' using install.packages('adehabitatHR') and try again.")
  if(!inherits(bandwidth, "numeric")) errors <- c(errors, "Error: 'bandwidth' must be a numeric value.")
  if(!inherits(contour.percent, "numeric")) errors <- c(errors, "Error: 'contour.percent' must be a numeric value.")
  if(any(contour.percent < 1 | contour.percent > 100)) errors <- c(errors, "Error: 'contour.percent' must be between 1 and 100.")
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # manage spatial objects
  coords <- sf::st_as_sf(data, coords=c(lon.col, lat.col))
  coords_crs <- .checkProjection(coords)
  spatial_data <- .processSpatial(coords, land.shape, epsg.code)
  coords <- spatial_data$coords
  land.shape <- spatial_data$spatial.layer
  epsg.code <- spatial_data$epsg.code

  # retrieve coords bounding box
  coords_bbox <- sf::st_bbox(coords)

  # convert spatial.grid from Raster to SpatialPixels if required
  if(!is.null(spatial.grid)){
    if(inherits(spatial.grid, "Raster")) spatial.grid <- methods::as(spatial.grid, 'SpatialPixels')
  }

  # convert subset variable to factor
  if(!is.null(subset)){
     if(!inherits(subset, "factor")){
       coords[[subset]] <- as.factor(coords[[subset]])
       warning(paste0("- The subset variable (", subset, ") has been converted to a factor."), call.=FALSE)
     }
  }

  # set up id.groups
  if(!is.null(id.groups)){
    id_groups_df <- data.frame("ID"=unlist(id.groups), "id.groups"=rep(names(id.groups), times=sapply(id.groups, length)), stringsAsFactors=FALSE)
    colnames(id_groups_df)[1] <- id.col
    coords <- merge(coords, id_groups_df, by=id.col, all.x=TRUE)
    if(is.null(subset)){
      subset <- "id.groups"
    }else{
      coords[[paste0(subset, ".by.ids")]] <- interaction(coords[[subset]], coords[["id.group"]], drop=TRUE)
      subset <- paste0(subset, ".by.ids")
    }
  }

  # check if KUDs need to be calculated for different groups
  multiple <- !is.null(subset)


  ##############################################################################
  ## Calculate KUDs ############################################################
  ##############################################################################

  #######################################################################
  # if no subset provided, calculate the KUDs for the whole dataset #####
  if(!multiple) {

    # print status message to console
    cat("Estimating kernel utilization distributions...\n")

    # compute KUDs for the entire dataset
    final_results <- .computeKUDs(coords, id.col, bandwidth, coords_bbox, contour.percent, spatial.grid, land.shape, epsg.code, verbose)

    # extract the 'spatial.grid' from the final results
    spatial.grid <- final_results$spatial.grid

    # remove the 'spatial.grid' element from the final results list
    final_results <- final_results[-length(final_results)]


  #######################################################################
  # else, split the data and calculate KUDs for each group separately ###
  } else{

    # print the grouping criteria (defined by 'subset') to the console
    cat(paste0("Grouping data by: ", paste(subset, collapse=" / "), "\n"))
    if(verbose) cat("\n")

    # split the data into groups based on the 'subset' variable
    group_coords <- split(coords, f=coords[[subset]], drop=TRUE)
    group_coords <- lapply(group_coords, function(x){x[[id.col]] <- droplevels(x[[id.col]]); return(x)})

    # check if any individuals (IDs) are repeated across the different groups
    repeated_ids <- list()
    for(i in 1:length(group_coords)){
      # get the IDs within the current group
      subset_ids <- levels(group_coords[[i]][[id.col]])
      # get the IDs present in all other groups
      remaining_ids <- levels(do.call("rbind", group_coords[-i])[[id.col]])
      # check if any ID from the current group is found in the other groups
      repeated_ids[[i]] <- any(subset_ids %in% remaining_ids)
    }
    repeated_ids <- any(unlist(repeated_ids))

    # compute KUDs for each group separately
    kud_results <- lapply(1:length(group_coords), function(i) {
      cat(paste0("Estimating kernel utilization distributions [", names(group_coords)[i], "]...\n"))
      group_results <- .computeKUDs(group_coords[[i]], id.col, bandwidth, coords_bbox, contour.percent, spatial.grid, land.shape, epsg.code, verbose)
      if(is.null(group_results)) return(NA)
      else return(group_results)
    })

    # assign names to each group result
    names(kud_results) <- names(group_coords)

    # remove empty elements from list
    kud_results <- kud_results[!unlist(lapply(kud_results, function(x) all(is.na(x))), recursive=FALSE)]

    #  extract and remove the 'spatial.grid' from the final results
    spatial.grid <- lapply(kud_results, function(x) x$spatial.grid)[[1]]
    kud_results <- lapply(kud_results, function(x) x[-length(x)])

    # combine the summary tables
    if(repeated_ids){
      summary_table <- lapply(kud_results, function(x) x$summary_table)
      summary_table <- mapply(function(table, group){colnames(table)[-1] <- paste0(colnames(table)[-1], " [", group, "]"); return(table)},
                              table=summary_table, group=names(kud_results), SIMPLIFY=FALSE)
      summary_table <- Reduce(function(x,y) merge(x, y, all=TRUE, by=id.col), summary_table)
      summary_table <- summary_table[match(levels(data[[id.col]]), summary_table[[id.col]]), ]
      summary_table <- summary_table[!is.na(summary_table[[id.col]]),]
      coa_cols <- grepl("COAs", colnames(summary_table), fixed=TRUE)
      summary_table[coa_cols][is.na(summary_table[coa_cols])] <- 0
      summary_table[is.na(summary_table)] <- "-"
    }else{
      summary_table <- lapply(kud_results, function(x) x$summary_table)
      summary_table <- do.call("rbind", summary_table)
      summary_table$group <- sub("(.*)\\.[^\\.]*$", "\\1", rownames(summary_table))
      rownames(summary_table) <- NULL
      summary_table <- summary_table[,c(ncol(summary_table), 1, 2:(ncol(summary_table)-1))]
      summary_table[["N COAs"]][is.na(summary_table[["N COAs"]])] <- 0
      summary_table[is.na(summary_table)] <- "-"
    }

    # initialize an empty list to store combined results for each contour percent
    contour_results <- list()

    # iterate through each group in 'kud_results'
    for (group in names(kud_results)) {
      group_results <- kud_results[[group]]
      # iterate through each contour percent defined in contour.percent
      for (k in names(group_results)[grepl("K", names(group_results))]) {
        # initialize a list for the current contour percent if it doesn't exist
        if (is.null(contour_results[[k]])) {
          contour_results[[k]] <- list()
        }
        # retrieve the kernel density result for the current group and contour percent
        contour_sf <- group_results[[k]]
        # add a new column with the group level to the sf object
        contour_sf$group <- group
        colnames(contour_sf)[colnames(contour_sf)=="group"] <- subset
        # append the group results to the corresponding contour result
        contour_results[[k]][[length(contour_results[[k]]) + 1]] <- contour_sf
      }
    }

    # combine all the sf objects for each contour percent into one
    for (k in names(contour_results)) {
      contour_results[[k]] <- do.call(rbind, contour_results[[k]])
    }

    # prepare the final results, including kernel density estimates, summary table and specified contours
    final_results <- list("kernel_density"=lapply(kud_results, function(x) x$kernel_density),
                         "summary_table"=summary_table)
    final_results <- c(final_results, contour_results)

  }

  ##############################################################################
  ## Wrap it up ################################################################
  ##############################################################################

  # print the total execution time to the console
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n"))

  # create attributes to save relevant metadata
  attr(final_results, 'id.groups') <- id.groups
  attr(final_results, 'bandwidth') <- bandwidth
  attr(final_results, 'contour.percent') <- contour.percent
  attr(final_results, 'subset') <- subset
  attr(final_results, 'land.shape') <- land_shape_name
  attr(final_results, 'epsg.code') <- epsg.code
  attr(final_results, 'grid.extent') <- raster::extent(spatial.grid)
  attr(final_results, 'grid.res') <- unique(spatial.grid@grid@cellsize)
  attr(final_results, 'grid.dims') <- spatial.grid@grid@cells.dim
  attr(final_results, 'processing.date') <- Sys.time()

  # return results
  return(final_results)
}


################################################################################
# Helper function I - remove individuals with less than 5 records ##############
# required to avoid errors in the kernel_UD function ###########################
################################################################################

#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.cleanData <- function(coords, id.col) {

  # remove unused factor levels from the specified ID column
  coords[[id.col]] <- droplevels(coords[[id.col]])

  # get unique IDs that have 5 or more occurrences
  valid_ids <- names(which(table(coords[[id.col]]) >= 5))

  # filter the data to keep only valid IDs
  coords <- coords[coords[[id.col]] %in% valid_ids, ]

  # drop unused levels again in the id column
  coords[[id.col]] <- droplevels(coords[[id.col]])

  # return the filtered data
  return(coords)

}



################################################################################
# Helper function II - subtract area on land ###################################
################################################################################

#' Subtracts kernel density values that overlap with land areas.
#'
#' This function iterates through kernel densities and sets density values to
#' zero for coordinates that overlap with a given land shape.
#'
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.subtractLand <- function(kernel.densities, land.shape, epsg.code, verbose) {

  # extract grid coordinates from the kernel density object
  grid_coords <- as.data.frame(kernel.densities[[1]]@coords)
  colnames(grid_coords) <- c("x","y")

  # create an sf object from kernel density coordinates
  grid_coords <- sf::st_as_sf(grid_coords, coords=c("x", "y"), crs=epsg.code)

  # check for overlaps with the land.shape sf object
  overlap_indexes <- as.logical(sf::st_intersects(grid_coords, sf::st_union(land.shape), sparse=FALSE))

  # initialize a counter for the number of individuals with overlapping areas
  n_corrected <- 0

  # iterate over each kernel density object
  for(i in 1:length(kernel.densities)){

    # access the density values from the current kernel density object
    density_values <- kernel.densities[[i]]$ud

    # total density before correction (including any overlap with land)
    total_density <- sum(density_values, na.rm=TRUE)

    # check if there are any density values overlapping with land
    if(any(density_values[overlap_indexes]>0)) n_corrected <- n_corrected + 1

    # set the overlapping density values to zero
    density_values[overlap_indexes] <- 0

    # total density after correction (land overlap set to zero)
    corrected_density <- sum(density_values, na.rm=TRUE)

    # standardize values back to the original total density
    if(corrected_density>0) {
      density_values <- density_values * (total_density / corrected_density)
    }

    # update the density values in the current kernel density object
    kernel.densities[[i]]$ud <- density_values
  }

  # print a message indicating how many kernel densities were corrected
  if(verbose){
    if(n_corrected==0){
      cat("No land overlap detected\n")
    }else{
      kud_label <- ifelse(n_corrected==1, "KUD", "KUDs")
      cat(paste0("Land overlap areas successfully subtracted from ", n_corrected, " ", kud_label, "\n"))
    }
  }

  # return the corrected kernel densities
  return(kernel.densities)
}


##############################################################################
## Helper function III - Set grid ############################################
##############################################################################

#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.createGrid <- function(coords.bbox, expand.factor=0.1, epsg.code){

  # format expand.factor
  expand.factor <- rep(expand.factor*c(-1,1), each=2)

  # expand bounding box by a given % in all directions
  dx <- coords.bbox["xmax"]-coords.bbox["xmin"]
  dy <- coords.bbox["ymax"]-coords.bbox["ymin"]
  grid_bbox <- coords.bbox + c(dx, dy, dx, dy) * expand.factor

  # get extent dimensions
  bbox_width <- grid_bbox["xmax"] - grid_bbox["xmin"]
  bbox_height <- grid_bbox["ymax"] - grid_bbox["ymin"]

  # calculate the total area of the bounding box
  grid_area <- bbox_width * bbox_height

  # calculate the approximate size of each cell (in square units)
  cell_size <- sqrt(grid_area/1000000)

  # round the cell size to the nearest 50
  rounded_cell_size <- round(cell_size/50) * 50
  if(rounded_cell_size==0) rounded_cell_size <- 5

  # check if rounded_cell_size results in at least 250 x 250 cells
  ncol <- ceiling(bbox_width / rounded_cell_size)
  nrow <- ceiling(bbox_height / rounded_cell_size)
  if (ncol < 250 | nrow < 250) { rounded_cell_size <- 1}

  # create raster grid
  spatial.grid <- raster::raster(raster::extent(grid_bbox), res=rounded_cell_size, crs=epsg.code$proj4string)

  # convert to  SpatialPixels class
  return(methods::as(spatial.grid, 'SpatialPixels'))
}


################################################################################
# Define relocate point function ###############################################
################################################################################

#' Internal function to relocate a point to the nearest water feature
#'
#' This internal function takes a point and relocates it to the nearest
#' water feature within a specified search region. If no water feature is
#' found within the defined limits, it returns an empty point.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.computeKUDs <- function(coords, id.col, bandwidth, coords.bbox, contour.percent,
                         spatial.grid, land.shape, epsg.code, verbose){

  # initialize the expand factor, which will be used to gradually increase the spatial grid extent
  expand_factor <- 0.05

  # check if a spatial.grid was supplied
  grid_supplied <- !is.null(spatial.grid)

   # filter out individuals with less than 5 detections
  filtered_coords <- .cleanData(coords, id.col)

  # if no coordinates remain after filtering, return an empty result
  if(nrow(filtered_coords)==0) return(NULL)

  # convert data to class SpatialPointsDataFrame (required input for the kernelUD function)
  filtered_coords <- sf::as_Spatial(filtered_coords)

  #######################################################################################
  # start the repeat loop to compute KUDs and gradually expand the spatial grid if needed
  repeat {

    # increase the grid extent incrementally by a factor of 0.05 (5%)
    expand_factor <- expand_factor + 0.05
    if(verbose && expand_factor>0.1) cat(paste0("Expanding grid bounding box by: ",round(expand_factor*100),  "%\n"))

    # expand bounding box by a given % in all directions
    # create a custom spatial grid based on the coordinate extent, if not supplied by the user
    if(!grid_supplied) spatial.grid <-.createGrid(coords.bbox, expand.factor=expand_factor, epsg.code)

    # estimate the kernel utilization distributions (KUD) for the provided coordinates
    kernel_density <- adehabitatHR::kernelUD(filtered_coords[,id.col], h=bandwidth, grid=spatial.grid)

    # clip out areas of the kernel density that overlap with landmasses
    if(!is.null(land.shape)){
      if(verbose) cat("Clipping out land masses...\n")
      kernel_density <- .subtractLand(kernel_density, land.shape, epsg.code, verbose)
    }

    # initialize a list to store kernel contours for each contour percentage (e.g., 50%, 95%)
    kernel_contours <- vector("list", length(contour.percent))
    names(kernel_contours) <- paste0("K", contour.percent)

    # loop through each contour percentage and calculate corresponding KUD contours
    for(c in 1:length(contour.percent)){
      # print progress message
      if(verbose) cat(paste0("Calculating ", contour.percent[c], "% contours...\n"))
      # attempt to calculate kernel contours for the specified contour percent
      kernel_contours[[c]] <- tryCatch({
        adehabitatHR::getverticeshr(kernel_density, percent=contour.percent[c], unin="m", unout="km2")
      }, error = function(e) {
        # check if grid extent needs to be increased
        if (grepl("grid is too small", e$message, fixed=TRUE)) {
          if(grid_supplied) stop(e$message)
          else {
            if(verbose) cat("Grid too small. Generating a new one...\n")
            return(NA)
          }
        } else {
          stop(e$message)
        }
      })
    }

    # if kernel contours were successfully generated
    if(all(unlist(lapply(kernel_contours, function(x) inherits(x, "SpatialPolygonsDataFrame"))))){
      # convert each contour to an `sf` object
      kernel_contours <- lapply(kernel_contours, sf::st_as_sf)
      # extract contour IDs and areas, removing geometry data
      kernel_areas <- lapply(kernel_contours, sf::st_drop_geometry)
      # break out of the repeat loop since contours were successfully calculated
      break
    }
  }

  #######################################################################################
  # merge the individual data frames for each contour percentage into a single table
  kud_table <- Reduce(function(x, y) merge(x, y, by="id", all=TRUE), kernel_areas)

  # format the area values to two decimal places
  kud_table[,-1] <- apply(kud_table[,-1], 2, function(x) sprintf("%.2f", x))

  # set column names of the KUD table
  colnames(kud_table) <- c("ID",  sprintf("KUD %d%% (Km2)", contour.percent))

  # add the number of COAs to the KUD table
  ncoas <- table(coords[[id.col]])
  ncoas <- as.data.frame(ncoas)
  ncoas[is.na(ncoas)] <- 0
  colnames(ncoas) <- c("ID", "N COAs")
  kud_table <- plyr::join(ncoas, kud_table, by="ID", type="left")
  colnames(kud_table)[1] <- id.col

  # create results list
  results <- list("kernel_density"=kernel_density, "summary_table"=kud_table)

  # append the kernel contours to the results list
  results <- c(results, kernel_contours, "spatial.grid"=spatial.grid)

  # print new line
  if(verbose) cat("\n")

  # return the results list, containing the kernel density, KUD table, and contours
  return(results)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
