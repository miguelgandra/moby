#######################################################################################################
# Calculate UDs - Kernel Utilization Distributions  ##################################################
#######################################################################################################

#' Calculate animals' kernel utilization areas

#' @description This function estimates animal utilization distributions (UDs) and home-range
#' areas. By default it uses **autocorrelated kernel density estimation (AKDE)** via the
#' `ctmm` package (`method = "akde"`), which fits a continuous-time movement model to account
#' for the serial autocorrelation inherent in tracking data and returns area estimates with
#' confidence intervals. A classic (IID) fixed-bandwidth kernel density estimator
#' (`method = "kde"`, via \code{\link[adehabitatHR]{kernelUD}}) is also available for speed or
#' backward compatibility. Utilization areas are returned for the specified isopleths
#' (default 50% and 95%); for display, isopleth polygons can be clipped to land, and UDs can
#' be estimated independently for different groups or time periods.
#'
#' @inheritParams as_moby
#' @param data A data frame with animal positions (e.g. COAs), containing longitude/latitude
#' (or projected x/y) columns and, for `method = "akde"`, a POSIXct time column.
#' @param method Estimation method: `"akde"` (default; autocorrelated KDE via `ctmm`, with
#' confidence intervals) or `"kde"` (classic fixed-bandwidth KDE via `adehabitatHR`).
#' @param bandwidth Numeric smoothing parameter (h) for `method = "kde"` only; ignored for AKDE
#' (which estimates smoothing from the fitted movement model). Larger values produce smoother,
#' broader distributions; smaller values give more localized, potentially fragmented estimates.
#' @param timebin.col Name of the POSIXct time column used by `method = "akde"` to model temporal
#' autocorrelation. The canonical input is a time-binned series (e.g. the `timebin` column of
#' \code{\link{calculateCOAs}} output), but a raw date-time column is also accepted. If `NULL`, the
#' function uses the `mobyData` time-bin (then date-time) column and, failing that, a `"timebin"` or
#' `"datetime"` column present in the data.
#' @param model.selection For `method = "akde"`: `"fit"` (default) fits a single movement model
#' from an automated guess (faster); `"select"` runs \code{ctmm::ctmm.select} to choose among
#' candidate models (more thorough, slower).
#' @param spatial.grid Optional. A `Raster` or `SpatialPixels` object representing the
#' grid over which the animal kernel utilization distributions (UDs) will be estimated
#' (see the `grid` argument in \code{\link[adehabitatHR]{kernelUD}}). If set to `NULL`,
#' the function will automatically generate an appropriate grid based on the spatial extent
#' of the supplied animal's positions.
#' @param subset Optional. A variable used to subset the data, allowing UDs to be calculated independently
#' for each level of this variable. This should be the name of a column in the provided dataset.
#' If left `NULL`, UDs are calculated for the whole monitoring period.
#' @param id.groups Optional. A named list where each element represents a group (e.g., species, sex, or age class),
#' containing a vector of IDs for that group. If supplied, UDs will be calculated
#' independently for each group. The names of the list correspond to the group labels.
#' @param land.shape Optional. A shapefile containing coastlines or landmasses,
#' provided either as an 'sf' object or as a 'SpatialPolygonsDataFrame'/'SpatialPolygons' object.
#' If the input is not in 'sf' format, the function will automatically convert it to 'sf'
#' to ensure compatibility with subsequent spatial operations. Used to clip and exclude
#' any portions of the estimated areas that overlap with landmasses.
#' @param contour.percent Numeric vector. The percentages for which isopleths (contour areas)
#' are calculated. Defaults to 50% and 95%, representing core and total areas of utilization.
#' @param verbose Logical. If TRUE, the function will print detailed processing
#' information. Defaults to \code{getOption("moby.verbose", TRUE)}.
#'
#' @return
#' A list containing:
#' \item{ud}{The estimated utilization distributions: an `adehabitatHR` "estUDm"
#' object for `method = "kde"`, or a named list of `ctmm` UD objects for `method = "akde"`.}
#' \item{summary_table}{A data frame summarizing the (point-estimate) area for each contour
#' (isopleth) per individual, in km2. Suitable for \code{\link{movementTable}}.}
#' \item{`K50`, `K95`, ...}{One 'sf' object per requested `contour.percent`, named
#' `paste0("K", contour.percent)`, holding the isopleth polygons.
#' For `method = "akde"` each individual contributes three polygons tagged by a `ci` column
#' (`"low"`, `"est"`, `"high"`) - the confidence envelope of the isopleth, which
#' \code{\link{plotMaps}} can draw. If a subset/group was provided, a column distinguishes the groups.}
#' \item{area_estimates}{(`method = "akde"` only) A tidy data frame of area estimates with
#' lower/upper confidence limits (in km2) and the effective sample size (DOF) for each
#' individual and contour.}
#'
#' The results list also contains multiple attributes to store relevant metadata,
#' such as function options and processing details. These attributes might be useful
#' for tracking parameters and ensuring reproducibility of the analysis.
#'
#'
#' @details
#'
#' This function estimates kernel utilization distributions (UDs) based on animal location data,
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
#' kernel utilization distribution (UD) analysis, representing the standard deviation of the kernel.
#' It defines the extent to which a location can influence the home range estimation and
#' the overall density estimate. The choice of bandwidth significantly impacts the results
#' of the analysis. The bandwidth can either be fixed (using a single value for all data points)
#' or variable (adapting based on point density). The `method = "kde"` pathway uses a fixed bandwidth;
#' the default `method = "akde"` instead estimates smoothing from a fitted continuous-time movement
#' model and is generally preferred for autocorrelated tracking data (Fleming et al. 2015).
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
#' topography in space use estimation from node-based telemetry data.
#' Methods in Ecology and Evolution, 11(12), 1733-1742.
#' https://doi.org/10.1111/2041-210X.13484
#'
#' Worton, B. J. (1989).
#' Kernel methods for estimating the utilization distribution in home-range studies.
#' Ecology, 70(1), 164-168.
#' https://doi.org/10.2307/1938423
#'
#' Fleming, C. H., Fagan, W. F., Mueller, T., Olson, K. A., Leimgruber, P., & Calabrese, J. M. (2015).
#' Rigorous home range estimation with movement data: a new autocorrelated kernel density estimator.
#' Ecology, 96(5), 1182-1188.
#' https://doi.org/10.1890/14-2010.1
#'
#' @seealso \code{\link[adehabitatHR]{kernelUD}}, \code{\link[ctmm]{akde}}
#'
#' @examples
#' \donttest{
#' data(rays)
#' if (requireNamespace("adehabitatHR", quietly = TRUE)) {
#'   # a coarse estimation grid keeps this example fast
#'   grid <- terra::rast(terra::ext(-9.05, -8.95, 38.43, 38.48),
#'                       ncol = 60, nrow = 60, crs = "EPSG:4326")
#'   terra::values(grid) <- 0
#'   grid <- terra::project(grid, "EPSG:32629")
#'   kud <- calculateUDs(rays, method = "kde", bandwidth = 500,
#'                        spatial.grid = grid)
#'   kud$summary_table
#' }
#' }
#'
#' @export


################################################################################
## Main function - Estimate Kernel Density #####################################

calculateUDs <- function(data,
                         id.col = NULL,
                         timebin.col = NULL,
                         lon.col = NULL,
                         lat.col = NULL,
                         id.groups = NULL,
                         subset = NULL,
                         land.shape = NULL,
                         epsg.code = NULL,
                         spatial.grid = NULL,
                         bandwidth = NULL,
                         method = c("akde", "kde"),
                         contour.percent = c(50,95),
                         model.selection = c("fit", "select"),
                         verbose = getOption("moby.verbose", TRUE)) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  method_missing <- missing(method)
  method <- match.arg(method)
  model.selection <- match.arg(model.selection)

  # graceful default: if the user did not explicitly request a method and 'ctmm' (needed for the
  # default AKDE) is unavailable, fall back to classic KDE with a clear warning rather than failing
  # on a first run. An explicitly requested method='akde' still errors below if ctmm is missing.
  if(method_missing && method=="akde" && !requireNamespace("ctmm", quietly=TRUE) &&
     requireNamespace("adehabitatHR", quietly=TRUE)){
    warning(paste(strwrap(paste("- The 'ctmm' package (required for the default method='akde') is not",
      "installed; falling back to method='kde'. Install 'ctmm' for the recommended autocorrelated KDE,",
      "or pass method='kde' explicitly to silence this warning."), width=getOption("width")),
      collapse="\n"), call.=FALSE)
    method <- "kde"
  }

  # capture mobyData metadata (before the data is coerced to a plain data.frame) so the
  # AKDE pathway can locate the time column even when not explicitly supplied
  .meta <- attr(data, "moby")
  akde_time <- timebin.col
  if (is.null(akde_time) && !is.null(.meta)) {
    akde_time <- if (!is.null(.meta$timebin.col)) .meta$timebin.col else .meta$datetime.col
  }

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
  if(method=="kde" && !requireNamespace("adehabitatHR", quietly=TRUE)) errors <- c(errors, "The 'adehabitatHR' package is required for method='kde' but is not installed. Please install it using install.packages('adehabitatHR') and try again.")
  if(method=="akde" && !requireNamespace("ctmm", quietly=TRUE)) errors <- c(errors, "The 'ctmm' package is required for method='akde' (the default) but is not installed. Please install it using install.packages('ctmm'), or use method='kde'.")
  if(method=="kde" && (is.null(bandwidth) || !is.numeric(bandwidth))) errors <- c(errors, "Error: method='kde' requires a numeric 'bandwidth' value.")
  if(!inherits(contour.percent, "numeric")) errors <- c(errors, "Error: 'contour.percent' must be a numeric value.")
  if(any(contour.percent < 1 | contour.percent > 100)) errors <- c(errors, "Error: 'contour.percent' must be between 1 and 100.")
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }
  if(method=="akde" && !is.null(bandwidth)){
    message("Note: 'bandwidth' is ignored when method='akde' (AKDE estimates smoothing from the fitted movement model).")
  }

  # dispatch to the autocorrelated KDE pipeline (default)
  if(method=="akde"){
    return(.calculateUDs_akde(data=data, id.col=id.col, time.col=akde_time,
                               lon.col=lon.col, lat.col=lat.col, epsg.code=epsg.code,
                               subset=subset, id.groups=id.groups, land.shape=land.shape,
                               land_shape_name=land_shape_name, contour.percent=contour.percent,
                               model.selection=model.selection, start.time=start.time, verbose=verbose))
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

  # convert a user-supplied raster grid to SpatialPixels if required (adehabitatHR needs sp)
  if(!is.null(spatial.grid)){
    if(inherits(spatial.grid, "SpatRaster")){
      spatial.grid <- methods::as(sf::as_Spatial(sf::st_as_sf(terra::as.points(spatial.grid[[1]]))), 'SpatialPixels')
    } else if(inherits(spatial.grid, "Raster")){
      spatial.grid <- methods::as(spatial.grid, 'SpatialPixels')
    }
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
      coords[[paste0(subset, ".by.ids")]] <- interaction(coords[[subset]], coords[["id.groups"]], drop=TRUE)
      subset <- paste0(subset, ".by.ids")
    }
  }

  # check if UDs need to be calculated for different groups
  multiple <- !is.null(subset)


  ##############################################################################
  ## Calculate UDs ############################################################
  ##############################################################################

  #######################################################################
  # if no subset provided, calculate the UDs for the whole dataset #####
  if(!multiple) {

    # print status message to console
    .mobyInform("Estimating kernel utilization distributions...", verbose = verbose)

    # compute UDs for the entire dataset
    final_results <- .computeUDs(coords, id.col, bandwidth, coords_bbox, contour.percent, spatial.grid, land.shape, epsg.code, verbose)

    # extract the 'spatial.grid' from the final results
    spatial.grid <- final_results$spatial.grid

    # remove the 'spatial.grid' element from the final results list
    final_results <- final_results[-length(final_results)]


  #######################################################################
  # else, split the data and calculate UDs for each group separately ###
  } else{

    # print the grouping criteria (defined by 'subset') to the console
    .mobyInform("Grouping data by: ", paste(subset, collapse=" / "), verbose = verbose)

    # split the data into groups based on the 'subset' variable
    group_coords <- split(coords, f=coords[[subset]], drop=TRUE)
    group_coords <- lapply(group_coords, function(x){x[[id.col]] <- droplevels(x[[id.col]]); return(x)})

    # check if any individuals (IDs) are repeated across the different groups
    repeated_ids <- list()
    for(i in seq_along(group_coords)){
      # get the IDs within the current group
      subset_ids <- levels(group_coords[[i]][[id.col]])
      # get the IDs present in all other groups
      remaining_ids <- levels(do.call("rbind", group_coords[-i])[[id.col]])
      # check if any ID from the current group is found in the other groups
      repeated_ids[[i]] <- any(subset_ids %in% remaining_ids)
    }
    repeated_ids <- any(unlist(repeated_ids))

    # compute UDs for each group separately
    kud_results <- lapply(seq_along(group_coords), function(i) {
      .mobyInform("Estimating kernel utilization distributions [", names(group_coords)[i], "]...", verbose = verbose)
      group_results <- .computeUDs(group_coords[[i]], id.col, bandwidth, coords_bbox, contour.percent, spatial.grid, land.shape, epsg.code, verbose)
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
    final_results <- list("ud"=lapply(kud_results, function(x) x$ud),
                         "summary_table"=summary_table)
    final_results <- c(final_results, contour_results)

  }

  ##############################################################################
  ## Wrap it up ################################################################
  ##############################################################################

  # print the total execution time to the console
  .reportRuntime(start.time, verbose)

  # create attributes to save relevant metadata
  attr(final_results, 'method') <- "kde"
  attr(final_results, 'id.groups') <- id.groups
  attr(final_results, 'bandwidth') <- bandwidth
  attr(final_results, 'contour.percent') <- contour.percent
  attr(final_results, 'subset') <- subset
  attr(final_results, 'land.shape') <- land_shape_name
  attr(final_results, 'epsg.code') <- epsg.code
  attr(final_results, 'grid.extent') <- spatial.grid@bbox
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

.subtractLand <- function(uds, land.shape, epsg.code, verbose) {

  # extract grid coordinates from the kernel density object
  grid_coords <- as.data.frame(uds[[1]]@coords)
  colnames(grid_coords) <- c("x","y")

  # create an sf object from kernel density coordinates
  grid_coords <- sf::st_as_sf(grid_coords, coords=c("x", "y"), crs=epsg.code)

  # check for overlaps with the land.shape sf object
  overlap_indexes <- as.logical(sf::st_intersects(grid_coords, sf::st_union(land.shape), sparse=FALSE))

  # initialize a counter for the number of individuals with overlapping areas
  n_corrected <- 0

  # iterate over each kernel density object
  for(i in seq_along(uds)){

    # access the density values from the current kernel density object
    density_values <- uds[[i]]$ud

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
    uds[[i]]$ud <- density_values
  }

  # print a message indicating how many kernel densities were corrected
  if(n_corrected==0){
    .mobyInform("No land overlap detected", verbose = verbose)
  }else{
    kud_label <- ifelse(n_corrected==1, "UD", "UDs")
    .mobyInform("Land overlap areas successfully subtracted from ", n_corrected, " ", kud_label, verbose = verbose)
  }

  # return the corrected kernel densities
  return(uds)
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

  # build the sp SpatialPixels grid directly (adehabitatHR::kernelUD requires sp SpatialPixels; terra
  # cannot produce them). Dims are recomputed from the FINAL cell size; anchored at the bottom-left
  # cell centre. No raster dependency.
  ncol <- as.integer(ceiling(bbox_width  / rounded_cell_size))
  nrow <- as.integer(ceiling(bbox_height / rounded_cell_size))
  gt <- sp::GridTopology(cellcentre.offset = c(as.numeric(grid_bbox["xmin"]) + rounded_cell_size/2,
                                               as.numeric(grid_bbox["ymin"]) + rounded_cell_size/2),
                         cellsize = c(rounded_cell_size, rounded_cell_size),
                         cells.dim = c(ncol, nrow))
  sg <- sp::SpatialGrid(gt, proj4string = sp::CRS(SRS_string = epsg.code$wkt))

  # convert to SpatialPixels class
  return(methods::as(sg, 'SpatialPixels'))
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

.computeUDs <- function(coords, id.col, bandwidth, coords.bbox, contour.percent,
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
  # start the repeat loop to compute UDs and gradually expand the spatial grid if needed
  repeat {

    # increase the grid extent incrementally by a factor of 0.05 (5%)
    expand_factor <- expand_factor + 0.05
    if(expand_factor>0.1) .mobyInform("Expanding grid bounding box by: ", round(expand_factor*100), "%", verbose = verbose)

    # expand bounding box by a given % in all directions
    # create a custom spatial grid based on the coordinate extent, if not supplied by the user
    if(!grid_supplied) spatial.grid <-.createGrid(coords.bbox, expand.factor=expand_factor, epsg.code)

    # estimate the kernel utilization distributions (UD) for the provided coordinates
    ud <- adehabitatHR::kernelUD(filtered_coords[,id.col], h=bandwidth, grid=spatial.grid)

    # clip out areas of the kernel density that overlap with landmasses
    if(!is.null(land.shape)){
      .mobyInform("Clipping out land masses...", verbose = verbose)
      ud <- .subtractLand(ud, land.shape, epsg.code, verbose)
    }

    # initialize a list to store kernel contours for each contour percentage (e.g., 50%, 95%)
    kernel_contours <- vector("list", length(contour.percent))
    names(kernel_contours) <- paste0("K", contour.percent)

    # loop through each contour percentage and calculate corresponding UD contours
    for(c in seq_along(contour.percent)){
      # print progress message
      .mobyInform("Calculating ", contour.percent[c], "% contours...", verbose = verbose)
      # attempt to calculate kernel contours for the specified contour percent
      kernel_contours[[c]] <- tryCatch({
        adehabitatHR::getverticeshr(ud, percent=contour.percent[c], unin="m", unout="km2")
      }, error = function(e) {
        # check if grid extent needs to be increased
        if (grepl("grid is too small", e$message, fixed=TRUE)) {
          if(grid_supplied) .mobyAbort("The supplied 'spatial.grid' is too small for the requested ",
                                       "contour; enlarge its extent or omit it to auto-generate one. (",
                                       e$message, ")")
          else {
            .mobyInform("Grid too small; generating a larger one...", verbose = verbose)
            return(NA)
          }
        } else {
          .mobyAbort("Kernel utilization distribution estimation failed: ", e$message)
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

  # set column names of the UD table
  colnames(kud_table) <- c("ID",  sprintf("UD %d%% (Km2)", contour.percent))

  # add the number of COAs to the UD table
  ncoas <- table(coords[[id.col]])
  ncoas <- as.data.frame(ncoas)
  ncoas[is.na(ncoas)] <- 0
  colnames(ncoas) <- c("ID", "N COAs")
  kud_table <- .joinKeep(ncoas, kud_table, by="ID", type="left")
  colnames(kud_table)[1] <- id.col

  # create results list
  results <- list("ud"=ud, "summary_table"=kud_table)

  # append the kernel contours to the results list
  results <- c(results, kernel_contours, "spatial.grid"=spatial.grid)

  # return the results list, containing the kernel density, UD table, and contours
  return(results)
}


################################################################################
# AKDE pipeline (autocorrelated kernel density estimation, via ctmm) ###########
################################################################################

#' Fit an autocorrelated kernel density estimate for a single individual
#'
#' @description Internal helper. Builds a `ctmm` telemetry object, fits a continuous-time
#' movement model (accounting for autocorrelation) and returns the AKDE utilization
#' distribution. Returns `NULL` if the fit fails.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.akdeModel <- function(times, lon, lat, id, proj4, tz, model.selection, verbose) {
  tel_df <- data.frame(timestamp = times, longitude = lon, latitude = lat,
                       individual.local.identifier = as.character(id))
  run <- function() {
    tel <- ctmm::as.telemetry(tel_df, projection = proj4, timezone = tz)
    guess <- ctmm::ctmm.guess(tel, interactive = FALSE)
    fit <- if (model.selection == "select") ctmm::ctmm.select(tel, guess) else ctmm::ctmm.fit(tel, guess)
    list(tel = tel, fit = fit)
  }
  # ctmm's own fitting chatter is an implementation detail, always suppressed; moby narrates progress
  # via the caller's progress bar instead (see the pass-1 loop in .calculateUDs_akde).
  tryCatch(
    suppressMessages(suppressWarnings(run())),
    error = function(e) {
      .mobyWarn("AKDE model fit failed for individual '", id, "': ", conditionMessage(e))
      NULL
    }
  )
}

#' AKDE utilization distribution from a fitted telemetry + model (single-individual / fallback path)
#' @keywords internal
#' @noRd
.akdeUD <- function(tel, fit, verbose) {
  ud <- tryCatch(
    suppressMessages(suppressWarnings(ctmm::akde(tel, fit))),
    error = function(e) NULL
  )
  # ctmm returns the fit object (not a UD) when it cannot form a usable distribution
  if (!is.null(ud) && !inherits(ud, "UD")) return(NULL)
  ud
}


#' Autocorrelated kernel density estimation pipeline
#'
#' @description Internal helper implementing the `method = "akde"` pathway of
#' \code{\link{calculateUDs}}: per-individual continuous-time movement-model fitting and
#' AKDE estimation (via the `ctmm` package), returning area estimates with confidence
#' intervals, isopleth polygons (optionally clipped to land for display) and a summary table
#' compatible with \code{\link{movementTable}}.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.calculateUDs_akde <- function(data, id.col, time.col, lon.col, lat.col, epsg.code,
                                subset, id.groups, land.shape, land_shape_name,
                                contour.percent, model.selection, start.time, verbose) {

  ## ---- resolve a POSIXct time column (AKDE requires timestamps) -------------
  candidates <- unique(c(time.col, "timebin", "datetime"))
  candidates <- candidates[candidates %in% colnames(data)]
  time_col <- NULL
  for (cc in candidates) if (inherits(data[[cc]], "POSIXct")) { time_col <- cc; break }
  if (is.null(time_col)) {
    stop(paste("method='akde' requires a POSIXct time column (it models temporal",
               "autocorrelation). Please supply 'timebin.col' (e.g. the time-bin column",
               "of your COAs), or use method='kde'."), call. = FALSE)
  }
  tz <- .dataTZ(data[[time_col]])

  ## ---- project coordinates --------------------------------------------------
  coords <- sf::st_as_sf(data, coords = c(lon.col, lat.col))
  spatial_data <- .processSpatial(coords, land.shape, epsg.code)
  coords <- spatial_data$coords
  land.shape <- spatial_data$spatial.layer
  epsg.code <- spatial_data$epsg.code
  proj4 <- sf::st_crs(epsg.code)$proj4string
  # geographic coordinates for ctmm::as.telemetry (which projects internally to proj4)
  coords_wgs <- sf::st_transform(coords, 4326)
  lonlat <- sf::st_coordinates(coords_wgs)
  data$.lon_wgs <- lonlat[, 1]
  data$.lat_wgs <- lonlat[, 2]
  data$.time <- data[[time_col]]

  ## ---- define grouping (id.groups / subset), mirroring the KDE pathway ------
  gcol <- NULL
  if (!is.null(id.groups)) {
    grp_map <- data.frame(.id = as.character(unlist(id.groups)),
                          .grp = rep(names(id.groups), lengths(id.groups)),
                          stringsAsFactors = FALSE)
    data$.group <- grp_map$.grp[match(as.character(data[[id.col]]), grp_map$.id)]
    if (is.null(subset)) { gcol <- ".group" }
    else { data$.gcomb <- interaction(data[[subset]], data$.group, drop = TRUE); gcol <- ".gcomb" }
  } else if (!is.null(subset)) {
    gcol <- subset
  }

  units <- if (is.null(gcol)) list(all = data) else split(data, droplevels(factor(data[[gcol]])))

  ## ---- per-unit, per-individual AKDE ----------------------------------------
  ud_list <- list()
  summary_rows <- list()
  area_rows <- list()
  contour_sf <- stats::setNames(vector("list", length(contour.percent)), paste0("K", contour.percent))

  .mobyInform("Fitting AKDE movement models (this may take a while for many individuals)...", verbose = verbose)

  for (u in names(units)) {
    udata <- units[[u]]
    ids <- unique(as.character(udata[[id.col]]))

    # pass 1: fit a continuous-time movement model per individual (>= 5 locations required).
    # A per-individual progress bar reassures the user through what is the package's slowest step.
    tels <- list(); fits <- list(); nloc <- list()
    pb <- .progressBar(length(ids), verbose)
    for (i in seq_along(ids)) {
      id <- ids[i]
      .progressSet(pb, i)
      d_id <- udata[udata[[id.col]] == id, ]
      n_loc <- nrow(d_id)
      if (n_loc < 5) {
        .mobyWarn("Skipping individual '", id, "'", if (!is.null(gcol)) paste0(" [", u, "]") else "",
                  ": at least 5 locations are required for AKDE (found ", n_loc, ").")
        next
      }
      mdl <- .akdeModel(d_id$.time, d_id$.lon_wgs, d_id$.lat_wgs, id, proj4, tz, model.selection, verbose)
      if (is.null(mdl)) next
      tels[[id]] <- mdl$tel; fits[[id]] <- mdl$fit; nloc[[id]] <- n_loc
    }
    .progressEnd(pb)
    if (length(tels) == 0) next

    # pass 2: estimate the AKDE UDs. With more than one individual, ctmm's list interface places all
    # UDs on a COMMON grid, so pairwise overlap (calculateUDOverlap) is well defined; a lone individual
    # uses akde() directly. Falls back to per-individual estimation if the joint call fails.
    if (length(tels) == 1) {
      uds <- stats::setNames(list(.akdeUD(tels[[1]], fits[[1]], verbose)), names(tels))
    } else {
      uds <- tryCatch(
        suppressMessages(suppressWarnings(ctmm::akde(tels, fits))),
        error = function(e) stats::setNames(lapply(seq_along(tels),
                              function(i) .akdeUD(tels[[i]], fits[[i]], verbose)), names(tels)))
      if (is.null(names(uds))) names(uds) <- names(tels)
    }

    # pass 3: per-UD extraction (area CIs, isopleths, summary)
    for (id in names(uds)) {
      ud <- uds[[id]]
      if (is.null(ud) || !inherits(ud, "UD")) {
        warning(paste0("- AKDE could not produce a reliable utilization distribution for individual '",
                       id, "' (degenerate model fit); individual skipped. Consider more data or method='kde'."),
                call. = FALSE)
        next
      }
      n_loc <- nloc[[id]]
      label <- if (is.null(gcol)) NA_character_ else u
      ud_list[[paste0(u, "::", id)]] <- ud

      # area estimates with confidence intervals, per contour
      areas_est <- stats::setNames(numeric(length(contour.percent)), paste0("UD ", contour.percent, "% (Km2)"))
      for (k in seq_along(contour.percent)) {
        p <- contour.percent[k]
        ci <- summary(ud, level.UD = p / 100, units = FALSE)$CI
        ci_km2 <- ci[1, c("low", "est", "high")] / 1e6
        areas_est[k] <- ci_km2[["est"]]
        area_rows[[length(area_rows) + 1]] <- data.frame(
          ID = id, group = label, contour = p,
          area_est = ci_km2[["est"]], area_low = ci_km2[["low"]], area_high = ci_km2[["high"]],
          units = "km2", DOF_area = unname(summary(ud, level.UD = p / 100)$DOF["area"]),
          stringsAsFactors = FALSE, check.names = FALSE)

        # isopleth polygons (low / estimate / high CI contours), optionally clipped to land for display.
        # ctmm::as.sf returns all three CI contours per level; keep them (tagged in a 'ci' column) so the
        # AKDE confidence envelope can be drawn by plotMaps, instead of discarding the low/high bounds.
        poly <- tryCatch(ctmm::as.sf(ud, level.UD = p / 100), error = function(e) NULL)
        if (!is.null(poly) && nrow(poly) > 0) {
          ci_lab <- sub(".*\\s", "", as.character(poly[[1]]))          # trailing token: "low"/"est"/"high"
          ord <- stats::na.omit(match(c("low", "est", "high"), ci_lab))
          if (length(ord) == 0) ord <- seq_len(nrow(poly))
          poly_sf <- sf::st_sf(ci = ci_lab[ord], geometry = sf::st_geometry(poly)[ord])
          suppressWarnings(sf::st_crs(poly_sf) <- sf::st_crs(epsg.code))
          if (!is.null(land.shape)) {
            poly_sf <- tryCatch(suppressWarnings(sf::st_difference(poly_sf, sf::st_union(sf::st_geometry(land.shape)))),
                                error = function(e) poly_sf)
          }
          if (nrow(poly_sf) > 0) {
            poly_sf[[id.col]] <- id
            if (!is.null(gcol)) poly_sf[[if (gcol %in% c(".group", ".gcomb")) "group" else gcol]] <- label
            kname <- paste0("K", p)
            contour_sf[[kname]][[length(contour_sf[[kname]]) + 1]] <- poly_sf
          }
        }
      }

      srow <- data.frame(ID = id, stringsAsFactors = FALSE, check.names = FALSE)
      if (!is.null(gcol)) srow$group <- label
      srow[["N COAs"]] <- n_loc
      for (k in seq_along(contour.percent)) srow[[names(areas_est)[k]]] <- sprintf("%.2f", areas_est[k])
      summary_rows[[length(summary_rows) + 1]] <- srow
    }
  }

  if (length(summary_rows) == 0) stop("AKDE could not be estimated for any individual (insufficient data or failed fits).", call. = FALSE)

  ## ---- assemble outputs -----------------------------------------------------
  summary_table <- do.call(.rbindFill, summary_rows)
  colnames(summary_table)[1] <- id.col
  area_estimates <- do.call(rbind, area_rows)
  colnames(area_estimates)[1] <- id.col

  contour_out <- lapply(contour_sf, function(x) if (length(x) > 0) do.call(rbind, x) else NULL)
  contour_out <- contour_out[!vapply(contour_out, is.null, logical(1))]

  final_results <- c(list(ud = ud_list, summary_table = summary_table),
                     contour_out, list(area_estimates = area_estimates))

  .reportRuntime(start.time, verbose)

  attr(final_results, 'method') <- "akde"
  attr(final_results, 'bandwidth') <- NA
  attr(final_results, 'model.selection') <- model.selection
  attr(final_results, 'id.groups') <- id.groups
  attr(final_results, 'subset') <- subset
  attr(final_results, 'contour.percent') <- contour.percent
  attr(final_results, 'land.shape') <- land_shape_name
  attr(final_results, 'epsg.code') <- epsg.code
  attr(final_results, 'processing.date') <- Sys.time()

  return(final_results)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
