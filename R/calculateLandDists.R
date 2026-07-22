################################################################################
# Function to calculate distance to nearest shore/land feature ################
################################################################################
#'
#' Calculate distance to nearest shore/land feature and related movement metrics.
#'
#' @description This function calculates the distance to the nearest coastline or land feature for each
#'  animal position. Instead of computing pairwise distances between all points and coastline vertices,
#' it uses a rasterised land mask and a pre-computed distance-to-land raster, providing major gains in
#' speed and memory efficiency for large datasets. If step lengths between consecutive positions are
#' supplied, the function also quantifies movement direction relative to land (inshore, offshore,
#' alongshore).
#'
#' @inheritParams as_moby
#' @param data A data frame containing animal positions.
#' @param land.shape An sf object representing landmass.
#' @param epsg.code Integer. The projected EPSG code (must use meters). If not supplied,
#' it is taken from the `mobyData` metadata (see \code{\link{as_moby}}); a projected
#' (metric) CRS is required.
#' @param mov.threshold Numeric (0-1). Proportion of movement perpendicular to
#' shore required for inshore/offshore classification.
#' @param dist.col Optional. Column name for step length (meters).
#' @param grid.resolution Numeric. Resolution of the distance raster in meters.
#'
#' @return A data frame with added spatial and movement columns.
#'
#' @examples
#' \donttest{
#' data(rays)
#' # a small land polygon spanning the study area
#' land <- sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
#'   c(-9.05, 38.49), c(-8.90, 38.49), c(-8.90, 38.43),
#'   c(-9.05, 38.49)))), crs = 4326))
#' rays_land <- calculateLandDists(rays, land.shape = land)
#' head(rays_land$land_dist)
#' }
#'
#' @export

calculateLandDists <- function(data,
                               land.shape,
                               epsg.code = NULL,
                               mov.threshold = 0.5,
                               id.col = NULL,
                               lon.col = NULL,
                               lat.col = NULL,
                               dist.col = NULL,
                               grid.resolution = 100) {

  # --- 1. Initial Setup and Validation ---
  start.time <- Sys.time()
  .printConsole("Calculating distances to nearest land (raster method)...")

  # resolve NULL column/CRS arguments from the mobyData metadata (or canonical defaults)
  .args <- .resolveArgs(data, list(id.col=id.col, lon.col=lon.col, lat.col=lat.col, epsg.code=epsg.code))
  id.col <- .args$id.col; lon.col <- .args$lon.col; lat.col <- .args$lat.col; epsg.code <- .args$epsg.code

  # coerce to a plain data.frame for consistent indexing (tibble/data.table safe)
  data <- as.data.frame(data)

  # a projected (metric) CRS is required for distance calculations
  if (is.null(epsg.code)) {
    stop(paste("No 'epsg.code' supplied. Please provide a projected EPSG code (in meters),",
               "either directly or via as_moby()."), call. = FALSE)
  }

  # required columns
  for (col in c(id.col, lon.col, lat.col)) {
    if (!col %in% names(data)) stop(paste0("Column '", col, "' not found in data."), call. = FALSE)
  }

  # Basic check for distance column if provided
  if (!is.null(dist.col) && !(dist.col %in% names(data))) {
    stop(paste("Column", dist.col, "not found in data."), call. = FALSE)
  }

  # --- 2. Spatial Projection ---
  # Convert data to sf and project to the metric system
  pts_sf <- sf::st_as_sf(data, coords = c(lon.col, lat.col), crs = 4326)
  pts_proj <- sf::st_transform(pts_sf, epsg.code)
  
  # Ensure land shape is sf and projected to the same system
  if (!inherits(land.shape, "sf")) land.shape <- sf::st_as_sf(land.shape)
  coast_proj <- sf::st_transform(land.shape, epsg.code)
  
  # --- 3. Raster Distance Calculation ---
  # Convert to terra objects
  coast_vect <- terra::vect(coast_proj)
  
  # Create a template raster based on the coastline extent
  r_template <- terra::rast(terra::ext(coast_vect), 
                            resolution = grid.resolution, 
                            crs = terra::crs(coast_vect))
  
  # Rasterize land: land = 1, water = NA
  # terra::distance calculates distance from NA cells to the nearest non-NA cell
  land_rast <- terra::rasterize(coast_vect, r_template, field = 1)
  dist_rast <- terra::distance(land_rast)
  
  # --- 4. Extraction ---
  pts_vect <- terra::vect(pts_proj)
  
  # terra::extract returns a data.frame [ID, value]
  # We extract the second column to get the actual distance values
  extracted_data <- terra::extract(dist_rast, pts_vect)
  data$land_dist <- round(extracted_data[, 2], 2)
  
  # --- 5. Movement Metrics (Conditional) ---
  if (!is.null(dist.col)) {
    message("Quantifying movement directionality...")
    
    # Split by ID to ensure calculations don't bleed between different animals
    data_list <- split(data, data[[id.col]])
    
    data_list <- lapply(data_list, function(d) {
      # land_diff: Change in distance to shore (+ is offshore, - is inshore)
      # We use lead() logic: distance at T+1 minus distance at T
      land_next <- c(d$land_dist[-1], NA)
      d$land_diff <- round(land_next - d$land_dist, 2)
      
      # land_ratio: Proportion of the total step length that was moving toward/away from shore
      d$land_ratio <- round(abs(d$land_diff) / d[[dist.col]], 2)
      
      # Handle cases where movement was 0 (avoid division by zero/NaN)
      d$land_ratio[d[[dist.col]] == 0] <- NA
      
      # Directional Classification
      d$direction <- "alongshore" # Default state
      d$direction[d$land_ratio >= mov.threshold & d$land_diff > 0] <- "offshore"
      d$direction[d$land_ratio >= mov.threshold & d$land_diff < 0] <- "inshore"
      
      # Cleanup NAs for static points
      d$direction[is.na(d$land_ratio)] <- NA
      
      return(d)
    })
    
    data <- do.call(rbind, data_list)
    rownames(data) <- NULL
    data$direction <- factor(data$direction, levels = c("inshore", "alongshore", "offshore"))
  }
  
  # --- 6. Execution Summary ---
  time.taken <- Sys.time() - start.time
  cat(paste("Done! Total execution time:", 
            round(as.numeric(time.taken), 2), units(time.taken), "\n"))
  
  return(data)
}