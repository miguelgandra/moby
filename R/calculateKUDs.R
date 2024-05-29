#######################################################################################################
# Calculate KUDs - Kernel Utilization Distributions  ##################################################
#######################################################################################################

#' Calculate animals' kernel utilization areas

#' @description Pipeline to calculate animals' utilization areas using kernel
#' density estimation, based on \code{\link[adehabitatHR]{kernelUD}} function. It estimates
#' areas corresponding to 50% and 95% of occurrence probability (automatically subtracting
#' portions over land masses) and allows estimation of KUDs independently by group.
#'
#' @param data A data frame with animal positions (COAs), containing 'longitude' and 'latitude' columns.
#' @param subset If defined, KUDs are calculated independently for each level of this variable.
#' If left NULL, single KUDs are calculated for the whole monitoring period.
#' @param bandwidth Kernel bandwidth (h) - smoothing parameter used to calculate KUDs.
#' @param grid SpatialPixels object containing the grid over which animal KUDs are going to be estimated
#' (see grid argument in the \code{\link[adehabitatHR]{kernelUD}} function).
#' @param land.shape A shape file containing coastlines.
#' @param id.col Name of the column containing animal IDs. Defaults to "ID".
#' @param lon.col Name of the column containing projected longitudes. Defaults to "lon".
#' @param lat.col Name of the column containing projected latitudes. Defaults to "lat".
#' @return A list containing the estimated UDs, areas corresponding to
#' 50% and 95% of occurrence probability, and a formatted data frame.
#' @seealso \code{\link{plotTemporalKUDs}}
#' @seealso \code{\link[adehabitatHR]{kernelUD}}
#' @export


################################################################################
## Main function - Estimate Kernel Density #####################################

calculateKUDs <- function(data, subset=NULL, bandwidth, grid, land.shape,
                          id.col="ID", lon.col="lon", lat.col="lat") {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # check if adehabitatHR package is installed
  if (!requireNamespace("adehabitatHR", quietly=TRUE)) {
    stop("The 'adehabitatHR' package is required for this function but is not installed. Please install 'adehabitatHR' using install.packages('adehabitatHR') and try again.")
  }

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # check if data contains lon.col and lat.col
  if(any(!c(lon.col, lat.col) %in% colnames(data))){
    stop("Longitude/latitude columns not found. Please assign the correct column names with 'lon.col' and 'lat.col'")
  }

  # check if dataset coordinates are projected
  geographic_coords <-  all(data[,lon.col]>=(-180) & data[,lon.col]<=180 & data[,lat.col]>=(-90) & data[,lat.col]<=90)
  if(is.na(geographic_coords)){
    stop("Longitudes/latitudes values are missing and/or in the wrong format")
  }else if(geographic_coords==T & !grepl("+proj=longlat +datum=WGS84", raster::projection(grid), fixed=T)){
    stop("Longitudes/latitudes and supplied UD grid appear to have different projections")
  }

  # check if data contains subset col
  if(!is.null(subset)){
    if(!subset %in% colnames(data)){
      stop(paste(subset, "column not found in the supplied data"))
    }
  }

  ##############################################################################
  ## Calculate KUDs ############################################################
  ##############################################################################

  # check if KUDs need to be calculated for different groups
  multiple <- !is.null(subset)

  ####################################################################
  # if not, calculate KUDS for the entire study duration #############
  if(multiple==F) {
    cat("Estimating kernel densities for the entire monitoring period...\n")
    # convert data to class SpatialPointsDataFrame
    coord_cols <- which(colnames(data) %in% c(lon.col, lat.col))
    data <- sp::SpatialPointsDataFrame(data[,coord_cols], data[,-coord_cols], proj4string=land.shape@proj4string)
    data <- cleanData(data, id.col)
    # estimate kernel density and calculate areas (km²)
    kernel_density <- adehabitatHR::kernelUD(data[,id.col], h=bandwidth, grid=grid)
    cat("Calculating 50% and 95% contours...\n")
    k50 <- adehabitatHR::getverticeshr(kernel_density, percent=50, unin = "m", unout = "km2")
    k95 <- adehabitatHR::getverticeshr(kernel_density, percent=95, unin = "m", unout = "km2")
    terra::crs(k50) <- terra::crs(land.shape)
    terra::crs(k95) <- terra::crs(land.shape)
    cat("Clipping out land masses...\n")
    k50 <- subtractLand(k50, land.shape=land.shape, out.prefix="KUD 50%")
    k95 <- subtractLand(k95, land.shape=land.shape, out.prefix="KUD 95%")
    # format summary table
    kud_table <- plyr::join(k50@data, k95@data, by="id", type="left")
    colnames(kud_table) <- c("ID", "KUD 50% (km2)", "KUD 95% (km2)")
    kud_table[,2] <- sprintf("%.2f", kud_table[,2])
    kud_table[,3] <- sprintf("%.2f", kud_table[,3])
    fish_ids <- kud_table$ID[!is.na(kud_table[,2]) & !is.na(kud_table[,3])]
    # calculate number of COAs
    ncoas <- table(data@data[,id.col])
    ncoas <- as.data.frame(ncoas)
    ncoas[is.na(ncoas)] <- 0
    colnames(ncoas) <- c("ID", "Nº COAs")
    kud_table <- plyr::join(ncoas, kud_table, by="ID", type="left")

    # create results list
    kud_results <- list("kernel_density"=kernel_density, "k50"=k50, "k95"=k95, "kud_table"=kud_table)
  }

  ####################################################################
  # else, split data and calculate KUDs for each group ###############
  if(multiple==T) {

    # check if subset groups are available in the supplied data
    if(any(!subset %in% colnames(data))) {stop("supplied grouping variables could not be found")}
    cat(paste0("Grouping data by ", paste(subset, collapse=" / "), "\n"))

    # split data
    data_group <- split(data, f=data[,subset], drop=T)
    data_group <- lapply(data_group, function(x){x[,id.col]<-droplevels(x[,id.col]); return(x)})


    # check if subset variable relates to id.groups
    repeated_ids <- list()
    for(i in 1:length(data_group)){
      subset_ids <- levels(data_group[[i]][,id.col])
      remaining_ids <- levels(do.call("rbind", data_group[-i])[,id.col])
      repeated_ids[[i]] <- any(subset_ids %in% remaining_ids)
    }
    repeated_ids <- any(unlist(repeated_ids)==T)


    # convert data to class SpatialPointsDataFrame
    coords_col <- match(c(lon.col, lat.col), colnames(data))
    data_group <- lapply(data_group, function(x) sp::SpatialPointsDataFrame(x[,coords_col], x[,-coords_col], proj4string=land.shape@proj4string))
    data_group <- lapply(data_group, function(x) {x@data[,id.col] <- droplevels(x@data[,id.col]); return(x)})

    # filter data and remove individuals with less than 5 records
    cat("Filtering out individuals with < 5 detections within each group\n")
    data_group <- lapply(data_group, function(x) cleanData(x, id.col=id.col))
    data_group <- data_group[sapply(data_group, nrow)>0]

    # estimate kernel density and calculate areas (km²)
    cat("Estimating kernel densities for each group...\n")
    kernel_density <- lapply(data_group, function(x) adehabitatHR::kernelUD(x[,id.col], h=bandwidth, grid=grid))
    cat("Calculating 50% and 95% contours...\n")
    k50 <- lapply(kernel_density, adehabitatHR::getverticeshr, percent=50, unin="m", unout="km2")
    k95 <- lapply(kernel_density, adehabitatHR::getverticeshr, percent=95, unin="m", unout="km2")
    k50 <- lapply(k50, function(x){terra::crs(x) <- terra::crs(land.shape); return(x)})
    k95 <- lapply(k95, function(x){terra::crs(x) <- terra::crs(land.shape); return(x)})
    cat("Clipping out land masses...\n")
    k50 <- mapply(function(data, name) {subtractLand(data, land.shape, out.prefix=paste("KUD 50%", name))},
                  data=k50, name=names(k50))
    k95 <- mapply(function(data, name) {subtractLand(data, land.shape, out.prefix=paste("KUD 95%", name))},
                  data=k95, name=names(k95))

    # create summary table
    k50_table <- lapply(k50, function(x) data.frame("ID"=x$id, "KUD 50% Area (km2)"=sprintf("%.2f", x$area), check.names=F))
    k95_table <- lapply(k95, function(x) data.frame("ID"=x$id, "KUD 95% Area (km2)"=sprintf("%.2f", x$area), check.names=F))

    # calculate number of COAs used in KUD estimation for each individual in each group
    ncoas <- lapply(data_group, function(x) table(x@data[,id.col]))
    ncoas <- lapply(ncoas, as.data.frame)
    ncoas <- lapply(ncoas, function(x) {colnames(x)<-c("ID", "Nº COAs"); return(x)})

    # merge and format
    if(repeated_ids==T){
      k50_table <- mapply(function(data,subset) {colnames(data)[2] <- paste(colnames(data)[2], "-", subset); return(data)},
                          k50_table, names(k50_table), SIMPLIFY=F)
      k95_table <- mapply(function(data,subset) {colnames(data)[2] <- paste(colnames(data)[2], "-", subset); return(data)},
                          k95_table, names(k95_table), SIMPLIFY=F)
      k50_table <- Reduce(function(x,y) merge(x, y, all=T, by="ID", sort=T), k50_table)
      k95_table <- Reduce(function(x,y) merge(x, y, all=T, by="ID", sort=T), k95_table)
      ncoas <- mapply(function(x, subset) {colnames(x)[2] <- paste(colnames(x)[2], "-", subset); return(x)},
                      x=ncoas, subset=names(ncoas), SIMPLIFY=F)
      ncoas <- Reduce(function(x,y) merge(x, y, all=T, by="ID", sort=T), ncoas)
      ncoas[is.na(ncoas)] <- 0

    }else{
      k50_table <- do.call("rbind", k50_table)
      k95_table <- do.call("rbind", k95_table)
      ncoas <- do.call("rbind", ncoas)
    }
    kud_table <- plyr::join(k50_table, k95_table, by="ID", type="full")
    kud_table <- plyr::join(kud_table, ncoas, by="ID", type="full")


    # compose final list
    kud_results <- list("kernel_density"=kernel_density, "k50"=k50, "k95"=k95,
                        "kud_table"=kud_table)
  }

  # return results
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(kud_results)
}


################################################################################
# Helper function I - remove individuals with less than 5 records ##############
# required to avoid errors in the kernel_UD function ###########################
################################################################################

cleanData <- function(data, id.col) {

  data@data[,id.col] <- droplevels(data@data[,id.col])
  for (i in levels(data@data[,id.col])) {
    coords <- data@coords[data@data[,id.col]==i,]
    if (nrow(coords) < 5 || is.null(nrow(coords))) {
      data <- subset(data, data@data[,id.col]!=i)
    }
  }
  data@data[,id.col] <- droplevels(data@data[,id.col])

  return(data)
}


################################################################################
# Helper function II - subtract area on land ###################################
################################################################################

subtractLand <- function(kud.polygons, land.shape, out.prefix="") {

  overlaps_number <- 0
  for (i in 1:length(kud.polygons)){

    # check if kernel utilization distribution overlaps with land
    polygon <- sf::st_as_sf(kud.polygons[i,])
    land.shape <-  sf::st_as_sf(land.shape)
    in_land <- lengths(sf::st_intersects(polygon, land.shape, sparse=T))>0

    # if not, ignore and go to next
    if(in_land==F){next}

    # else subtract overlapping land area
    overlap <- sf::st_intersection(sf::st_geometry(polygon), sf::st_geometry(land.shape))
    kud.polygons$area[i] <- kud.polygons$area[i] - sum(as.numeric(sf::st_area(overlap)))/1000000
    overlaps_number <- overlaps_number+1
  }

  if(overlaps_number==0){cat(paste0(out.prefix, " - No overlapping land areas were detected\n"))}
  if(overlaps_number>0){cat(paste0(out.prefix, " - Overlapping land areas subtracted from ", overlaps_number, " KUD(s)\n"))}
  return(kud.polygons)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
