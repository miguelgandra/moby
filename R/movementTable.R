#######################################################################################################
## Create movement table ###############################################################################
#######################################################################################################

#' Create movement stats table
#'
#' @description Creates a table containing  movement metrics for each animal
#' (total distances traveled, rate of movements, KUDs, etc.).
#'
#' @param data A data frame containing binned animal detections and distances traveled,
#' as returned by \code{\link{calculateDistances}}.
#' @param kud.results Output of \code{\link{calculateKUDs}}.
#' @param land.shape A projected shape file containing coastlines.
#' @param epsg.code Coordinate reference system used to project positions (class 'CRS').
#' If not supplied, CRS is assumed to be the same as in land.shape.
#' @param transition.layer Transition layer used to calculate shortest in-water paths, as
#' retrieved by the \code{\link{calculateDistances}} function.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param id.col Name of the column containing animal IDs Defaults to 'ID'.
#' @param lon.col Name of the column containing longitude values (unprojected). Defaults to 'lon'.
#' @param lat.col Name of the column containing latitude values (unprojected). Defaults to 'lat'.
#' @param dist.col Name of the column containing distance values (in meters). Defaults to 'dist_m'.
#' @param discard.missing If true, only individuals with detections are included.

#' @export


movementTable <- function(data, kud.results, land.shape, epsg.code=NULL, transition.layer, id.groups=NULL, id.col="ID",
                          lon.col="lon", lat.col="lat", dist.col="dist_m", discard.missing=T) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # check if data contains id.col
  if(!c("timebin") %in% colnames(data)){
    stop("'timebin' column not found in the supplied data")
  }

  # convert ids to factor
  if(class(data[,id.col])!="factor"){
    cat(paste("Warning: converting", id.col, "column to factor\n"))
    data[,id.col] <- as.factor(data[,id.col])
  }

  # check if data contains lon.col and lat.col
  if(!lon.col %in% colnames(data)){
    stop("Longitude column not found. Please assign the correct column name with 'lon.col'")
  }
  if(!lat.col %in% colnames(data)){
    stop("Latitude column not found. Please assign the correct column name with 'lat.col'")
  }

  # check if dataset coordinates are in geographic format (unprojected)
   geographic_coords <-  all(data[,lon.col]>=(-180) & data[,lon.col]<=180 & data[,lat.col]>=(-90) & data[,lat.col]<=90)
   if(is.na(geographic_coords)){
     stop("Longitudes/latitudes values are missing and/or in the wrong format\n")
   }

  # retrieve epsg.code if not provided
  if(is.null(epsg.code)){
    if(!grepl("+units=m", land.shape@proj4string, fixed=T)){
      stop("Please supply a projected land.shape (in metres)")
    }else{
      epsg.code <- land.shape@proj4string
      warning(paste0("Assuming CRS projection '", epsg.code, "'\n"))
    }
  } else {
    land.shape <- sp::spTransform(land.shape, epsg.code)
  }

  # check if data contains dist.col
  if(!dist.col %in% colnames(data)){
    stop("Distance column not found. Please assign the correct column name with 'dist.col'")
  }

  # set id groups
  if(is.null(id.groups)) {
    id.groups <- list(levels(data[,id.col]))
  } else {
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {"Some of the ID(s) in id.groups are not present in the data"}
  }

   # get time bins interval (in minutes) and interpolate distances if required
   individual_data <- split(data, f=data[,id.col])
   interval <- lapply(individual_data, function(x) difftime(x$timebin, data.table::shift(x$timebin), units="min"))
   interval <- unlist(lapply(interval, unique))
   interval <- unique(interval[!is.na(interval)])
   if(length(unique(interval))>2) {
     data <- interpolateDistances(data, keep.intermediate=T)
     individual_data <- split(data, f=data[,id.col])
     interval <- lapply(individual_data, function(x) difftime(x$timebin, data.table::shift(x$timebin), units="min"))
     interval <- unlist(lapply(interval, unique))
     interval <- unique(interval[!is.na(interval)])
   }

  # subset KUD results
  kud.results <- lapply(id.groups, function(x) kud.results$kud_table[kud.results$kud_table$ID %in% x,])
  data_list <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])



  #####################################################################
  ## Calculate stats ##################################################

  movement_table <- list()

  for(i in 1:length(data_list)){

    # grab data
    data <- data_list[[i]]
    data[,id.col] <- droplevels(data[,id.col])

    #total distance traveled (km)
    total_dist <- as.numeric(aggregate(data[,dist.col], by=list(data[,id.col]), sum, na.rm=T, drop=F)$x)
    total_distance <- sprintf("%.1f", total_dist/1000)

    # rate of movement (hourly distance - m/h)
    mean_rom <- as.numeric(aggregate(data[,dist.col], by=list(data[,id.col]), mean, na.rm=T, drop=F)$x)
    mean_rom <- mean_rom*60/interval
    max_rom <- as.numeric(aggregate(data[,dist.col], by=list(data[,id.col]), max, na.rm=T, drop=F)$x)
    max_rom <- max_rom*60/interval
    if(mean(mean_rom, na.rm=T)>1000 & mean(max_rom, na.rm=T)>1000){
      mean_rom <- mean_rom/1000
      max_rom <- max_rom/1000
      rom_units <- "km/h"
    }else{
      rom_units <- "m/h"
    }
    mean_rom <- sprintf("%.1f", mean_rom)
    max_rom <- sprintf("%.1f", max_rom)

    ## linearity index
    first_coas <- by(data, data[,id.col], function(x) x[which.min(x$timebin),c(lon.col, lat.col)], simplify=F)
    last_coas <- by(data, data[,id.col], function(x) x[which.max(x$timebin),c(lon.col, lat.col)], simplify=F)
    dist_diff <- mapply(calculatePairDistance, coord1=first_coas, coord2=last_coas,
                        MoreArgs=list(trCost=transition.layer, land.shape=land.shape, epsg.code=epsg.code, geographic_coords=geographic_coords))
    li_index <- sprintf("%.2f", dist_diff/total_dist)
    li_index[li_index=="NaN"] <- "NA"


    # overall movement stats
    movement_stats <- data.frame("ID"=levels(data[,id.col]), "Distance (km)"=total_distance,
                                 "ROM"=mean_rom, "Max ROM"=max_rom,
                                 "LI"=li_index, check.names=F, row.names=NULL)
    colnames(movement_stats)[3] <- paste0(colnames(movement_stats)[3], " (", rom_units,")")
    colnames(movement_stats)[4] <- paste0(colnames(movement_stats)[4], " (", rom_units,")")
    movement_stats <- plyr::join(movement_stats, kud.results[[i]], by="ID", type="left")


    # calculate means ± se and format missing values
    values <- suppressWarnings(as.matrix(sapply(movement_stats[,-1], as.numeric)))
    decimal_digits <- apply(values, 2, moby:::decimalPlaces)
    decimal_digits <- apply(decimal_digits, 2, max, na.rm=T)
    movement_stats$ID <- as.character(movement_stats$ID)
    movement_stats[nrow(movement_stats)+1,] <- NA
    movement_stats$ID[nrow(movement_stats)] <- "mean"
    movement_stats[nrow(movement_stats), -1] <- sprintf(paste0("%.", decimal_digits, "f"), colMeans(values, na.rm=T))
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(values, 2, plotrix::std.error)))
    movement_stats[nrow(movement_stats), -1] <- paste(movement_stats[nrow(movement_stats), -1], "±", errors)
    movement_stats[is.na(movement_stats)] <- "-"
    movement_stats[movement_stats=="NA"] <- "-"

    # discard missing IDs if required
    if(discard.missing==T){
      missing_ids <- names(table(data[,id.col])[table(data[,id.col])==0])
      movement_stats <- movement_stats[!movement_stats$ID %in% missing_ids,]
    }

    # add title
    if(length(id.groups)>1){
      table_title <- movement_stats[0,]
      table_title[1,] <- ""
      table_title$ID <- names(id.groups)[i]
      movement_stats <- rbind(table_title, movement_stats)
    }

    # save table
    movement_table[[i]] <- movement_stats
  }


  # aggregate group tables
  movement_table <- do.call("rbind", movement_table)

  # return table
  return(movement_table)
}



################################################################################
# Helper function I - calculate shortest path between two points ###############

calculatePairDistance <- function(coord1, coord2, trCost, land.shape, epsg.code, geographic_coords){
  if(is.null(coord1) | is.null(coord2)) {return(NA)}
  if(class(coord1)[1]!="matrix" | class(coord2)[1]!="matrix"){
    coord1 <- as.matrix(coord1)
    coord2 <- as.matrix(coord2)
  }

  # project coords if required
  coords <- data.frame(rbind(coord1, coord2))
  colnames(coords) <- c("lon", "lat")
  if(geographic_coords==T){
    coords <- st_as_sf(coords, coords=c("lon", "lat"), crs=sf::st_crs("+proj=longlat +datum=WGS84"), agr="constant")
    coords <- sf::st_transform(coords, crs=sf::st_crs(epsg.code))
  }else{
    coords <- st_as_sf(coords, coords=c("lon", "lat"), crs=sf::st_crs(epsg.code), agr="constant")
  }

  # create segment
  coords <- sf::st_coordinates(coords)
  segment <- sf::st_linestring(coords)
  segment <- sf::st_sfc(segment, crs=sf::st_crs(epsg.code))
  is_pt <- all(coord1 == coord2)
  if(is_pt==T){return(0)}
  in_land <- lengths(sf::st_intersects(segment, sf::st_as_sf(land.shape), sparse=T))>0
  if(in_land==T){
    pts <- sf::st_coordinates(segment)[,1:2]
    segment <- gdistance::shortestPath(trCost, pts[1,], pts[2,], output="SpatialLines")
    segment <- sf::st_as_sf(segment)
    sf::st_crs(segment) <- sf::st_crs(epsg.code)
  }

  # convert back to WGS84 and estimate distance
  segment <- sf::st_transform(segment, crs=sf::st_crs("+proj=longlat +datum=WGS84"))
  dist <- sum(geosphere::distVincentyEllipsoid(sf::st_coordinates(segment)[,1:2]))
  return(dist)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
