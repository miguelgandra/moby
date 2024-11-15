#######################################################################################################
## Create movement table ###############################################################################
#######################################################################################################

#' Create movement stats table
#'
#' @description Creates a table containing  movement metrics for each animal
#' (total distances traveled, rate of movements, KUDs, etc.).
#'
#' @inheritParams setDefaults
#' @param data A data frame containing binned animal detections and distances traveled,
#' as returned by \code{\link{calculateTracks}}.
#' @param kud.results Output of \code{\link{calculateKUDs}}.
#' @param land.shape A projected shape file containing coastlines.
#' @param epsg.code Coordinate reference system used to project positions (class 'CRS').
#' If not supplied, CRS is assumed to be the same as in land.shape.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param dist.col Name of the column containing distance values (in meters). Defaults to 'dist_m'.
#' @param discard.missing If true, only individuals with detections are included.
#' @param ... Further arguments passed to \code{\link{calculateTracks}} function,
#' in order to calculate distances between the first and last recorded detection, for each individual.
#' @export


movementTable <- function(data,
                          kud.results,
                          land.shape,
                          epsg.code = getDefaults("epsg"),
                          id.groups = NULL,
                          id.col = getDefaults("ID"),
                          timebin.col = getDefaults("timebin"),
                          lon.col = getDefaults("lon"),
                          lat.col = getDefaults("lat"),
                          dist.col = "dist_m",
                          discard.missing = TRUE,
                          ...) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # check kud.results
  if(!c("bandwidth") %in% names(attributes(kud.results))) stop("The supplied kud.results do not seem to be in the right format. Please use the output of the 'calculateKUDs' function.")

  # check if dataset coordinates are in geographic format (unprojected)
   geographic_coords <-  all(data[,lon.col]>=(-180) & data[,lon.col]<=180 & data[,lat.col]>=(-90) & data[,lat.col]<=90)
   if(is.na(geographic_coords)){
     stop("Longitudes/latitudes values are missing and/or in the wrong format\n")
   }

  # retrieve epsg.code if not provided
  if(is.null(epsg.code)){
    if(!grepl("+units=m", land.shape@proj4string, fixed=TRUE)){
      stop("Please supply a projected land.shape (in metres)")
    }else{
      epsg.code <- land.shape@proj4string
      warning(paste0("Assuming CRS projection '", epsg.code, "'\n"))
    }
  } else {
    land.shape <- sp::spTransform(land.shape, epsg.code)
  }

   # get time bins interval (in minutes) and interpolate distances if required
   individual_data <- split(data, f=data[,id.col])
   interval <- lapply(individual_data, function(x) difftime(x[,timebin.col], dplyr::lag(x[,timebin.col]), units="min"))
   interval <- unlist(lapply(interval, unique))
   interval <- unique(interval[!is.na(interval)])
   if(length(unique(interval))>2) {
     data <- interpolateDistances(data, keep.intermediate=TRUE)
     individual_data <- split(data, f=data[,id.col])
     interval <- lapply(individual_data, function(x) difftime(x[,timebin.col], dplyr::lag(x[,timebin.col]), units="min"))
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
    total_dist <- as.numeric(stats::aggregate(data[,dist.col], by=list(data[,id.col]), sum, na.rm=TRUE, drop=FALSE)$x)
    total_distance <- sprintf("%.1f", total_dist/1000)

    # rate of movement (hourly distance - m/h)
    mean_rom <- as.numeric(stats::aggregate(data[,dist.col], by=list(data[,id.col]), mean, na.rm=TRUE, drop=FALSE)$x)
    mean_rom <- mean_rom*60/interval
    max_rom <- as.numeric(stats::aggregate(data[,dist.col], by=list(data[,id.col]), max, na.rm=TRUE, drop=FALSE)$x)
    max_rom <- max_rom*60/interval
    if(mean(mean_rom, na.rm=TRUE)>1000 & mean(max_rom, na.rm=TRUE)>1000){
      mean_rom <- mean_rom/1000
      max_rom <- max_rom/1000
      rom_units <- "km/h"
    }else{
      rom_units <- "m/h"
    }
    mean_rom <- sprintf("%.1f", mean_rom)
    max_rom <- sprintf("%.1f", max_rom)

    ## linearity index
    first_coas <- by(data, data[,id.col], function(x) x[which.min(x[,timebin.col]),c(lon.col, lat.col)], simplify=FALSE)
    first_coas <- do.call("rbind", first_coas)
    first_coas$id <- rownames(first_coas)
    first_coas$type <- "start"
    last_coas <- by(data, data[,id.col], function(x) x[which.max(x[,timebin.col]),c(lon.col, lat.col)], simplify=FALSE)
    last_coas <- do.call("rbind", last_coas)
    last_coas$id <- rownames(last_coas)
    last_coas$type <- "end"
    start_end_coas <- rbind(first_coas, last_coas)
    rownames(start_end_coas) <- NULL
    start_end_coas$type <- factor(start_end_coas$type, levels=c("start","end"))
    start_end_coas <- start_end_coas[order(start_end_coas$id, start_end_coas$type),]
    start_end_dists <- calculateTracks(start_end_coas, land.shape=land.shape,
                                       epsg.code=epsg.code, id.col="id", verbose=FALSE,
                                       lon.col=lon.col, lat.col=lat.col, ...)$data
    start_end_dists <- start_end_dists$dist_m[!is.na(start_end_dists$dist_m)]
    li_index <- sprintf("%.2f", start_end_dists/total_dist)
    li_index[li_index=="NaN"] <- "NA"


    # overall movement stats
    movement_stats <- data.frame("ID"=levels(data[,id.col]), "Distance (km)"=total_distance,
                                 "ROM"=mean_rom, "Max ROM"=max_rom,
                                 "LI"=li_index, check.names=FALSE, row.names=NULL)
    colnames(movement_stats)[3] <- paste0(colnames(movement_stats)[3], " (", rom_units,")")
    colnames(movement_stats)[4] <- paste0(colnames(movement_stats)[4], " (", rom_units,")")
    movement_stats <- plyr::join(movement_stats, kud.results[[i]], by="ID", type="left")


    # calculate means ± se and format missing values
    values <- suppressWarnings(as.matrix(sapply(movement_stats[,-1], as.numeric)))
    decimal_digits <- apply(values, 2, .decimalPlaces)
    decimal_digits <- apply(decimal_digits, 2, max, na.rm=TRUE)
    movement_stats$ID <- as.character(movement_stats$ID)
    movement_stats[nrow(movement_stats)+1,] <- NA
    movement_stats$ID[nrow(movement_stats)] <- "mean"
    movement_stats[nrow(movement_stats), -1] <- sprintf(paste0("%.", decimal_digits, "f"), colMeans(values, na.rm=TRUE))
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(values, 2, plotrix::std.error)))
    movement_stats[nrow(movement_stats), -1] <- paste(movement_stats[nrow(movement_stats), -1], "\u00b1", errors)
    movement_stats[is.na(movement_stats)] <- "-"
    movement_stats[movement_stats=="NA"] <- "-"

    # discard missing IDs if required
    if(discard.missing){
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


#######################################################################################################
#######################################################################################################
#######################################################################################################
