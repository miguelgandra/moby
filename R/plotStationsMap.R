#######################################################################################################
# Plot stations map  ##################################################################################
#######################################################################################################

#' Plot map with receiver locations

#' @description Function to plot all receivers locations on a map. Pie charts representing the
#' average of detection frequencies per receiver are displayed above. Additionally,
#' detection frequencies can be calculated separately by a grouping variable (e.g. diel phase).
#'
#' @param data A data frame containing animal detections with corresponding time-bins.
#' @param split.by If defined, detection frequencies on each receiver are calculated
#' separately for each level of this variable. If NULL, frequencies are calculated
#' for the entire study duration. Defaults to "timeofday".
#' @param lon.name Name of the column containing longitude values.
#' @param lat.name Name of the column containing latitude values.
#' @param pie.color Color of the pie charts. Single value if 'split.by' is not defined, otherwise
#' a vector containing the same number of colors as the number of levels in the grouping variable.
#' @param land.shape A shape file containing coastlines.
#' @param land.color Color of land areas.
#' @param background.col Background color if a background layer is not supplied.
#' @param background.layer A projected raster containing a variable to be displayed in the background
#' (e.g. bathymetry, temperature, etc.).
#' @param background.pal Color palette for the background layer.
#' @param scale.meters Distance covered by the scale bar, in meters. If null, it is
#' automatically defined as 20% of the plot region.
#' @param scale.pos Position of the scale bar, specified by keyword.
#' @param scale.inset Inset distance(s) from the margins as a fraction of the
#' plot region, relatively to the scale bar position.
#' @export


#######################################################################################
## Main function - plots map with the receivers #######################################

plotStationsMap <- function(data, split.by="timeofday", lon.name="longitude", lat.name="latitude",
                            pie.color=NULL, land.shape, land.color="gray50", background.col="#F3F7F7",
                            background.layer=NULL,  background.pal=NULL, scale.meters=NULL,
                            scale.pos="bottomright", scale.inset=0.1) {

  ############################################################################
  ## Initial checks ##########################################################

  if(!is.null(split.by) & !split.by %in% colnames(data)){
    stop("split.by variable not found in the supplied data")
  }

  if(!lon.name %in% colnames(data) | !lat.name %in% colnames(data)){
    stop("longitude/latitude columns not found in the supplied data")
  }

  if(missing(land.shape)){stop("land.shape argument missing")}

  cat("Plotting stations map\n")

  ############################################################################
  ## Calculate average proportion of detections per receiver #################

  # split COAs from different individuals
  data_individual <- split(data, f=data$ID)
  data_plot <- list()

  # if split.by is not defined, calculate average of overall detection frequencies
  if(is.null(split.by)){
    for(f in 1:length(data_individual)) {
      detections <- data_individual[[f]]
      data_stations <- aggregate(detections$timebin, by=list(detections$ID, detections$station), length)
      colnames(data_stations) <- c("ID", "station", "detections")
      data_stations$detections[is.na(data_stations$detections)] <- 0
      data_stations$freq <- data_stations$detections / sum(data_stations$detections)
      data_stations <- data_stations[order(data_stations$station),]
      data_plot[[f]] <- data_stations
    }
  }

  # else calculate average of detection frequencies by each grouping level
  if(!is.null(split.by)){
    for(f in 1:length(data_individual)) {
      data[,split.by] <- as.factor(data[,split.by])
      groups <- levels(data[,split.by])
      detections <- data_individual[[f]]
      if(nrow(detections)==0 | length(unique(detections[,split.by]))==1){data_plot[[f]]<-NULL; next}
      data_stations <- aggregate(detections$timebin, by=list(detections$ID, detections$station, detections[,split.by]), length)
      colnames(data_stations) <- c("ID", "station", split.by, "detections")
      data_stations <- reshape2::dcast(data_stations, as.formula(paste0("ID+station~", split.by)), value.var="detections", fill=0)
      data_stations$total <- rowSums(data_stations[,-c(1:2)])
      total_detections <- sum(data_stations$total)
      data_freqs <- data_stations[,-c(1:2)] / total_detections
      colnames(data_freqs) <- paste0(colnames(data_stations[,-c(1:2)]), "_freq")
      data_stations <- cbind(data_stations, data_freqs)
      data_stations <- data_stations[order(data_stations$station),]
      data_plot[[f]] <- data_stations
    }
  }

  data_plot <- data_plot[!sapply(data_plot,is.null)]
  data_plot <- do.call("rbind", data_plot)
  freq_cols <- which(grepl("freq", colnames(data_plot), fixed=T))
  data_plot <- aggregate(data_plot[,freq_cols], by=list(data_plot$station), mean)
  colnames(data_plot)[1] <- "station"
  sector_cols <- which(grepl("freq", colnames(data_plot), fixed=T) & colnames(data_plot)!=c("total_freq"))
  data_sector <- data_plot[,sector_cols]
  data_sector[] <- lapply(data_sector, function(x) as.numeric(x))
  for(i in 1:nrow(data_sector)){data_sector[i,] <- data_sector[i,]/sum(data_sector[i,])}
  full_pies <- which(apply(data_sector, 1, function(r) any(r==1)))
  sector_pies <- which(apply(data_sector, 1, function(r) all(r!=1)))
  full_col <- apply(data_sector, 1, function(r) which(r==1))
  full_col <- unlist(lapply(full_col, function(x) ifelse(length(x)==0, NA, x)))


  ############################################################################
  ## Grab station coordinates ################################################
  stations_list <- aggregate(cbind(data[,lon.name], data[,lat.name]), by=list(data$station), mean)
  colnames(stations_list) <- c("station", "longitude", "latitude")
  coordinates <- sp::SpatialPoints(cbind(stations_list$longitude, stations_list$latitude))
  raster::projection(coordinates) <- CRS("+proj=longlat +datum=WGS84")
  coordinates <- sp::spTransform(coordinates, land.shape@proj4string)
  stations_list$longitude <- coordinates@coords[,1]
  stations_list$latitude <- coordinates@coords[,2]
  data_plot <- plyr::join(data_plot, stations_list, by="station", type="left")
  template_raster <- raster::raster(raster::extent(coordinates)*1.2, res=100)
  template_raster[] <- 0

  # set pie chart radius based on calculated detection freqs
  # (between 1% and 5% of the longitudinal axis)
  pie_min <- (extent(coordinates)[2] - extent(coordinates)[1]) * 0.01
  pie_max <- (extent(coordinates)[2] - extent(coordinates)[1])* 0.05
  data_plot$total_freq <- moby:::rescale(data_plot$total_freq, c(pie_min,pie_max))

  # set color palette
  if(is.null(pie.color)){
    if(!is.null(split.by)){pie.color <- colorRampPalette(colors=c("white", "gray40"))(length(groups))}
    if(is.null(split.by)){pie.color <- adjustcolor("orange", alpha.f=0.6)}
  }

  # repel pie charts
  circle_coords <- data_plot[,c("longitude", "latitude","total_freq")]
  circle_coords$total_freq <- circle_coords$total_freq+pie_min
  long_lims <- extent(template_raster)[1:2]
  lat_lims <- extent(template_raster)[3:4]
  pie_coords <- packcircles::circleRepelLayout(circle_coords, xlim=long_lims, ylim=lat_lims,
                                               xysizecols=c(1,2,3), sizetype = c("radius"),
                                               maxiter=1000, wrap=F, weights=1)$layout

  ############################################################################
  ## Plot map ################################################################

  par(mar=c(.5,.5,.5,1), oma=c(0,0,0,0))
  plot(template_raster, col=adjustcolor("white", alpha.f=0), axes=F, legend=F)
  if(!is.null(background.layer)){plot(background.layer, col=background.pal, legend=F, axes=F, add=T)}
  if(is.null(background.layer)){rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.col, border=NA)}
  plot(land.shape, col=land.color, border=NA, add=T)
  if(is.null(scale.meters)){scale.meters <- min(pretty((par("usr")[2]-par("usr")[1])*0.15))}
  scale_xy <- moby:::getPosition(scale.pos, inset=scale.inset)
  scale_km <- scale.meters/1000
  moby:::scalebar(d=scale.meters, xy=scale_xy, type="bar", divs=2, below="km", label=c(0, scale_km/2, scale_km), lwd=0.2, cex=0.5, bar.lwd=0.2)
  points(coordinates, pch=16, bg="black", cex=0.2)
  for(i in 1:nrow(data_plot)){segments(x0=coordinates@coords[i,1], y0=coordinates@coords[i,2], x1=pie_coords$x[i], y1=pie_coords$y[i], lwd=0.2, lty=2)}
  for(i in sector_pies){
    plotrix::floating.pie(xpos=pie_coords$x[i], ypos=pie_coords$y[i], x=as.numeric(data_sector[i,]),
                          lwd=0.2, radius=data_plot$total_freq[i], col=pie.color)
  }
  for(i in full_pies){
    plotrix::draw.circle(x=pie_coords$x[i], y=pie_coords$y[i], radius=data_plot$total_freq[i],
                         border="black", lwd=0.2, col=pie.color[full_col[i]])
  }

  par(lwd=0.4)
  legend("topright", legend=groups, fill=pie.color, bty="n", cex=0.8, inset=c(-0.14, 0.1), xpd=T)
}


##################################################################################################
##################################################################################################
##################################################################################################
