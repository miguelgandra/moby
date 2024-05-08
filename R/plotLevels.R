#######################################################################################################
# Function to generate level plots  #################################################################
#######################################################################################################

#' Generate level plots
#'
#' @description Function to represent a given metric (nº detections, presences, individuals or co-occurrences)
#' per time bin. It generates 2-dimensional plots (hour x date) highlighting the chosen variable
#' either through color-coded cells or size-variable points, and illustrating the
#' variation in diel phases' hours across the study duration (annual variation of daylight time).
#'
#' @param data A data frame containing binned animal detections. Needs to include a 'timebin' column.
#' @param tagging.dates Optional. A POSIXct vector containing the tag/release date of each animal.
#' If supplied, plots start is set the earliest release date, instead of the first detection.
#' @param split.by Optional. If defined, plots are generated individually for each level
#' of this variable (e.g. species, ontogeny or habitat).
#' @param variables The type(s) of metric to plot. Accepted types: "detections", "individuals" and "co-occurrences".
#' @param style Style of the plot. Either "raster" (nº represented by color) or "points" (nº represented by size).
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param station.col Name of the column containing station/receiver IDs Defaults to 'station'.
#' @param color.by Variable defining the color group of the plotted metric, when style = "points".
#' Can be used for example to display detections by receiver, animal trait or temporal category.
#' @param color.pal Color palette for the level plot.
#' @param date.format Date-time format (as used in \code{\link[base]{strptime}}),
#' defining the x-axis labels. Defaults to month ("%b").
#' @param date.interval Number defining the interval between each
#' displayed date (x-axis label). Defaults to 4.
#' @param date.start Integer defining the first displayed date (can be used in combination
#'  with 'date.interval" to better control the x-axis labels). Defaults to 1.
#' @param diel.lines Number indicating the number of diel phase lines (boundaries)
#' to display. Either 2 (corresponding to sunrise and sunset) or 4 (depicting
#' dawn, sunrise, sunset, dusk).
#' @param polygons Either "diel", "season". If style is set to points, can be used to highlight
#' diel or seasonal phases using color-coded polygons. Defaults to F (no polygons).
#' @param sunriset.coords A SpatialPoints or matrix object containing longitude and
#' latitude coordinates (in that order) at which to estimate sunrise and sunset times.
#' @param solar.depth Angle of the sun below the horizon (in degrees). Passed the
#' solarDep argument in \code{\link[maptools]{crepuscule}} function.
#' @param lunar.info Include top mural containing lunar illumination?
#' @param background.col If style is set to points, defines the background color of the plot.
#' @param grid If true, a grid is plotted (horizontal and vertical guides across the entire plot).
#' @param grid.color Color of the grid lines, if 'grid' is set to true.
#' @param highlight.isolated If a 'color.by' variable is defined and style="points, points are ordered
#' and plotted according to their "density", with isolated levels being brought forward
#' to prevent them from being hidden behind denser point clouds.
#' @param pt.cex Defines the min and max point radius when style="point. Defaults to c(0,3).
#' @param uniformize.scale Use the same color scale for all plots when 'split.by' is defined.
#' @param uniformize.dates Use the same date range for all plots when 'split.by' is defined.
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @param ... Further arguments passed to the \code{\link[graphics]{image}} function (if style="raster"),
#' or to the \code{\link[graphics]{plot}} function (if style="points").
#' @export


plotLevels <- function(data, tagging.dates=NULL, variables="detections", split.by=NULL, style="raster",  id.col="ID",
                       station.col="station", color.by=NULL, color.pal=NULL, date.format="%d/%b",
                       date.interval=4, date.start=1, diel.lines=2, sunriset.coords, solar.depth=18, polygons=F, lunar.info=F,
                       background.col="white", grid=T, grid.color="white", highlight.isolated=F,
                       pt.cex=c(0.5, 3), uniformize.scale=F, uniformize.dates=T, cols=NULL, ...) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  # print to console
  cat("Generating level plot(s)\n")

  if(any(!variables %in% c("detections", "individuals", "co-occurrences"))){
    stop("Wrong variable(s) specified. Accepted values: 'detections','individuals' or 'co-occurrences'")
  }

  # check if style contains one of the accepted types
  if(!style %in% c("raster", "points")){
    stop("Style should be either 'raster' or 'points'")
  }

  # check if data contains id.col
  if(!id.col %in% colnames(data)) {
    stop("'id.col' variable not found in the supplied data")
  }

  # check if data contains timebin column
  if(!c("timebin") %in% colnames(data)) {
    stop("'timebin' column not found in the supplied data")
  }

  # check if data contains the station variable
  if(!station.col %in% colnames(data)){
    stop("'station.col' variable not found in the supplied data")
  }

  # check if data contains the split.by variable
  if(!is.null(split.by) && !split.by %in% colnames(data)){
    stop("'split.by' variable not found in the supplied data")
  }

  # check if data contains the color.by variable
  if(!is.null(color.by)){

    if(!color.by %in% colnames(data)) {
      stop("'color.by' variable not found in the supplied data")}

    if(any(is.na(data[,color.by]))) {
      stop("Missing values in 'color.by' variable")}

    if(class(data[,color.by])!="factor"){
      data[,color.by] <- as.factor(data[,color.by])
      cat("Warning: 'color.by' variable converted to factor\n")}

    if(!is.null(color.pal)){
      n_groups <- length(unique(data[,color.by])[!is.na(unique(data[,color.by]))])
      if(!class(data[,color.by])[1] %in% c("integer", "numeric") & length(color.pal)!=n_groups){
        cat("Warning: The nº of colors doesn't match the number of levels in 'color.by' variable")
      }}

  }

  if(!polygons %in% c(F, "diel","season")) {
    stop("polygons variable can only be set to 'diel', 'season' or 'F' (disabled)")
  }

  if(length(pt.cex)!=2) {
    stop("pt.cex should be a vector of length 2 (min and max)")
  }

  if(style=="raster"){
    if(!is.null(color.by)){
      cat("Warning: 'color.by' argument is discarded when plot if of type raster")
    }

    if(!is.null(color.pal) & class(color.pal)!="function"){
      stop("Please supply a color palette *function* when plot if of type raster")
    }

  }


  #####################################################################################
  # Aggregate data and calculate stats ################################################
  #####################################################################################

  # aggregate data in a list
  if(!is.null(split.by)){
    if(class(data[,split.by])!="factor"){
      data[,split.by] <- as.factor(data[,split.by] )
    }
    data_list <- split(data, f=data[,split.by])
    data_list <- lapply(data_list, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
    groups <- levels(data[,split.by])
  }else{
    data_list <- list(data)
    groups <- "complete"
  }

  # calculate all variable/group combinations
  plot_combs <- expand.grid(variables, groups)
  colnames(plot_combs) <- c("variable", "group")
  plot_combs$variable <- as.character(plot_combs$variable)

  # initialize list variable
  plot_template <- list()
  plot_data_list <- list()

  # iterate through each data group/variable
  for(i in 1:nrow(plot_combs)){

    group <- plot_combs$group[i]
    var <-  plot_combs$variable[i]
    data_group <- data_list[[group]]
    group_ids <- which(levels(data[,id.col]) %in% levels(data_group[,id.col]))

    # set dates scale
    if(uniformize.dates==T){
      if(!is.null(tagging.dates)){first_dates <- min(tagging.dates)}
      if(is.null(tagging.dates)){first_dates <- min(data$timebin, na.rm=T)}
      last_dates <- max(data$timebin, na.rm=T)
    }else{
      if(!is.null(tagging.dates)){first_dates <- tagging.dates[group_ids]}
      if(is.null(tagging.dates)){first_dates <- min(data_group$timebin, na.rm=T)}
      last_dates <- max(data_group$timebin, na.rm=T)
    }

    # create table with complete seq of timebins
    all_timebins <- createWideTable(data_group, start.dates=first_dates, end.dates=last_dates, value.col="detections", round.dates=T)

    # create day x time matrix
    all_timebins$day <- strftime(all_timebins$timebin,"%Y-%m-%d", tz="UTC")
    all_timebins$time <- strftime(all_timebins$timebin,"%H:%M:%S", tz="UTC")
    days <- unique(all_timebins$day)
    bins <- unique(all_timebins$time)
    plot_matrix <- t(matrix(NA, nrow=length(bins), ncol=length(days)))
    rownames(plot_matrix) <- days
    colnames(plot_matrix) <- bins

    # aggregate data by timebin, station and ID
    data_group$row <- 1:nrow(data_group)
    detections_table <- aggregate(data_group$row, by=list(data_group$timebin, data_group[,id.col], data_group[,station.col]), length)
    colnames(detections_table) <- c("timebin", "id", "station", "detections")

    # calculate metrics
    if(var=="detections"){
      plot_data <- aggregate(detections_table$detections, by=list(detections_table$timebin, detections_table$station), sum)
      colnames(plot_data) <- c("timebin", "station", "var")
    }else{
      plot_data <- aggregate(detections_table$id, by=list(detections_table$timebin, detections_table$station), function(x) length(unique(x)))
      colnames(plot_data) <- c("timebin", "station", "var")
      if(var=="co-occurrences"){
        plot_data <- plot_data[plot_data$var>1,]
      }
    }
    plot_data <- plot_data[order(plot_data$timebin),]
    plot_data$day <- strftime(plot_data$timebin, "%Y-%m-%d", tz="UTC")
    plot_data$hour <- strftime(plot_data$timebin, "%H:%M:%S", tz="UTC")

    # calculate stats for color-by variable if required
    if(!is.null(color.by)){
      color_data <- aggregate(data_group[,color.by], by=list(data_group$timebin, data_group[,id.col], data_group[,station.col]), aggFun)
      colnames(color_data) <- c("timebin", "id", "station", "level")
      color_data <-  aggregate(color_data$level, by=list(color_data$timebin, color_data$station), aggFun)
      colnames(color_data) <- c("timebin", "station", "level")
      if(class(color_data$level)[1] %in% c("numeric","integer")){
        color_data$color <- round(scales::rescale(color_data$level, to=c(1,256)))
      }else{
        color_data$level <- factor(color_data$level, levels=levels(data_group[,color.by]))
        color_data$color <- color_data$level
        levels(color_data$color) <- 1:nlevels(color_data$color)
        color_data$color <- as.integer(color_data$color)
      }
      plot_data <- plyr::join(plot_data, color_data, by=c("timebin", "station"), type="left")
    }else{
      plot_data$color <- 1
    }

    # save formatted data to list
    plot_template[[i]] <- plot_matrix
    plot_data_list[[i]] <- plot_data
    names(plot_data_list)[i] <- paste0(group, "<->", var)
  }



  #####################################################################################
  # Get variables range  ##############################################################
  #####################################################################################

  var_range_list <- list()
  for(i in 1:length(variables)){
    var <- variables[i]
    var_indexes <- which(grepl(var, names(plot_data_list), fixed=T))
    var_range <- plot_data_list[var_indexes]
    var_range <- unlist(lapply(var_range, function(x) range(x$var, na.rm=T)))
    var_range_list[[i]] <- range(var_range, na.rm=T)
    names(var_range_list)[i] <- var
  }

  dates_range <- list()
  for(i in 1:length(plot_data_list)){
    dates_range[[i]] <- lapply(plot_data_list[i], function(x) range(as.character(x$timebin)))
    dates_range[[i]] <- paste(unlist(dates_range[[i]]), collapse=" / ")
  }
  same_dates <- length(unique(unlist(dates_range)))==1



  #####################################################################################
  # Set plot variables ################################################################
  #####################################################################################

  # set nº cols based on nº variables and groups
  if(is.null(cols)){
    if(length(groups)>1 & length(variables)>1){
      cols <- length(variables)
    }else{
      cols <- 1
    }
  }

  # calculate plot grid and set layout
  if(lunar.info==F){
    nplots <- nrow(plot_combs)
    rows <- ceiling(nplots/cols)
    plot_grid <- matrix(1:nplots, ncol=cols, byrow=T)
    layout(plot_grid)
  }else{
    nplots <- nrow(plot_combs)+cols
    rows <- ceiling(nplots/cols)
    fig_height <- dev.size("cm")[2]
    lunar_height <- fig_height*0.2
    plots_heights <- (fig_height-lunar_height)/(rows-1)
    layout_heights <- c(lunar_height, rep(plots_heights, rows-1))
    plot_grid <- matrix(1:nplots, ncol=cols, byrow=T)
    layout(plot_grid, heights=lcm(layout_heights))
  }

  # get indexes of bottom and right plots (to draw axes)
  #bottom_plots <- apply(plot_grid, 2, max, na.rm=T)
  #if(lunar.info==T){bottom_plots <- bottom_plots-cols}
  bottom_plots <- as.numeric(plot_grid)
  if(lunar.info==T){bottom_plots <- bottom_plots-cols}
  right_plots <- apply(plot_grid, 1, min, na.rm=T)

  # set color palette
  if(is.null(color.pal)){
    if(style=="points"){
      if(is.null(color.by)){
        color.pal <- "grey25"
      } else{
        if(class(data[,color.by])[1] %in% c("numeric","integer")){
          color.pal <- viridis::viridis(256)
        }else{
          ngroups <- length(unique(data[,color.by]))
          if(ngroups==3) {color.pal <- c("#326FA5","#D73134","#1C8E43")
          }else if (ngroups>3 & ngroups<10){color.pal <- ggthemes::economist_pal()(ngroups)
          }else {color.pal <- rainbow(ngroups)}
        }
      }
    }
  }


  #####################################################################################
  # Plot lunar graphics if lunar.illumination is set true #############################
  #####################################################################################

  if(lunar.info==T){
    timebins <- as.POSIXct(days, "%Y-%m-%d", tz="UTC")
    lunar_illumination <- data.frame("timebin"=timebins, "fraction"=lunar::lunar.illumination(timebins))
    lunar_illumination$phase <- as.character(lunar::lunar.phase(timebins, name=T))

    for(l in 1:cols){
      par(mar=c(2,5,7,6), xpd=T)
      xcoords <- c(1, seq(1, nrow(lunar_illumination)), nrow(lunar_illumination))
      ycoords <- c(0, lunar_illumination$fraction, 0)
      plot(lunar_illumination$fraction, type="l", lwd=0.5, axes=F, ylim=c(0,1), xlim=c(1, nrow(lunar_illumination)),
           xlab="", ylab="Lunar\nIllumination", xaxs="i", yaxs="i")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="#3C3C3C", border=NULL)
      polygon(x=xcoords, y=ycoords, col=adjustcolor("#FFFBEB", alpha.f=0.8), lwd=0.6)
      axis(2, at=seq(0,1,by=0.5), labels=sprintf("%.1f", seq(0,1,by=0.5)), cex.axis=0.8, las=1)
      box()
      moon_cycles <- rle(lunar_illumination$phase)
      lunar_illumination$cycle <- rep(1:length(moon_cycles$lengths), moon_cycles$lengths)
      lunar_illumination$cycle <- paste0(lunar_illumination$phase, "_", lunar_illumination$cycle)
      lunar_illumination$index <- 1:nrow(lunar_illumination)
      moon_indexes <- split(lunar_illumination, f=lunar_illumination$cycle)
      moon_lookup <- data.frame(phase=c("New", "Waxing", "Full", "Waning"), fraction1=c(0, 0.5, 1, 0.5))
      moon_indexes <- unlist(lapply(moon_indexes, function(x) x$index[which.min(abs(x$fraction - moon_lookup$fraction1[moon_lookup$phase==unique(x$phase)]))]))
      moon_indexes <- moon_indexes[order(as.numeric(sub('.*_', '', names(moon_indexes))))]
      moon_indexes <- data.frame("phase"=sub("\\_.*", "", names(moon_indexes)), "index"=moon_indexes, row.names=NULL)
      moon_indexes <- join(moon_indexes, moon_lookup, by="phase", type="left")
      moon_indexes$fraction2 <- abs(1- moon_indexes$fraction1)
      for(i in 1:nrow(moon_indexes)){
        if(moon_indexes$phase[i]!="Waning"){moon_color <- c("#FFFBEB", "#3C3C3C")}
        if(moon_indexes$phase[i]=="Waxing"){moon_color <- c("#3C3C3C", "#FFFBEB")}
        plotrix::floating.pie(xpos=moon_indexes$index[i], ypos=1.6, x=as.numeric(moon_indexes[i,c("fraction1", "fraction2")]),
                              startpos=90*pi/180, col=moon_color, lwd=0.5, border=NULL, xpd=T)
      }
    }
  }

  #####################################################################################
  # Generate level plots ##############################################################
  #####################################################################################

  ############################################################################
  # iterate through each data group/variable #################################

  for(i in 1:nrow(plot_combs)){

    # set variables
    group <- plot_combs$group[i]
    var <-  plot_combs$variable[i]
    plot_data <- plot_data_list[[i]]
    plot_title <- switch(var, detections="Nº of detections", individuals="Nº of individuals",
                         "co-occurrences"="Nº of co-occurring animals")
    if(!is.null(split.by)){plot_title <- paste(group, "-", plot_title)}

    par(mar=c(5,5,3,6), mgp=c(3,0.7,0))

    # get diel phase timelines
    days <- as.POSIXct(rownames(plot_template[[i]]), "%Y-%m-%d", tz="UTC")
    bins <- colnames(plot_template[[i]])
    interval <- as.numeric(difftime(data$timebin, data.table::shift(data$timebin), units="mins"))
    interval <- min(interval[interval>0], na.rm=T)
    daytimes_table <- getSunTimes(sunriset.coords, min(days), max(days), by="%Y-%m-%d", solar.depth)
    daytimes_table[,-1] <- daytimes_table[,-1] * (60/interval)
    daytimes_table[,-1] <- daytimes_table[,-1] + 1

    # prepare date variables
    all_dates <- strftime(days, date.format)
    consec_dates <- rle(all_dates)
    consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
    unique_dates <- unique(consec_dates)
    disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
    indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
    disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
    disp_dates <- sub("\\_.*", "", disp_dates)

    # prepare hour variables
    hours <- lubridate::hms(bins)
    hours <- lubridate::hour(hours)
    hour_indexes <- sapply(0:23, function(x) min(which(hours==x)))
    disp_hours <- paste0(unique(hours), "h")


    # set scale
    if(uniformize.scale==T){
      var_range <- var_range_list[[var]]
    }else{
      var_range <- range(plot_data$var, na.rm=T)
    }


    ##########################################################
    # if style is set to raster use "image"   ################
    if(style=="raster"){
      if(is.null(color.pal)){raster_pal <- grey.colors(max(var_range)+1, start=1, end=0)}
      if(!is.null(color.pal)){raster_pal <- color.pal(max(var_range)+1)}
      image(y=1:length(bins), x=1:length(days), zlim=var_range, z=plot_template[[i]],  xlab="", ylab="",
            main="", col=raster_pal, axes=F, cex.lab=0.9, ...)
    }

    ##########################################################
    # else, if style is set to points use "plot"
    if(style=="points"){

      # create empty plot and fill background
      plot(0,0, type="n", xlim=c(1, nrow(plot_template[[i]])), ylim=c(1, ncol(plot_template[[i]])), pch=16,
           xlab="", ylab="", axes=F, xaxs="i", ...)
      if(polygons=="season"){
        seasons_table <- shadeSeasons(min(data$timebin, na.rm=T), max(data$timebin, na.rm=T), interval)
        seasons_table$start <- strftime(seasons_table$start, "%Y-%m-%d", tz="UTC")
        seasons_table$end <- strftime(seasons_table$end, "%Y-%m-%d", tz="UTC")
        seasons_table$start <- sapply(seasons_table$start, function(x) min(which(as.character(days)==x)))
        seasons_table$end <- sapply(seasons_table$end, function(x) min(which(as.character(days)==x)))
        seasons_table$start[1] <- as.numeric(par("usr")[1])
        seasons_table$end[nrow(seasons_table)] <- par("usr")[2]
        rect(xleft=seasons_table$start, xright=seasons_table$end, ybottom=par("usr")[3], ytop=par("usr")[4], col=seasons_table$color, border=NA)
        seasons_legend <- seasons_table[!duplicated(seasons_table$season), c("season","color")]
        seasons_legend <- seasons_legend[order(match(seasons_legend$season, c("spring", "summer", "autumn", "winter"))),]
        moby:::legend2("right", legend=seasons_legend$season, fill=seasons_legend$color, bty="n", border="black",
                           inset=c(-0.18, 0), xpd=T, y.intersp=1.6, box.cex=c(1.6, 1.2), cex=0.8, horiz=F)
      } else if(polygons=="diel"){
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="grey98", border=NULL)
        x <- 1:length(daytimes_table$sunrises)
        polygon(c(x, rev(x)), c(daytimes_table$sunrises, rep(par("usr")[3], length(x))), col=adjustcolor("#CAD9EA", alpha.f=0.8), border=F)
        polygon(c(x, rev(x)), c(daytimes_table$sunsets, rep(par("usr")[4], length(x))), col=adjustcolor("#CAD9EA", alpha.f=0.8), border=F)
        polygon(c(x, rev(x)), c(daytimes_table$sunrises, rev(daytimes_table$sunsets)), col=adjustcolor("#FFFBE1", alpha.f=0.8), border=F)
      } else {
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.col, border=NULL)
      }

      # draw points
      plot_data$cex <- scales::rescale(plot_data$var, from=var_range, to=pt.cex)
      if(highlight.isolated==T){
        runs <- rle(plot_data$color)
        runs_length <-  rep(runs$lengths, runs$lengths)
        runs_value <- rep(runs$values, runs$lengths)
        runs_table <- data.frame(runs_length, runs_value)
        sorted_rows <- order(runs_table$runs_length, decreasing=T)
        plot_data <- plot_data[sorted_rows,]
      }
      x_indexes <- unlist(sapply(plot_data$day, function(x) which(as.character(days)==x), simplify=F))
      y_indexes <- unlist(sapply(plot_data$hour, function(x) which(colnames(plot_template[[i]])==x), simplify=F))
      points(x=x_indexes, y=y_indexes, pch=16, cex=plot_data$cex, col=adjustcolor(color.pal[plot_data$color], alpha.f=0.7))

      # draw legend
      size_labs <- pretty(var_range)
      size_labs <- unique(round(size_labs))
      size_labs <- size_labs[size_labs>=min(var_range) & size_labs<=max(var_range) & size_labs>0]
      size_scale <- scales::rescale(size_labs, from=var_range, to=pt.cex)
      legend("topright", inset=c(-0.12, 0), legend=size_labs, pch=16, pt.cex=size_scale,
             bty="n", col=adjustcolor("grey25", alpha.f=0.7), cex=0.7, y.intersp=1.6, xpd=T)
      if(!is.null(color.by)){
        if(!class(data[,color.by])[1] %in% c("numeric", "integer")){
          legend("bottomright", inset=c(-0.12, 0), pch=16, pt.cex=1.4, legend=levels(data[,color.by]), col=color.pal, bty="n",
                 cex=0.8, y.intersp=1.2, xpd=T)
        }
      }
    }

    # add title
    title(main=plot_title, cex.main=1.2, font=2, line=1)

    # add x axis to the bottom plots
    if(same_dates==T){
      if(i %in% bottom_plots){
        title(xlab="Date", cex.lab=1)
        axis(side=1, labels=disp_dates, at=disp_indexes, cex.axis=1)
        axis(side=1, labels=F, at=indexes, tck=-0.015, lwd.ticks=0.5)
      }
    }else{
      title(xlab="Date", cex.lab=1)
      axis(side=1, labels=disp_dates, at=disp_indexes, cex.axis=1)
      axis(side=1, labels=F, at=indexes, tck=-0.015, lwd.ticks=0.5)
    }

    # add x axis to the bottom plots
    if(i %in% right_plots){
      title(ylab="Hour", cex.lab=1)
      axis(2, at=hour_indexes[c(F,T)], labels=F, tck=-0.015, lwd.ticks=0.5)
      axis(2, at=hour_indexes[c(T,F)], labels=disp_hours[c(T,F)], cex.axis=1, las=1)
    }

    # draw grid
    if(grid==T) {
      usr <- par("usr")
      grid_x <- seq.int(usr[1], usr[2], length.out=length(days)+1)[-c(1, length(days) + 1)]
      grid_y <- seq.int(usr[3], usr[4], length.out=length(bins)+1)[-c(1, length(bins) + 1)]
      segments(x0=usr[1], x1=usr[2], y0=grid_y, lty="solid", lwd=0.05, col=grid.color)
      segments(x0=grid_x, y0=usr[3], y1=usr[4], lty="solid", lwd=0.05, col=grid.color)
    }
    # draw diel lines
    lines(1:length(daytimes_table$sunrises), daytimes_table$sunrises, lty=2)
    if(diel.lines==4){lines(1:length(daytimes_table$dawn), daytimes_table$dawn, lty=2)}
    lines(1:length(daytimes_table$sunsets), daytimes_table$sunsets, lty=2)
    if(diel.lines==4){lines(1:length(daytimes_table$dusks), daytimes_table$dusks, lty=2)}
    box()

    # draw color legend
    if(style=="raster"){
      fields::image.plot(y=1:length(bins), x=1:length(days), zlim=var_range, z=plot_template[[i]], legend.only=T,
                         col=raster_pal, add=T, graphics.reset=T, axis.args=list(cex.axis=0.7), legend.mar=3)
    } else if (!is.null(color.by)){
      if(class(data[,color.by])[1] %in% c("numeric", "integer")){
        color_labs <- pretty(c(min(plot_data$level), max(plot_data$level)), min.n=4)
        color_labs <- color_labs[color_labs>=min(plot_data$level) & color_labs<=max(plot_data$level)]
        shape::colorlegend(col=color.pal, zlim=range(plot_data$level), zval=color_labs,
                           posx=c(0.935, 0.95), posy=c(0, 0.45), main=color.by, main.cex=0.6, digit=1, cex=0.7)
      }
    }

  }

}



#######################################################################################################
#######################################################################################################
#######################################################################################################

aggFun <- function(data){
  if(class(data)[1] %in% c("numeric","integer")){
    return(mean(data, na.rm=T))
  }else{
    return(names(table(data))[which.max(table(data))])
  }
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
