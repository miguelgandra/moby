#######################################################################################################
# Function to generate chronograms (heatmaps) #########################################################
#######################################################################################################

#' Generate chronograms (heatmaps)
#'
#' @description Function to represent a given metric (nº detections, presences, individuals or co-occurrences)
#' per time bin. It generates 2-dimensional plots (hour x date) highlighting the chosen variable
#' either through color-coded cells or size-variable points, and illustrating the
#' variation in diel phases' hours across the study duration (annual variation of daylight time).
#'
#' @inheritParams setDefaults
#' @param data A data frame containing binned animal detections.
#' @param split.by Optional. If defined, plots are generated individually for each level
#' of this variable (e.g. species, ontogeny or habitat).
#' @param variables The type(s) of metric to plot. Accepted types: "detections", "individuals" and "co-occurrences".
#' @param style Style of the plot. Either "raster" (nº represented by color) or "points" (nº represented by size).
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
#' to display. Either 0 (no lines), 2 (corresponding to sunrise and sunset) or 4
#' (depicting dawn, sunrise, sunset, dusk). Defaults to 4.
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
#' @param cex.axis Determines the size of the text labels on the axes. Defaults to 1.
#' @param cex.legend Determines the size of the color legend. Defaults to 0.8.
#' @param cex.moon Determines the size of the moon phase symbols, when lunar.info is set to TRUE. Defaults to 2.5.
#' @param legend.intersp Vertical distances between legend elements
#' (in lines of text shared above/below each legend entry). Defaults to 1.2.
#' @param legend.cols Integer. The number of columns in which to set the legend items.
#' If NULL, it is set automatically based on the number of levels. Defaults to NULL.
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @param ... Further arguments passed to the \code{\link[graphics]{image}} function (if style="raster"),
#' or to the \code{\link[graphics]{plot}} function (if style="points").
#' @export


plotChronogram <- function(data, tagging.dates=getDefaults("tagging.dates"), variables="detections",
                           split.by=NULL, style="points", id.col=getDefaults("id"),
                           timebin.col=getDefaults("timebin"), station.col=getDefaults("station"),
                           color.by=NULL, color.pal=NULL, date.format="%d/%b", date.interval=4, date.start=1,
                           diel.lines=4, sunriset.coords, solar.depth=18, polygons=F, lunar.info=F,
                           background.col="white", grid=TRUE, grid.color="white", highlight.isolated=FALSE,
                           pt.cex=c(0.5, 3), uniformize.scale=FALSE, uniformize.dates=TRUE,
                           cex.axis=1, cex.legend=0.8, cex.moon=2.5, legend.intersp=1.2,
                           legend.cols=NULL, cols=NULL, ...) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates

  # validate additonal variables
  if(length(variables) == 0 || any(!variables %in% c("detections", "individuals", "co-occurrences"))) stop("Invalid variable(s) specified. Accepted values are: 'detections', 'individuals', or 'co-occurrences'.", call.=FALSE)
  if(length(pt.cex)!=2) stop("pt.cex should be a vector of length 2 (min and max)", call.=FALSE)

  # check if the lunar package is installed
  if(lunar.info==T & !requireNamespace("lunar", quietly=TRUE)){
    stop("The 'lunar' package is required but is not installed. Please install 'lunar' using install.packages('lunar') and try again.", call.=FALSE)
  }

  # print to console
  .printConsole("Generating chronogram(s)")


  #####################################################################################
  # Validate color-coding and color palette ###########################################
  #####################################################################################

  #######################################################
  # set color palette if a color.by variable was supplied (style=="points")
  if(!is.null(color.by) && style=="points"){

    # check if data contains the color.by variable
    if(!color.by %in% colnames(data)) stop("'color.by' variable not found in the supplied data", call.=FALSE)
    if(any(is.na(data[,color.by]))) stop("Missing values in 'color.by' variable",  call.=FALSE)

    # if class character, convert to factor
    if(inherits(data[,color.by], "character")) {
      data[,color.by] <- as.factor(data[,color.by])
      warning("'color.by' variable converted to factor", call.=FALSE)}

    # set color palette for categorical levels (style = "points")
    if(inherits(data[,color.by], "factor")){
      ngroups <- nlevels(data[,color.by])
      if(!is.null(color.pal)){
        if(style=="points" && inherits(color.pal, "function")) color.pal <- color.pal(ngroups)
        if(length(color.pal)!=ngroups)  warning("The number of colors doesn't match the number of levels in 'color.by' variable", call.=FALSE)
      }else{
        if(ngroups==3) color.pal <- c("#326FA5","#D73134","#1C8E43")
        else if (ngroups>3 & ngroups<10) color.pal <- .economist_pal(ngroups)
        else color.pal <- rainbow(ngroups)
      }
    # set color palette for numeric scale (style=="points")
    }else if(inherits(data[, color.by], "numeric")) {
      if(!is.null(color.pal)) {
        if(style == "points" && inherits(color.pal, "function")) {
          color.pal <- color.pal(100)
        }else{
          color.pal <- colorRampPalette(color.pal)(100)
          warning("The color palette was generated using the provided colors", call.=FALSE)
        }
      }else{
        color.pal <- .viridis_pal(100)
      }
    }
  #######################################################
  # set color palette if no color.by variable was supplied (style=="points")
  }else if(is.null(color.by) && style=="points"){
    if(!is.null(color.pal)) color.pal <- color.pal(1)
    else color.pal <- "grey25"
  }

  #######################################################
  # set color palette for style=="raster"
  if(style=="raster"){
    if(!is.null(color.pal) && !inherits(color.pal, "function")) color.pal <- colorRampPalette(color.pal)
    if(is.null(color.pal)) color.pal <- .viridis_pal
    if(!is.null(color.by)) warning("'color.by' argument is discarded when plot if of type raster", call.=FALSE)
  }


  #######################################################
  # set nº of columns in legend
  if(is.null(legend.cols)){
    legend.cols <- 1
    if(!is.null(color.by) && inherits(data[,color.by], "factor")){
      if(nlevels(data[,color.by])>=15) legend.cols<-2
    }
  }


  #####################################################################################
  # Aggregate data and calculate stats ################################################
  #####################################################################################

  # aggregate data in a list
  if(!is.null(split.by)){
    if(!inherits(data[,split.by], "factor")){
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
      if(is.null(tagging.dates)){first_dates <- min(data[,timebin.col], na.rm=T)}
      last_dates <- max(data[,timebin.col], na.rm=T)
    }else{
      if(!is.null(tagging.dates)){first_dates <- tagging.dates[group_ids]}
      if(is.null(tagging.dates)){first_dates <- min(data_group[,timebin.col], na.rm=T)}
      last_dates <- max(data_group[,timebin.col], na.rm=T)
    }

    # get time bins interval (in minutes)
    interval <- difftime(data[,timebin.col], dplyr::lag(data[,timebin.col]), units="min")
    interval <- as.numeric(min(interval[interval>0], na.rm=T))

    #round dates
    first_dates <- lubridate::floor_date(first_dates, unit="day")
    last_dates <- lubridate::ceiling_date(last_dates, unit="day")-60*60*(interval/60)

    # create complete seq of time-bins
    all_timebins <- seq.POSIXt(first_dates, last_dates, by=paste(interval, "mins"))

    # create day x time matrix
    all_timebins <- data.frame("timebin"=all_timebins)
    all_timebins$day <- strftime(all_timebins$timebin,"%Y-%m-%d", tz="UTC")
    all_timebins$time <- strftime(all_timebins$timebin,"%H:%M:%S", tz="UTC")
    days <- unique(all_timebins$day)
    bins <- unique(all_timebins$time)
    plot_matrix <- t(matrix(NA, nrow=length(bins), ncol=length(days)))
    rownames(plot_matrix) <- days
    colnames(plot_matrix) <- bins

    # aggregate data by timebin, station and ID
    if("detections" %in% colnames(data)){
      data_group$row <- data_group$detections
    }else{
      data_group$row <- 1:nrow(data_group)
    }
    detections_table <- stats::aggregate(data_group$row, by=list(data_group[,timebin.col], data_group[,id.col], data_group[,station.col]), length)
    colnames(detections_table) <- c("timebin", "id", "station", "detections")

    # calculate metrics
    if(var=="detections"){
      plot_data <- stats::aggregate(detections_table$detections, by=list(detections_table$timebin, detections_table$station), sum)
      colnames(plot_data) <- c("timebin", "station", "var")
    }else{
      plot_data <- stats::aggregate(detections_table$id, by=list(detections_table$timebin, detections_table$station), function(x) length(unique(x)))
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
      color_data <- stats::aggregate(data_group[,color.by], by=list(data_group[,timebin.col], data_group[,id.col], data_group[,station.col]), aggFun)
      colnames(color_data) <- c("timebin", "id", "station", "level")
      color_data <-  stats::aggregate(color_data$level, by=list(color_data$timebin, color_data$station), aggFun)
      colnames(color_data) <- c("timebin", "station", "level")
      color_data <- color_data[order(color_data$timebin),]
      if(inherits(color_data$level, "numeric")){
        color_data$color <- round(.rescale(color_data$level, to=c(1,100)))
      }else{
        color_data$level <- factor(color_data$level, levels=levels(data_group[,color.by]))
        color_data$color <- as.integer(color_data$level)
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
  # Extend dates (add left and right 4% margin) #######################################
  #####################################################################################

  # Define the extension process as a function
  extend_range <- function(template) {
    dates <- as.Date(rownames(template))
    date_range <- as.numeric(max(dates) - min(dates))
    extend_by <- floor(date_range*0.04)
    new_start <- min(dates) - extend_by
    new_end <- max(dates) + extend_by
    new_dates <- seq(new_start, new_end, by="day")
    new_template <- data.frame(matrix(NA, nrow=length(new_dates), ncol=ncol(template)))
    rownames(new_template) <- new_dates
    colnames(new_template) <- colnames(template)
    common_dates <- intersect(rownames(new_template), rownames(template))
    new_template[common_dates, ] <- template[common_dates, ]
    new_template <- new_template[order(new_dates),]
    return(new_template)
  }

  # Apply the function to each data frame in the list using lapply
  plot_template <- lapply(plot_template, extend_range)


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


  # set par margins
  if(style=="raster"){
    mars <- c(4,4,3,4)
  }else if (style=="points"){
    if(!is.null(color.by)){
      if(inherits(data[,color.by], "numeric")) mars <- c(4,4,3,12)
      else
        if (legend.cols==1) mars <- c(4,4,3,10)
        else mars <- c(4,4,3,14)
        labels_size <- max(unlist(sapply(levels(data[,color.by]), nchar)))
        if(labels_size>8) mars <- c(4,4,3,14)
    }else{
      mars <- c(4,4,3,4)
    }
  }
  if(lunar.info==T){
    mars[2] <- 6
  }


  #####################################################################################
  # Plot lunar graphics if lunar.illumination is set true #############################
  #####################################################################################

  if(lunar.info==T){
    timebins <- as.POSIXct(days, "%Y-%m-%d", tz="UTC")
    lunar_illumination <- data.frame("timebin"=timebins, "fraction"=lunar::lunar.illumination(timebins))
    lunar_illumination$phase <- as.character(lunar::lunar.phase(timebins, name=T))
    monitoring_period <- round(as.numeric(difftime(max(timebins), min(timebins), units="days")))
    if(monitoring_period>10*30){
      warning(paste0("Time range (", monitoring_period, " days) may be too large to display lunar illumination accurately. Consider splitting the data into smaller intervals to better visually inspect lunar patterns."), call.=FALSE)
    }

    for(l in 1:cols){
      par(mar=c(3,6,6,mars[4]), xpd=T)
      xcoords <- c(1, seq(1, nrow(lunar_illumination)), nrow(lunar_illumination))
      ycoords <- c(0, lunar_illumination$fraction, 0)
      plot(lunar_illumination$fraction, type="l", lwd=0.5, axes=F, ylim=c(0,1), xlim=c(1, nrow(lunar_illumination)),
           xlab="", ylab="Lunar\nIllumination", xaxs="i", yaxs="i")
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="#3C3C3C", border=NULL)
      polygon(x=xcoords, y=ycoords, col=adjustcolor("#FFFBEB", alpha.f=0.8), lwd=0.6)
      axis(2, at=seq(0,1,by=0.5), labels=sprintf("%.1f", seq(0,1,by=0.5)), cex.axis=cex.axis, las=1)
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
      moon_indexes <- plyr::join(moon_indexes, moon_lookup, by="phase", type="left")
      moon_indexes$fraction2 <- abs(1- moon_indexes$fraction1)
      moon_indexes <- moon_indexes[moon_indexes$phase %in% c("Full","New"),]
      for(i in 1:nrow(moon_indexes)){
        if(moon_indexes$phase[i]!="Waning"){moon_color <- c("#FFFBEB", "#3C3C3C")}
        if(moon_indexes$phase[i]=="Waxing"){moon_color <- c("#3C3C3C", "#FFFBEB")}
        plotrix::floating.pie(xpos=moon_indexes$index[i], ypos=1.4, x=as.numeric(moon_indexes[i,c("fraction1", "fraction2")]),
                              startpos=90*pi/180, col=moon_color, lwd=0.2, radius=cex.moon, xpd=T)
        plotrix::floating.pie(xpos=moon_indexes$index[i], ypos=1.4, x=as.numeric(moon_indexes[i,c("fraction1", "fraction2")]),
                              startpos=90*pi/180, col=moon_color, lwd=0.2, border=NA, radius=cex.moon, xpd=T)
      }
    }
  }

  #####################################################################################
  # Generate level plots ##############################################################
  #####################################################################################

  ############################################################################
  # iterate through each data group/variable #################################

  par(mar=mars, mgp=c(3,0.7,0))

  for(i in 1:nrow(plot_combs)){

    # set variables
    group <- plot_combs$group[i]
    var <-  plot_combs$variable[i]
    plot_data <- plot_data_list[[i]]
    plot_title <- switch(var, detections="n\u00ba of detections", individuals="n\u00ba of individuals",
                         "co-occurrences"="n\u00ba of co-occurring animals")
    if(!is.null(split.by)){plot_title <- paste(paste0(toupper(substring(group, 1, 1)), substring(group, 2)), "-", plot_title)}

    # get unique days and hours
    days <- as.POSIXct(rownames(plot_template[[i]]), "%Y-%m-%d", tz="UTC")
    bins <- colnames(plot_template[[i]])

    # get diel phase timelines
    daytimes_table <- getSunTimes(sunriset.coords, min(days), max(days), by="%Y-%m-%d", solar.depth)
    daytimes_table[,-1] <- daytimes_table[,-1] * (60/interval)
    daytimes_table[,-1] <- daytimes_table[,-1] + 1

    # prepare date variables
    all_dates <- strftime(days, date.format)
    consec_dates <- rle(all_dates)
    consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
    unique_dates <- unique(consec_dates)
    disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
    indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))[-1]
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
      if(is.null(color.pal)){raster_pal <- .viridis_pal(max(var_range))}
      if(!is.null(color.pal)){raster_pal <- color.pal(max(var_range))}
      plot_matrix <- reshape2::dcast(plot_data, formula="day~hour", value.var="var", fill=0, fun.aggregate=max, drop=F)
      missing_hours <- base::setdiff(colnames(plot_template[[i]]), colnames(plot_matrix))
      missing_days <- base::setdiff(rownames(plot_template[[i]]), plot_matrix$day)
      if(length(missing_hours)>0) {
        plot_matrix[missing_hours] <- NA
        plot_matrix <- plot_matrix[,order(colnames(plot_matrix))]
      }
      if(length(missing_days)>0) {
        missing_days <- data.frame("day"=missing_days)
        plot_matrix <- dplyr::bind_rows(plot_matrix, missing_days)
      }
      rownames(plot_matrix) <- plot_matrix$day
      plot_matrix <- plot_matrix[,-which(colnames(plot_matrix)=="day")]
      plot_matrix <- plot_matrix[order(rownames(plot_matrix)),]
      plot_matrix <- as.matrix(plot_matrix)
      image(y=1:length(bins), x=1:length(days), z=plot_matrix,
            zlim=var_range, xlim=xlim, xlab="", main="", ylab="",
            col=raster_pal, axes=F, cex.lab=0.9, ...)
    }

    ##########################################################
    # else, if style is set to points use "plot" ############
    if(style=="points"){

      # create empty plot and fill background
      plot(0,0, type="n", xlim=c(1, nrow(plot_template[[i]])), ylim=c(1, ncol(plot_template[[i]])),
           pch=16, xlab="", ylab="", axes=F, xaxs="i")
      coords_list <- list()
      if(polygons=="season"){
        seasons_table <- shadeSeasons(min(days), max(days), interval)
        seasons_table$start <- strftime(seasons_table$start, "%Y-%m-%d", tz="UTC")
        seasons_table$end <- strftime(seasons_table$end, "%Y-%m-%d", tz="UTC")
        seasons_table$start <- sapply(seasons_table$start, function(x) min(which(as.character(days)==x)))
        seasons_table$end <- sapply(seasons_table$end, function(x) min(which(as.character(days)==x)))
        seasons_table$start[1] <- as.numeric(par("usr")[1])
        seasons_table$end[nrow(seasons_table)] <- par("usr")[2]
        rect(xleft=seasons_table$start, xright=seasons_table$end, ybottom=par("usr")[3], ytop=par("usr")[4], col=seasons_table$color, border=NA)
        seasons_legend <- seasons_table[!duplicated(seasons_table$season), c("season","color")]
        seasons_legend <- seasons_legend[order(match(seasons_legend$season, c("spring", "summer", "autumn", "winter"))),]
        coords <- .legend(x=par("usr")[2], y=par("usr")[4], legend=seasons_legend$season, fill=seasons_legend$color, bty="n", border="black",
                                 xpd=T, y.intersp=legend.intersp+0.2, box.cex=c(1.6, 1.2), cex=cex.legend, horiz=F)
        coords_list <- c(coords_list, list(coords))
      } else if(polygons=="diel"){
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="grey98", border=NULL)
        x <- 1:length(daytimes_table$sunrises)
        if(diel.lines==2){
          polygon(c(x, rev(x)), c(daytimes_table$sunrises, rep(par("usr")[3], length(x))), col=adjustcolor("#CAD9EA", alpha.f=0.8), border=F)
          polygon(c(x, rev(x)), c(daytimes_table$sunsets, rep(par("usr")[4], length(x))), col=adjustcolor("#CAD9EA", alpha.f=0.8), border=F)
          polygon(c(x, rev(x)), c(daytimes_table$sunrises, rev(daytimes_table$sunsets)), col=adjustcolor("#FFFBE1", alpha.f=0.8), border=F)
          coords <- .legend(x=par("usr")[2], y=par("usr")[4], legend=c("day","night"), fill=adjustcolor(c("#FFFBE1","#CAD9EA"), alpha.f=0.8),
                            bty="n", xpd=T, y.intersp=legend.intersp+0.2, box.cex=c(1.6, 1.2), cex=cex.legend, horiz=F)
        }
         if(diel.lines==4){
           polygon(c(x, rev(x)), c(daytimes_table$dawns, rep(par("usr")[3], length(x))), col=adjustcolor("#CAD9EA", alpha.f=0.8), border=F)
           polygon(c(x, rev(x)), c(daytimes_table$dusks, rep(par("usr")[4], length(x))), col=adjustcolor("#CAD9EA", alpha.f=0.8), border=F)
           polygon(c(x, rev(x)), c(daytimes_table$sunrises, rev(daytimes_table$sunsets)), col=adjustcolor("#FFFBE1", alpha.f=0.8), border=F)
           polygon(c(x, rev(x)), c(daytimes_table$dawns, rev(daytimes_table$sunrises)), col=adjustcolor("#EDE8EC", alpha.f=0.8), border=F)
           polygon(c(x, rev(x)), c(daytimes_table$sunsets, rev(daytimes_table$dusks)), col=adjustcolor("#EDE8EC", alpha.f=0.8), border=F)
           coords <- .legend(x=par("usr")[2], y=par("usr")[4], legend=c("day","night", "crepuscule"), fill=adjustcolor(c("#FFFBE1","#CAD9EA","#EDE8EC"), alpha.f=0.8),
                             bty="n", xpd=T, y.intersp=legend.intersp+0.2, box.cex=c(1.6, 1.2), cex=cex.legend, horiz=F)
         }
        coords_list <- c(coords_list, list(coords))
      } else {
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.col, border=NULL)
      }

      # draw points
      plot_data$cex <- .rescale(plot_data$var, from=var_range, to=pt.cex)
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

      # draw legend for point sizes
      size_labs <- pretty(var_range)
      size_labs <- unique(round(size_labs))
      size_labs <- size_labs[size_labs>=min(var_range) & size_labs<=max(var_range) & size_labs>0]
      size_scale <- .rescale(size_labs, from=var_range, to=pt.cex)
      if(length(coords_list)==0){legend.y<-par("usr")[4]
      }else{legend.y<-abs(coords_list[[1]]$rect$top-coords$rect$h)-0.5}
      legend.x <- par("usr")[2] + (par("usr")[2]-par("usr")[1])*0.01
      coords <- legend(x=legend.x, y=legend.y, legend=size_labs, pch=16, pt.cex=size_scale,
                       bty="n", col=adjustcolor("grey25", alpha.f=0.7), cex=cex.legend, y.intersp=legend.intersp+0.2, xpd=T)
      coords_list <- c(coords_list, list(coords))

      # plot color.by legend
      if(!is.null(color.by)){
        if(inherits(data[,color.by], "factor")) {
          # plot legend
          legend.x <- max(unlist(lapply(coords_list, function(x) x$rect$left + x$rect$w + 5)))
          legend.y <- par("usr")[4]
          legend(x=legend.x, y=legend.y, pch=16, pt.cex=1.4, legend=levels(data[,color.by]), col=color.pal, bty="n",
                 cex=cex.legend, y.intersp=legend.intersp, ncol=legend.cols, xpd=T)
        }
      }
    }

    ##########################################################
    # add additional elements ################################

    # add title
    title(main=plot_title, cex.main=1.2, font=2, line=1)

    # add x axis to the bottom plots
    if(same_dates==T){
      if(i %in% bottom_plots){
        title(xlab="Date", cex.lab=1, line=2.5)
        axis(side=1, labels=disp_dates, at=disp_indexes, cex.axis=cex.axis)
        axis(side=1, labels=F, at=indexes, tck=-0.015, lwd.ticks=0.5)
      }
    }else{
      title(xlab="Date", cex.lab=1, line=2.5)
      axis(side=1, labels=disp_dates, at=disp_indexes, cex.axis=cex.axis)
      axis(side=1, labels=F, at=indexes, tck=-0.015, lwd.ticks=0.5)
    }

    # add x axis to the bottom plots
    if(i %in% right_plots){
      title(ylab="Hour", cex.lab=1, line=3)
      axis(2, at=hour_indexes[c(F,T)], labels=F, tck=-0.015, lwd.ticks=0.5)
      axis(2, at=hour_indexes[c(T,F)], labels=disp_hours[c(T,F)], cex.axis=cex.axis, las=1)
    }

    # draw grid
    if(grid==T) {
      segments(x0=par("usr")[1], x1=par("usr")[2], y0=hour_indexes, lty="solid", lwd=0.05, col=grid.color)
      segments(x0=indexes, y0=par("usr")[3], y1=par("usr")[4], lty="solid", lwd=0.05, col=grid.color)
    }

    # draw diel lines
    if(diel.lines>0){
      lines(1:length(daytimes_table$sunrises), daytimes_table$sunrises, lty=2)
      lines(1:length(daytimes_table$sunsets), daytimes_table$sunsets, lty=2)
    }
    if(diel.lines==4){
      lines(1:length(daytimes_table$dawn), daytimes_table$dawn, lty=2)
      lines(1:length(daytimes_table$dusks), daytimes_table$dusks, lty=2)
    }

    # draw box
    box()

    # draw color legend
    if(style=="raster"){
      scale_labs <- pretty(var_range, min.n=4)
      scale_labs <- unique(round(scale_labs))
      fields::image.plot(y=1:length(bins), x=1:length(days), zlim=var_range, z=plot_matrix, legend.only=T, legend.shrink=0.7, legend.width=1,
                        col=raster_pal, add=T, graphics.reset=T, axis.args=list(scale_labs, cex.axis=cex.axis), legend.mar=3, legend.cex=cex.legend)

    }else if(!is.null(color.by) && inherits(data[,color.by], "numeric")){
      # plot color scale
      dx <- par("usr")[2]-par("usr")[1]
      legend.x <- max(unlist(lapply(coords_list, function(x) x$rect$left+x$rect$w)))
      legend.x <- legend.x + dx*0.025
      legend.x <- c(legend.x, legend.x + dx*0.016)
      legend.x <- graphics::grconvertX(legend.x, from="user", to="ndc")
      scale_labs <- pretty(data[,color.by], min.n=4)
      scale_labs <- scale_labs[scale_labs>=min(data[,color.by]) & scale_labs<=max(data[,color.by])]
      digits <- max(.decimalPlaces(scale_labs))
      .colorlegend(col=color.pal, zlim=range(data[,color.by], na.rm=T), zval=scale_labs, digit=digits, xpd=T,
                         posx=c(legend.x[1], legend.x[2]), posy=c(0.4,0.9), main=color.by, main.cex=cex.legend, cex=cex.legend-0.1)
    }


  }

}



#######################################################################################################
#######################################################################################################
#######################################################################################################

aggFun <- function(data){
  if(inherits(data, "numeric")){
    return(mean(data, na.rm=T))
  }else{
    return(names(table(data))[which.max(table(data))])
  }
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
