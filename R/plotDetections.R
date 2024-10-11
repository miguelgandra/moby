##################################################################################################
## Plot detections ###############################################################################
##################################################################################################

#' Detection plots
#'
#' @description Plots color-coded detections over time and date, independently for each individual.
#'
#' @inheritParams setDefaults
#' @param data  A data frame containing animal detections.
#' @param tag.durations Optional. A numeric vector containing the estimated battery
#' duration of the deployed tags (in days). The length of this vector should match the number of
#' unique animal IDs, and the values must be in the same order as the ID levels.
#' Alternatively, if a single value is provided, it will be applied to all IDs.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param color.by Optional. Variable defining the color group of individual detections.
#' Can be used for example to depict detections of different receivers or in different habitats.
#' @param color.pal Color palette for the points.
#' @param date.format Date-time format (as used in \code{\link[base]{strptime}}),
#' defining the x-axis labels. Defaults to month-year ("%b/%y").
#' @param date.interval Number defining the interval between each
#' displayed date (x-axis label). Defaults to 4.
#' @param date.start Integer defining the first displayed date (can be used in combination
#'  with 'date.interval" to better control the x-axis labels). Defaults to 1.
#' @param sunriset.coords A SpatialPoints or matrix object containing longitude and
#' latitude coordinates (in that order) at which to estimate sunrise and sunset times.
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param discard.missing If true, only individuals with detections are included.
#' If false, empty plots are drawn for missing individuals.
#' @param background.color Color for the plot background. Defaults to "gray96".
#' @param highlight.isolated If a 'color.by' variable is defined, detections are ordered
#' and plotted according to their "density", with isolated detections being brought forward
#' to prevent them from being hidden behind points' clouds. Defaults to True.
#' @param pt.cex Expansion factor(s) for the points (detections). Defaults to 1.6.
#' @param diel.lines Number indicating the number of diel phase lines (boundaries)
#' to display. Either 0 (no lines), 2 (corresponding to sunrise and sunset) or 4
#' (depicting dawn, sunrise, sunset, dusk). Defaults to 4.
#' @param grid If true, a grid is plotted (horizontal and vertical guides across the entire plot). Defaults to False.
#' @param grid.color Color of the grid lines. Defaults to white.
#' @param cex.title Determines the size of the plot title (animal ID). Defaults to 2.2.
#' @param cex.lab Determines the size of the axes titles. Defaults to 1.8.
#' @param cex.axis Determines the size of the tick mark labels on the axes. Defaults to 1.4.
#' @param cex.legend Determines the size of the color legend. Defaults to 1.6.
#' @param cols Number of columns in the final panel. Defaults to 2.
#' @param legend.cols Number of columns in the legend (when 'color.by' is defined).  Defaults to 3.
#' @export


plotDetections <- function(data, id.col=getDefaults("id"), datetime.col=getDefaults("datetime"),
                           tagging.dates=getDefaults("tagging.dates"), tag.durations=NULL,
                           id.groups=NULL, discard.missing=T, color.by=NULL, color.pal=NULL,
                           date.format="%b/%y", date.interval=4, date.start=1, sunriset.coords,
                           diel.lines=2, background.color="gray96", highlight.isolated=T,
                           pt.cex=1.6, grid=F, grid.color="white", cex.title=2.2,
                           cex.lab=1.8, cex.axis=1.4, cex.legend=1.6, cols=2, legend.cols=3){

  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  tag.durations <- reviewed_params$tag.durations

  # create temporary color.by column
  if(is.null(color.by)){
    data$temp_group <- as.factor(1)
    color.by <- "temp_group"
  }

  # check nº of colors in supplied color palette
  if(!is.null(color.pal) & length(color.pal)<nlevels(data[,color.by])){
    warning("The number of supplied colors doesn't match number of group levels", call.=FALSE)
  }

  # print to console
  .printConsole("Plotting detections")


  ##############################################################################
  # Prepare data ###############################################################
  ##############################################################################

  # define single id.group if needed
  if(is.null(id.groups)){
    id.groups <- list(levels(data[,id.col]))
  }

  # discard missing animals if required
  missing_IDs <- which(table(data[,id.col])==0)
  if(discard.missing==T){
    if(length(missing_IDs)>0){
      tagging.dates <- tagging.dates[-missing_IDs]
      tag.durations <- tag.durations[-missing_IDs]
      data[,id.col] <- droplevels(data[,id.col])
      id.groups <- lapply(id.groups, function(x) x[!x %in% names(missing_IDs)])
    }
  }

  # get nº of animals
  nindividuals <- nlevels(data[,id.col])

  # assign hour and day
  data <- data[order(data[,datetime.col]),]
  data$day <- lubridate::floor_date(data[,datetime.col], "day")
  data$hour <- strftime(data[,datetime.col], "%H:%M:%S", tz="UTC")
  data$hour <- as.numeric(lubridate::hms(data$hour), "hours")
  data_individual <- split(data, f=data[,id.col])

  # estimate tag lifetime dates
  if(!is.null(tag.durations)){
    end.dates <- NA
    for(e in 1:length(tagging.dates)){
      end.dates[e] <- tagging.dates[e] + tag.durations[e]*60*60*24
    }
  }

  ##############################################################################
  # Set layout variables #######################################################
  ##############################################################################

  # estimate sunrise/sunset times
  start <- min(data$day)
  end <- max(data$day)
  daytimes_table <- getSunTimes(sunriset.coords, start, end, by="%Y-%m-%d")
  daytimes_table$interval <- as.POSIXct(daytimes_table$interval, "%Y-%m-%d", tz="UTC")

  # prepare date variables
  complete_dates <- seq.POSIXt(start, end, by="day", tz="UTC")
  all_dates <- strftime(complete_dates, date.format)
  consec_dates <- rle(all_dates)
  consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
  unique_dates <- unique(consec_dates)
  indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
  detec_dates <- strftime(seq.POSIXt(min(tagging.dates, na.rm=T), max(data[,datetime.col], na.rm=T), "day"), date.format)
  start <- min(which(sub("\\_.*", "", unique_dates)==detec_dates[1]))
  end <- max(which(sub("\\_.*", "", unique_dates)==detec_dates[length(detec_dates)]))
  indexes <- indexes[start:end]
  unique_dates <- unique_dates[start:end]
  disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
  disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
  disp_dates <- sub("\\_.*", "", disp_dates)

  # define hour labels
  hour_labels <- paste0(formatC(0:24, 1, flag=0),"h")[c(TRUE, FALSE)]

  # set color palette if required
  if(is.null(color.pal)){color.pal <- rainbow(nlevels(data[,color.by]))}

  # set layout grid
  layout_params <- .setLayout(cols, id.groups, plots.height=6, dividers.height=1, legend=TRUE)
  nplots <- max(layout_params$matrix, na.rm=T)
  layout_params$matrix[is.na(layout_params$matrix)] <- nplots + 1
  layout(mat=layout_params$matrix, heights=layout_params$heights)

  # set outer margins
  oma <- if (length(id.groups) > 1) c(1,8,1,2) else c(1,4,1,2)
  par(mar=c(3,2,3,0), oma=oma)


  ##############################################################################
  # Draw plots #################################################################
  ##############################################################################

  # iterate through each individual
  for (i in 1:(nplots-1)) {

    selected_id <- levels(data[,id.col])[i]
    data_plot <- data_individual[[selected_id]]

    #################################################################
    # plot detections ###############################################

    # set plot
    plot(x=data$day, y=data$hour, type="n", axes=F, xlim=range(complete_dates), xaxs="i",
         ylim=c(0,24), xlab="", ylab="", main=selected_id, cex.main=cex.title)
    # add background
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col=background.color)
    # add grid
    if(grid==TRUE) abline(h=0:24, v=complete_dates[disp_indexes], lwd=0.05, col=grid.color)
    # add horizontal x-axis
    axis(side=1, labels=disp_dates, at=complete_dates[disp_indexes], cex.axis=cex.axis)
    axis(side=1, labels=F, at=complete_dates[indexes], tck=-0.015, lwd.ticks=0.5)
    # add vertical y-axis
    if(i %in% layout_params$matrix[,1]){
      #mtext(text="Hour", side=2, line=4.5, cex=cex.lab)
      title(ylab="Hour", cex.lab=cex.lab, line=4, xpd=NA)
      axis(2, labels=hour_labels, at=seq(0, 24, by=2), tck=-0.03, cex.axis=cex.axis, las=1)
      axis(2, labels=F, at=0:24, tck=-0.015, lwd.ticks=0.5)
    }
    # draw points
    if(nrow(data_plot)>0){
      # if highlight.isolated, reorder detections using run length encoding
      if(highlight.isolated==T){
        data_plot <- data_plot[!is.na(data_plot[,color.by]),]
        runs <- rle(as.character(data_plot[,color.by]))
        data_plot$run_length <- rep(runs$lengths, runs$lengths)
        data_plot$run_value <- rep(runs$values, runs$lengths)
        data_plot <- data_plot[!is.na(data_plot$run_value),]
        data_plot <- data_plot[order(data_plot$run_length, decreasing=T),]
      }
      points(x=data_plot$day, y=data_plot$hour, col=color.pal[data_plot[,color.by]], pch=16, cex=pt.cex)
    }
    # draw tagging date line
    abline(v=tagging.dates[i], lwd=1.4)
    # draw tag end line
    if(!is.null(tag.durations)) abline(v=end.dates[i], lwd=1.4)
    # draw diel lines
    if(diel.lines>0){
      lines(x=daytimes_table$interval, y=daytimes_table$sunrises, lty=2)
      lines(x=daytimes_table$interval, y=daytimes_table$sunsets, lty=2)
    }
    if(diel.lines==4){
      lines(x=daytimes_table$interval, y=daytimes_table$dawns, lty=2)
      lines(x=daytimes_table$interval, y=daytimes_table$dusks, lty=2)
    }
    #draw box
    box()
  }

  #################################################################
  # add 'color.by' legend #########################################
  par(mar=c(1,1,3,1))
  plot.new()
  legend("center", legend=levels(data[,color.by]), pch=16, col=color.pal, ncol=legend.cols,
         cex=cex.legend, bg=background.color, pt.cex=cex.legend+0.6, y.intersp=1.4)

  #################################################################
  # add id.group labels ###########################################
  if (length(id.groups)>1) {
    label_pos <- layout_params$group_positions
    layout_height <- sum(layout_params$heights)
    label_pos <- grconvertY(1 - (label_pos / layout_height), "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels = names(id.groups),
         srt = 90, cex = cex.title + 0.2, font = 2, xpd = NA, adj = c(0.5, 0.5))
  }

}

##################################################################################################
##################################################################################################
##################################################################################################
