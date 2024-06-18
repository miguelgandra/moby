##################################################################################################
## Plot detections ###############################################################################
##################################################################################################

#' Detection plots
#'
#' @description Plots color-coded detections over time and date, independently for each individual.
#'
#' @inheritParams setDefaults
#' @param data  A data frame containing animal detections.
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
#' @param pt.cex Expansion factor(s) for the points (detections). Defaults to 1.5.
#' @param diel.lines Number indicating the number of diel phase lines (boundaries)
#' to display. Either 0 (no lines), 2 (corresponding to sunrise and sunset) or 4
#' (depicting dawn, sunrise, sunset, dusk).
#' @param grid If true, a grid is plotted (horizontal and vertical guides across the entire plot). Defaults to False.
#' @param grid.color Color of the grid lines. Defaults to white.
#' @param cols Number of columns in the final panel (passed to the mfrow argument). Defaults to 2.
#' @param legend.cols Number of columns in the legend (when 'color.by' is defined).  Defaults to 3.
#' @export


plotDetections <- function(data, id.col=getDefaults("id"), datetime.col=getDefaults("datetime"),
                           tagging.dates=getDefaults("tagging.dates"), tag.durations=NULL,
                           id.groups=NULL, discard.missing=T, color.by=NULL, color.pal=NULL,
                           date.format="%b/%y", date.interval=4, date.start=1, sunriset.coords,
                           diel.lines=2, background.color="gray96", highlight.isolated=T, pt.cex=1.5,
                           grid=F, grid.color="white", cols=2, legend.cols=3){

  ##################################################################################
  # Initial checks #################################################################

  # perform argument checks
  data <- .validateArguments()

  # check if data contains color.by
  if(is.null(color.by)){
    data$temp_group <- as.factor(1)
    color.by <- "temp_group"
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {warning("Some of the ID(s) in id.groups don't match the IDs in the data")}
    data <- data[data[,id.col] %in% unlist(id.groups),]
    tagging.dates <- tagging.dates[match(unlist(id.groups), levels(data[,id.col]))]
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
  }

  if(!is.null(color.pal) & length(color.pal)<nlevels(data[,color.by])){
    warning("The number of supplied colors doesn't match number of group levels")
  }

  # check tag durations
  if(!is.null(tag.durations)){
    if(length(tag.durations)>1 && length(tag.durations)!=nlevels(data[,id.col])){
      stop("Incorrect number of tag.durations. Must be either a single value or
           a vector containing the estimated tag duration for each individual")
    }
    if(length(tag.durations)==1){
      tag.durations <- rep(tag.durations, nlevels(data[,id.col]))
    }
  }

  # print to console
  .printConsole("Plotting detections")


  ################################################################################
  # Prepare data #################################################################

  if(discard.missing==T){
    missing_IDs <- which(table(data[,id.col])==0)
    if(length(missing_IDs)>0){
      tagging.dates <- tagging.dates[-missing_IDs]
      data[,id.col] <- droplevels(data[,id.col])
    }
  }
  nindividuals <- nlevels(data[,id.col])

  data_table <- createWideTable(data, timebin.col=timebin.col, start.dates=tagging.dates, value.col=color.by)
  data_flat <- reshape2::melt(data_table, id.vars=timebin.col)
  colnames(data_flat)[1:3] <- c("timebin", "id", color.by)
  data_flat$hour <- as.numeric(format(data_flat$timebin, "%H")) + as.numeric(format(data_flat$timebin, "%M"))/60
  data_flat$date <-strftime(data_flat$timebin, "%Y-%m-%d", tz="UTC")
  data_flat$date <- as.POSIXct(data_flat$date, "%Y-%m-%d" , tz="UTC")
  data_flat[,color.by] <- factor(data_flat[,color.by], levels=levels(data[,color.by]))
  data_individual <- split(data_flat, f=data_flat$id)
  id_order <- match(levels(data[,id.col]), names(data_individual))
  data_individual <- data_individual[id_order]


  ################################################################################
  # Set layout variables #########################################################

  # estimate sunrise/sunset times
  start <- data_table$timebin[1]
  end <- data_table$timebin[nrow(data_table)]
  daytimes_table <- getSunTimes(sunriset.coords, start, end, by="%Y-%m-%d")
  daytimes_table$interval <- as.POSIXct(daytimes_table$interval, "%Y-%m-%d", tz="UTC")

  # set date labels
  all_dates <- strftime(data_table$timebin, date.format)
  consec_dates <- rle(all_dates)
  consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
  unique_dates <- unique(consec_dates)
  disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
  indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
  disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
  disp_dates <- sub("\\_.*", "", disp_dates)

  # define hour labels
  hour_labels <- paste0(formatC(0:24, 1, flag=0),"h")[c(TRUE, FALSE)]

  # set color palette if required
  if(is.null(color.pal)){color.pal <- rainbow(nlevels(data$station))}

  # set layout variables
  if(!is.null(id.groups)){
    group_ids_selected <- lapply(id.groups, function(x) x[x %in% levels(data[,id.col])])
    group_numbers <- lapply(group_ids_selected,  length)
    group_rows <- lapply(group_numbers, function(x) ceiling(x/cols))
    rows <- do.call("sum", group_rows)
    group_plots <- lapply(group_rows, function(x) x*cols)
    animal_indexes <- mapply(function(nids, nplots) {if(nids<nplots){c(1:nids, rep(NA, nplots-nids))}else{1:nids}},
                             nids=group_numbers, nplots=group_plots, SIMPLIFY=F)
    for(i in 2:length(animal_indexes)){animal_indexes[[i]]<-animal_indexes[[i]]+max(animal_indexes[[i-1]], na.rm=T)}
    animal_indexes <- unlist(animal_indexes, use.names=F)
  } else{
    rows <- ceiling(nindividuals/cols)
    nplots <- rows*cols
    if(nindividuals<nplots){
      animal_indexes <- c(1:nindividuals, rep(NA, nplots-nindividuals))
    }else{
      animal_indexes <- 1:nindividuals
    }
  }

  plot_layout <- matrix(animal_indexes, nrow=rows, ncol=cols, byrow=T)
  if(!is.na(plot_layout[length(plot_layout)])){rows<-rows+1}
  plot_ids <- as.integer(apply(plot_layout, 1, function(x) x))
  par(mfrow=c(rows, cols), mar=c(3,2,3,0), oma=c(6,10,1,1))

  #if(same.scale==T){par(mfrow=c(rows, cols), mar=c(2,1,0,0), oma=c(6,6,1,1))}
  #if(same.scale==F){par(mfrow=c(rows, cols), mar=c(2,4,0,4), oma=c(6,6,1,1))}


  ################################################################################
  # Generate plots ###############################################################

  # iterate through each individual
  for (i in plot_ids) {

    # if NA, §add empty plot and go to next
    if(is.na(i)) {
      plot.new()
      next
    }

    selected_id <- levels(data[,id.col])[i]
    data_plot <- data_individual[[selected_id]]
    plot(x=data_plot$date, y=data_plot$hour, type="n", axes=F,
         ylim=c(0,24), xlab="", ylab="", main=unique(data_plot$id), cex.main=2.5)
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col=background.color)
    if(grid==T){
      abline(h=0:24, v=data_table$timebin[disp_indexes], lwd=0.05, col=grid.color)
      #grid(nx=length(disp_dates), ny=24, lty="solid", lwd=0.05, col=grid.color)
    }
    # if highlight.isolated, reorder detections using run length encoding,
    # plotting isolated detections in front
    if(highlight.isolated==T){
      data_plot <- data_plot[!is.na(data_plot[,color.by]),]
      runs <- rle(as.character(data_plot[,color.by]))
      data_plot$run_length <- rep(runs$lengths, runs$lengths)
      data_plot$run_value <- rep(runs$values, runs$lengths)
      data_plot <- data_plot[!is.na(data_plot$run_value),]
      data_plot <- data_plot[order(data_plot$run_length, decreasing=T),]
    }
    points(x=data_plot$date, y=data_plot$hour, col=color.pal[data_plot[,color.by]], pch=16, cex=pt.cex)
    axis(side=1, labels=disp_dates, at=data_table$timebin[disp_indexes], tck=-0.03, cex.axis=1.6)
    axis(side=1, labels=F, at=data_table$timebin[indexes], tck=-0.015, lwd.ticks=0.5)

    # y axis
    if(i %in% plot_layout[,1]){
      mtext(text="Hour", side=2, line=4.5, cex=1.3)
      axis(2, labels=hour_labels, at=seq(0, 24, by=2), tck=-0.03, cex.axis=1.6, las=1)
      axis(2, labels=F, at=0:24, tck=-0.015, lwd.ticks=0.5)
    }

    abline(v=tagging.dates[i], lwd=1.4)
    lines(x=daytimes_table$interval, y=daytimes_table$sunrises, lty=2, col="black")
    lines(x=daytimes_table$interval, y=daytimes_table$sunsets, lty=2, col="black")
    if(diel.lines==T){
      lines(x=daytimes_table$interval, y=daytimes_table$dawns, lty=2, col="black")
      lines(x=daytimes_table$interval, y=daytimes_table$dusks, lty=2, col="black")
    }
    box()
  }

  # create empty plot with 'color.by' legend
  par(mar=c(1,1,1,1))
  if(!is.na(i)){plot.new()}
  legend("center", legend=levels(data[,color.by]), pch=16, col=color.pal, ncol=legend.cols,
         cex=1.8, bg=background.color, pt.cex=2.4)

  # add id.group names if required
  if(!is.null(id.groups)){
    label_pos <- rev(unlist(lapply(group_rows, function(x) x/2)))
    label_pos <- cumsum(rev(group_rows)) - label_pos
    label_pos <- grconvertY(label_pos/rows, "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels=rev(names(id.groups)),
         srt=90, cex=2.8, font=2, xpd=NA, adj=c(0.5, 0.5))
  }

}

##################################################################################################
##################################################################################################
##################################################################################################
