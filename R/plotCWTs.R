#######################################################################################################
# Check for detections periodicity through CWT (Continuous wavelet transform) #########################
#######################################################################################################

#' Continuous wavelet transform (CWT)
#' Decompose a signal into its constituent frequency components over time.
#'
#' @description Function that analyzes and plots periodic patterns in a time series using a
#' Continuous Wavelet Transform (CWT) framework. CWT analysis provides an alternative method
#' to the Fast Fourier Transform (FFT) or other time-frequency decomposition techniques,
#' enabling the examination of periodicities over time. This function serves mostly as a
#' wrapper for the \code{\link[wavScalogram]{cwt_wst}} function.
#'
#'
#' @param data A data frame containing binned animal detections (must contain a "timebin" column).
#' @param variable Name of the column containing the numeric variable to be analyzed.
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param period.range The range of period scales (y-axis limits) to be considered, specified in hours.
#' @param axis.periods Periods (in hours) to include/highlight in the y-axis.
#' @param color.pal Color palette. Defaults to \code{\link[pals]{jet}}
#' @param date.format Date-time format (as used in \code{\link[base]{strptime}}),
#' defining the x-axis labels. Defaults to month ("%d/%b").
#' @param date.interval Number defining the interval between each
#' displayed date (x-axis label). Defaults to 4.
#' @param date.start Integer defining the first displayed date (can be used in combination
#'  with 'date.interval" to better control the x-axis labels). Defaults to 1.
#' @param min.days Discard individuals that were detected in less than x days.
#' @param detrend Detrend time series using differences (\code{\link[base]{diff}})
#' rather than the actual values. Defaults to False.
#' @param cex.main Determines the size of the title(s). Defaults to 1.2.
#' @param cex.lab Determines the size of the y-axis and y-axis labels. Defaults to 1.1.
#' @param cex.axis Determines the size of the text labels on the axes. Defaults to 1.
#' @param legend.xpos Relative position of left and right edge of color bar on first axis (0-1).
#' Passed to \code{\link[shape]{colorlegend}} function. Defaults to c(0.90, 0.915).
#' @param legend.ypos Relative position of left and right edge of color bar on second axis (0-1).
#' Passed to \code{\link[shape]{colorlegend}} function. Defaults to c(0.15, 0.85).
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @param same.scale Forces same spectral scale (zlims) across all plots,
#' allowing for density comparison between individuals.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param ... Further arguments passed to the \code{\link[wavScalogram]{cwt_wst}} function
#' (e.g. wname, border_effects, waverad).
#' @export


plotCWTs <- function(data, variable, id.col="ID", period.range=c(3, 48), axis.periods=c(6,12,16,24,48),
                     color.pal=NULL, date.format="%d/%b", date.interval=4, date.start=1,
                     min.days=NULL, detrend=F, cex.main=1.2, cex.lab=1.1, cex.axis=1,
                     legend.xpos=c(0.90, 0.915), legend.ypos=c(0.15, 0.85), cols=2,
                     same.scale=F, id.groups=NULL, ...) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################


  # check if data contains id.col
  if(!variable %in% colnames(data)) {
    stop("'variable' variable not found in the supplied data")
  }

  if(!class(data[,variable]) %in% c("numeric", "integer")){
    stop("Please convert signal to class numeric")
  }

  # check if data contains id.col
  if(!id.col %in% colnames(data)) {
    stop("'id.col' variable not found in the supplied data")
  }

  # convert ids to factor
  if(class(data[,id.col])!="factor"){
    cat(paste("Warning: converting", id.col, "column to factor\n"))
    data[,id.col] <- as.factor(data[,id.col])
  }

  # check if data contains timebin column
  if(!c("timebin") %in% colnames(data)) {
    stop("'timebin' column not found in the supplied data")
  }

  if(length(period.range)!=2){
    stop("Please supply two values (min and max) in the period.range argument")
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {"Some of the ID(s) in id.groups don't match the IDs in the data"}
    data <- data[data[,id.col] %in% unlist(id.groups),]
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
  }



  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # create data frame with signal value
  cwt_table <- createWideTable(data, value.col=variable, id.col=id.col, agg.fun=mean)

  # convert hour periods to minutes
  period.range <- period.range*60

  # offset min period by 30 mins to avoid truncation after CWT
  period.range[1] <-  period.range[1] - 30


  # subset individuals based on minimum number of days with data
  if(!is.null(min.days)){
    data$day <- strftime(data$timebin, "%Y-%m-%d", tz="UTC")
    days_detected <- by(data$day, data[,id.col], function(x) length(unique(x)))
    selected_individuals <- which(as.numeric(days_detected) >= min.days)
    nindividuals <- length(selected_individuals)
    cat(paste(nindividuals, "individuals with >", min.days, "logged days\n"))
    cat(paste(nlevels(data[,id.col])-nindividuals, "individual(s) excluded\n"))
  }else{
    selected_individuals <- 1:nlevels(data[,id.col])
  }

  # split data by individual
  data_individual <- sapply(selected_individuals, function(i) cwt_table[,i+1], simplify=F)
  data_individual <- lapply(data_individual, function(x) x[!is.na(x)])
  names(data_individual) <- levels(data[,id.col])[selected_individuals]
  data_ts <- lapply(data_individual, ts)


  # get time bins interval (in minutes)
  interval <- difftime(data$timebin, dplyr::lag(data$timebin), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=T))


  ##############################################################################
  ## Set layout variables ######################################################
  ##############################################################################

  # set layout variables
  if(!is.null(id.groups)){
    group_ids_selected <- lapply(id.groups, function(x) x[x %in% levels(data[,id.col])[selected_individuals]])
    group_numbers <- lapply(group_ids_selected,  length)
    group_rows <- lapply(group_numbers, function(x) ceiling(x/cols))
    rows <- do.call("sum", group_rows)
    group_plots <- lapply(group_rows, function(x) x*cols)
    animal_indexes <- mapply(function(nids, nplots) {if(nids<nplots){c(1:nids, rep(NA, nplots-nids))}else{1:nids}},
                             nids=group_numbers, nplots=group_plots, SIMPLIFY=F)
    for(i in 2:length(animal_indexes)){animal_indexes[[i]]<-animal_indexes[[i]]+max(animal_indexes[[i-1]], na.rm=T)}
    animal_indexes <- unlist(animal_indexes, use.names=F)
    background_pal <- grey.colors(length(id.groups), start=0.97, end=0.93)
  } else{
    rows <- ceiling(nindividuals/cols)
    nplots <- rows*cols
    if(nindividuals<nplots){
      animal_indexes <- c(1:nindividuals, rep(NA, nplots-nindividuals))
    }else{
      animal_indexes <- 1:nindividuals
    }
    background_col <- "gray96"
  }

  plot_layout <- matrix(animal_indexes, nrow=rows, ncol=cols, byrow=T)
  bottom_plots <- apply(plot_layout, 2, max, na.rm=T)
  plot_ids <- as.integer(apply(plot_layout, 1, function(x) x))

  # set par
  par(mfrow=c(rows, cols), mar=c(4,4,4,6), oma=c(5,5,5,5))

  # set color pallete
  if(is.null(color.pal)){
    color.pal <- pals::jet(100)
    color.pal <- colorRampPalette(color.pal[5:100])(100)
  }

  ##############################################################################
  ## Calculate CWTs - Morlet wavelet spectrum ##################################
  ##############################################################################

  # compute CWTs
  cat("Calculating wavelet periodograms...\n")
  if(detrend==T) {data_ts <-  lapply(data_ts, diff)}
  cwts <- lapply(data_ts, function(x) wavScalogram::cwt_wst(x, dt=interval, scales=c(period.range, 20), powerscales=T, wname="MORLET",
                                                            border_effects="BE", makefigure=F, energy_density=T, figureperiod=T, ...))

  density_range <- lapply(cwts, function(x) range(abs(x$coefs)^2))
  density_range <- range(unlist(density_range))


  #############################################################################################
  # Plot CWTs #################################################################################

  for (i in plot_ids) {

    # if NA, add empty plot and go to next
    if(is.na(i)) {
      plot.new()
      next
    }

    # else generate CWTs
    id <- selected_individuals[i]
    cwt <- cwts[[i]]
    Z <- abs(cwt$coefs)^2

    # plot CWT spectrum
    if(same.scale==F){
      graphics::image(x=1:nrow(cwt$coefs), y=1:ncol(cwt$coefs), z=Z, col=color.pal, axes=F, useRaster=T, main="",
                      xlab="Date", ylab="Period (h)", frame.plot=T, cex.lab=cex.lab, xaxs="i", yaxs="i")
    }else{
      graphics::image(x=1:nrow(cwt$coefs), y=1:ncol(cwt$coefs), z=Z, zlim=density_range, col=color.pal, axes=F, useRaster=T, main="",
                      xlab="Date", ylab="Period (h)", frame.plot=T, cex.lab=cex.lab, xaxs="i", yaxs="i")
    }

    # add individual id
    title(main=levels(data[,id.col])[id], cex.main=cex.main, line=1)

    # plot date axis
    all_dates <- cwt_table$timebin[!is.na(cwt_table[,id+1])]
    all_dates <- strftime(all_dates, date.format)
    consec_dates <- rle(all_dates)
    consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
    unique_dates <- unique(consec_dates)
    disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
    indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
    disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
    disp_dates <- sub("\\_.*", "", disp_dates)
    axis(1, labels=disp_dates, at=disp_indexes, cex.axis=cex.axis)
    axis(1, labels=F, at=indexes, tck=-0.02, lwd.ticks=0.5)

    # plot scale axis
    period_indexes <- moby:::rescale(log2(axis.periods*60), from=log2(range(cwt$scales*cwt$fourierfactor)), to=c(1, length(cwt$scales)))
    axis(2, at=period_indexes, labels=axis.periods, cex.axis=cex.axis, las=1)

    # add guide lines
    abline(h=period_indexes, lty=5, col="grey70", lwd=0.7)

    # plot cone of influence
    x <- 1:nrow(cwt$coefs)
    coi_indexes <- moby:::rescale(cwt$coi_maxscale, from=range(cwt$scales), to=c(1, length(cwt$scales)))
    segments(x0=x[-nrow(cwt$coefs)], y0=coi_indexes[-nrow(cwt$coefs)], x1=x[-1], y1=coi_indexes[-1])
    polygon(c(x, rev(x)), c(coi_indexes, rep(par("usr")[4], length(x))), col=adjustcolor("white", alpha.f=0.5), border=F)

    # add color scale
    if(same.scale==F){
      density_labs <- pretty(c(min(Z), max(Z)), min.n=4)
      density_labs <- density_labs[density_labs>=min(Z) & density_labs<=max(Z)]
      shape::colorlegend(col=color.pal, zlim=range(Z), zval=density_labs,
                         posx=legend.xpos, posy=legend.ypos, main="", main.cex=1, digit=1, cex=0.9)
    }else{
      density_labs <- pretty(density_range)
      density_labs <- density_labs[density_labs>=min(density_range) & density_labs<=max(density_range)]
      shape::colorlegend(col=color.pal, zlim=density_range, zval=density_labs,
                         posx=legend.xpos, posy=legend.ypos, main="", main.cex=1, digit=0, cex=0.9)
    }

  }

  # add top title
  mtext(text=paste("Wavelet Power Spectrum -", tools::toTitleCase(variable)), side=3, line=1, outer=T, cex=cex.main, font=2)

  # if id.groups defined, add group legend
  if(!is.null(id.groups)){
    label_pos <- rev(unlist(lapply(group_rows, function(x) x/2)))
    for(i in 2:length(group_rows)){label_pos[i]<-label_pos[i]+rev(group_rows)[[i-1]]}
    label_pos <- grconvertY(label_pos/rows, "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels=rev(names(id.groups)),
         srt=90, cex=cex.main, font=2, xpd=NA)
  }

  #reset par
  par(mar=c(5, 4, 4, 2) + 0.1)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
