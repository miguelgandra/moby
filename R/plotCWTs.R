#######################################################################################################
# Check for detections periodicity through CWT (Continuous wavelet transform) #########################
#######################################################################################################

#' Continuous wavelet transform (CWT)
#'
#' @description Function that analyzes and plots periodic patterns in a time series using a
#' Continuous Wavelet Transform (CWT) framework. CWT analysis provides an alternative method
#' to the Fast Fourier Transform (FFT) or other time-frequency decomposition techniques,
#' enabling the examination of periodicities over time. This function serves mostly as a
#' wrapper for the \code{\link[wavScalogram]{cwt_wst}} function from the `wavScalogram` package
#' (Bolos & Benitez, 2022).
#'
#' @inheritParams setDefaults
#' @param data A data frame containing binned animal detections or other time-based measurements.
#' The time series does not need to be regular; gaps between time bins are allowed.
#' If there are gaps between measurements (i.e., missing time bins), the function will
#' automatically assume a value of zero for those missing time steps.
#' @param variable Name of the column containing the numeric variable to be analyzed.
#' @param plot.title A string specifying the title of the plot. By default, this title is
#' automatically generated as "Wavelet Power Spectrum" followed by the `variable` name.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param wavelet.type A character string specifying the wavelet type to be used in the continuous wavelet transform.
#' This is passed to the `wname` argument in the \code{\link[wavScalogram]{cwt_wst}} function.
#' Possible values are: "MORLET", "DOG", "PAUL", "HAAR", or "HAAR2". The default is "MORLET".
#' @param same.scale Forces same spectral scale (zlims) across all plots,
#' allowing for density comparison between individuals.
#' @param period.range The range of period scales (y-axis limits) to be considered,
#' specified in the units defined by \code{time.unit}. Defaults to c(3, 48) in hours.
#' @param axis.periods Periods to include/highlight on the y-axis, specified in the units
#' defined by \code{time.unit}. Defaults to c(6,12,16,24,48) in hours.
#' @param time.unit Time unit for `period.range` and `axis.periods`.
#' Options are "mins", "hours" and "days". Defaults to "hours".
#' @param color.pal Color palette. Defaults to the Jet colormap.
#' @param min.days Discard individuals that were detected in less than x days.
#' @param detrend Detrend time series using differences (\code{\link[base]{diff}})
#' rather than the actual values. Defaults to False.
#' @param date.format Date-time format (as used in \code{\link[base]{strptime}}),
#' defining the x-axis labels. Defaults to month ("%d/%b").
#' @param date.interval Number defining the interval between each
#' displayed date (x-axis label). Defaults to 4.
#' @param date.start Integer defining the first displayed date (can be used in combination
#'  with 'date.interval" to better control the x-axis labels). Defaults to 1.
#' @param cex.main Determines the size of the title(s). Defaults to 1.2.
#' @param cex.lab Determines the size of the y-axis and y-axis labels. Defaults to 1.1.
#' @param cex.axis Determines the size of the text labels on the axes. Defaults to 1.
#' @param cex.legend Determines the size of the text labels on the color legend. Defaults to 0.9.
#' @param legend.xpos Relative position of left and right edge of color bar on first axis (0-1).
#' Defaults to c(0.90, 0.915).
#' @param legend.ypos Relative position of left and right edge of color bar on second axis (0-1).
#' Defaults to c(0.15, 0.85).
#' @param cols Number of columns to arrange plots in when multiple individuals are displayed.
#' Ignored if \code{par.args$mfrow} is provided, in which case the layout is controlled
#' entirely by the user via \code{par.args}.
#' @param par.args Optional named list of graphical parameters to pass to [graphics::par()].
#' This allows fine control over the multi-panel layout and plot spacing.
#' Any arguments supplied here will override the internal defaults used by the function
#' (such as \code{mfrow}, \code{mar}, \code{oma}, or \code{mgp}).
#' For example, to reduce margins and adjust outer spacing:
#' \preformatted{
#' par.args = list(
#'   mar = c(3, 4, 2, 1),
#'   oma = c(1, 1, 2, 1)
#' )
#' }
#' By default, an empty list (\code{list()}) is used, which applies the function’s built-in layout settings.
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations.
#' Run \code{parallel::detectCores()} to check the number of available cores.
#' @param ... Further arguments passed to the \code{\link[wavScalogram]{cwt_wst}} function.
#' (e.g. border_effects, waverad).
#'
#' @references
#' Bolos, V. J., & Benitez, R. (2022). wavScalogram: an R package with scalogram wavelet tools for time series analysis. The R Journal, 14(2), 164-185.
#'
#' @export


plotCWTs <- function(data,
                     variable,
                     plot.title = paste("Wavelet Power Spectrum -", tools::toTitleCase(variable)),
                     id.col = getDefaults("ID"),
                     timebin.col = getDefaults("timebin"),
                     id.groups = NULL,
                     wavelet.type = "MORLET",
                     same.scale = FALSE,
                     period.range = c(3, 48),
                     axis.periods = c(6, 12, 16, 24, 48),
                     time.unit = "hours",
                     color.pal = NULL,
                     min.days = NULL,
                     detrend = FALSE,
                     date.format = "%d/%b",
                     date.interval = 4,
                     date.start = 1,
                     cex.main = 1.2,
                     cex.lab = 1.1,
                     cex.axis = 1,
                     cex.legend = 0.9,
                     legend.xpos = c(0.895, 0.91),
                     legend.ypos = c(0.15, 0.85),
                     cols = 1,
                     par.args = list(),
                     cores = 1,
                     ...) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # validate additional parameters
  errors <- c()
  if(!requireNamespace("wavScalogram", quietly=TRUE)) errors <- c(errors, "The 'wavScalogram' package is required for this function but is not installed. Please install 'wavScalogram' using install.packages('wavScalogram') and try again.")
  if(!class(data[[variable]]) %in% c("numeric", "integer")) errors <- c(errors, "Please convert signal to class numeric")
  if (!is.numeric(period.range) || length(period.range) != 2) errors <- c(errors, "`period.range` must be a numeric vector of length 2 (min and max).")
  if(!time.unit %in% c("mins", "hours", "days")) errors <- c(errors, "Invalid time unit specified. Use 'mins', 'hours', or 'days'.")
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # save the current par settings and ensure they are restored upon function exit
  original_par <- par(no.readonly=TRUE)
  original_par$pin <- NULL
  original_par$cin <- NULL
  original_par$cra <- NULL
  original_par$csi <- NULL
  original_par$cxy <- NULL
  original_par$din <- NULL
  original_par$page <- NULL
  on.exit(par(original_par), add = TRUE)

  # drop missing ID levels
  data[[id.col]] <- droplevels( data[[id.col]])

  # convert time.unit periods
  if (time.unit=="mins") {
    period.range <- period.range
    time_factor <- 1
    unit_abbrev <- "mins"
  } else if (time.unit=="hours") {
    period.range <- period.range * 60
    time_factor <- 60
    unit_abbrev <- "h"
  } else if (time.unit=="days") {
    period.range <- period.range * 1440
    time_factor <- 1440
    unit_abbrev <- "days"
  }


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # create data frame with signal value
  cwt_table <- createWideTable(data, id.col=id.col, timebin.col=timebin.col,
                               value.col=variable,  agg.fun=mean)

  # offset min period by 30 mins to avoid truncation after CWT
  period.range[1] <-  period.range[1] - 30

  # subset individuals based on minimum number of days with data
  if(!is.null(min.days)){
    data$day <- strftime(data[[timebin.col]], "%Y-%m-%d", tz="UTC")
    days_detected <- by(data$day, data[[id.col]], function(x) length(unique(x)))
    selected_individuals <- which(as.numeric(days_detected) >= min.days)
    nindividuals <- length(selected_individuals)
    cat(paste(nindividuals, "individuals with >", min.days, "logged days\n"))
    cat(paste(nlevels(data[[id.col]])-nindividuals, "individual(s) excluded\n"))
  }else{
    selected_individuals <- 1:nlevels(data[[id.col]])
    nindividuals <- length(selected_individuals)
  }

  # split data by individual
  data_individual <- lapply(selected_individuals, function(i) cwt_table[,i+1])
  data_individual <- lapply(data_individual, function(x) x[!is.na(x)])
  names(data_individual) <- levels(data[[id.col]])[selected_individuals]
  data_ts <- lapply(data_individual, ts)


  # get time bins interval (in minutes)
  interval <- difftime(data[[timebin.col]], dplyr::lag(data[[timebin.col]]), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=TRUE))


  ##############################################################################
  ## Set layout variables ######################################################
  ##############################################################################

  # set layout variables
  if(!is.null(id.groups)){
    group_ids_selected <- lapply(id.groups, function(x) x[x %in% levels(data[[id.col]])[selected_individuals]])
    group_numbers <- lapply(group_ids_selected,  length)
    group_rows <- lapply(group_numbers, function(x) ceiling(x/cols))
    rows <- do.call("sum", group_rows)
    group_plots <- lapply(group_rows, function(x) x*cols)
    animal_indexes <- mapply(function(nids, nplots) {if(nids<nplots){c(1:nids, rep(NA, nplots-nids))}else{1:nids}},
                             nids=group_numbers, nplots=group_plots, SIMPLIFY=FALSE)
    for(i in 2:length(animal_indexes)){animal_indexes[[i]]<-animal_indexes[[i]]+max(animal_indexes[[i-1]], na.rm=TRUE)}
    animal_indexes <- unlist(animal_indexes, use.names=FALSE)
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

  plot_layout <- matrix(animal_indexes, nrow=rows, ncol=cols, byrow=TRUE)
  bottom_plots <- apply(plot_layout, 2, max, na.rm=TRUE)
  plot_ids <- as.integer(apply(plot_layout, 1, function(x) x))
  plot_ids <- plot_ids[!is.na(plot_ids)]

  # set par layout with user-overridable graphical parameters
  par.defaults <- list(mfrow = c(rows, cols),
                       mar   = c(4, 5, 4, 6),
                       oma   = c(2, 2, 3, 2),
                       mgp   = c(3, 0.8, 0))

  # allow user to override any of these
  par.defaults[names(par.args)] <- par.args
  # apply combined settings
  par(par.defaults)

  # set color pallete
  if(is.null(color.pal)){
    color.pal <- .jet_pal(100)
    color.pal <- colorRampPalette(color.pal[5:100])(100)
  }

  ##############################################################################
  ## Calculate CWTs - Morlet wavelet spectrum ##################################
  ##############################################################################

  # print message to console
  cat("Calculating wavelet periodograms...\n")

  # if detrending is requested, apply a difference operation to the time series
  if(detrend) data_ts <- lapply(data_ts, diff)

  # initialize progress bar
  pb <- txtProgressBar(min=1, max=length(data_individual), initial=0, style=3)


  ########################################################################
  # use parallel computing ###############################################
  if(cores>1) {

    # print information to console
    cat(paste0("Starting parallel computation: ", cores, " cores\n"))

    # register parallel backend with the specified number of cores
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)

    # ensure the cluster is properly stopped when the function exits
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # define the `%dopar%` operator locally for parallel execution
    `%dopar%` <- foreach::`%dopar%`

    # set progress bar option
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))

    # perform parallel computation over each individual's data using foreach
    cwts <- foreach::foreach(i=1:length(data_ts), .options.snow=opts, .packages="wavScalogram") %dopar% {
      wavScalogram::cwt_wst(data_ts[[i]], dt=interval, scales=c(period.range, 20), powerscales=TRUE, wname=wavelet.type,
                            border_effects="BE", makefigure=FALSE, energy_density=TRUE, figureperiod=TRUE, ...)
    }

  ########################################################################
  # fallback to sequential processing if cores == 1   ####################
  } else {

    # initialize list to store the results
    cwts <- vector("list", length(data_ts))

    # loop through each individual
    for(i in 1:length(data_ts)){
      # calculate wavelet periodograms
      cwts[[i]] <- wavScalogram::cwt_wst(data_ts[[i]], dt=interval, scales=c(period.range, 20), powerscales=TRUE, wname=wavelet.type,
                                         border_effects="BE", makefigure=FALSE, energy_density=TRUE, figureperiod=TRUE, ...)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
  }

  # close the progress bar
  close(pb)

  # calculate the range of the absolute squared coefficients across all individuals
  density_range <- lapply(cwts, function(x) t(t(abs(x$coefs) ^ 2) / x$scales))
  density_range <- range(unlist(density_range))



  ##############################################################################
  # Plot CWTs ##################################################################
  ##############################################################################

  # loop through each selected individual to plot the CWTs
  for (i in plot_ids) {

    # if the ID is NA, skip this iteration and create an empty plot
    if(is.na(i)) {
      plot.new()
      next
    }

    # otherwise, generate the CWT plot for the current individual
    id <- selected_individuals[i]
    cwt <- cwts[[i]]

    # compute the wavelet power spectrum (Z), scaling the coefficients by the scales (energy density)
    Z <- t(t(abs(cwt$coefs) ^ 2) / cwt$scales)

    # set density limits for the plot
    if(same.scale) zlim <- density_range
    else zlim <- range(Z)

    # plot CWT spectrum
    graphics::image(x=1:nrow(cwt$coefs), y=1:ncol(cwt$coefs), z=Z, zlim = zlim,
                    col=color.pal, axes=FALSE, useRaster=TRUE, main="", xlab="Date",
                    ylab = paste0("Period (", unit_abbrev, ")"),
                    frame.plot = TRUE, cex.lab = cex.lab, xaxs = "i", yaxs = "i")

    # add the individual ID as a title
    title(main=levels(data[[id.col]])[id], cex.main=cex.main, line=1)

    # add date axis
    all_dates <- cwt_table[,timebin.col][!is.na(cwt_table[,id+1])]
    all_dates <- strftime(all_dates, date.format)
    consec_dates <- rle(all_dates)
    consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
    unique_dates <- unique(consec_dates)
    disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
    indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
    disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
    disp_dates <- sub("\\_.*", "", disp_dates)
    axis(1, labels=disp_dates, at=disp_indexes, cex.axis=cex.axis)
    axis(1, labels=FALSE, at=indexes, tck=-0.02, lwd.ticks=0.5)

    # add scale (period) axis
    period_indexes <- .rescale(log2(axis.periods*time_factor), from=log2(range(cwt$scales*cwt$fourierfactor)), to=c(1, length(cwt$scales)))
    axis(2, at=period_indexes, labels=axis.periods, cex.axis=cex.axis, las=1)

    # add guide lines for each period
    abline(h=period_indexes, lty=5, col="grey70", lwd=0.7)

    ############################################################################
    # plot the cone of influence (COI) to indicate regions with reduced reliability
    x <- 1:nrow(cwt$coefs)
    coi_values <- .rescale(log2(cwt$coi_maxscale*cwt$fourierfactor+1e-20),
                           from=log2(range(cwt$scales*cwt$fourierfactor)),
                           to=c(1, length(cwt$scales)))
    # find valid COI indices (within plot bounds)
    coi_indices <- which(coi_values>=1 & coi_values<=par("usr")[4])
    coi_values_filtered <- coi_values[coi_indices]
    x_filtered <- x[coi_indices]

    # extend polygon to plot boundaries
    if(length(x_filtered) > 0) {
      # add left boundary point
      x_polygon <- c(par("usr")[1], x_filtered, par("usr")[2])
      # extend COI values to boundaries (use edge values)
      coi_polygon <- c(coi_values_filtered[1], coi_values_filtered, coi_values_filtered[length(coi_values_filtered)])
      # create polygon from COI line to top of plot
      polygon(c(x_polygon, rev(x_polygon)),
              c(coi_polygon, rep(par("usr")[4], length(x_polygon))),
              col=adjustcolor("white", alpha.f=0.5),
              border=TRUE)
    }

    ############################################################################
    # add a color legend to the plot (only once if same.scale = TRUE)
    if (!same.scale) {
      # current behaviour: draw one legend per plot
      density_labs <- pretty(zlim)
      density_labs <- density_labs[density_labs >= min(zlim) & density_labs <= max(zlim)]
      digits <- max(.decimalPlaces(density_labs))
      .colorlegend(col=color.pal, zlim=zlim, zval=density_labs, posx=legend.xpos,
                   posy=legend.ypos, main="", main.cex=1, digit=digits, cex=cex.legend, xpd = NA)
    } else if (same.scale && i == max(plot_ids)) {
      # only draw the legend once (after the last plot)
      # current behaviour: draw one legend per plot
      density_labs <- pretty(zlim)
      density_labs <- density_labs[density_labs >= min(zlim) & density_labs <= max(zlim)]
      digits <- max(.decimalPlaces(density_labs))
      .colorlegend(col=color.pal, zlim=zlim, zval=density_labs, posx=legend.xpos,
                   posy=legend.ypos, main="", main.cex=1, digit=digits, cex=cex.legend, xpd = NA)
    }
  }

  # add a top title for the whole plot
  mtext(text=plot.title, side=3, line=1.6, outer=T, cex=cex.main, font=2)

  # if groupings are defined, add a legend for the groups
  if(!is.null(id.groups)){
    label_pos <- rev(unlist(lapply(group_rows, function(x) x/2)))
    for(i in 2:length(group_rows)){label_pos[i]<-label_pos[i]+rev(group_rows)[[i-1]]}
    label_pos <- grconvertY(label_pos/rows, "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels=rev(names(id.groups)),
         srt=90, cex=cex.main, font=2, xpd=NA)
  }

}


#######################################################################################################
#######################################################################################################
#######################################################################################################
