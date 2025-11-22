#######################################################################################################
# Check for signal periodicity through CWT (Continuous wavelet transform) #############################
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
#' If there are gaps between measurements (i.e., missing time bins), the handling method is
#' controlled by the \code{gap.handling} parameter.
#' @param variable Name of the column containing the numeric variable to be analyzed.
#' @param plot.title A string specifying the title of the plot. By default, this title is
#' automatically generated as "Wavelet Power Spectrum" followed by the `variable` name.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param wavelet.type A character string specifying the wavelet type to be used in the continuous wavelet transform.
#' This is passed to the `wname` argument in the \code{\link[wavScalogram]{cwt_wst}} function.
#' Possible values are: "MORLET", "DOG", "PAUL", "HAAR", or "HAAR2". The default is "MORLET".
#' @param gap.handling Method for handling missing time bins in the time series. Options are:
#' \itemize{
#'     \item "zero" (default): Fill gaps with zeros (appropriate for detection/count data)
#'     \item "mean": Fill gaps with the mean of the non-missing values for that individual
#'     \item "locf": Last observation carried forward (use the last valid observation)
#'     \item "interpolate": Linear interpolation between valid observations (recommended for environmental data)
#' }
#' @param same.scale Forces same spectral scale (zlims) across all plots,
#' allowing for density comparison between individuals.
#' @param mask.coi Logical. If TRUE, regions outside the cone of influence (COI)
#' are excluded from color scale calculations and plotted in black. This prevents
#' edge effects from dominating the color scale. Defaults to FALSE.
#' @param period.range The range of period scales (y-axis limits) to be considered,
#' specified in the units defined by \code{time.unit}. Defaults to c(3, 48) in hours.
#' @param axis.periods Periods to include/highlight on the y-axis, specified in the units
#' defined by \code{time.unit}. Defaults to c(6,12,16,24,48) in hours.
#' @param time.unit Time unit for `period.range` and `axis.periods`.
#' Options are "mins", "hours" and "days". Defaults to "hours".
#' @param color.pal Color palette. Defaults to the Jet colormap.
#' @param min.days Discard individuals that were detected in less than x days.
#' @param detrend.method Method for detrending the time series prior to CWT analysis.
#' Options are:
#' \itemize{
#'   \item "none" (default): No detrending is applied.
#'   \item "diff": Applies first-order differencing (\code{\link[base]{diff}}) to remove
#'   linear trends. Note that this analyzes the *rate of change* of the signal.
#'   \item "loess": Fits a LOESS regression (using \code{\link[stats]{loess}}) and
#'   performs the CWT on the residuals. This is effective for removing complex,
#'   non-linear trends. The smoothness is controlled by \code{loess.span}.
#' }
#' @param loess.span The span (alpha) parameter for LOESS detrending, used only if
#' \code{detrend.method = "loess"}. Controls the degree of smoothing; larger values
#' fit a smoother, less "wiggly" trend line. Defaults to 0.75.
#' @param standardize Logical. If TRUE, the time series for each individual is
#' standardized to have a mean of 0 and a standard deviation of 1 (i.e.,
#' Z-scored using \code{\link[base]{scale}}) prior to CWT analysis.
#' This is applied *after* detrending. This can help reveal periodicities
#' in signals with low amplitude relative to their mean. Defaults to FALSE.
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
#'    mar = c(3, 4, 2, 1),
#'    oma = c(1, 1, 2, 1)
#' )
#' }
#' By default, an empty list (\code{list()}) is used, which applies the function's built-in layout settings.
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
                     gap.handling = "zero",
                     same.scale = FALSE,
                     mask.coi = FALSE,
                     period.range = c(3, 48),
                     axis.periods = c(6, 12, 16, 24, 48),
                     time.unit = "hours",
                     color.pal = NULL,
                     min.days = NULL,
                     detrend.method = "none",
                     loess.span = 0.75,
                     standardize = FALSE,
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
  if(!gap.handling %in% c("zero", "mean", "locf", "interpolate")) errors <- c(errors, "Invalid gap.handling method. Use 'zero', 'mean', 'locf', or 'interpolate'.")
  if(!detrend.method %in% c("none", "diff", "loess")) errors <- c(errors, "Invalid detrend.method. Use 'none', 'diff', or 'loess'.")
  if(detrend.method == "loess" && (!is.numeric(loess.span) || loess.span <= 0)) errors <- c(errors, "`loess.span` must be a positive numeric value.")
  if(!is.logical(standardize)) errors <- c(errors, "`standardize` must be logical (TRUE or FALSE).")
  if(gap.handling %in% c("locf", "interpolate") && !requireNamespace("zoo", quietly=TRUE)) {
    errors <- c(errors, paste0("The 'zoo' package is required for gap.handling='", gap.handling, "' but is not installed. Please install 'zoo' using install.packages('zoo') and try again."))
  }
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

  # print message to console
  .printConsole("Calculating Wavelet Periodograms")

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

  # get time bins interval (in minutes) - THIS IS A ROBUST WAY
  sorted_times <- sort(unique(data[[timebin.col]]))
  if (length(sorted_times) < 2) {
    stop("Need at least two unique time bins to determine interval.", call. = FALSE)
  }
  interval <- difftime(sorted_times, dplyr::lag(sorted_times), units = "min")
  interval <- as.numeric(min(interval[interval > 0], na.rm = TRUE))

  if (!is.finite(interval) || interval == 0) {
    stop(paste("Could not determine a valid time interval (dt > 0) from 'timebin.col'."), call. = FALSE)
  }
  # print info to console
  cat(paste0("Inferred sampling interval (dt): ", interval, " minutes\n"))


  # create complete sequence of time bins
  min_time <- min(data[[timebin.col]])
  max_time <- max(data[[timebin.col]])
  all_timebins <- seq(min_time, max_time, by = paste(interval, "min"))

  # initialize wide table with timebin column
  cwt_table <- data.frame(all_timebins)
  colnames(cwt_table) <- timebin.col

  # loop over individuals
  for(id in levels(data[[id.col]])) {

    # extract individual data
    id_data <- data[data[[id.col]] == id, c(timebin.col, variable)]

    # aggregate if multiple values per timebin
    if(any(duplicated(id_data[[timebin.col]]))) {
      id_data <- aggregate(id_data[[variable]],
                           by = list(id_data[[timebin.col]]),
                           FUN = mean, na.rm = TRUE)
      colnames(id_data) <- c(timebin.col, variable)
    }

    # initialize full-length vector with NAs
    values <- rep(NA_real_, length(all_timebins))

    # fill available data
    matched_idx <- match(id_data[[timebin.col]], all_timebins)
    values[matched_idx] <- id_data[[variable]]

    # identify first and last non-NA positions
    non_na_idx <- which(!is.na(values))
    if(length(non_na_idx) > 0) {
      first_idx <- non_na_idx[1]
      last_idx  <- non_na_idx[length(non_na_idx)]

      # replace leading/trailing NAs with NaN
      if(first_idx > 1) values[1:(first_idx-1)] <- NaN
      if(last_idx < length(values)) values[(last_idx+1):length(values)] <- NaN
    } else {
      # all values NA → mark all as NaN
      values[] <- NaN
    }

    # add to wide table
    cwt_table[[id]] <- values
  }

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

  # remove any remaining NAs at the edges (before first/after last valid observation)
  data_individual <- lapply(data_individual, function(x) x[!is.nan(x)])
  names(data_individual) <- levels(data[[id.col]])[selected_individuals]

  # handle gaps according to user specification
  cat(paste0("Handling gaps using method: '", gap.handling, "'\n"))

  data_individual <- lapply(data_individual, function(x) {
    if(gap.handling == "zero") {
      # Fill NAs with zeros (original behavior for detection data)
      x[is.na(x)] <- 0
    } else if(gap.handling == "mean") {
      # Fill with mean of non-missing values
      x[is.na(x)] <- mean(x, na.rm=TRUE)
    } else if(gap.handling == "locf") {
      # Last observation carried forward
      x <- zoo::na.locf(x, na.rm=FALSE)
    } else if(gap.handling == "interpolate") {
      # Linear interpolation between valid observations
      x <- zoo::na.approx(x, na.rm=FALSE)
    }
    return(x)
  })

  # convert to time-series
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
                       mar    = c(4, 5, 4, 6),
                       oma    = c(2, 2, 3, 2),
                       mgp    = c(3, 0.8, 0))

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

  # apply detrending based on the chosen method
  if (detrend.method == "diff") {
    cat("Detrending time series using first-differencing (diff)\n")
    data_ts <- lapply(data_ts, diff)
  } else if (detrend.method == "loess") {
    cat(paste0("Detrending time series using loess (span = ", loess.span, ")\n"))
    # Iterate over names to provide better warnings
    for (id_name in names(data_ts)) {
      ts_data <- data_ts[[id_name]]
      # Loess needs a minimum number of points (e.g., > 3)
      if (length(ts_data) < 4) {
        warning(paste("Time series for", id_name, "is too short (< 4) for loess detrending, skipping."), call. = FALSE)
        next # Keep original data in data_ts
      }
      time_index <- 1:length(ts_data)
      # Use tryCatch for robustness, as loess can fail
      loess_fit <- try(stats::loess(ts_data ~ time_index, span = loess.span), silent = TRUE)
      if (inherits(loess_fit, "try-error")) {
        warning(paste("Loess detrending failed for", id_name, ", proceeding with original data."), call. = FALSE)
        next # Keep original data
      }
      data_ts[[id_name]] <- stats::residuals(loess_fit) # Replace with residuals
    }
  } else {
    cat("Proceeding without detrending\n")
  }

  # standardize the time series (Z-score) if requested
  if(standardize) {
    cat("Standardizing (Z-scoring) time series\n")
    data_ts <- lapply(data_ts, scale)
  }

  # initialize progress bar
  pb <- txtProgressBar(min=0, max=length(data_individual), initial=0, style=3)


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
    # fallback to sequential processing if cores == 1    ####################
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
  # optionally masking regions outside the cone of influence
  if(mask.coi) {
    density_range <- lapply(1:length(cwts), function(i) {
      cwt <- cwts[[i]]
      Z <- t(t(abs(cwt$coefs) ^ 2) / cwt$scales)

      # create COI mask
      coi_mask <- matrix(TRUE, nrow = nrow(Z), ncol = ncol(Z))
      for(j in 1:nrow(Z)) {
        coi_scale <- cwt$coi_maxscale[j] * cwt$fourierfactor
        # mark scales larger than COI as outside (to be masked)
        coi_mask[j, cwt$scales * cwt$fourierfactor > coi_scale] <- FALSE
      }

      # return only values inside COI
      Z[coi_mask]
    })
    density_range <- range(unlist(density_range))
  } else {
    density_range <- lapply(cwts, function(x) t(t(abs(x$coefs) ^ 2) / x$scales))
    density_range <- range(unlist(density_range))
  }



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
    if(same.scale) {
      zlim <- density_range
    } else {
      # calculate zlim from valid (inside COI) regions only if mask.coi is TRUE
      if(mask.coi) {
        # create COI mask
        coi_mask <- matrix(TRUE, nrow = nrow(Z), ncol = ncol(Z))
        for(j in 1:nrow(Z)) {
          coi_scale <- cwt$coi_maxscale[j] * cwt$fourierfactor
          # mark scales larger than COI as outside (to be masked)
          coi_mask[j, cwt$scales * cwt$fourierfactor > coi_scale] <- FALSE
        }
        zlim <- range(Z[coi_mask])
      } else {
        zlim <- range(Z)
      }
    }

    # plot CWT spectrum
    graphics::image(x=1:nrow(cwt$coefs), y=1:ncol(cwt$coefs), z=Z, zlim = zlim,
                    col=color.pal, axes=FALSE, useRaster=TRUE, main="", xlab="Date",
                    ylab = paste0("Period (", unit_abbrev, ")"),
                    frame.plot = TRUE, cex.lab = cex.lab, xaxs = "i", yaxs = "i")

    # add the individual ID as a title
    title(main=levels(data[[id.col]])[id], cex.main=cex.main, line=1)

    # add date axis
    all_dates <- cwt_table[,timebin.col][!is.nan(cwt_table[,id+1])]
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
    # 1 - plot the cone of influence (COI) to indicate regions with reduced reliability (edge effects)

    # use ALL x-points
    x_polygon <- 1:nrow(cwt$coefs)

    # get the rescaled y-values for the COI line
    coi_y_values <- .rescale(log2(cwt$coi_maxscale * cwt$fourierfactor + 1e-20),
                             from = log2(range(cwt$scales * cwt$fourierfactor)),
                             to = c(1, length(cwt$scales)))

    # get plot boundaries from 'par'
    y_bottom <- par("usr")[3] # This is 1
    y_top    <- par("usr")[4] # This is length(cwt$scales)

    # "clamp" the COI y-values to be *within* the plot's y-boundaries
    coi_y_clamped <- pmin(pmax(coi_y_values, y_bottom), y_top)

    # define the mask color based on user parameter
    mask_col <- if(mask.coi) "black" else adjustcolor("white", alpha.f = 0.5)

    # create polygon from the clamped COI line to the top of the plot
    polygon(c(x_polygon, rev(x_polygon)),
            c(coi_y_clamped, rep(y_top, length(x_polygon))),
            col = mask_col, border = TRUE)


    ############################################################################
    # 2 - plot the cone of influence (COI) to indicate regions with reduced reliability (minimum resolvable period)

    # get the minimum COI scale log-value
    min_coi_log_value <- log2(cwt$coi_minscale * cwt$fourierfactor + 1e-20)

    # rescale it to the plot's y-axis coordinates
    min_coi_y_coord <- .rescale(min_coi_log_value,
                                from = log2(range(cwt$scales * cwt$fourierfactor)),
                                to = c(1, length(cwt$scales)))

    # draw the masking polygon at the bottom
    if (length(min_coi_y_coord) == 1 && is.numeric(min_coi_y_coord) && min_coi_y_coord > y_bottom) {
      x_min <- par("usr")[1]
      x_max <- par("usr")[2]
      # draw a rectangle from the bottom (y_bottom) up to the min_coi_y_coord
      rect(x_min, y_bottom, x_max, min_coi_y_coord, col = mask_col, border = TRUE)
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
