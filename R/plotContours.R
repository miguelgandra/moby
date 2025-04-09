#######################################################################################################
# Function to generate contour plots  #################################################################
#######################################################################################################

#' Generate contour plot
#'
#' @description Function to generate a contour plot (color-coded 2-dimensional plots)
#' of a given variable, by customizable time and date intervals, either for all individuals
#' combined or independently for each ID. It also plots lines illustrating the
#' variation in diel phases' hours across the study duration (annual variation of daylight time).
#'
#' @inheritParams setDefaults
#' @param data A data frame containing the variable(s) to plot.
#' @param variables Variables to plot. Should correspond columns in the data.
#' Can be either a vector or a single variable.
#' @param var.titles Display name of the variables (e.g. "Depth (m)"). If NULL, defaults to supplied variable names.
#' @param plot.title Optional. Main title displayed on top of the plots/panel.
#' @param datetime.col Name of the column containing datetimes in POSICX format. Defaults to 'datetime'.
#' @param split.by Optional. If defined, plots are generated individually for each level (or combination of levels)
#' of this variable(s) (e.g. species, ontogeny, sex, habitat, etc).
#' @param aggregate.fun Function used to aggregate values by hour and month. Defaults to mean.
#' @param color.pal Color palette for the contour plot. If NULL, defaults to the viridis color palette.
#' @param time.interval Time interval for binning data along the y-axis. Can be specified as a
#' lubridate-style period string (e.g., "30 mins", "1 hour", "2 hours"). Defaults to "hour".
#' @param date.interval Date interval for binning data along the x-axis. Can be specified as:
#' \itemize{
#'   \item A lubridate-style period string (e.g., "1 day", "15 days", "1 week", "1 month")
#'   \item "month/year" or "month-year" for month-year combinations (e.g., "Jan-2023", "Feb-2023")
#' }
#' Defaults to "month".
#' @param diel.lines Number indicating the number of diel phase lines (boundaries)
#' to display. Either 0 (no lines), 2 (corresponding to sunrise and sunset) or 4
#' (depicting dawn, sunrise, sunset, dusk). Defaults to 4.
#' @param diel.lines.col Color of the diel lines. Defaults to black.
#' @param sunriset.coords A SpatialPoints or matrix object containing longitude and
#' latitude coordinates (in that order) at which to estimate sunrise and sunset times
#' (if `diel.lines` > 0).
#' @param solar.depth Angle of the sun below the horizon (in degrees). Passed to the
#' `solarDep` argument in \code{\link[suntools]{crepuscule}} function. Defaults to 18 (astronomical twilight).
#' @param cex.main Determines the size of the title(s). Defaults to 1.1.
#' @param cex.lab Determines the size of the y-axis and y-axis labels. Defaults to 1.
#' @param cex.axis Determines the size of the text labels on the axes. Defaults to 0.8.
#' @param cex.legend Determines the size of the text labels on the color legend. Defaults to 0.7.
#' @param grid Add horizontal and vertical guiding lines? Defaults to True.
#' @param zlim Optional. A numeric vector of length 2 specifying the range for the color scale (e.g., `c(-10, 10)`).
#' If provided, this overrides the default range calculated from the data.
#' Useful, for example, for ensuring symmetric color scales centered on zero. Defaults to NULL.
#' @param invert.scale Invert color scale legend? Defaults to False.
#' @param uniformize.scale Use the same color scale for all plots when 'split.by' is defined.
#' @param legend.xpos Relative position of left and right edge of color bar on first axis (0-1).
#' Defaults to c(0.89, 0.92).
#' @param legend.ypos Relative position of left and right edge of color bar on second axis (0-1).
#' Defaults to c(0.15, 0.85).
#' @param tz A character string specifying the time zone to be used for the conversion. Defaults to "UTC".
#' @param cols Number of columns in the final panel (passed to the mfrow argument). Only used when
#' by.id is set to true or the split.by argument is defined.
#' @param disable.par Boolean. Disables internal par settings, so that multi-figure/panel layouts can be
#' generated if required. If set to True, 'cols' argument will have no effect in the output. Defaults to False.
#' @param ... Further arguments passed to the \code{\link{filled.contour}} function
#' @export


plotContours <- function(data,
                         variables,
                         var.titles = NULL,
                         plot.title = NULL,
                         split.by = NULL,
                         aggregate.fun = function(x) mean(x, na.rm=TRUE),
                         id.col = getDefaults("ID"),
                         datetime.col = getDefaults("datetime"),
                         color.pal = NULL,
                         time.interval = "hour",
                         date.interval = "month",
                         diel.lines = 4,
                         diel.lines.col = "black",
                         sunriset.coords = NULL,
                         solar.depth = 18,
                         cex.main = 1.1,
                         cex.lab = 1,
                         cex.axis = 0.8,
                         cex.legend = 0.7,
                         grid = TRUE,
                         zlim = NULL,
                         invert.scale = FALSE,
                         uniformize.scale = FALSE,
                         legend.xpos = c(0.89, 0.92),
                         legend.ypos = c(0.15, 0.85),
                         tz = "UTC",
                         cols = 1,
                         disable.par = FALSE,
                         ...) {


  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # validate additional parameters
  errors <- c()
  if(!is.null(var.titles) & length(variables)!=length(var.titles)) errors <- c(errors, "Number of variables and variables' titles do not match.")
  if(diel.lines>0 && is.null(sunriset.coords))  errors <- c(errors, "Sunrise/sunset coordinates must be provided if 'diel.lines' is greater than 0.")

  # validate time.interval and date.interval
  time_interval_valid <- tryCatch({
    lubridate::as.period(lubridate::period(time.interval))
    TRUE
  }, error = function(e) FALSE)
  date_interval_valid <- tryCatch({
    lubridate::as.period(lubridate::period(date.interval))
    TRUE
  }, error = function(e) FALSE)
  if(!time_interval_valid) errors <- c(errors, "time.interval must be a valid lubridate period string (e.g., '30 mins', '1 hour')")
  if(!date_interval_valid) errors <- c(errors, "date.interval must be a valid lubridate period string (e.g., '1 day', '15 days', '1 month')")

  # return errors (if any)
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }


  # set the color palette if none provided
  if(is.null(color.pal)) color.pal <- .viridis_pal
  if(!inherits(color.pal, "function")) color.pal <- colorRampPalette(color.pal)

  # print to console
  .printConsole("Generating contour plot(s)")


  #####################################################################################
  # Process grouping variables ########################################################
  #####################################################################################

  # process grouping variable(s) if provided
  if(!is.null(split.by)){
    # ensure grouping columns are factors for consistent plotting
    for(s in split.by){
      if(!inherits(data[,s], "factor")){
        warning("Converting split.by column to factor", call.=FALSE)
        data[,s] <- as.factor(data[,s])
      }
    }
    # if multiple grouping variables, create a single interaction column to track combinations
    if(length(split.by)>1){
      data$dummy_group <- apply(data[split.by], 1, function(x) paste(x, collapse=" "))
      data$dummy_group <- as.factor(data$dummy_group)
      split.by <- "dummy_group"
    }

  # if no grouping variable, assign all data to a single artificial group
  }else{
    data$dummy_group <- "1"
    data$dummy_group <- as.factor(data$dummy_group)
    split.by <- "dummy_group"
  }


  #####################################################################################
  # Aggregate data and calculate stats ################################################
  #####################################################################################

  # define start and end of the day
  start_time <- as.POSIXct("00:00", format = "%H:%M", tz = tz)
  end_time <- as.POSIXct("23:59", format = "%H:%M", tz = tz)
  # generate sequence of time breaks at specified interval
  time_breaks <- seq(start_time, end_time, by = time.interval)
  time_breaks <- format(time_breaks, "%H:%M")
  # assign time bins using lubridate
  data$time_bin <- lubridate::floor_date(data[[datetime.col]], unit = time.interval)
  data$time_bin <- format(data$time_bin, "%H:%M")
  data$time_bin <- factor(data$time_bin, levels=time_breaks)

  # generate custom date intervals
  if (grepl("month[/-]year", date.interval, ignore.case = TRUE)) {
    # for month-year grouping
    data$date_bin <- format(data[[datetime.col]], "%b-%Y")
    # create all possible month-year levels for proper ordering
    date_range <- range(data[[datetime.col]], na.rm = TRUE)
    all_months <- seq(floor_date(date_range[1], "month"), floor_date(date_range[2], "month"),  by = "1 month")
    date_breaks <- format(all_months, "%b-%Y")
    data$date_bin <- factor(data$date_bin, levels = date_breaks)
  } else {
    # original date binning logic for other interval types
    start_date <- as.Date("2023-01-01")
    end_date <- as.Date("2023-12-31")
    date_breaks <- seq.Date(from = start_date, to = end_date, by = date.interval)
    date_breaks <- unique(lubridate::floor_date(date_breaks, unit = date.interval))
    date_breaks <- format(date_breaks, "%d/%b")
    data$date_bin <- lubridate::floor_date(data[[datetime.col]], unit = date.interval)
    data$date_bin <- format(data$date_bin, "%d/%b")
    data$date_bin <- factor(data$date_bin, levels=date_breaks)
  }

  # aggregate data by ID, date interval, time interval, and group
  aggregated_data <- stats::aggregate(data[, variables],  by = list(data[, id.col], data$date_bin, data$time_bin, data[, split.by]),
                                      aggregate.fun, simplify = TRUE, drop = FALSE)
  colnames(aggregated_data)[1:4] <- c("id", "date", "time", "group")
  colnames(aggregated_data)[5:ncol(aggregated_data)] <- variables

  # calculate number of individuals with data per variable and group
  nids <- list()
  for(v in 1:length(variables)){
    var <- variables[v]
    var_data <- aggregated_data[,c("id","group",var)]
    var_data <- var_data[!is.na(var_data[,var]),]
    nids[[v]] <- stats::aggregate(var_data$id, by=list(var_data$group), function(x) length(unique(x)))
    colnames(nids[[v]]) <- c("group", "nids")
    nids[[v]]$var <- var
  }
  nids <- do.call("rbind", nids)

  # calculate monthly and hourly averages for each group
  aggregated_data <- stats::aggregate(aggregated_data[,variables], by=list(aggregated_data$group, aggregated_data$date, aggregated_data$time),
                                      mean, na.rm=TRUE, drop = FALSE)
  colnames(aggregated_data)[1:3] <- c("group", "date", "time")
  colnames(aggregated_data)[4:ncol(aggregated_data)] <- variables
  aggregated_data <- aggregated_data[order(aggregated_data$group, aggregated_data$date, aggregated_data$time),]
  aggregated_data <- split(aggregated_data, f=aggregated_data$group, drop=TRUE)


  # join all data (to uniformize scales, if required)
  complete_data <- do.call("rbind", aggregated_data)


  #####################################################################################
  # Prepare plot variables ############################################################
  #####################################################################################

  # calculate sunrise and sunset times for the specified coordinates and date range
  if(diel.lines>0){
    sunriset_start <- min(data[,datetime.col])
    sunriset_end <- max(data[,datetime.col])
    # determine the appropriate time interval for sunrise/sunset calculations
    date_format <- ifelse(grepl("month[/-]year", date.interval, ignore.case = TRUE), "%b-%Y", "%d/%b")
    # calculate diel phase boundary times
    daytimes_table <- getSunTimes(sunriset.coords, sunriset_start, sunriset_end, solar.depth, by=date_format)
    # convert interval to numeric position matching the plot's x-axis
    daytimes_table$interval <- match(daytimes_table$interval, date_breaks)
    daytimes_table <- daytimes_table[!is.na(daytimes_table$interval),]
    daytimes_table <- daytimes_table[order(daytimes_table$interval),]
    # convert sunrise/sunset times to match time.interval bins
    time_bin_size <- as.numeric(lubridate::duration(time.interval)) / 3600
    # Convert decimal hours to binned time positions
    time_cols <- c("dawns", "sunrises", "sunsets", "dusks")
    daytimes_table[time_cols] <- lapply(daytimes_table[time_cols], function(x) {x / time_bin_size})
}

  # precompute x-axis parameters based on date interval type
  axis_params <- if (grepl("month[/-]year", date.interval, ignore.case = TRUE)) {
    list(at = seq_along(date_breaks), labels = date_breaks, xlab = "Month-Year")
  } else if (date.interval %in% c("1 month", "month")) {
    list(at = seq_along(date_breaks), labels = substr(date_breaks, 4, 6), xlab = "Month")
  } else {
    date_breaks_dates <- as.Date(date_breaks, format = "%d/%b")
    month_numbers <- as.numeric(format(date_breaks_dates, "%m"))
    month_centers <- sapply(1:12, function(m) which(month_numbers == m)[1])
    list(at = month_centers, labels = month.abb, xlab = "Month")
  }

  # precompute y-axis parameters based on time interval type
  time_hours <- as.numeric(substr(time_breaks, 1, 2)) + as.numeric(substr(time_breaks, 4, 5))/60
  even_hours <- seq(0, 22, by = 2)
  display_hour_pos <- time_hours[sapply(even_hours, function(h) {which.min(abs(time_hours - h))})]
  hour_labels <- sprintf("%02dh", even_hours)
  hour_pos <- time_hours[sapply(0:23, function(h) which.min(abs(time_hours - h)))]

  # set up the graphical layout for multi-panel plots if necessary
  if(!disable.par){
    nplots <- length(aggregated_data)*length(variables)
    rows <- ceiling(nplots/cols)
    par(mfrow=c(rows, cols), mgp=c(2.2,0.6,3), mar=c(4,4,4,6))
    if(!is.null(plot.title)){
      par(oma=c(0,0,2.5,0))
    }else{
      par(oma=c(0,0,0,0))
    }
  }


  #####################################################################################
  # Plot contours #####################################################################
  #####################################################################################

  # iterate through each group in the aggregated data
  for(i in 1:length(aggregated_data)){

    # iterate through each specified variable
    for(v in 1:length(variables)){

      # generate a 2D matrix for the current variable across hours and months
      var <- variables[v]
      contour_matrix <- matrix(aggregated_data[[i]][,var], nrow=length(time_breaks), ncol=length(date_breaks))
      colnames(contour_matrix) <- date_breaks
      rownames(contour_matrix) <- time_breaks
      contour_matrix <- rbind(contour_matrix, contour_matrix[1,])

      # set variable scale
      if(uniformize.scale){
        var_min <- min(complete_data[,var], na.rm=TRUE)
        var_max <- max(complete_data[,var], na.rm=TRUE)
        scale <- c(min(var_min), max(var_max))
      } else {
        scale <- range(contour_matrix, na.rm=TRUE)
      }

      # override scale if zlim is provided
      if(!is.null(zlim)){
        scale <- zlim
      }

      # define the plot title based on group name and variable title if available
      if(length(aggregated_data)>1){
        if(!is.null(var.titles)){plot_title <- paste(names(aggregated_data)[i], "-", var.titles[v])
        }else{plot_title <- paste(names(aggregated_data)[i], "-", var)}
      }else{
        if(!is.null(var.titles)){plot_title <- var.titles[v]
        }else{plot_title <- var}
      }

      # add a subtitle for sample size in the plot
      plot_sub <- nids$nids[nids$group==names(aggregated_data)[i] & nids$var==var]
      plot_sub <- paste0("(n=", plot_sub, ")")


      ##########################################################################
      # generate the contour plot ##############################################
      .filled.contour(x=1:length(date_breaks), y=0:length(time_breaks), z=t(contour_matrix), main=plot_title, cex.main=cex.main, invert.scale=invert.scale,
                      xlab="Month", ylab="Hour", nlevels=100, color.palette=color.pal, cex.lab=cex.lab, zlim=scale,
                      plot.axes = {
                        # configure x-axis - month labels
                        axis(1, at = axis_params$at, labels = axis_params$labels, pos = 0, cex.axis = cex.axis)
                        # configure y-axis - hour labels
                        axis(2, at=display_hour_pos, labels=hour_labels, pos=1, cex.axis=cex.axis, las=1)
                        axis(2, at = hour_pos,  labels=FALSE, tck=-0.015, lwd.ticks=0.5, pos=1)
                        # plot lines for different daylight periods (dawn, sunrise, sunset, dusk)
                        if(diel.lines>0){
                          lines(daytimes_table$interval, daytimes_table$sunrises, lty=2, col=diel.lines.col)
                          lines(daytimes_table$interval, daytimes_table$sunsets, lty=2, col=diel.lines.col)
                        }
                        if(diel.lines==4){
                          lines(daytimes_table$interval, daytimes_table$dawns, lty=2, col=diel.lines.col)
                          lines(daytimes_table$interval, daytimes_table$dusks, lty=2, col=diel.lines.col)
                        }
                        # add grid lines if enabled
                        if(grid){
                          abline(v=axis_params$at, lwd=0.03)
                          abline(h=hour_pos, lwd=0.03)
                        }}, ...)
      # add the sample size as a title below the main title
      legend("topright", legend=plot_sub, bty="n", cex=cex.legend)
      # set up the color scale legend
      scale_labs <- pretty(scale, min.n=4)
      scale_labs <- scale_labs[scale_labs>=min(scale) & scale_labs<=max(scale)]
      color_scale <- color.pal(100)
      if(invert.scale){scale<-rev(scale); color_scale<-rev(color_scale)}
      digits <- max(.decimalPlaces(scale_labs))
      .colorlegend(col=color_scale, zlim=scale, zval=scale_labs, digit=digits, xpd=TRUE,
                         posx=legend.xpos, posy=legend.ypos, main="", main.cex=cex.axis, cex=cex.legend)
    }
  }

  # add an overall title for the entire plot, if specified
  if(!is.null(plot.title)){
    mtext(text=plot.title, side=3, cex=1.25, font=2, line=-0.5, outer=TRUE)
  }

}

#######################################################################################################
#######################################################################################################
#######################################################################################################
