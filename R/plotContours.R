#######################################################################################################
# Function to generate contour plots  #################################################################
#######################################################################################################

#' Generate contour plot
#'
#' @description Function to generate a contour plot (color-coded 2-dimensional plots)
#' of a given variable, by hour and month, either for all individuals combined (average) or
#' independently for each ID. It also plots lines illustrating the variation in
#' diel phases' hours across the study duration (annual variation of daylight time).
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
#' @param grid Add horizontal and vertical guiding lines? Defaults to True.
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
                         id.col = getDefaults("id"),
                         datetime.col = getDefaults("datetime"),
                         color.pal = NULL,
                         diel.lines = 4,
                         diel.lines.col = "black",
                         sunriset.coords = NULL,
                         solar.depth = 18,
                         cex.main = 1.1,
                         cex.lab = 1,
                         cex.axis = 0.8,
                         grid = TRUE,
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

  # extract month and hour from datetime column
  data$month_tmp <- strftime(data[,datetime.col], "%m", tz=tz)
  data$hour_tmp <- strftime(data[,datetime.col], "%H", tz=tz)
  all_months <- sprintf("%02d", 1:12)

  # aggregate data by ID, month, hour, and group
  aggregated_data <- stats::aggregate(data[,variables], by=list(data[,id.col], data$month_tmp, data$hour_tmp, data[,split.by]), aggregate.fun, simplify=TRUE, drop=TRUE)
  colnames(aggregated_data)[1:4] <- c("id", "month", "hour", "group")
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
  aggregated_data <- stats::aggregate(aggregated_data[,variables], by=list(aggregated_data$group, aggregated_data$month, aggregated_data$hour), mean, na.rm=TRUE)
  colnames(aggregated_data)[1:3] <- c("group", "month", "hour")
  colnames(aggregated_data)[4:ncol(aggregated_data)] <- variables

  # fill in missing month/hour combinations for all groups
  complete_seqs <- expand.grid("group"=unique(data[,split.by]), "hour"=formatC(0:23, 1, flag=0), "month"=all_months)
  missing_seqs <- dplyr::anti_join(complete_seqs, aggregated_data, by=c("group", "hour","month"))
  if(nrow(missing_seqs)>0){
    aggregated_data <- plyr::rbind.fill(aggregated_data, missing_seqs)
  }
  aggregated_data$month <- factor(aggregated_data$month, levels=all_months)
  aggregated_data <- aggregated_data[order(aggregated_data$group, aggregated_data$month, aggregated_data$hour),]
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
    daytimes_table <- getSunTimes(sunriset.coords, sunriset_start, sunriset_end, solar.depth, by="%m")
    daytimes_table <- daytimes_table[order(daytimes_table$interval),]
    daytimes_table$interval <- as.numeric(daytimes_table$interval)
  }

  # define hour labels to mark every two hours on the x-axis
  hour_labels <- paste0(formatC(0:23, 1, flag=0),"h")[c(TRUE, FALSE)]

  # set a default color palette if none is provided
  if(is.null(color.pal)){
    blues <- colorRampPalette(colors=c("#06405C", "#00537B", "#0985C2", "white"))(60)
    reds <- colorRampPalette(colors=c("white", "#B23E42", "#5C1315"))(40)
    color.pal <- colorRampPalette(colors=c(blues, reds))
  }

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
      contour_matrix <- matrix(aggregated_data[[i]][,var], nrow=24, ncol=12)
      colnames(contour_matrix) <- month.abb
      rownames(contour_matrix) <- formatC(0:23, 1, flag=0)
      contour_matrix <- rbind(contour_matrix, contour_matrix[1,])

      # set variable scale
      if(uniformize.scale){
        var_min <- min(complete_data[,var], na.rm=TRUE)
        var_max <- max(complete_data[,var], na.rm=TRUE)
        scale <- c(min(var_min), max(var_max))
      } else {
        scale <- range(contour_matrix, na.rm=TRUE)
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
      .filled.contour(x=1:12, y=0:24, z=t(contour_matrix), main=plot_title, cex.main=cex.main, invert.scale=invert.scale,
                      xlab="Months", ylab="Hours", nlevels=100, color.palette=color.pal, cex.lab=cex.lab, zlim=scale,
                      plot.axes = {
                        # configure x-axis for months and y-axis for hours
                        axis(1, labels=month.abb, at=1:12, cex.axis=cex.axis, pos=0)
                        axis(2, labels=hour_labels, at=seq(0, 23, by=2), cex.axis=cex.axis, las=1, pos=1)
                        axis(2, labels=FALSE, at=seq(0, 23, by=1), tck=-0.015, lwd.ticks=0.5, pos=1)
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
                          abline(v=1:12, lwd=0.03)
                          abline(h=seq(0, 23, by=1), lwd=0.03)
                        }}, ...)
      # add the sample size as a title below the main title
      title(main=plot_sub, line=0.7, font.main=2, cex.main=0.9)
      # set up the color scale legend
      scale_labs <- pretty(scale, min.n=4)
      scale_labs <- scale_labs[scale_labs>=min(scale) & scale_labs<=max(scale)]
      color_scale <- color.pal(100)
      if(invert.scale){scale<-rev(scale); color_scale<-rev(color_scale)}
      digits <- max(.decimalPlaces(scale_labs))
      .colorlegend(col=color_scale, zlim=scale, zval=scale_labs, digit=digits, xpd=TRUE,
                         posx=legend.xpos, posy=legend.ypos, main="", main.cex=0.8, cex=cex.axis)
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
