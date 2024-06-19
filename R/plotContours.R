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
#' @param data A data frame containing the variable(s) to plot.
#' @param variables Variables to plot. Should correspond columns in the data.
#' Can be either a vector or a single variable.
#' @param var.titles Display name of the variables (e.g. "Depth (m)"). If NULL, defaults to supplied variable names.
#' @param plot.title Optional. Main title displayed on top of the plots/panel.
#' @param id.col Name of the column containing animal IDs Defaults to 'ID'.
#' @param date.col Name of the column containing datetimes in POSICX format. Defaults to 'datetime'.
#' @param split.by Optional. If defined, plots are generated individually for each level (or combination of levels)
#' of this variable(s) (e.g. species, ontogeny, sex, habitat, etc).
#' @param aggregate.fun Function used to aggregate values by hour and month. Defaults to mean.
#' @param color.pal Color palette for the contour plot.
#' @param diel.lines Number indicating the number of diel phase lines (boundaries)
#' to display. Either 2 (corresponding to sunrise and sunset) or 4 (depicting
#' dawn, sunrise, sunset, dusk).
#' @param diel.lines.col Color of the diel lines. Defaults to black.
#' @param sunriset.coords A SpatialPoints or matrix object containing longitude and
#' latitude coordinates (in that order) at which to estimate sunrise and sunset times.
#' @param solar.depth Angle of the sun below the horizon (in degrees). Passed to the
#' solarDep argument in \code{\link[maptools]{crepuscule}} function. Defaults to 18 (astronomical twilight).
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


plotContours <- function(data, variables, var.titles=NULL, plot.title=NULL, split.by=NULL, aggregate.fun=function(x) mean(x, na.rm=T),
                         id.col="ID", date.col="datetime", color.pal=NULL, diel.lines=4, diel.lines.col="black",
                         sunriset.coords, sunriset.start=NULL, sunriset.end=NULL, solar.depth=18, cex.main=1.1, cex.lab=1, cex.axis=0.8,
                         grid=T, invert.scale=F, uniformize.scale=F, legend.xpos=c(0.89, 0.92), legend.ypos=c(0.15, 0.85),
                         tz="UTC", cols=1, disable.par=F, ...) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  if(any(!variables %in% colnames(data))){
    stop("Some of the variables not found in the supplied data")}

  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please specify the correct column using 'id.col'")}

  if(class(data[,id.col])!="factor"){
    cat("Converting ids to factor\n")
    data[,id.col] <- as.factor(data[,id.col])}

  if(!date.col %in% colnames(data)){
    stop("Datetime column not found. Please specify the correct column using 'date.col'")}

  if(any(class(data[,date.col])=="POSIXc")){
    stop("Please supply datetimes in POSIXct format")}

  # check grouping variable(s)
  if(!is.null(split.by)){
    if(!all(split.by %in% colnames(data))){
      stop("Groupping column not found. Please specify the correct column using 'split.by'")
    }
    # convert all grouping columns to factor
    for(s in split.by){
      if(class(data[,s])!="factor"){
        cat("Converting split.by column to factor\n")
        data[,s] <- as.factor(data[,s])
      }
    }
    # if more than one grouping variable, create a new column using interactions
    if(length(split.by)>1){
      data$dummy_group <- apply(data[split.by], 1, function(x) paste(x, collapse=" "))
      data$dummy_group <- as.factor(data$dummy_group)
      split.by <- "dummy_group"
    }
  # else create a single artificial group
  }else{
    data$dummy_group <- "1"
    data$dummy_group <- as.factor(data$dummy_group)
    split.by <- "dummy_group"
  }

  if(!is.null(var.titles) & length(variables)!=length(var.titles)){
    stop("Number of variables and variables' titles do not match")
  }

  if(is.null(color.pal)){
    color.pal <- .viridis_pal
  }

  if(class(color.pal)!="function"){
    color.pal <- colorRampPalette(color.pal)
  }


  # print to console
  cat(paste0("Generating contour plot(s)\n"))


  #####################################################################################
  # Aggregate data and calculate stats ################################################
  #####################################################################################

  # get month and hour
  data$month_tmp <- strftime(data[,date.col], "%m", tz=tz)
  data$hour_tmp <- strftime(data[,date.col], "%H", tz=tz)
  all_months <- sprintf("%02d", 1:12)

  # aggregate data
  aggregated_data <- stats::aggregate(data[,variables], by=list(data[,id.col], data$month_tmp, data$hour_tmp, data[,split.by]), aggregate.fun, simplify=T, drop=T)
  colnames(aggregated_data)[1:4] <- c("id", "month", "hour", "group")
  colnames(aggregated_data)[5:ncol(aggregated_data)] <- variables

  # calculate nº individuals with data for each variable and group
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

  # calculate variables' averages for each group (per month and hour)
  aggregated_data <- stats::aggregate(aggregated_data[,variables], by=list(aggregated_data$group, aggregated_data$month, aggregated_data$hour), mean, na.rm=T)
  colnames(aggregated_data)[1:3] <- c("group", "month", "hour")
  colnames(aggregated_data)[4:ncol(aggregated_data)] <- variables
  complete_seqs <- expand.grid("group"=unique(data[,split.by]), "hour"=formatC(0:23, 1, flag=0), "month"=all_months)
  missing_seqs <- dplyr::anti_join(complete_seqs, aggregated_data, by=c("group", "hour","month"))
  if(nrow(missing_seqs)>0){
    aggregated_data <- plyr::rbind.fill(aggregated_data, missing_seqs)
  }
  aggregated_data$month <- factor(aggregated_data$month, levels=all_months)
  aggregated_data <- aggregated_data[order(aggregated_data$group, aggregated_data$month, aggregated_data$hour),]
  aggregated_data <- split(aggregated_data, f=aggregated_data$group, drop=T)


  # join all data (to uniformize scales, if required)
  complete_data <- do.call("rbind", aggregated_data)


  #####################################################################################
  # Prepare plot variables ############################################################
  #####################################################################################

  sunriset_start <- min(data[,date.col])
  sunriset_end <- max(data[,date.col])
  daytimes_table <- getSunTimes(sunriset.coords, sunriset_start, sunriset_end, solar.depth, by="%m")
  daytimes_table <- daytimes_table[order(daytimes_table$interval),]
  daytimes_table$interval <- as.numeric(daytimes_table$interval)
  hour_labels <- paste0(formatC(0:23, 1, flag=0),"h")[c(TRUE, FALSE)]
  if(is.null(color.pal)){
    blues <- colorRampPalette(colors=c("#06405C", "#00537B", "#0985C2", "white"))(60)
    reds <- colorRampPalette(colors=c("white", "#B23E42", "#5C1315"))(40)
    color.pal <- colorRampPalette(colors=c(blues, reds))
  }

  # set layout
  if(disable.par==F){
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

  # iterate through each group
  for(i in 1:length(aggregated_data)){

    # iterate through each variable
    for(v in 1:length(variables)){

      # generate 2D matrix
      var <- variables[v]
      contour_matrix <- matrix(aggregated_data[[i]][,var], nrow=24, ncol=12)
      colnames(contour_matrix) <- month.abb
      rownames(contour_matrix) <- formatC(0:23, 1, flag=0)
      contour_matrix <- rbind(contour_matrix, contour_matrix[1,])

      # set variable scale
      if(uniformize.scale==T){
        var_min <- min(complete_data[,var], na.rm=T)
        var_max <- max(complete_data[,var], na.rm=T)
        scale <- c(min(var_min), max(var_max))
      } else {
        scale <- range(contour_matrix, na.rm=T)
      }

      # set plot title
      if(length(aggregated_data)>1){
        if(!is.null(var.titles)){plot_title <- paste(names(aggregated_data)[i], "-", var.titles[v])
        }else{plot_title <- paste(names(aggregated_data)[i], "-", var)}
      }else{
        if(!is.null(var.titles)){plot_title <- var.titles[v]
        }else{plot_title <- var}
      }

      plot_sub <- nids$nids[nids$group==names(aggregated_data)[i] & nids$var==var]
      plot_sub <- paste0("(n=", plot_sub, ")")

      # generate plot
      .filled.contour(x=1:12, y=0:24, z=t(contour_matrix), main=plot_title, cex.main=cex.main, invert.scale=invert.scale,
                      xlab="Months", ylab="Hours", nlevels=100, color.palette=color.pal, cex.lab=cex.lab, zlim=scale,
                      plot.axes = {
                        axis(1, labels=month.abb, at=1:12, cex.axis=cex.axis, pos=0)
                        axis(2, labels=hour_labels, at=seq(0, 23, by=2), cex.axis=cex.axis, las=1, pos=1)
                        axis(2, labels=F, at=seq(0, 23, by=1), tck=-0.015, lwd.ticks=0.5, pos=1)
                        if(diel.lines==4){lines(daytimes_table$interval, daytimes_table$dawns, lty=2, col=diel.lines.col)}
                        lines(daytimes_table$interval, daytimes_table$sunrises, lty=2, col=diel.lines.col)
                        lines(daytimes_table$interval, daytimes_table$sunsets, lty=2, col=diel.lines.col)
                        if(diel.lines==4){lines(daytimes_table$interval, daytimes_table$dusks, lty=2, col=diel.lines.col)}
                        if(grid==T){
                          abline(v=1:12, lwd=0.03)
                          abline(h=seq(0, 23, by=1), lwd=0.03)
                        }}, ...)
      title(main=plot_sub, line=0.7, font.main=2, cex.main=0.9)
      scale_labs <- pretty(scale, min.n=4)
      scale_labs <- scale_labs[scale_labs>=min(scale) & scale_labs<=max(scale)]
      color_scale <- color.pal(100)
      if(invert.scale==T){scale<-rev(scale); color_scale<-rev(color_scale)}
      digits <- max(.decimalPlaces(scale_labs))
      .colorlegend(col=color_scale, zlim=scale, zval=scale_labs, digit=digits, xpd=T,
                         posx=legend.xpos, posy=legend.ypos, main="", main.cex=0.8, cex=cex.axis)

    }
  }

  if(!is.null(plot.title)){
    mtext(text=plot.title, side=3, cex=1.25, font=2, line=-0.5, outer=T)
  }

}

#######################################################################################################
#######################################################################################################
#######################################################################################################
