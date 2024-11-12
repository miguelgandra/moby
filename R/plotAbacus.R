#######################################################################################################
## Function to generate abacus plots  #################################################################
#######################################################################################################

#' Abacus plot
#'
#' @description Produces an abacus plot, displaying color-coded detections of
#' tagged individuals over the monitoring period. The function allows for various
#' customization options to cater to specific visualization needs, including
#' different color schemes, date formats, and handling of missing data.

#' @inheritParams setDefaults
#' @param data A data frame containing animal detections.
#' @param color.by Name of the column used to color-code individual detections.
#' This parameter can be used for example to differentiate detections by station,
#' animal trait, or temporal category. If NULL, all detections will be plotted in a single color.
#' @param tag.durations Optional. A numeric vector containing the estimated battery
#' duration of the deployed tags (in days). The length of this vector should match the number of
#' unique animal IDs, and the values must be in the same order as the ID levels.
#' Alternatively, if a single value is provided, it will be applied to all IDs.
#' @param id.groups Optional. A named list containing ID groups used to visually aggregate
#' animals belonging to the same class (e.g., different species, sexes or ages). Each element of
#' the list should be a vector of IDs that belong to the same group.
#' @param discard.missing Logical. If TRUE, animals without detections will not be
#' included in the plot. Defaults to FALSE.
#' @param color.pal A vector of colors defining the color palette used to plot detections per group level
#' (as defined by the 'color.by' argument).
#' @param date.format A string defining the date-time format for the x-axis labels
#' (as used in \code{\link[base]{strptime}}). Defaults to month ("%b").
#' @param date.interval A numeric value defining the interval between each
#' displayed date on the x-axis. Defaults to 4.
#' @param date.start An integer defining the first displayed date (can be used in
#' combination with 'date.interval' to better control the x-axis labels). Defaults to 1.
#' @param top.mural A string defining the date-time format for the optional top mural
#' containing additional time labels (as used in \code{\link[base]{strptime}}).
#' If set to FALSE, the top mural is not included. Defaults to year ("%Y").
#' @param season.shade Logical. If TRUE, the background is shaded based on
#' annual seasons. See \code{\link{shadeSeasons}} for details. Defaults to TRUE.
#' @param background.col A string specifying the background color of the plot.
#' Used if season.shade is set to FALSE. Defaults to "grey96".
#' @param pch Plotting character, i.e., the symbol to use for detections. Defaults to 16 (filled circle).
#' @param release.pch The symbol to use for the release/tagging date. Defaults to 8 (star).
#' @param end.pch The symbol to use for the estimated end date of the tag's lifetime. Defaults to 4 (cross).
#' @param pt.cex Numeric. The expansion factor for the points (detections). Defaults to 1.
#' @param transparency Numeric. The transparency level of the detections (points),
#' ranging from 0 (fully opaque) to 1 (fully transparent). Defaults to 0.
#' @param highlight.isolated Logical. If TRUE and a 'color.by' variable is defined,
#' detections are sorted and plotted based on their "density" distribution, with isolated detections
#' brought forward to prevent them from being hidden behind denser point clouds
#' (i.e., more densely clustered points of the same 'color.by' class). Defaults to TRUE.
#' @param cex.lab Determines the size of the axes titles. Defaults to 0.8.
#' @param cex.axis Determines the size of the tick mark labels on the axes. Defaults to 0.7.
#' @param cex.legend Determines the size of the color legend. Defaults to 0.7.
#' @param cex.mural Determines the size of the labels in the optional top mural. Defaults to 0.7.
#' @param legend.intersp Vertical distances between legend elements
#' (in lines of text shared above/below each legend entry). Defaults to 1.2.
#' @param legend.cols Integer. The number of columns in which to set the legend items.
#' If NULL, it is set automatically based on the number of levels. Defaults to NULL.
#'
#' @return Generates an abacus plot
#' @export


plotAbacus <- function(data,
                       id.col = getDefaults("id"),
                       datetime.col = getDefaults("datetime"),
                       color.by = NULL,
                       tagging.dates = getDefaults("tagging.dates"),
                       tag.durations = NULL,
                       id.groups=NULL,
                       discard.missing = FALSE,
                       color.pal = NULL,
                       date.format = "%b",
                       date.interval = 4,
                       date.start = 1,
                       top.mural = "%Y",
                       season.shade = TRUE,
                       background.col = "grey96",
                       pch = 16,
                       pt.cex = 1,
                       release.pch = 8,
                       end.pch = 4,
                       transparency = 0,
                       highlight.isolated = TRUE,
                       cex.lab = 0.8,
                       cex.axis = 0.7,
                       cex.legend = 0.7,
                       cex.mural = 0.7,
                       legend.intersp = 1.2,
                       legend.cols = NULL) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  tag.durations <- reviewed_params$tag.durations

  # check color.by variable
  if(!is.null(color.by) && !is.null(color.pal)) {
    if (length(color.pal) < nlevels(data[, color.by])) stop("The number of supplied colors needs to be greater than or equal to the number of color.by levels", call.=FALSE)
    if (length(color.pal) > nlevels(data[, color.by])) warning("The number of specified colors exceeds the number of levels in color.by", call.=FALSE)
  }

  # save the current par settings and ensure they are restored upon function exit
  original_par <- par(no.readonly=TRUE)
  on.exit(par(original_par))

  # print to console
  .printConsole("Generating abacus plot")


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # define single id.group if needed
  if(is.null(id.groups)){
    id.groups <- list(levels(data[,id.col]))
  }

  # discard missing animals if required
  missing_IDs <- which(table(data[,id.col])==0)
  if(discard.missing){
    if(length(missing_IDs)>0){
      tagging.dates <- tagging.dates[-missing_IDs]
      tag.durations <- tag.durations[-missing_IDs]
      data[,id.col] <- droplevels(data[,id.col])
      id.groups <- lapply(id.groups, function(x) x[!x %in% names(missing_IDs)])
    }
  }else{
    dummy_data <- data.frame("id"=names(missing_IDs))
    colnames(dummy_data) <- id.col
    data <- plyr::rbind.fill(data, dummy_data)
    data[,id.col] <- as.factor(data[,id.col])
  }

  # estimate tag lifetime dates
  if(!is.null(tag.durations)){
    end.dates <- NA
    for(e in 1:length(tagging.dates)){
      end.dates[e] <- tagging.dates[e] + tag.durations[e]*60*60*24
    }
  }

  # calculate nº of rows in plot
  row.groups <- id.groups
  for(i in 1:length(id.groups)){
    if(i==length(id.groups)){break}
    row.groups[[i]] <- c(id.groups[[i]], paste0("blank",i))
  }
  data$row_index <- data[,id.col]
  data$row_index <- factor(data$row_index, levels=unlist(row.groups))
  total_rows <- nlevels(data$row_index)
  data$row_index <- as.numeric(data$row_index)
  data$id_index <- as.numeric(data[,id.col])


  # color data by variable
  if(!is.null(color.by)) {
    groups <- levels(data[,color.by])
    ngroups <- nlevels(data[,color.by])
    if(is.null(color.pal)) {
      if(ngroups==3) {color.pal <- c("#326FA5","#D73134","#1C8E43")
      }else if (ngroups>3 & ngroups<=10){color.pal <- .economist_pal(ngroups)
      }else {color.pal <- topo.colors(ngroups)}
    }
    data$plot_color <- color.pal[data[,color.by]]

    # set nº of columns in legend if required
    if(is.null(legend.cols)){
      if(ngroups<(total_rows-5)){legend.cols<-1
      }else{legend.cols<-2}
    }

  }else{
    if(is.null(color.pal)) {
      data$plot_color <- "black"
    }else{
      data$plot_color <- color.pal[1]
    }
  }

  # set margins automatically
  mar <- c(4,4,2,10)
  if(!is.null(color.by)){
    if(legend.cols==2){mar <- c(4,4,2,14)}
    if(legend.cols>2){mar <- c(4,4,2,18)}
  }

  ##############################################################################
  ## Generate plot #############################################################
  ##############################################################################

  # set margins
  par(mar=mar, mgp=c(2.5,0.6,0), xpd=TRUE)
  xmin <- min(tagging.dates, na.rm=TRUE)
  xmax <- max(data[,datetime.col], na.rm=TRUE)
  # create empty plot
  plot(x=data[,datetime.col], y=data$id_index, ylim=c(total_rows+1, 0), axes=FALSE, type="n", xlab="", ylab="",
       xlim=c(xmin, xmax))
  date_lim1 <- lubridate::ceiling_date(as.POSIXct(par("usr")[1], origin='1970-01-01', tz="UTC"), "day")
  date_lim2 <- lubridate::floor_date(as.POSIXct(par("usr")[2], origin='1970-01-01', tz="UTC"), "day")
  complete_dates <- seq.POSIXt(date_lim1, date_lim2, by="day", tz="UTC")
  coords <- c()
  # add axis labels
  title(xlab="Date", cex.lab=0.8, line=1)
  if(length(id.groups)==1){
    title(ylab="Animal ID", cex.lab=cex.lab, line=2.4)
  }else{
    label_pos <- lapply(id.groups, function(x) mean(unique(data$row_index[data[,id.col] %in% x])))
    mtext(names(id.groups), side=2, at=label_pos, line=2.4, cex=cex.lab)
  }
  # shade seasons or add backgound
  if(season.shade) {
    seasons_table <- moby::shadeSeasons(date_lim1, date_lim2, interval=60*24)
    seasons_table$start <- as.numeric(seasons_table$start)
    seasons_table$end <- as.numeric(seasons_table$end)
    seasons_table$start[1] <- as.numeric(par("usr")[1])
    seasons_table$end[nrow(seasons_table)] <- par("usr")[2]
    rect(xleft=seasons_table$start, xright=seasons_table$end, ybottom=0, ytop=total_rows+1, col=seasons_table$color, border=NA)
    seasons_legend <- seasons_table[!duplicated(seasons_table$season), c("season","color")]
    seasons_legend <- seasons_legend[order(match(seasons_legend$season, c("spring", "summer", "autumn", "winter"))),]
    coords <- .legend(x=par("usr")[2], y=0, legend=seasons_legend$season, fill=seasons_legend$color,
                      bty="n", border="black", xpd=TRUE, y.intersp=1.6, box.cex=c(1.6, 1.2), cex=cex.legend)
  } else {
    rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=0, ytop=total_rows+1, col=background.col, border="black")
  }
  # plot detections of each individual
  for(i in 1:total_rows) {
     # subset data
    pts <- data[data$row_index==i,]
    if(nrow(pts)==0){next}
    id_index <- unique(pts$id_index)
    # add horizontal guide
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=i, lty=2, lwd=0.5)
    # if highlight.isolated, reorder detections using run length encoding,
    # plotting isolated detections in front
    if(highlight.isolated){
      pts <- pts[order(pts[,datetime.col]),]
      runs <- rle(pts$plot_color)
      runs_length <-  rep(runs$lengths, runs$lengths)
      runs_value <- rep(runs$values, runs$lengths)
      runs_table <- data.frame(runs_length, runs_value)
      sorted_rows <- order(runs_table$runs_length, decreasing=TRUE)
      pts <- pts[sorted_rows,]
    }
    # add points
    if(nrow(pts)>0){
      points(x=pts[,datetime.col], y=pts$row_index, pch=pch, col=adjustcolor(pts$plot_color, alpha.f=(1-transparency)), cex=pt.cex)
    }
    # add tagging symbols
    points(x=tagging.dates[id_index], y=i, pch=release.pch, cex=1.2)
    # add estimated tag end dates, if provided
    if(!is.null(tag.durations)){
      if(end.dates[id_index]<=date_lim2){
        points(x=end.dates[id_index], y=i, pch=end.pch, cex=1.2)
      }
    }
  }
  # add legend for release and tag end dates symbols
  if(length(coords)==0){legend.y<-0
  }else{legend.y<-abs(coords$rect$top-coords$rect$h)+1}
  legend.x <- par("usr")[2] + (par("usr")[2]-par("usr")[1])*0.01
  if(!is.null(tag.durations)){
    coords <- legend(x=legend.x, y=legend.y, legend=c("release date", "tag lifetime"),
                     pch=c(release.pch, end.pch), pt.cex=1.2, bty="n", y.intersp=1.4, cex=cex.legend)
  }else{
    coords <- legend(x=legend.x, y=legend.y, legend=c("release date"), pch=release.pch, pt.cex=1.2, bty="n", cex=cex.legend)
  }
  # prepare date variables
  all_dates <- strftime(complete_dates, date.format)
  consec_dates <- rle(all_dates)
  consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
  unique_dates <- unique(consec_dates)
  indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
  detec_dates <- strftime(seq.POSIXt(min(tagging.dates, na.rm=TRUE), max(data[,datetime.col], na.rm=TRUE), "day"), date.format)
  start <- min(which(sub("\\_.*", "", unique_dates)==detec_dates[1]))
  end <- max(which(sub("\\_.*", "", unique_dates)==detec_dates[length(detec_dates)]))
  indexes <- indexes[start:end]
  unique_dates <- unique_dates[start:end]
  disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
  disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
  disp_dates <- sub("\\_.*", "", disp_dates)
  # add axes
  axis(side=1, labels=disp_dates, at=complete_dates[disp_indexes], cex.axis=cex.axis, pos=total_rows+1)
  axis(side=1, labels=FALSE, at=complete_dates[indexes], pos=total_rows+1, tck=-0.006, lwd.ticks=0.5)
  id_rows <- stats::aggregate(data$row_index, by=list(data[,id.col]), unique)$x
  axis(side=2, labels=levels(data[,id.col]), at=id_rows, las=1, cex.axis=cex.axis)
  # add color legend if necessary
  if(!is.null(color.by)){
    legend.y <- abs(coords$rect$top-coords$rect$h)+1
    legend.x <- par("usr")[2] + (par("usr")[2]-par("usr")[1])*0.01
    legend(x=legend.x, y=legend.y, legend=groups, pch=pch, col=adjustcolor(color.pal, alpha.f=(1-transparency)),
           bty="n", border=NA, xpd=TRUE, y.intersp=legend.intersp, ncol=legend.cols, pt.cex=pt.cex*1.2, cex=cex.legend)
  }
  # add box
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=0, ytop=total_rows+1, col=NA, border="black", xpd=TRUE)
  if(top.mural!=FALSE){
    mural_height <- 1.3
    mural_gap <- (total_rows+1)/60
    mural_bottom <- 0 - mural_gap
    mural_center <-  0 - mural_gap - mural_height/2
    mural_top <- 0 - mural_gap - mural_height
    mural_vals <- strftime(complete_dates, top.mural)
    consec_vals <- rle(mural_vals)
    consec_vals <- paste0(mural_vals, "_", rep(1:length(consec_vals$lengths), consec_vals$lengths))
    mural_vals <- unique(consec_vals)
    mural_divs <- unlist(lapply(mural_vals, function(x) min(which(consec_vals==x))))
    mural_divs <- c(mural_divs, length(complete_dates))
    mural_vals <- sub("\\_.*", "", mural_vals)
    mural_disp <- zoo::rollapply(mural_divs, width=2, FUN=mean)
    # do not display 1st div if not enough space
    if(mural_divs[2]/length(complete_dates)<0.15){
      mural_vals <- mural_vals[-1]; mural_divs <- mural_divs[-1]; mural_disp <- mural_disp[-1]
    }
    # do not display last div if not enough space
    if(mural_divs[length(mural_divs)]/length(complete_dates)<0.15){
      mural_vals <- mural_vals[-length(mural_divs)]; mural_divs <- mural_divs[-length(mural_divs)]; mural_disp <- mural_disp[-length(mural_divs)]
    }
    # discard divs spanning further than the data
    valid_vals<- which(mural_vals %in% unique(strftime(data[,datetime.col], top.mural, tz="UTC")))
    mural_vals <- mural_vals[valid_vals]; mural_divs <- mural_divs[valid_vals]; mural_disp <- mural_disp[valid_vals]
    # generate rectangle
    rect(xleft=date_lim1, xright=date_lim2, ybottom=mural_bottom, ytop=mural_top, col="black", border="black")
    segments(x0=complete_dates[mural_divs], y0=mural_bottom, y1=mural_top, col="white", lwd=1.5)
    text(x=complete_dates[mural_disp], y=mural_center, labels=mural_vals, col="white", cex=cex.mural, adj=0.5)
  }
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
