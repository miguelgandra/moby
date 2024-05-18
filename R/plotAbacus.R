#######################################################################################################
## Function to generate abacus plots  #################################################################
#######################################################################################################

#' Abacus plot
#'
#' @description Produces an abacus plot, displaying color-coded detections of
#' tagged individuals over the monitoring period.
#'
#' @param data A data frame containing animal detections.
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param datetime.col Name of the column containing time bins in POSIXct format.
#' Defaults to 'datetime'.
#' @param color.by Variable defining the color group of individual detections.
#' Can be used for example to display detections by receiver, animal trait or temporal category.
#' @param tagging.dates A POSIXct vector containing the tag/release date of each animal.
#' @param tag.durations Optional. A  numeric vector containing the estimated battery
#' duration of the deployed tags (in days). If a single value is provided, it will be used to all IDs.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param discard.missing If true, animals without detections are not included. Defaults to false.
#' @param color.pal Color palette used to plot detections per group level
#' (as defined by the color.by argument). Defaults to ggthemes::economist_pal or
#' rainbow, depending on the number of 'color.by' levels.
#' @param date.format Date-time format (as used in \code{\link[base]{strptime}}),
#' defining the x-axis labels. Defaults to month ("%b").
#' @param date.interval Number defining the interval between each
#' displayed date (x-axis label). Defaults to 4.
#' @param date.start Integer defining the first displayed date (can be used in combination
#'  with 'date.interval" to better control the x-axis labels). Defaults to 1.
#' @param top.mural Include top mural containing additional time labels?
#' If so, defines the date-time format of the labels (as used by \code{\link[base]{strptime}}).
#' @param season.shade Boolean indicating if the background should be shaded
#' according with the respective season. See \code{\link{shadeSeasons}}.
#' @param background.col If season.shade is set to false, defines the background
#' color of the plot.
#' @param pch Plotting ‘character’, i.e., symbol to use. Defaults to 16.
#' @param pt.cex Expansion factor(s) for the points (detections).
#' @param transparency Transparency level of the detections (points).
#' @param highlight.isolated If a 'color.by' variable is defined, detections are ordered
#' and plotted according to their "density", with isolated detections being brought forward
#' to prevent them from being hidden behind denser point clouds.
#' @param transparency Transparency level of the detections (points).
#' @param cex.lab Determines the size of the y-axis and y-axis labels. Defaults to 0.8.
#' @param cex.axis Determines the size of the text labels on the axes. Defaults to 0.7.
#' @param cex.legend Determines the size of the color legend. Defaults to 0.7.
#' @param cex.mural Determines the size of the labels in the optional date top mural. Defaults to 0.7.
#' @param legend.intersp Vertical distances between legend elements
#' (in lines of text shared above/below each legend entry). Defaults to 1.2.
#' @param legend.cols Determines the  number of columns in which to set the legend items.
#' If left null, it is set automatically based on the number of levels.
#' @export

plotAbacus <- function(data, id.col="ID", datetime.col="datetime", color.by=NULL,
                       tagging.dates, tag.durations=NULL, id.groups=NULL, discard.missing=F,
                       color.pal=NULL, date.format="%b", date.interval=4, date.start=1, top.mural=F,
                       season.shade=F, background.col="grey96", pch=16, pt.cex=1, transparency=0,
                       highlight.isolated=T, cex.lab=0.8, cex.axis=0.7, cex.legend=0.7, cex.mural=0.7,
                       legend.intersp=1.2, legend.cols=NULL) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # print to console
  cat("Generating detections chronogram\n")

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please specify the correct column using 'id.col'")
  }

  # check if data contains datetime.col
  if(!datetime.col %in% colnames(data)){
    stop("Datetime column not found. Please specify the correct column using 'datetime.col'")
  }

  # check if timebins are in the right format
  if(!grepl("POSIXct", paste(class(data[,datetime.col]), collapse=" "))){
    stop("Datetimes must be provided in POSIXct format")
  }

  # check if tagging.dates are in the right format
  if(!grepl("POSIXct", paste(class(tagging.dates), collapse=" "))){
    stop("'tagging.dates' must be provided in POSIXct format")
  }

  # check color.by variable
  if(!is.null(color.by)){

    if(!color.by %in% colnames(data)) {
      stop("'color.by' variable not found in the supplied data")}

    if(any(is.na(data[,color.by]))){
      stop("Missing values in color.by variable")}

    if(class(data[,color.by])!="factor"){
      data[,color.by] <- as.factor(data[,color.by])
      cat("Warning: 'color.by' variable converted to factor\n")}

    if(length(color.pal)<nlevels(data[,color.by])){
      stop("The number of supplied colors needs to be greater than or equal to the number of color.by levels")}

    if(length(color.pal)>nlevels(data[,color.by])){
      cat("Warning: The number of specified colors exceeds the number of levels in color.by.")}
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {cat("Warning: Some of the ID(s) in id.groups don't match the IDs in the data")}
    data <- data[data[,id.col] %in% unlist(id.groups),]
    tagging.dates <- tagging.dates[match(unlist(id.groups), levels(data[,id.col]))]
    if(!is.null(tag.durations)){
      tag.durations <- tag.durations[match(unlist(id.groups), levels(data[,id.col]))]
    }
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
  }else{
    id.groups <- list(levels(data[,id.col]))
  }

  # check tag durations
  if(!is.null(tag.durations)){
    if(length(tag.durations)>1 & length(tag.durations)!=nlevels(data[,id.col])){
      stop("Incorrect number of tag.durations. Must be either a single value or
           a vector containing the estimated tag duration for each individual")
    }
    if(length(tag.durations)==1){
      tag.durations <- rep(tag.durations, nlevels(data[,id.col]))
    }
  }


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # discard missing animals if required
  if(discard.missing==T){
    missing_IDs <- which(table(data[,id.col])==0)
    if(length(missing_IDs)>0){
      tagging.dates <- tagging.dates[-missing_IDs]
      tag.durations <- tag.durations[-missing_IDs]
      data[,id.col] <- droplevels(data[,id.col])
      id.groups <- lapply(id.groups, function(x) x[!x %in% names(missing_IDs)])
    }
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
      }else if (ngroups>3 & ngroups<10){color.pal <- ggthemes::economist_pal()(ngroups)
      }else {color.pal <- rainbow(ngroups)}
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

  # prepare date variables
  start_date <- lubridate::floor_date(min(tagging.dates, na.rm=T), unit="day")
  end_date <- lubridate::ceiling_date(max(data[,datetime.col]), unit="day")
  complete_dates <- seq.POSIXt(start_date, end_date, by="day")
  all_dates <- strftime(complete_dates, date.format)
  consec_dates <- rle(all_dates)
  consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
  unique_dates <- unique(consec_dates)
  disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
  indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
  disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
  disp_dates <- sub("\\_.*", "", disp_dates)


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
  par(mar=mar, mgp=c(2.5,0.6,0), xpd=T)
  # create empty plot
  plot(x=data[,datetime.col], y=data$id_index, ylim=c(total_rows+1, 0), axes=F, type="n", xlab="", ylab="")
  date_lim1 <- as.POSIXct(par("usr")[1], origin='1970-01-01', tz="UTC")
  date_lim2 <- as.POSIXct(par("usr")[2], origin='1970-01-01', tz="UTC")
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
  if(season.shade==T) {
    seasons_table <- moby::shadeSeasons(date_lim1, date_lim2, interval=60*24)
    seasons_table$start <- as.numeric(seasons_table$start)
    seasons_table$end <- as.numeric(seasons_table$end)
    seasons_table$start[1] <- as.numeric(par("usr")[1])
    seasons_table$end[nrow(seasons_table)] <- par("usr")[2]
    rect(xleft=seasons_table$start, xright=seasons_table$end, ybottom=0, ytop=total_rows+1, col=seasons_table$color, border=NA)
    seasons_legend <- seasons_table[!duplicated(seasons_table$season), c("season","color")]
    seasons_legend <- seasons_legend[order(match(seasons_legend$season, c("spring", "summer", "autumn", "winter"))),]
    coords <- moby:::legend2(x=par("usr")[2], y=0, legend=seasons_legend$season, fill=seasons_legend$color,
                             bty="n", border="black", xpd=T, y.intersp=1.6, box.cex=c(1.6, 1.2), cex=cex.legend)
  } else {
    rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=0, ytop=total_rows+1, col=background.col, border="black")
  }
  # plot detections of each individual
  for(i in 1:total_rows) {
    # add horizontal guide
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=i, lty=2, lwd=0.5)
    # subset data
    pts <- data[data$row_index==i,]
    if(nrow(pts)==0){next}
    id_index <- unique(pts$id_index)
    # if highlight.isolated, reorder detections using run length encoding,
    # plotting isolated detections in front
    if(highlight.isolated==T){
      runs <- rle(pts$plot_color)
      runs_length <-  rep(runs$lengths, runs$lengths)
      runs_value <- rep(runs$values, runs$lengths)
      runs_table <- data.frame(runs_length, runs_value)
      sorted_rows <- order(runs_table$runs_length, decreasing=T)
      pts <- pts[sorted_rows,]
    }
    # add points
    if(nrow(pts)>0){
      points(x=pts[,datetime.col], y=pts$row_index, pch=pch, col=adjustcolor(pts$plot_color, alpha.f=(1-transparency)), cex=pt.cex)
    }
    # add tagging symbols
    points(x=tagging.dates[id_index], y=i, pch=4, cex=1.2)
    # add estimated tag end dates, if provided
    if(!is.null(tag.durations)){
      if(end.dates[id_index]<=date_lim2){
        points(x=end.dates[id_index], y=i, pch=8, cex=1.2)
      }
    }
  }
  # add legend for release and tag end dates symbols
  if(length(coords)==0){legend.y<-0
  }else{legend.y<-abs(coords$rect$top-coords$rect$h)+1}
  legend.x <- par("usr")[2] + (par("usr")[2]-par("usr")[1])*0.01
  if(!is.null(tag.durations)){
    coords <- legend(x=legend.x, y=legend.y, legend=c("release date", "tag lifetime"),
                     pch=c(4,8), pt.cex=1.2, bty="n", y.intersp=1.4, cex=cex.legend)
  }else{
    coords <- legend(x=legend.x, y=legend.y, legend=c("release date"), pch=4, pt.cex=1.2, bty="n", cex=cex.legend)
  }
  # add axes
  axis(side=1, labels=disp_dates, at=complete_dates[disp_indexes], cex.axis=cex.axis, pos=total_rows+1)
  axis(side=1, labels=F, at=complete_dates[indexes], pos=total_rows+1, tck=-0.006, lwd.ticks=0.5)
  id_rows <- aggregate(data$row_index, by=list(data[,id.col]), unique)$x
  axis(side=2, labels=levels(data[,id.col]), at=id_rows, las=1, cex.axis=cex.axis)
  # add color legend if necessary
  if(!is.null(color.by)){
    legend.y <- abs(coords$rect$top-coords$rect$h)+1
    legend.x <- par("usr")[2] + (par("usr")[2]-par("usr")[1])*0.01
    legend(x=legend.x, y=legend.y, legend=groups, pch=pch, col=adjustcolor(color.pal, alpha.f=(1-transparency)),
           bty="n", border=NA, xpd=T, y.intersp=legend.intersp, ncol=legend.cols, pt.cex=pt.cex*1.2, cex=cex.legend)
  }
  # add box
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=0, ytop=total_rows+1, col=NA, border="black", xpd=T)
  if(top.mural!=F){
    mural_height <- 1.3
    mural_gap <- (total_rows+1)/60
    mural_bottom <- 0 - mural_gap
    mural_center <-  0 - mural_gap - mural_height/2
    mural_top <- 0 - mural_gap - mural_height
    complete_dates <- seq.POSIXt(date_lim1, date_lim2, by="day", tz="UTC")
    #mural_vals <- unique(strftime(complete_dates, top.mural))
    mural_vals <- strftime(complete_dates, top.mural)
    consec_vals <- rle(mural_vals)
    consec_vals <- paste0(mural_vals, "_", rep(1:length(consec_vals$lengths), consec_vals$lengths))
    mural_vals <- unique(consec_vals)
    mural_divs <- unlist(lapply(mural_vals, function(x) min(which(consec_vals==x))))
    mural_vals <- sub("\\_.*", "", mural_vals)
    #mural_divs <- unlist(lapply(mural_vals, function(x) min(which(strftime(complete_dates, top.mural)==x))))
    mural_disp <- zoo::rollapply(mural_divs, width=2, FUN=mean)
    # do not display 1st div if not enough space
    if(mural_divs[2]/length(complete_dates)<0.06){
      mural_vals <- mural_vals[-1]; mural_divs <- mural_divs[-1]; mural_disp <- mural_disp[-1]
    }
    # discard divs spanning further than the data
    valid_vals<- which(mural_vals %in% unique(strftime(data[,datetime.col], top.mural, tz="UTC")))
    mural_vals <- mural_vals[valid_vals]; mural_divs <- mural_divs[valid_vals]; mural_disp <- mural_disp[valid_vals]
    # generate rectangle
    rect(xleft=date_lim1, xright=date_lim2, ybottom=mural_bottom, ytop=mural_top, col="black", border="black")
    segments(x0=complete_dates[mural_divs], y0=mural_bottom, y1=mural_top, col="white", lwd=1.5)
    text(x=complete_dates[mural_disp], y=mural_center, labels=mural_vals, col="white", cex=cex.mural, adj=0.5)
  }
  #reset par
  par(mar=c(5, 4, 4, 2) + 0.1)

}



#######################################################################################################
#######################################################################################################
#######################################################################################################
