#######################################################################################################
## Function to generate abacus plots  #################################################################
#######################################################################################################

#' Abacus plot
#'
#' @description Produces an abacus plot, displaying color-coded detections of
#' tagged individuals over the monitoring period.
#'
#' @param data A data frame containing binned animal detections.
#' @param tagging.dates A POSIXct vector containing the tag/release date of each animal.
#' @param date.format Date-time format (as used in \code{\link[base]{strptime}}),
#' defining the x-axis labels. Defaults to month ("%b").
#' @param date.interval Number defining the interval between each
#' displayed date (x-axis label). Defaults to 4.
#' @param date.start Integer defining the first displayed date (can be used in combination
#'  with 'date.interval" to better control the x-axis labels). Defaults to 1.
#' @param discard.missing If true, animals without detections are not included. Defaults to false.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param top.mural Include top mural containing additional time labels?
#' If so, defines the date-time format of the labels (as used by \code{\link[base]{strptime}}).
#' @param season.shade Boolean indicating if the background should be shaded
#' according with the respective season. See \code{\link{shadeSeasons}}.
#' @param background.col If season.shade is set to false, defines the background
#' color of the plot.
#' @param color.by Variable defining the color group of individual detections.
#' Can be used for example to display detections by receiver, animal trait or temporal category.
#' @param color.pal Color palette used to plot detections per group level
#' (as defined by the color.by argument). Defaults to ggthemes::economist_pal or
#' rainbow, depending on the number of 'color.by' levels.
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param pch Plotting ‘character’, i.e., symbol to use. Defaults to 16.
#' @param pt.cex Expansion factor(s) for the points (detections).
#' @param transparency Transparency level of the detections (points).
#' @param highlight.isolated If a 'color.by' variable is defined, detections are ordered
#' and plotted according to their "density", with isolated detections being brought forward
#' to prevent them from being hidden behind denser point clouds.
#' @export

plotAbacus <- function(data, tagging.dates, date.format="%b", date.interval=4, date.start=1, discard.missing=F,
                       id.groups=NULL, top.mural=F, season.shade=F, background.col="grey96", color.by=NULL,
                       color.pal=NULL, id.col="ID", pch=16, pt.cex=1, transparency=0, highlight.isolated=T) {

  #####################################################################################
  # Prepare data ######################################################################

  # print to console
  cat("Generating detections chronogram\n")

  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please specify the correct column using 'id.col'")}

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
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
  }

  # discard missing animals if required
  if(discard.missing==T){
    missing_IDs <- which(table(data[,id.col])==0)
    if(length(missing_IDs)>0){
      tagging.dates <- tagging.dates[-missing_IDs]
      data[,id.col] <- droplevels(data[,id.col])
    }
  }

  # create data table
  table <- createWideTable(data, start.dates=tagging.dates, value.col="receiver", id.col=id.col, verbose=F)

  # retrieve last detections dates
  last_detections <- tapply(X=data$timebin, INDEX=data[,id.col], FUN=max)
  last_detections <- as.POSIXct(last_detections, origin='1970-01-01', tz="UTC")

  # get time bins interval (in minutes)
  interval <- difftime(data$timebin, data.table::shift(data$timebin), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=T))

  # add empty lines if IDs groups were supplied
  if(!is.null(id.groups)){
    n_ids <-  unlist(lapply(id.groups, function(x) length(x[x %in% colnames(table)])))+1
    sep_indexes <- n_ids[-length(n_ids)]
    table <- tibble::add_column(table, "SPACE"=NA, .after=sep_indexes)
  }

  # format data
  timebins <- table$timebin
  data_plot <- as.matrix(table[,-1])
  replace_values <- which(!is.na(data_plot), arr.ind=T)
  data_plot[replace_values] <- replace_values[,2]
  #data_plot[replace_values] <- (1:ncol(data_plot))[replace_values[,2]]
  color_matrix <- as.matrix(table[,-1])
  color_matrix[!is.na(color_matrix)] <- "black"
  id_cols <- which(colnames(data_plot) %in% levels(data[,id.col]))
  total_cols <- ncol(data_plot)
  label_pos <- unlist(lapply(id.groups, function(x) mean(which(colnames(table)[-1] %in% x))))

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
    replace_timebins <- timebins[replace_values[,1]]
    replace_ids <- colnames(data_plot)[replace_values[,2]]
    data_aggr <- stats::aggregate(data$plot_color, by=list(data$timebin, data[,id.col]), function(x) names(which.max(table(x))))
    colnames(data_aggr) <- c("timebin", "ID", "plot_color")
    color_matrix[replace_values] <- data_aggr$plot_color[data_aggr$timebin==replace_timebins & data_aggr$ID==replace_ids]
  }

  # prepare date variables
  all_dates <- strftime(timebins, date.format)
  consec_dates <- rle(all_dates)
  consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
  unique_dates <- unique(consec_dates)
  disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
  indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
  disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
  disp_dates <- sub("\\_.*", "", disp_dates)
  if(top.mural!=F) {
    mural_vals <- unique(strftime(timebins, top.mural))
    mural_divs <- unlist(lapply(mural_vals, function(x) min(which(strftime(timebins, top.mural)==x))))
    mural_divs <- mural_divs[-1]
  }

  #####################################################################################
  # Generate chronogram ######################################################################
  par(mar=c(4,4,2,6), mgp=c(2.5,0.6,0), xpd=T)
  plot(x=timebins, y=data_plot[,1], ylim=c(total_cols+1, 0), axes=F, type="n", xlab="", ylab="")
  title(xlab="Date", cex.lab=0.8, line=1.6)
  if(is.null(id.groups)){title(ylab="Animal ID", cex.lab=0.8, line=2.2)}
  if(!is.null(id.groups)){mtext(names(id.groups), side=2, at=label_pos, line=2.2, cex=0.8)}
  if(season.shade==T) {
    seasons_table <- moby::shadeSeasons(min(tagging.dates, na.rm=T), max(last_detections, na.rm=T), interval)
    seasons_table$start <- as.numeric(seasons_table$start)
    seasons_table$end <- as.numeric(seasons_table$end)
    seasons_table$start[1] <- as.numeric(par("usr")[1])
    seasons_table$end[nrow(seasons_table)] <- par("usr")[2]
    rect(xleft=seasons_table$start, xright=seasons_table$end, ybottom=0, ytop=total_cols+1, col=seasons_table$color, border=NA)
    seasons_legend <- seasons_table[!duplicated(seasons_table$season), c("season","color")]
    seasons_legend <- seasons_legend[order(match(seasons_legend$season, c("spring", "summer", "autumn", "winter"))),]
    moby:::legend2("topright", legend=seasons_legend$season, fill=seasons_legend$color, bty="n", border="black",
           inset=c(-0.1, 0.2), xpd=T, y.intersp=1.8, box.cex=c(1.8, 1.4), cex=0.6)
  } else {
    rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=0, ytop=total_cols+1, col=background.col, border="black")
  }
  for(i in 1:nlevels(data[,id.col])) {
    index <- id_cols[i]
    segments(x0=par("usr")[1], x1=par("usr")[2], y0=index, lty=2, lwd=0.5)
    pts <- data.frame("timebins"=table$timebin, "detections"=data_plot[,index],
                      "colors"=color_matrix[,index])
    pts <- pts[!is.na(pts$detections),]
    # if highlight.isolated, reorder detections using run length encoding,
    # plotting isolated detections in front
    if(highlight.isolated==T){
      runs <- rle(pts$colors)
      runs_length <-  rep(runs$lengths, runs$lengths)
      runs_value <- rep(runs$values, runs$lengths)
      runs_table <- data.frame(runs_length, runs_value)
      sorted_rows <- order(runs_table$runs_length, decreasing=T)
      pts <- pts[sorted_rows,]
    }
    # add points and tagging symbols
    if(nrow(pts)>0){
      points(x=pts$timebins, y=pts$detections, pch=pch, col=adjustcolor(pts$colors, alpha.f=(1-transparency)), cex=pt.cex)
    }
    points(x=tagging.dates[i], y=index, pch=4, cex=1.2)

  }
  axis(side=1, labels=disp_dates, at=timebins[disp_indexes], cex.axis=0.65, pos=total_cols+1)
  axis(side=1, labels=F, at=timebins[indexes], pos=total_cols+1, tck=-0.006, lwd.ticks=0.5)
  axis(side=2, labels=levels(data[,id.col]), at=id_cols, las=1, cex.axis=0.75)
  if(!is.null(color.by)){
    legend("right", legend=groups, pch=pch, col=adjustcolor(color.pal, alpha.f=(1-transparency)), bty="n", border=NA,
           inset=c(-0.1, 0), xpd=T, y.intersp=1.4, pt.cex=pt.cex*1.2, cex=0.6)
  }
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=0, ytop=total_cols+1, col=NA, border="black")
  if(top.mural!=F){
    usr <- par("usr")
    mural_height <- 1.3
    mural_gap <- (total_cols+1)/60
    mural_bottom <- 0 - mural_gap
    mural_center <-  0 - mural_gap - mural_height/2
    mural_top <- 0 - mural_gap - mural_height
    rect(xleft=usr[1], xright=usr[2], ybottom=mural_bottom, ytop=mural_top, col="black", border="black")
    mural_disp <- zoo::rollapply(c(usr[1], mural_divs, usr[2]), width=2, FUN=mean)
    segments(x0=timebins[mural_divs], y0=mural_bottom, y1=mural_top, col="white", lwd=1.5)
    text(x=timebins[mural_disp], y=mural_center, labels=mural_vals, col="white", cex=0.7, adj=0.5)
  }
  #reset par
  par(mar=c(5, 4, 4, 2) + 0.1)

}

#######################################################################################################
#######################################################################################################
#######################################################################################################
