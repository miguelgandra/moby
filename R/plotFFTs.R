#######################################################################################################
# Check for detections periodicity through FFTs #######################################################
#######################################################################################################

#' Analyze fine-scale rhythmic patterns in detections
#'
#' @description Function that analyzes and plots periodic patterns in detections using a
#' FFT (Fast Fourier Transform) framework. By decomposing data series into the
#' frequency domain, FFT allows the identification of dominant spectral peaks
#' that may reflect, for example, tidal (6 h–12 h) or diel (24 h) cyclic patterns in habitat use.
#'
#' @param data A data frame containing binned animal detections (must contain a "timebin" column).
#' @param tagging.dates A POSIXct vector containing the tag/release date of each animal.
#' @param type Type of response used to calculate FFTs. One of  "detections" or "presences".
#' Presences are estimated by time bin (0 - 1), independently of the number of detections.
#' @param axis.periods Periods (in hours) to display below the plots.
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param min.days Discard individuals that were detected in less than x days.
#' @param detrend Detrend time series using differences (\code{\link[base]{diff}})
#' rather than the actual detections.
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @param same.scale Forces same spectral scale (ylims) across all plots,
#' allowing for density comparison between individuals. If set to false, all y-axis are displayed.
#' Otherwise they are only displayed in the left-most plots to save space.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @export


plotFFTs <- function(data, tagging.dates, type="detections", axis.periods=c(48,24,12,8,6),
                     id.col="ID", min.days=10, detrend=T, cols=2, same.scale=T, id.groups=NULL) {

  # check if data contains id.col
  if(!id.col %in% colnames(data)) {
    stop("'id.col' variable not found in the supplied data")
  }

  # check if data contains timebin column
  if(!c("timebin") %in% colnames(data)) {
    stop("'timebin' column not found in the supplied data")
  }

  # check if type contains one of the accepted terms
  if(!type %in% c("detections", "presences")){
    stop("Style should be either 'detections' or 'presences'")
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {"Some of the ID(s) in id.groups don't match the IDs in the data"}
    data <- data[data[,id.col] %in% unlist(id.groups),]
    tagging.dates <- tagging.dates[match(unlist(id.groups), levels(data[,id.col]))]
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
  }

  # get time bin interval
  interval <- difftime(data$timebin, data.table::shift(data$timebin), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=T))
  time_fraction <- 60/interval

  # create data frame with number of detections
  fft_table <- createWideTable(data, start.dates=tagging.dates, id.col=id.col, value.type="detections")

  # if type = presences, convert detections to binary (0-1)
  if(type=="presences"){
    fft_table[,-1] <- apply(fft_table[,-1], 2, function(x) {x[x>1]<-1; return(x)})
  }

  # subset individuals based on minimum number of days with detections
  data$day <- strftime(data$timebin, "%Y-%m-%d", tz="UTC")
  days_detected <- by(data$day, data[,id.col], function(x) length(unique(x)))
  selected_individuals <- which(as.numeric(days_detected) >= min.days)
  nindividuals <- length(selected_individuals)
  cat(paste(nindividuals, "individuals with >", min.days, "logged days\n"))
  cat(paste(nlevels(data[,id.col])-nindividuals, "individual(s) excluded\n"))

  # split data by individual
  data_individual <- sapply(selected_individuals, function(i) fft_table[,i+1], simplify=F)
  data_individual <- lapply(data_individual, function(x) x[!is.na(x)])
  names(data_individual) <- levels(data[,id.col])[selected_individuals]
  data_ts <- lapply(data_individual, ts)

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

  if(same.scale==T){par(mfrow=c(rows, cols), mar=c(2,1,0,0), oma=c(6,6,1,1))}
  if(same.scale==F){par(mfrow=c(rows, cols), mar=c(2,4,0,4), oma=c(6,6,1,1))}

  # compute periodograms
  ffts <- list()
  if(detrend==T) {ffts <- lapply(data_ts, function(x) TSA::periodogram(diff(x), plot=F))}
  if(detrend==F) {ffts <- lapply(data_ts, function(x) TSA::periodogram(x, plot=F))}
  max_freq <- 1/(max(axis.periods*time_fraction)*1.25)
  max_density <- lapply(ffts, function(x) max(x$spec[x$freq>max_freq & x$freq<0.25]))
  max_density <- max(unlist(max_density))


  for (i in plot_ids) {

    # if NA, add empty plot and go to next
    if(is.na(i)) {
      plot.new()
      next
    }

    # else calculate and plot FFTs
    fft <- ffts[[i]]
    id <- selected_individuals[i]
    periods <- 1/fft$freq[fft$freq>max_freq & fft$freq<0.25]
    densities <- fft$spec[which(fft$freq>max_freq & fft$freq<0.25)]
    if(!is.null(id.groups)){
      group <- which(unlist(lapply(id.groups, function(x) levels(data[,id.col])[id] %in% x)))
      background_col<-background_pal[group]
    }

    if(same.scale==T) {
      plot(x=periods, y=densities, type="n", axes=F, ylab="", xlab="", ylim=c(0, max_density))
      rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4], col=background_col, border=NULL)
      points(x=periods, y=densities, type="h")
      # only add y axis to the left plots
      if(i %in% plot_layout[,1]){
        mtext(text="spectral density", side=2, line=4.5, cex=1.3)
        vals <- pretty(c(0, max_density))
        disp_vals <- vals
        if(any(vals%%1!=0)) {
          digits <- max(moby:::decimalPlaces(vals))
          disp_vals <- sprintf(paste0("%.", digits, "f"), disp_vals)
        }
        axis(2, at=vals, lab=disp_vals, las=1, cex.axis=1.8, col=NA, col.ticks=1)
      }}

    if(same.scale==F) {
      plot(x=periods, y=densities, type="n", axes=F, ylab="", xlab="")
      rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4], col=background_col, border=NULL)
      points(x=periods, y=densities, type="h")
      # add y axis to all plots
      mtext(text="spectral density", side=2, line=4.5, cex=1.3)
      vals <- pretty(densities)
      disp_vals <- vals
      if(any(vals%%1!=0)) {
        digits <- max(moby:::decimalPlaces(vals))
        disp_vals <- sprintf(paste0("%.", digits, "f"), disp_vals)
      }
      axis(2, at=vals, lab=disp_vals, las=1, cex.axis=1.8, col=NA, col.ticks=1)
     }

    # add x axis to the bottom plots
    if(i %in% bottom_plots){
      mtext(text="period (h)", side=1, line=3, cex=1.3)
      axis(1, at=axis.periods*time_fraction, labels=axis.periods, cex.axis=1.8, col=NA, col.ticks=1)
    }

    legend("topright", legend=levels(data[,id.col])[id], cex=2.2, bty="n")
    box()
  }

  if(!is.null(id.groups)){
    label_pos <- rev(unlist(lapply(group_rows, function(x) x/2)))
    for(i in 2:length(group_rows)){label_pos[i]<-label_pos[i]+rev(group_rows)[[i-1]]}
    label_pos <- grconvertY(label_pos/rows, "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels=rev(names(id.groups)),
         srt=90, cex=2.8, font=2, xpd=NA)
  }

  #reset par
  par(mar=c(5, 4, 4, 2) + 0.1)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
