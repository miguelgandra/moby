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
#' @inheritParams setDefaults
#' @param data A data frame containing binned animal detections.
#' @param type Type of response used to calculate FFTs. One of  "detections" or "presences".
#' Presences are estimated by time bin (0 - 1), independently of the number of detections.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param axis.periods Periods (in hours) to display below the plots.
#' @param min.days Discard individuals that were detected in less than x days.
#' @param detrend Detrend time series using differences (\code{\link[base]{diff}})
#' rather than the actual detections.
#' @param same.scale Forces same spectral scale (ylims) across all plots,
#' allowing for density comparison between individuals. If set to false, all y-axis are displayed.
#' Otherwise they are only displayed in the left-most plots to save space.
#' @param background.col A string specifying the background color of the plot. Defaults to "grey94".
#' @param cex.title Determines the size of the plot title (animal ID). Defaults to 2.
#' @param cex.lab Determines the size of the axes titles. Defaults to 1.8.
#' @param cex.axis Determines the size of the tick mark labels on the axes. Defaults to 1.4.
#' @param cols Number of columns in the final panel.
#' @export


plotFFTs <- function(data,
                     tagging.dates = getDefaults("tagging.dates"),
                     id.col = getDefaults("ID"),
                     timebin.col = getDefaults("timebin"),
                     type = "detections",
                     id.groups = NULL,
                     axis.periods = c(48,24,12,8,6),
                     min.days = 10,
                     detrend = TRUE,
                     same.scale = TRUE,
                     background.col = "grey94",
                     cex.title = 2,
                     cex.lab = 1.8,
                     cex.axis = 1.4,
                     cols = 2) {


  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates

  # check if type contains one of the accepted terms
  if(!type %in% c("detections", "presences")) stop("Style should be either 'detections' or 'presences'", call.=FALSE)

  # check if TSA package is installed
  if (!requireNamespace("TSA", quietly=TRUE)) {
    stop("The 'TSA' package is required for this function but is not installed. Please install 'TSA' using install.packages('TSA') and try again.")
  }

  # save the current par settings and ensure they are restored upon function exit
  original_par <- par(no.readonly=TRUE)
  on.exit(par(original_par))


  ##############################################################################
  # Prepare data ###############################################################
  ##############################################################################

  # define single id.group if needed
  if(is.null(id.groups)){
    id.groups <- list(levels(data[,id.col]))
  }

  # get time bin interval
  interval <- difftime(data[,timebin.col], dplyr::lag(data[,timebin.col]), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=TRUE))
  time_fraction <- 60/interval

  # create data frame with number of detections
  fft_table <- suppressWarnings(createWideTable(data, id.col=id.col, timebin.col=timebin.col,
                                                value.col="detections", start.dates=tagging.dates))

  # if type = presences, convert detections to binary (0-1)
  if(type=="presences"){
    fft_table[,-1] <- apply(fft_table[,-1], 2, function(x) {x[x>1]<-1; return(x)})
  }

  # subset individuals based on minimum number of days with detections
  nfish <- nlevels(data[,id.col])
  data$day <- strftime(data[,timebin.col], "%Y-%m-%d", tz="UTC")
  days_detected <- by(data$day, data[,id.col], function(x) length(unique(x)))
  missing_ids <- which(table(data[,id.col])==0)
  missing_ids <- names(missing_ids)
  discarded_ids <- names(days_detected)[as.numeric(days_detected)<min.days]
  discarded_ids <- c(missing_ids, discarded_ids[!is.na(discarded_ids)])
  if(length(discarded_ids)>0){
    data <- data[!data[,id.col] %in% discarded_ids,]
    data[,id.col] <- droplevels(data[,id.col])
    id.groups <- lapply(id.groups, function(x) x[!x %in% discarded_ids])
    selected_individuals <- levels(data[,id.col])
    nindividuals <- length(selected_individuals)
    cat(paste(nindividuals, "individuals with >", min.days, "logged days\n"))
    cat(paste(nfish-nindividuals, "individual(s) excluded\n"))
  }

  # split data by individual
  data_individual <- sapply(selected_individuals, function(i) fft_table[,i], simplify=FALSE)
  data_individual <- lapply(data_individual, function(x) x[!is.na(x)])
  names(data_individual) <- selected_individuals
  data_ts <- lapply(data_individual, ts)


  ##############################################################################
  # Set layout variables #######################################################
  ##############################################################################

  # set layout grid
  layout_params <- .setLayout(cols, id.groups, plots.height=6, dividers.height=2, legend=FALSE)
  bottom_indices <- apply(layout_params$matrix, 2, function(x) {
    plots_rle <- rle(!is.na(x))
    return(x[cumsum(plots_rle$lengths)[plots_rle$values]])
    }, simplify=FALSE)
  bottom_indices <- unlist(bottom_indices)
  bottom_indices <- bottom_indices[order(bottom_indices)]
  nplots <- max(layout_params$matrix, na.rm=TRUE)
  layout_params$matrix[is.na(layout_params$matrix)] <- nplots + 1
  layout(mat=layout_params$matrix, heights=layout_params$heights)

  # set margins
  if(length(id.groups)>1 && !same.scale){
    oma <- c(4, 3, 1, 1)
    mar <- c(1.5, 8, 0, 0)
  }else if(length(id.groups)>1 && same.scale){
    oma <- c(4, 10, 1, 1)
    mar <- c(1.5, 1, 0, 0)
  }else if(length(id.groups)==1 && !same.scale){
    oma <- c(4, 1, 1, 1)
    mar <- c(1.5, 8, 0, 0)
  }else{
    oma <- c(4, 6, 1, 1)
    mar <- c(1.5, 1, 0, 0)
  }
  # oma <- if(length(id.groups)>1) c(4,3,1,1) else c(4,1,1,1)
  # mar <- if(same.scale==TRUE) c(1.5,1,0,0) else c(1.5,8,0,0)
  par(mar=mar, oma=oma, mgp=c(3,0.8,0))

  # compute periodograms
  ffts <- list()
  if(detrend) ffts <- lapply(data_ts, function(x) TSA::periodogram(diff(x), plot=FALSE))
  else ffts <- lapply(data_ts, function(x) TSA::periodogram(x, plot=FALSE))
  max_freq <- 1/(max(axis.periods*time_fraction)*1.25)
  max_density <- lapply(ffts, function(x) max(x$spec[x$freq>max_freq & x$freq<0.25]))
  max_density <- max(unlist(max_density))
  all_densities <- lapply(ffts, function(x) x$spec[x$freq>max_freq & x$freq<0.25])
  digits <- lapply(all_densities, function(x) ifelse(any(pretty(x)%%1!=0), max(.decimalPlaces(pretty(x))), 0))
  max_digits <- max(unlist(digits))


  ##############################################################################
  # Draw plots #################################################################
  ##############################################################################

  # iterate through each individual
  for (i in 1:nplots) {

    # retrieve FFT results
    fft <- ffts[[i]]
    id <- selected_individuals[i]
    periods <- 1/fft$freq[fft$freq>max_freq & fft$freq<0.25]
    densities <- fft$spec[which(fft$freq>max_freq & fft$freq<0.25)]

    # define y-axis range and labels
    if(same.scale){
      yrange <- c(0, max_density)
      ylabels <- pretty(yrange)
    }else{
      yrange <- NULL
      ylabels <- pretty(densities)
    }

    # draw plot
    plot(x=periods, y=densities, type="n", axes=FALSE, ylab="", xlab="", ylim=yrange)
    rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4], col=background.col, border=NULL)
    points(x=periods, y=densities, type="h")

    # draw y axis
    if(!same.scale || i %in% layout_params$matrix[,1]){
      title(ylab="Spectral density", line=4.8, cex.lab=cex.lab, xpd=NA)
      disp_vals <- sprintf(paste0("%.", max_digits, "f"), ylabels)
      axis(2, at=ylabels, labels=disp_vals, las=1, cex.axis=cex.axis, col=NA, col.ticks=1)
    }

    # draw x-axis in the bottom plots
    if(i %in% bottom_indices){
      title(xlab="Period (h)", line=3, cex.lab=cex.lab, xpd=NA)
      axis(1, at=axis.periods*time_fraction, labels=axis.periods, cex.axis=cex.axis, col=NA, col.ticks=1)
    }

    # add animal ID
    legend("topright", legend=id, text.font=2, cex=cex.title, bty="n")

    # draw box
    box()

  }

  #################################################################
  # add id.group labels ###########################################
  if(!is.null(id.groups)){
    label_pos <- layout_params$group_positions
    layout_height <- sum(layout_params$heights) + oma[1] + oma[3]
    label_pos <- 1 - ((label_pos + oma[3]) / layout_height)
    label_pos <- grconvertY(label_pos, "ndc", "user")
    x_pos <- grconvertX(0.01, "ndc", "user")
    text(x=x_pos, y=label_pos, labels=names(id.groups),
         srt=90, cex = cex.title+0.2, font=2, xpd=NA, adj=c(0.5, 0.5))
  }


}

#######################################################################################################
#######################################################################################################
#######################################################################################################
