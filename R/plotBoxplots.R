##################################################################################################
## Plot boxplots #################################################################################
##################################################################################################

#' Generates boxplots of multiple responses per grouping variable.
#'
#' @description Produces a panel with several boxplots (one per response), for
#' comparing multiple behavioural metrics across different grouping levels (usually time frames).
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections and distances traveled,
#' as returned by \code{\link{calculateDistances}}.
#' @param subsetted.overlap List containing pairwise overlaps estimated by a given
#' time-frame (e.g. diel phase), as returned by \code{\link{calculateOverlap}} using "subset".
#' Only required if overlap is included as a response.
#' @param subsetted.kuds List containing kernel utilization areas estimated by a given
#' time-frame (e.g. diel phase), as returned by \code{\link{calculateKUDs}} using "subset".
#' Only required if kud50 or kud95 are included as responses.
#' @param responses Metrics to include in the boxplots panel. Any of "detections",
#' "overlap", "total_dist, "kud50", "kud95" or any other included in the supplied data.
#' @param by Time-frame corresponding to the grouping category used to subset the data.
#' @param discard.incomplete Boolean to indicate if individuals with missing values (for any of the
#' grouping categories) should be discarded. If "by" contains only two grouping levels
#' (e.g. "day" and "night"), paired tests are computed whenever the individuals with data match between groups.
#' @param display.n If true, the number of individuals used to calculate each boxplot is
#' displayed together with the grouping labels.
#' @param box.colors Color of the box plot bodies.
#' @param outliers Logical. Specifies whether outliers should be included. Defaults to FALSE.
#' @param background.color Color for the plot background. Defaults to "gray96".
#' @param box.pattern Boolean to indicate if boxplot bodies should be filled with
#' different patterns.
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @export


################################################################################
# Main function - generate boxplots ####################################################

plotBoxplots <- function(data, subsetted.overlap=NULL, subsetted.kuds=NULL, responses, by,
                         id.col=getDefaults("id"), discard.incomplete=T, display.n=F, box.colors=NULL,
                         outliers=FALSE, background.color="gray96", box.pattern=F, cols=2){

  ################################################################################
  # Initial checks ###############################################################

  responses <- tolower(responses)
  if(any(responses=="overlap") & is.null(subsetted.overlap)){
    stop("Subsetted overlaps need to be supplied")
  }
  if(any(grepl("kud", responses, fixed=T)) & is.null(subsetted.kuds)){
    stop("Subsetted KUDs need to be supplied")
  }

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  data_responses <- responses[!responses %in% c("detections", "presences", "overlap", "total_dist", "kud50", "kud95")]
  missing_responses <- data_responses[!data_responses %in% colnames(data)]
  if(length(missing_responses)>0){
    missing_responses <- knitr::combine_words(shQuote(missing_responses))
    stop(paste(missing_responses, "response(s) not found in the supplied data"))
  }

  if(class(data[,by])!="factor"){
    cat("Warning: converting grouping variable to factor\n")
    data[,by] <- as.factor(data[,by])
  }
  groups <- levels(data[,by])

  if(!is.null(box.colors) & length(groups)!=length(box.colors)){
    stop("Number of box colors should match the number of group levels")
  }
  cat(paste0("Generating boxplots (", by, ")\n"))


  ################################################################################
  # Prepare data #################################################################

  # get time bins interval (in minutes)
  interval <- difftime(data$timebin, dplyr::lag(data$timebin), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=T))

  rows <- ceiling(length(responses)/cols)
  par(mfrow=c(rows, cols),  mar=c(4,4,4,2))
  for (i in 1:length(responses)){
    subset_vars <- c(by, id.col)
    subset_vars <- lapply(subset_vars, function(x) {return(data[,x])})
    if(responses[i]=="detections"){
      plot_data <- stats::aggregate(data$detections, by=subset_vars, function(x) {sum(x, na.rm=T)/(length(x)*interval/60)})
      plot_title <- expression(bold(paste("Detections (h"^"-1",")")))
    }
    else if(responses[i]=="presences"){
      data$presences <- data$detections
      data$presences[data$presences>=1] <- 1
      plot_data <- stats::aggregate(data$presences, by=subset_vars, function(x) {sum(x, na.rm=T)/(length(x)*interval/60)})
      plot_title <- expression(bold(paste("Presences (h"^"-1",")")))
    }
    else if(responses[i]=="overlap"){
      subsetted.overlap <- lapply(subsetted.overlap, function(x) as.numeric(x$overlap))
      plot_data <- reshape2::melt(subsetted.overlap)
      colnames(plot_data) <- c("x", "Group.1")
      plot_data$Group.2 <- rep(1:length(subsetted.overlap$day),2)
      plot_title <- "Overlap (%)"
    }
    else if(responses[i]=="total_dist"){
      plot_data <- stats::aggregate(data$dist_m, by=subset_vars, sum, na.rm=T)
      plot_data$x <-  plot_data$x/1000
      plot_title <- expression(bold("Total Distance (km)"))
    }
    else if(responses[i]=="kud50"){
      plot_data <- reshape2::melt(subsetted.kuds$k50_table[,c("id",groups)], id.vars="id")
      colnames(plot_data) <- c("Group.2", "Group.1", "x")
      plot_title <- expression(bold(paste("KUD 50% (km"^"2",")")))
    }
    else if(responses[i]=="kud95"){
      plot_data <- reshape2::melt(subsetted.kuds$k95_table[,c("id",groups)], id.vars="id")
      colnames(plot_data) <- c("Group.2", "Group.1", "x")
      plot_title <- expression(bold(paste("KUD 95% (km"^"2",")")))
    } else {
      form <- as.formula(paste(responses[i], "~", paste(c(by, id.col), collapse="+")))
      plot_data <- stats::aggregate(form, data=data, mean, na.rm=T)
      colnames(plot_data) <- c("Group.1", "Group.2", "x")
      plot_title <- tools::toTitleCase(responses[i])
    }
    plot_data$Group.1 <- as.factor(plot_data$Group.1)
    plot_data$Group.2 <- as.factor(plot_data$Group.2)
    plot_data$x <- suppressWarnings(as.numeric(plot_data$x))
    # discard individuals with missing values or simply truncate data?
    if(discard.incomplete==T){
      plot_data <- plot_data[!is.na(plot_data$x),]
      incomplete_ids <- names(table(plot_data$Group.2))[table(plot_data$Group.2)<length(groups)]
      plot_data <- plot_data[!plot_data$Group.2 %in% incomplete_ids,]
      plot_data$Group.2 <- droplevels(plot_data$Group.2)
    } else {
      plot_data <- plot_data[!is.na(plot_data$x),]
    }
    ids <- stats::aggregate(as.character(plot_data$Group.2), by=list(plot_data$Group.1),
                     function(x) paste(unique(x), collapse=", "), simplify=T)
    same_ids <- ids[1,"x"] == ids[2,"x"]

    # set box colors
    if(is.null(box.colors)){
      if(nlevels(plot_data$Group.1)==2){boxcolors<-c("white", "grey50")}
      if(nlevels(plot_data$Group.1)==4){boxcolors<-c("white","grey80","grey50","grey95")}
    }

    ################################################################################
    # Generate boxplots #############################################################

    p<-boxplot(x~Group.1, data=plot_data, plot=F)

    if(outliers==F){
      ymin <- min(p$stats)
      ymax <- max(p$stats)
      yrange <- ymax-ymin
    } else {
      ymin<-min(c(p$stats, p$out))
      ymax<-max(c(p$stats, p$out))
      yrange<-ymax-ymin
    }

    if(ymax<max(p$stats[nrow(p$stats),]+yrange*0.05)){ymax<-max(p$stats[nrow(p$stats),]+yrange*0.05)}
    boxplot(x~Group.1, data=plot_data, axes=F, main="", xlab="", ylab="", boxwex=0.45, ylim=c(ymin-yrange*0.04, ymax+yrange*0.04), outline=outliers)
    rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col=background.color, border=NULL)
    p<-boxplot(x~Group.1, data=plot_data, axes=F, col=boxcolors, xlab="", ylab="",
               boxwex=0.45, ylim=c(ymin-yrange*0.04, ymax+yrange*0.04), outline=outliers, add=T)
    title(main=plot_title, cex.main=1.6, font=2)
    if(box.pattern==T){
      rect(xleft=c(1:length(groups)-0.225), xright=c(1:length(groups)+0.225), ybottom=p$stats[2,], ytop=p$stats[4,],
           density=12, angle=c(45,-45))
    }
    if(display.n==F) {group_labs <- levels(plot_data$Group.1)}
    if(display.n==T) {group_labs <- paste0(levels(plot_data$Group.1), " (", p$n, ")")}
    text(x=1:nlevels(plot_data$Group.1), y=p$stats[nrow(p$stats),]+yrange*0.05, labels=group_labs, cex=1.2)
    yvals <- pretty(c(ymin, ymax))
    digits <- max(.decimalPlaces(yvals))
    ylabs <- sprintf(paste0("%.", digits, "f"), yvals)
    axis(2, at=yvals, labels=ylabs, cex.axis=1.5, las=1)
    legend("topleft", legend=LETTERS[i], bty="n", cex=2.4, inset=c(-0.14,-0.04))
    if(nlevels(plot_data$Group.1)==2) {
      grouped_data <- lapply(split(plot_data, f=plot_data$Group.1), function(x) return(x$x))
      if(same_ids==F){pair<-F; stat<-"W="}
      if(same_ids==T){pair<-T; stat<-"V="}
      test <- wilcox.test(grouped_data[[1]], grouped_data[[2]], paired=pair)
      test_out <-paste0(stat, sprintf("%.2f", test$statistic), "; p=", sprintf("%.3f",test$p.value))}
    if(nlevels(plot_data$Group.1)>2) {
      test <-kruskal.test(x~Group.1, data=plot_data)
      test_out <-paste0("X=", sprintf("%.2f", test$statistic), "; p=", sprintf("%.3f",test$p.value))}
    legend("bottom", legend=test_out, bty="n", cex=1.5, xpd=T, inset=c(0,-0.14))
    if(!is.nan(test$p.value) & test$p.value<0.05){box(col="darkred", lwd=1.5)}else{box(lwd=1.5)}

  }

}

##################################################################################################
##################################################################################################
##################################################################################################
