#######################################################################################################
## Plot KUD time series ###############################################################################
#######################################################################################################

#' Plot KUD time series
#'
#' @description Function to plot the variation of 95% and 50% kernel utilization areas
#' across a given time-frame (usually month).
#'
#' @param subsetted.kuds List containing kernel utilization areas estimated by
#' time-frame (e.g. month), as returned by \code{\link{calculateKUDs} using "subset"}
#' @param label Time-frame used to subset KUDs, displayed on the x-axis (e.g. "Month").
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param plot.by Either "group" or "kud". If "group", the function plots a graph by id.group (overlaying KUD 50% and KUD95%).
#' If by "kud" the function plots a graph by KUD percentage (overlaying group series).
#' @param shade.season Boolean to indicate if seasons should be depicted in the background
#' (useful when plotting KUDs by month). Defaults to false.
#' @seealso \code{\link{calculateKUDs}}
#' @export


plotTemporalKUDs <- function(subsetted.kuds, label, id.groups=NULL, plot.by="group",
                             shade.season=F, same.scale=T) {

  #####################################################################################
  ## Prepare data #####################################################################

  cat(paste0("Plotting KUDs by ", label, "\n"))

  # check if data contains dist.col
  if(!plot.by %in% c("group", "kud")){
    stop("plot.by argument can only be set to  'group' or 'kud'")
  }

  # set id groups
  if(is.null(id.groups)) {
    id.groups <- list(subsetted.kuds$kud_table$ID)
  } else {
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% subsetted.kuds$kud_table$ID)) {"Some of the ID(s) in id.groups are not present in the supplied KUDs"}
  }

  k50_cols <- which(grepl("50%", colnames(subsetted.kuds$kud_table), fixed=T))
  k95_cols <- which(grepl("95%", colnames(subsetted.kuds$kud_table), fixed=T))
  k50_table <- subsetted.kuds$kud_table[,c(1, k50_cols)]
  k95_table <- subsetted.kuds$kud_table[,c(1, k95_cols)]
  colnames(k50_table) <- gsub("KUD 50% Area (km2) - ", "", colnames(k50_table), fixed=T)
  colnames(k95_table) <- gsub("KUD 95% Area (km2) - ", "", colnames(k95_table), fixed=T)
  suppressWarnings(k50_table[,-1] <- apply(k50_table[,-1], 2, as.numeric))
  suppressWarnings(k95_table[,-1] <- apply(k95_table[,-1], 2, as.numeric))
  k50_list <- lapply(id.groups, function(x) k50_table[k50_table$ID %in% x,])
  k95_list <- lapply(id.groups, function(x) k95_table[k95_table$ID %in% x,])


  # calculate stats
  k50_means <- lapply(k50_list, function(x) apply(x[,-1], 2, mean, na.rm=T))
  k50_se <- lapply(k50_list, function(x) apply(x[,-1], 2, plotrix::std.error))
  k95_means <- lapply(k95_list, function(x) apply(x[,-1], 2, mean, na.rm=T))
  k95_se <- lapply(k95_list, function(x) apply(x[,-1], 2, plotrix::std.error))
  ngroups <- length(k50_means[[1]])


  #####################################################################################
  ## Calculate layout variables #######################################################

  k50_min <- mapply(function(a, b) {min(a-b, na.rm=T)}, a=k50_means, b=k50_se)
  k50_max <- mapply(function(a, b) {max(a+b, na.rm=T)}, a=k50_means, b=k50_se)
  k95_min <- mapply(function(a, b) {min(a-b, na.rm=T)}, a=k95_means, b=k95_se)
  k95_max <- mapply(function(a, b) {max(a+b, na.rm=T)}, a=k95_means, b=k95_se)
  k50_total_range <- range(c(k50_min, k50_max))
  k95_total_range <- range(c(k95_min, k95_max))



  #####################################################################################
  ## Plot graphs ######################################################################

  par(mar=c(5,5,2,5))

  ##############################################################################
  ## Plot a graph by id.group (overlaying KUD 50% and KUD95%) ##################
  if(plot.by=="group"){
    par(mfrow=c(length(id.groups),1))
    for(i in 1:length(id.groups)){

      ## KUD 50% ##################################
      if(same.scale==T){
        yrange <- k50_total_range
      } else {
        yrange <- range(c(k50_min[i], k50_max[i]))
      }
      ydiff <- yrange[2] - yrange[1]
      yvals <- pretty(yrange)
      digits <- max(.decimalPlaces(yvals))
      ylabs <- sprintf(paste0("%.", digits, "f"), yvals)
      plot(0, axes=F, ylab="", xlab=label, xlim=c(1, ngroups), ylim=c(yrange[1]-ydiff*0.05, yrange[2]+ydiff*0.05), cex.lab=1.2)
      title(main=names(id.groups)[i], cex.main=1.2)
      mtext(expression(paste("KUD 50% (km"^"2"*")")), side=2, cex=1.2, line=3.5)
      if(shade.season==T) {
        rect(xleft=c(1,3,6,9), xright=c(3,6,9,12), ybottom=par("usr")[3], ytop=par("usr")[4],
             col=c("white","grey88", "grey72", "grey80"), border=NA)
      }
      if(shade.season==F) {
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="grey96", border=NA)
      }
      axis(1, at=1:ngroups, labels=colnames(k50_table)[-1], cex.axis=1.1)
      axis(2, at=yvals, labels=ylabs, cex.axis=1.1, las=1)
      points(k50_means[[i]], type="o", pch=16, cex=1)
      suppressWarnings(
        arrows(1:ngroups, k50_means[[i]]+k50_se[[i]], 1:ngroups, k50_means[[i]]-k50_se[[i]],
               angle=90, code=3, length=0.03, col="black", lwd=0.5))


      ## KUD 95% ##################################
      par(new=T)
      if(same.scale==T){
        yrange <- k95_total_range
      } else{
        yrange <- range(c(k95_min[i], k95_max[i]))
      }
      yvals <- pretty(yrange)
      ydiff <- yrange[2] - yrange[1]
      ylabs <- sprintf("%.1f", yvals)
      plot(k95_means[[i]], axes=F, type="o", lty=2, pch=15, cex=1, xlim=c(1,ngroups),
           ylim=c(yrange[1]-ydiff*0.05, yrange[2]+ydiff*0.05), cex.lab=1.2, ylab="", xlab="")
      suppressWarnings(
        arrows(1:ngroups, k95_means[[i]]+k95_se[[i]], 1:ngroups, k95_means[[i]]-k95_se[[i]],
             angle=90, code=3, length=0.05, col="black", lwd=0.5))
      axis(4, at=yvals, labels=ylabs, cex.axis=1.1, las=1)
      mtext(expression(paste("KUD 95% (km"^"2"*")")), side=4, cex=1.2, line=4)
      legend("top", legend=c("KUD 50%", "KUD 95%"), lty=c(1,2), pch=c(16,15), bty="n", horiz=T)
      box()
    }

    ############################################################################
    ## Plot a graph by KUD percentage (overlaying group series) ################
  } else if (plot.by=="kud"){
    par(mfrow=c(2,1))

    ## KUD 50% ##################################
    yrange <- k50_total_range
    ydiff <- yrange[2] - yrange[1]
    yvals <- pretty(c(yrange[1], yrange[2]))
    ylabs <- sprintf("%.1f", yvals)
    plot(0, axes=F, ylab="", xlab=label, xlim=c(1, ngroups),
         ylim=c(yrange[1]-ydiff*0.05, yrange[2]+ydiff*0.05), cex.lab=1.2)
    mtext(expression(paste("KUD 50% (km"^"2"*")")), side=2, cex=1.2, line=3.5)
    if(shade.season==T) {
      rect(xleft=c(1,3,6,9), xright=c(3,6,9,12), ybottom=par("usr")[3], ytop=par("usr")[4],
           col=c("white","grey88", "grey72", "grey80"), border=NA)
    }
    if(shade.season==F) {
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="grey96", border=NA)
    }
    axis(1, at=1:ngroups, labels=colnames(k50_table)[-1], cex.axis=1.1)
    axis(2, at=yvals, labels=ylabs, cex.axis=1.1, las=1)

    for(i in 1:length(id.groups)){
      points(k50_means[[i]], type="o", pch=15+i, lty=i, cex=1)
      suppressWarnings(
        arrows(1:ngroups, k50_means[[i]]+k50_se[[i]], 1:ngroups, k50_means[[i]]-k50_se[[i]],
                              angle=90, code=3, length=0.03, col="black", lwd=0.5))

    }
    legend("top", legend=names(id.groups), pch=15+(1:length(k50_means)), lty=1:length(k50_means),
           horiz=T, bty="n", cex=0.8)
    box()


    ## KUD 95% ##################################
    yrange <- k95_total_range
    ydiff <- yrange[2] - yrange[1]
    yvals <- pretty(c(yrange[1], yrange[2]))
    ylabs <- sprintf("%.1f", yvals)
    plot(0, axes=F, ylab="", xlab=label, xlim=c(1, ngroups),
         ylim=c(yrange[1]-ydiff*0.05, yrange[2]+ydiff*0.05), cex.lab=1.2)
    mtext(expression(paste("KUD 95% (km"^"2"*")")), side=2, cex=1.2, line=3.5)
    if(shade.season==T) {
      rect(xleft=c(1,3,6,9), xright=c(3,6,9,12), ybottom=par("usr")[3], ytop=par("usr")[4],
           col=c("white","grey88", "grey72", "grey80"), border=NA)
    }
    if(shade.season==F) {
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="grey96", border=NA)
    }
    axis(1, at=1:ngroups, labels=colnames(k50_table)[-1], cex.axis=1.1)
    axis(2, at=yvals, labels=ylabs, cex.axis=1.1, las=1)

    for(i in 1:length(id.groups)){
      points(k95_means[[i]], type="o", pch=15+i, lty=i, cex=1)
      suppressWarnings(
        arrows(1:ngroups, k95_means[[i]]+k95_se[[i]], 1:ngroups, k95_means[[i]]-k95_se[[i]],
               angle=90, code=3, length=0.03, col="black", lwd=0.5))

    }
    legend("top", legend=names(id.groups), pch=15+(1:length(k95_means)), lty=1:length(k95_means),
           horiz=T, bty="n", cex=0.8)
    box()
  }


  # reset par
  par(mar=c(5,4,4,2)+0.1)

}




#######################################################################################################
#######################################################################################################
#######################################################################################################
