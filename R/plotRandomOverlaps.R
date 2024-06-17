#######################################################################################################
# Plot overlap randomization results ##################################################################
#######################################################################################################

#'  Plot overlap randomization results

#' @description Plots frequency distributions of mean/median overlaps
#' obtained through the Monte Carlo permutation tests. Null hypotheses are tested
#' by comparing the observed mean/median overlaps (Obs. overlap, indicated by an arrow)
#' to the upper and lower tails of the null distributions, using a chosen confidence level.
#'
#' @param overlaps Similarity matrix containing pairwise overlaps, as returned by \code{\link{calculateOverlap}}.
#' @param random.results Overlap simulations, as returned by \code{\link{randomizeOverlaps}}.
#' @param include.table Boolean to indicate if a summary table should be drawn below the plot.
#' @param type Type of statistic used to plot the histogram and calculate overlap significance.
#' One of "mean" or "median". Defaults to "mean".
#' @param conf.level Confidence level. Defaults to 95%.
#' @export


plotRandomOverlaps <- function(overlaps, random.results, include.table=F, type="mean", conf.level=0.95) {


  ######################################################################
  # Calculate stats ####################################################

  # extract overlap matrix (if not directly supplied)
  if(class(overlaps)=="list"){overlaps<-overlaps$overlap}

  # calculate confidence intervals
  conf_intervals <- (1-conf.level)/2

  # calculate means/medians
  simulated_results <- random.results$randomized_overlaps
  if(type=="mean") {
    random_overlap <- unlist(lapply(simulated_results, function(x) mean(x, na.rm=T)))
    observed_overlap <- mean(overlaps, na.rm=T)
  } else if (type=="median") {
    random_overlap <- unlist(lapply(simulated_results, function(x) median(x, na.rm=T)))
    observed_overlap <- median(overlaps, na.rm=T)
  }

  bounds <- quantile(random_overlap, probs=c(conf_intervals, 1-conf_intervals))
  pval <- ifelse (observed_overlap %between% bounds, "Non-significant", "p < 0.05")

  if(include.table==T){
    pairwise_sign <- reshape2::melt(table(random.results$pairwise_significance$p.val))
    all_sign <- data.frame("Var1"=c("+", "ns","-"))
    pairwise_sign <- plyr::join(all_sign, pairwise_sign, by="Var1", type="left")
    pairwise_sign[is.na(pairwise_sign)] <- 0
    overlap_error <- sprintf("%.2f", plotrix::std.error(as.numeric(overlaps)))
    disp_overlap <- paste(sprintf("%.2f", mean(overlaps, na.rm=T)), "±", overlap_error)
    sign_table <- data.frame("Var1"="Overlap", "value"=disp_overlap)
    sign_table <- rbind(sign_table, pairwise_sign)
    sign_table$Var1 <- c("Overlap (%)", "No pairs with overlap > random",
                         "No pairs with non-signif. overlap",
                         "No pairs with overlap < random")
    colnames(sign_table) <- NULL
  }

  ######################################################################
  # Generate plot ######################################################

  xmin <- min(c(random_overlap, observed_overlap))
  xmax <- max(c(random_overlap, observed_overlap))
  xlabs <- pretty(c(xmin, xmax))
  digits <- max(moby:::decimalPlaces(xlabs))
  max_freq <- max(as.numeric(hist(random_overlap, breaks=25, plot=F)$counts))
  if(include.table==T){layout(matrix(c(1,1,1,1,2), ncol=1))}

  par(mar=c(4,5,4,2), lwd=0.6, mgp=c(3,0.8,0))
  plot_title <- paste0("Null Model (", type, ")")
  hist(random_overlap, xlim=c(xmin, xmax), axes=F, xlab="", ylab="", main=plot_title,
       breaks=25, ylim=c(0,max_freq), col=NA, border=NA, cex.main=1.4)
  rect(par("usr")[1], 0, par("usr")[2], par("usr")[4], col="gray96", border=NULL)
  hist(random_overlap, xlim=c(xmin, xmax), axes=F, xlab="", ylab="", main="",
       breaks=25, ylim=c(0,max_freq), col="gray75", add=T)
  par(lwd=1)
  rect(par("usr")[1], 0, par("usr")[2], par("usr")[4], col=NA, border=T, lwd=2)
  title(xlab="Overlap (%)", cex.lab=1.4, line=2)
  title(ylab="Frequency", cex.lab=1.4, line=3.5)
  arrows(x0=observed_overlap, y0=max_freq/6, y1=0, length=0.10, lwd=1.6)
  axis(1, at=xlabs, labels=sprintf(paste0("%.",digits,"f"), xlabs), pos=0, cex.axis=1.2)
  axis(2, at=pretty(c(0,max_freq)), labels=pretty(c(0, max_freq)), las=1, cex.axis=1.2)
  permut_info <- paste("Simulations =", length(random_overlap))
  real_stats <- paste0("Obs. overlap = ", round(observed_overlap,2), "%")
  random_stats <- paste0("Sim. overlap = ", round(mean(random_overlap),2), "%")
  legend("topright", legend = c(permut_info, real_stats, random_stats, pval), bty="n",
         cex=1.1, y.intersp=1.1, inset=c(0.025, 0.015))

  if(include.table==T){
    par(mar=c(3,1,3,1))
    plot.new()
    plotrix::addtable2plot("center", table=sign_table, bty='o', hlines=F, vlines=F, lwd=1,
                  display.rownames=F, display.colnames=F, ypad=1.5, xpad=0.5, bg="gray96", cex=1.2)
  }

}
##################################################################################################
##################################################################################################
##################################################################################################
