#######################################################################################################
## Plot summary table for pairwise overlap results ####################################################
#######################################################################################################

#' Plots pairwise stats in a contingency table

#' @description Plots a contingency table illustrating the pairwise distribution of overlap
#' scores/significance across individuals. IDs can be sorted by length or any other quantitative metric.
#'
#' @param randomized.overlaps Similarity matrix containing pairwise overlaps, as returned by \code{\link{calculateOverlap}}.
#' @param id.metadata Optional. A data frame containing individual's metadata. It should contain
#' at least animal IDs  and a numeric variable (e.g. length), defined by
#' @param id.groups Optional. A list containing ID groups, used to calculate stats independently
#' within each group and/or to compare relationships between ids of different groups.
#' @param type Type of metric to be plotted, one of "significance" or "mean overlap".
#' Defaults to "significance".
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' One of "within", "between" or "both". Useful to visualize comparisons between individuals
#' belonging to the same id.group or to only compare individuals in different id.groups
#' (e.g. intraspecific vs interspecific). Defaults to "both".
#' @param color.pal Color palette. If type is set to significance, it should contain 3 colors
#' ("+", "-", "ns").
#' @param discard.empty Boolean. Discard IDs with missing interactions.
#' @param id.col Name of the column in id.metadata containing animal IDs. Defaults to 'ID'.
#' @param full.scale Boolean. If type is set to 'mean overlap', sets the overlap scale from 0 to 100%.
#' @param sort.by Optional. Sort IDs based on this variable (e.g., length).
#' @seealso \code{\link{calculateOverlap}} and \code{\link{randomizeOverlaps}}
#' @export


################################################################################
# Main function - plot frequencies of co-occurring group sizes #################
################################################################################

plotOverlapTable <- function(randomized.overlaps,
                             id.metadata = NULL,
                             id.groups = NULL,
                             type = "significance",
                             group.comparisons = "both",
                             color.pal = NULL,
                             discard.empty = TRUE,
                             id.col = getDefaults("id"),
                             full.scale = FALSE,
                             sort.by = NULL) {


  ##############################################################################
  # initial checks #############################################################

  if(!type %in% c("mean overlap", "significance")){
    stop("Wrong type argument, please select either 'mean overlap', or 'significance'")
  }

  if(!is.null(sort.by) & is.null(id.metadata)){
    stop(paste("Please supply an id.metadata data frame containing a", sort.by, "column"))
  }

  if(!is.null(sort.by) & !is.null(id.metadata)){
    if(!sort.by %in% colnames(id.metadata)){
      stop("sort.by variable not found within the supplied id.metadata")
    }
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){

    # check for duplicated IDs
    if(any(duplicated(unlist(id.groups)))) {
      stop("Repeated ID(s) in id.groups")
    }

    # check if any ID is missing from the supplied metadata
    if(any(!unlist(id.groups) %in% levels(id.metadata[,id.col]))){
      stop("Some of the ID(s) in id.groups were not found in the supplied metadata")
    }

    if(!group.comparisons %in% c("both", "within", "between")){
      stop("Wrong group.comparisons argument, please select one of: 'within', 'between' or 'both'")
    }
  }

  # check color palette
  if(!is.null(color.pal)){
    if(type=="significance"){
      if(length(color.pal)!=3){stop("Color palette should contain 3 colors when type is set to 'significance'")}
    }else if (type=="mean overlap"){
      if(length(color.pal)<10){warning("Consider increasing the number of colors in color.pal")}
    }
  }else{
    if(type=="significance"){
      color.pal <- c("lightseagreen", "darkred", "grey85")
      color.pal <- adjustcolor(color.pal, alpha.f=0.3)
    }else if (type=="mean overlap"){
      color.pal <- .viridis_pal(100)
      color.pal <- c(adjustcolor(color.pal[1], alpha.f=0.65), adjustcolor(color.pal[2:100], alpha.f=0.7))
    }
  }

  ##############################################################################
  # create contigency table ####################################################

  pairwise_results <- randomized.overlaps$pairwise_significance

  if(!is.null(id.groups)){
    ids <- unique(unlist(id.groups))
  }else{
    ids <- pairwise_results$ids
    ids <- strsplit(ids, split="-", fixed=TRUE)
    ids <- unique(unlist(ids))
  }

  nids <- length(ids)
  pairwise_results$id1 <- unlist(lapply(strsplit(pairwise_results$ids, split="-", fixed=TRUE), function(x) x[[1]]))
  pairwise_results$id2 <- unlist(lapply(strsplit(pairwise_results$ids, split="-", fixed=TRUE), function(x) x[[2]]))

  # discard comparisons between individuals belonging to the same group or
  # skip comparisons between different groups, if required
  if(!is.null(id.groups)){
    discard_rows <- c()
    for(r in 1:nrow(pairwise_results)){
      group1 <- which(unlist(lapply(id.groups, function(x) pairwise_results$id1[r] %in% x)))
      group2 <- which(unlist(lapply(id.groups, function(x) pairwise_results$id2[r] %in% x)))
      if(group.comparisons=="within" & group1!=group2){discard_rows<-c(discard_rows, r)}
      if(group.comparisons=="between" & group1==group2){discard_rows<-c(discard_rows, r)}
    }
  }

  if(length(discard_rows)>0){
    pairwise_results <- pairwise_results[-discard_rows,]
  }

  # create contingency table
  if(type=="mean overlap"){
    pairwise_results <- reshape2::dcast(pairwise_results, id1 ~ id2, value.var="overlap")
    rownames(pairwise_results) <- pairwise_results$id1
    pairwise_results <- pairwise_results[,-1]
    if(discard.empty){
      pairwise_results <- pairwise_results[rowSums(is.na(pairwise_results)) != ncol(pairwise_results), ]
      pairwise_results <- pairwise_results[,colSums(is.na(pairwise_results)) != nrow(pairwise_results)]
    }
    overlap_range <- range(pairwise_results, na.rm=TRUE)
  }else if(type=="significance"){
    pairwise_results <- reshape2::dcast(pairwise_results, id1 ~ id2, value.var="p.val")
    rownames(pairwise_results) <- pairwise_results$id1
    pairwise_results <- pairwise_results[,-1]
    if(discard.empty){
      pairwise_results <- pairwise_results[rowSums(is.na(pairwise_results)) != ncol(pairwise_results), ]
      pairwise_results <- pairwise_results[,colSums(is.na(pairwise_results)) != nrow(pairwise_results)]
    }
  }

  ids1 <- rownames(pairwise_results)
  ids2 <- colnames(pairwise_results)
  pairwise_results <- cbind(ids1, pairwise_results)
  pairwise_results <- rbind(c(NA, ids2), pairwise_results)
  rownames(pairwise_results) <- NULL
  colnames(pairwise_results) <- NULL
  pairwise_results[1,1] <- "IDs"

  # sort IDs if required
  if(!is.null(sort.by) & !is.null(id.metadata)){
    vars1 <- unlist(sapply(ids1, function(x) id.metadata[,sort.by][id.metadata[,id.col]==x], simplify=FALSE))
    vars2 <- unlist(sapply(ids2, function(x) id.metadata[,sort.by][id.metadata[,id.col]==x], simplify=FALSE))
    pairwise_results <- cbind(pairwise_results, c(NA, vars1))
    pairwise_results <- rbind(pairwise_results, c(NA, vars2, NA))
    cols <- ncol(pairwise_results)
    rows <- nrow(pairwise_results)
    pairwise_results <- pairwise_results[,c(1, cols, 2:(cols-1))]
    pairwise_results <- pairwise_results[c(1, rows, 2:(rows-1)),]
    pairwise_results[-c(1:2),] <- pairwise_results[-c(1:2),][order(vars1),]
    pairwise_results[,-c(1:2)] <- pairwise_results[,-c(1:2)][,order(vars2)]
    pairwise_results[2,2] <- sort.by
    rownames(pairwise_results) <- NULL
    colnames(pairwise_results) <- NULL
    digits <- max(.decimalPlaces(c(vars1, vars2)), na.rm=TRUE)
    pairwise_results[-c(1:2), 2] <- sprintf(paste0("%.", digits, "f"), as.numeric(unlist(pairwise_results[-c(1:2), 2])))
    pairwise_results[2, -c(1:2)] <- sprintf(paste0("%.", digits, "f"), as.numeric(unlist(pairwise_results[2, -c(1:2)])))
  }

  ##############################################################################
  # plot table #################################################################

  n_cols <- ncol(pairwise_results)
  n_rows <- nrow(pairwise_results)

  par(mar=c(3,1,3,7))
  plot(0, 0, xlim=c(0.5, n_cols+0.5), ylim=c(n_rows+0.5, 0.5),
       type="n", axes=FALSE, main="", ylab="", xlab="")
  values <- expand.grid(1:n_rows, 1:n_cols)
  colnames(values) <- c("row", "col")
  values$value <- unlist(apply(values, 1, function(x) pairwise_results[x[1], x[2]]))
  values$font <- 1
  values$font[values$row==1 | values$col==1] <- 2
  values$color <- NA
  if(type=="significance"){
    values$color[values$value=="+"] <- color.pal[1]
    values$color[values$value=="-"] <- color.pal[2]
    values$color[values$value=="ns"] <- color.pal[3]
    values$color[is.na(values$value) & values$row>1 & values$col>1] <- color.pal[3]
  }else if (type=="mean overlap"){
    if(!is.null(sort.by)){overlap_cells <- which(values$row>2 & values$col>2)
    }else{overlap_cells <- which(values$row>1 & values$col>1)}
    overlaps <- as.numeric(values$value[overlap_cells])
    if(full.scale){
      overlaps <- round(.rescale(overlaps, from=c(0,100), to=c(1, length(color.pal)-1)))
    }else{
      overlaps <- round(.rescale(overlaps, from=range(overlaps, na.rm=TRUE), to=c(1, length(color.pal)-1)))
    }
    overlaps[is.na(overlaps)] <- 0
    values$color[overlap_cells] <- color.pal[overlaps+1]
    values$value[overlap_cells] <- sprintf("%.1f", as.numeric(values$value[overlap_cells]))
    values$value[values$value=="NA"] <- ""
  }
  rect(xleft=values$col-0.5 , xright=values$col+0.5 , ybottom=values$row-0.5,
       ytop=values$row+0.5, col=values$color, border=NA)
  text(x=values$col, y=values$row, labels=values$value, font=values$font, cex=0.7, xpd=TRUE, adj=0.5)
  if(!is.null(sort.by)){
    segments(x0=2.5, x1=ncol(pairwise_results)+0.5, y0=2.5)
    segments(y0=2.5, y1=nrow(pairwise_results)+0.5, x0=2.5)
  }else{
    segments(x0=1.5, x1=ncol(pairwise_results)+0.5, y0=1.5)
    segments(y0=1.5, y1=nrow(pairwise_results)+0.5, x0=1.5)
  }

  # add bounding box
  segments(x0=0.5, x1=n_cols+0.5, y0=0.5, xpd=TRUE)
  segments(x0=0.5, x1=n_cols+0.5, y0=n_rows+0.5, xpd=TRUE)
  segments(y0=0.5, y1=n_rows+0.5, x0=0.5, xpd=TRUE)
  segments(y0=0.5, y1=n_rows+0.5, x0=n_cols+0.5, xpd=TRUE)

  # add legend
  if(type=="significance"){
    plus <- "overlap > random\n(joint resource utilization)"
    minus <- "overlap < random\n(spatiotemporal segregation)"
    ns <- "non-significant overlap"
    .legend("right", inset=c(-0.3, 0), legend=c(plus, minus, ns), box.lwd=0.1,
                       fill=color.pal, box.cex=c(1.4, 0.9), bty="n", cex=0.6, y.intersp=1.4, xpd=TRUE)
  }else if(type=="mean overlap"){
    if(full.scale) {
      overlap_range <- c(0,100)
      overlap_labs <- pretty(0:100, min.n=4)
    }else{
      overlap_labs <- pretty(overlap_range, min.n=4)
      overlap_labs <- overlap_labs[overlap_labs>=min(overlap_range) & overlap_labs<=max(overlap_range)]
    }
    overlap_digits <- max(.decimalPlaces(overlap_labs), na.rm=TRUE)
    .colorlegend(col=color.pal, zlim=overlap_range, zval=overlap_labs,
                            posx=c(0.88, 0.895), posy=c(0.2, 0.65), main="Overlap (%)",
                            main.cex=0.8, digit=overlap_digits, main.adj=0, cex=0.8)
  }

}



##################################################################################################
##################################################################################################
##################################################################################################
