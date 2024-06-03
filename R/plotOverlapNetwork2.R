#######################################################################################################
# Plot overlap network  ###############################################################################
#######################################################################################################

#' Network representation of pairwise overlaps.

#' @description Function to plot a network of estimated pairwise overlaps.
#' Each node represents a single individual and can be dimensioned proportionally
#' to its size. Edges represent the extent of spatial overlap between each pair
#' (i.e. greater thickness indicates a higher overlap). If animals sizes are supplied,
#' two additional plots are produced, one representing the relationship between
#' animals' size/length and average overlap and the other representing the
#' relationship between pairwise differences in size and overlap scores.
#'
#' @param overlaps Similarity matrix containing pairwise overlaps, as returned by \code{\link{calculateOverlap}}.
#' @param id.metadata A data frame containing animal IDs as well as an optional 'color.by' and 'size.by' variables.
#' @param id.col Name of the column containing animal identifications. Defaults to "ID".
#' @param color.by Optional. Factor variable in id.metadata defining the color of each node (e.g. species).
#' @param scale.by Optional. Numeric variable in id.metadata defining the size of each node (e.g. animal length).
#' @param model Type of regression in case animal.sizes have been supplied ("linear" or "quadratic").
#' A vector with two elements should be supplied, the first to be used in the "average overlap ~ size" graph,
#' and the second in the "pairwise overlap ~ size difference" graph.
#' @param min.val Edges with overlaps lower than this value are not shown.
#' @param cut.val Edges with overlaps larger than this value will have the strongest color intensity
#' and become wider the stronger they are, and edges with overlaps lower than this value
#' will have the smallest width and become vaguer the weaker the weight.
#' @param color.pal Color of the nodes.
#' @param ... Further arguments passed to \code{\link[qgraph]{qgraph}}.
#' @seealso \code{\link{calculateOverlap}}
#' @seealso \code{\link{plotOverlapFrequences}}
#' @export


plotOverlapNetwork <- function(overlaps, id.metadata, id.col="ID", color.by=NULL, scale.by=NULL, min.size=0.7, max.size=1.4,
                               model=c("linear", "linear"), min.val=NULL, cut.val=NULL, color.pal=NULL, ...) {

  #######################################################################################################
  # Initial checks ###############################################################################


  # check if qgraph package is installed
  if (!requireNamespace("qgraph", quietly=TRUE)) {
    stop("The 'qgraph' package is required for this function but is not installed. Please install 'qgraph' using install.packages('qgraph') and try again.")
  }


  # check if id.metadata contains id.col
  if(!id.col %in% colnames(id.metadata)){
    stop("ID column not found in id.metadata. Please assign the correct column name with 'id.col'")
  }

  # check if id.metadata contains id.col
  if(!is.null(color.by) & !color.by %in% colnames(id.metadata)){
    stop("color.by column not found in id.metadata")
  }

  # check id.col format
  if(class(id.metadata[,id.col])!="factor"){
    cat("Converting metadata ids to factor\n")
    id.metadata[,id.col] <- as.factor(id.metadata[,id.col])
  }else{
    id.metadata[,id.col] <- droplevels(id.metadata[,id.col])
  }

  # check if id.metadata contains color.by
  if(!is.null(color.by)){
    if(!color.by %in% colnames(id.metadata)){
      stop("color.by column not found in id.metadata'")
    }
    if(class(id.metadata[,color.by])!="factor"){
      cat("Converting color.by column to factor\n")
      id.metadata[,color.by] <- as.factor(id.metadata[,color.by])
    }
  }

  # check if id.metadata contains scale.by
  if(!is.null(scale.by)){
    if(!scale.by %in% colnames(id.metadata)){
      stop("scale.by column not found in id.metadata'")
    }
    if(!class(id.metadata[,scale.by]) %in% c("numeric", "integer")){
      stop("scale.by should be of class numeric'")
    }
  }

  if(!is.null(color.pal)){
    if(length(color.pal)!=nlevels(id.metadata[,color.by])){
      cat("Warning: Number of colors in color.pal don't match the number of levels in color.by\n")
    }
  }else{
    if(!is.null(color.by)){
      color.pal <- grey.colors(nlevels(id.metadata[,color.by]))
    }else{
      color.pal <- "gray60"
    }
  }

  if(any(!model %in% c("linear", "quadratic"))){
    stop("Wrong 'model' argument. Please choose between 'linear' or 'quadratic'.")
  }

  cat("Generating overlap network\n")


  #######################################################################################################
  # Create overlap matrix ###############################################################################

  if(any(class(overlaps)=="list")){overlaps <- overlaps$overlap}

  overlaps <- as.matrix(overlaps)
  network_matrix <- round(overlaps, 1)
  network_matrix[is.na(network_matrix)] <- 0

  network_ids <- colnames(network_matrix)

  if(length(network_ids)>nlevels(id.metadata[,id.col])){
    stop("IDs appear to be missing from id.metadata")
  }else if(length(network_ids)<nlevels(id.metadata[,id.col])){
    id.metadata <- id.metadata[id.metadata[,id.col] %in% network_ids,]
    id.metadata[,id.col] <- droplevels(id.metadata[,id.col])
  }


  #######################################################################################################
  # Plot networks ###############################################################################

  # set nodes' sizes
  if(!is.null(scale.by)) {
    layout(mat=matrix(c(1,1,1,1,2,3), byrow=T, nrow=3, ncol=2))
    if(length(discard_ids)>0){animal.sizes <- animal.sizes[-discard_ids]}
    node_sizes <- moby:::rescale(id.metadata[,scale.by], c(min.size, max.size))
  } else {
    node_sizes <- 1
  }

  # set nodes' color
  if(!is.null(color.by)) {
   node_colors <- color.pal[id.metadata[,color.by]]
  }else{
   node_colors <- color.pal[1]
  }

  vals <- overlaps[!is.na(overlaps)]
  vals <- vals[order(vals, decreasing=T)]
  if(is.null(min.val)){min.val <- quantile(vals, 0.5)}
  if(is.null(cut.val)){cut.val <-  quantile(vals, 0.9)}

  # overlap network
  par(mar=c(1,1,1,1))
  qgraph::qgraph(network_matrix, edge.labels=T, color=node_colors, edge.color="black", edge.width=0.5,
         cut=cut.val, colFactor=0.5, minimum=min.val, label.color="white", node.width=node_sizes,
         edge.label.cex=1.3, edge.label.margin=0.01, label.cex=0.8, label.prop=1, ...)

  if(!is.null(color.by)){
    legend("topleft", legend=levels(id.metadata[,color.by]), pch=21, pt.bg=color.pal, pt.lwd=0.4, pt.cex=2.4,
             bty="n", cex=1.2, inset=c(0.1, 0.05))
  }


  if(!is.null(animal.sizes)) {
    # calculate size-related metrics
    size_diffs <- as.numeric(dist(animal.sizes))
    overlap_table <- reshape2::melt(overlaps)
    colnames(overlap_table) <- c("ID_1", "ID_2", "overlap")
    mean_overlaps <- aggregate(overlap_table$overlap, by=list(overlap_table$ID_1), mean, na.rm=T)
    colnames(mean_overlaps) <- c("ID", "mean_overlap")
    mean_overlaps$size <- animal.sizes
    par(mar=c(5,7,0,2), mgp=c(3,0.8,0))
    if(length(model==1)){model <- c(model, model)}
    model1 <- model[1]
    model2 <- model[2]

    # plot mean overlap ~ individual size
    yrange <- range(mean_overlaps$mean_overlap)*1.1
    plot(x=mean_overlaps$size, y=mean_overlaps$mean_overlap, type="n", ylim=yrange, xlab="", ylab="", axes=F)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NULL)
    points(x=mean_overlaps$size, y=mean_overlaps$mean_overlap, pch=16, cex=0.8)
    axis(1, at=pretty(mean_overlaps$size), labels=pretty(mean_overlaps$size))
    axis(2, at=pretty(mean_overlaps$mean_overlap), labels=sprintf(pretty(mean_overlaps$mean_overlap), fmt="%.1f"), las=2)
    title(xlab="Length", cex.lab=1.3, line=2.5)
    title(ylab="Mean overlap (%)", cex.lab=1.3, line=3)
    if(model1=="linear"){
      abline(lm(mean_overlaps$mean_overlap ~ mean_overlaps$size), col="red")
      correlation <-  cor.test(mean_overlaps$mean_overlap, mean_overlaps$size)
      legend("topright", c(paste("r =", sprintf("%.2f",correlation$estimate)),
                           paste("p =", sprintf("%.3f",correlation$p.value))), bty="n", cex=1.1)
    }
    if(model1=="quadratic"){
      res <- nls(mean_overlap ~ k*exp(-1/2*(size-mu)^2/sigma^2), start=c(mu=100,sigma=10, k=1) , data=mean_overlaps)
      v <- summary(res)$parameters[,"Estimate"]
      plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2),col=2,add=T, xlim=range(mean_overlaps$size))
    }
      legend("topleft", legend="B", bty="n", cex=1.8, inset=c(-0.08,-0.04))
    box()

    # plot overlaps ~ size difference
    total_overlaps <- as.numeric(proxy::as.simil(overlaps))
    yrange <- range(total_overlaps, na.rm=T)*1.1
    plot(x=size_diffs, y=total_overlaps, type="n", ylim=yrange, ylab="", xlab="", axes=F)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NULL)
    points(x=size_diffs, y=total_overlaps, pch=16, cex=0.8)
    axis(1, at=pretty(size_diffs), labels=pretty(size_diffs))
    axis(2, at=pretty(total_overlaps), labels=sprintf(pretty(total_overlaps), fmt="%.1f"), las=2)
    title(xlab="Size diff.", cex.lab=1.3, line=2.5)
    title(ylab="Overlap (%)", cex.lab=1.3, line=3)
    if(model2=="linear"){
      abline(lm(total_overlaps ~ size_diffs), col="red")
      correlation <-  cor.test(total_overlaps, size_diffs)
      legend("topright", c(paste("r =", sprintf("%.2f",correlation$estimate)),
                           paste("p =", sprintf("%.3f",correlation$p.value))), bty="n", cex=1.1)
    }
    if(model2=="quadratic"){
      res <- nls(total_overlaps ~ k*exp(-1/2*(size_diffs-mu)^2/sigma^2), start=c(mu=100,sigma=10, k=1))
      v <- summary(res)$parameters[,"Estimate"]
      plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2),col=2,add=T, xlim=range(mean_overlaps$size))
    }
    legend("topleft", legend="C", bty="n", cex=1.8, inset=c(-0.08,-0.04))
    box()
  }
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
