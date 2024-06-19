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
#'
#' @param min.value Defaults to the 50% quantile
#' @param cut.value Defaults to the 90% quantile
#' @param title.cex  Size of the plot title. Defaults to 1.
#' @param color.nodes.by Variable used to color nodes (locations). It can be set to any variable in the dataset.
#' Alternatively, it can be set to 'detections' (default) to distinguish between nodes with and without detections,
#' or to 'group' to assign different colors to different ID groups (one color per map). Defaults to "detection".
#' @param nodes.color Color(s) for the nodes. Defaults to c("darkblue", "black").
#' @param nodes.alpha Numeric. Transparency level for the nodes (0 to 1): 0=fully transparent, 1=fully opaque. Defaults to 0.8.
#' @param nodes.size A numeric vector of length 2. Represents the desired min and max vertex sizes relative
#' to the x-axis in terms of percentage (see details). See \code{\link[netdiffuseR]{rescale_vertex_igraph}}. Defaults to c(0.04, 0.08).
#' @param nodes.label.wrap Logical. If TRUE, splits node labels into multiple lines (to better fit within the nodes). Defaults to FALSE.
#' @param nodes.label.cex The font size for vertex labels. Defaults to 0.5.
#' @param nodes.label.color The font color for vertex labels. Defaults to white.
#' @param repel.nodes Logical. If TRUE, nodes are plotted using a repulsion algorithm to avoid overlap. Defaults to FALSE.
#' @param repel.buffer Controls the amount of space between nodes if repel.nodes is set to TRUE. Defaults to 1.1.
#' @param edge.type A character string indicating the metric to be used to calculate the network edges.
#' It can be either "movements" to represent the number of movements between sites or "individuals"
#' to represent the number of individual animals moving between sites. Defaults to "movements".
#' @param edge.color Color(s) for the edges. Defaults to "darkblue".
#' @param edge.curved Specifies whether to draw curved edges, or not. This can be a logical or a numeric vector or scalar.
#' A numeric value specifies the curvature of the edge; zero curvature means straight edges, negative values mean the edge bends
#' clockwise, positive values the opposite. TRUE means curvature 0.5, FALSE means curvature zero. Defaults to 0.5.
#' @param edge.width A numeric vector of length 2. Represents the desired min and max
#' width/thickness of the edges. Defaults to c(0.4, 3.5).
#' @param edge.arrow.size The size of the arrows. Defaults to 0.5.
#' @param edge.arrow.width The width of the arrows. Defaults to 1.5.
#' @param edge.label.cex The font size for edge labels. Defaults to 0.6.
#' @param edge.label.color The color of the edge labels. Defaults to black.
#' @param edge.label.font The font for the edge labels. It is interpreted the same way
#' as the font graphical parameter: 1 is plain text, 2 is bold face, 3 is italic,
#' 4 is bold and italic and 5 specifies the symbol font. Defaults to 1.


#' @seealso \code{\link{calculateOverlap}}
#' @seealso \code{\link{randomizeOverlaps}}
#' @export


plotOverlapNetwork <- function(overlaps,
                               random.results=NULL,
                               min.val=NULL,
                               cut.val=NULL,
                               group.order=NULL,
                               title.cex = 1,
                               color.nodes.by = "detection",
                               network.layout="spring",
                               nodes.color = NULL,
                               nodes.alpha = 0.8,
                               nodes.size = c(0.7, 1.4),
                               nodes.scale.by= NULL,
                               nodes.label.wrap = FALSE,
                               nodes.label.cex = 0.5,
                               nodes.label.color = "white",
                               repel.nodes = FALSE,
                               repel.buffer = 1.1,
                               edge.color = "darkblue",
                               edge.curved = 0.5,
                               edge.width = c(0.4, 3.5),
                               edge.arrow.size = 0.5,
                               edge.arrow.width = 1.5,
                               edge.label.cex = 0.6,
                               edge.label.color = "black",
                               edge.label.font = 1,
                               ...){


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # validate parameters
  errors <- c()
  if (!is.data.frame(overlaps) || !c("ids") %in% names(attributes(overlaps))) errors <- c(errors, "'overlaps' format not recognized. Please make sure to use the output from the 'calculateOverlap' function.")
  if(!inherits(random.results, "list"))  errors <- c(errors, "random.results must be a list, as expected from the output of the 'randomizeOverlaps' function.")
  if(!c("metric") %in% names(attributes(random.results))) errors <- c(errors, "random.results format not recognized. Please make sure to use the output from the 'randomizeOverlaps' function.")
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # check if qgraph package is installed
  if (!requireNamespace("qgraph", quietly=TRUE)) {
    stop("The 'qgraph' package is required for this function but is not installed. Please install 'qgraph' using install.packages('qgraph') and try again", call.=FALSE)
  }

  # print to console
  printConsole("Generating overlap network")


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # split overlaps by type (if required)
  if(!any(colnames(overlaps)=="type")) overlaps$type <- "All"
  types <- unique(overlaps$type)

  # split overlaps by group
  if(is.null(random.results)){
    pairwise_overlaps <- split(overlaps, f=overlaps$type)
  }else{
    pairwise_overlaps <- split(random.results$pairwise_results, f=random.results$pairwise_results$type)
  }
  unique_ids <- lapply(pairwise_overlaps, function(x) unique(c(x$id1, x$id2)))
  pairwise_overlaps <- mapply(function(x,ids){x$id1<-factor(x$id1, levels=ids); return(x)}, x=pairwise_overlaps, ids=unique_ids, SIMPLIFY=FALSE)
  pairwise_overlaps <- mapply(function(x,ids){x$id2<-factor(x$id2, levels=ids); return(x)}, x=pairwise_overlaps, ids=unique_ids, SIMPLIFY=FALSE)
  network_matrices <- lapply(pairwise_overlaps, function(x) reshape2::dcast(x, formula="id1~id2", value.var="overlap", drop=F))
  network_matrices <- lapply(network_matrices, function(x) {rownames(x)<-x$id1; return(x[,-1])})
  network_matrices <- lapply(network_matrices, function(x) {x[is.na(x)] <- 0; return(x)})
  network_matrices <- lapply(network_matrices, as.matrix)
  if(!is.null(group.order)) network_matrices <- network_matrices[group.order]

  # fetch null distribution
  if(!is.null(random.results)){
    type_indices <- lapply(types, function(x) which(overlaps$type==x))
    iterations <- attributes(random.results)$iterations
    null_dist <- lapply(1:iterations, function(i) lapply(type_indexes, function(p) random.results$randomized_overlaps[[i]][p]))
  }


  # filter out null comparisons
  valid_overlaps <- overlaps[!is.na(overlaps$overlap),]

  # calculate average shared monitoring period
  shared_period <- aggregate(valid_overlaps$shared_monit_days, by=list(valid_overlaps$type), mean, na.rm=T)
  colnames(shared_period) <- c("type", "days")

  # get nº of dyads
  n_dyads <- aggregate(valid_overlaps$overlap, by=list(valid_overlaps$type), length)
  colnames(shared_period) <- c("type", "n_dyads")

  # get binary degree (number of edges each node has)
  binary_degree <- network_matrices
  binary_degree <- lapply(binary_degree, function(x) unlist(apply(x, 1, function(y) length(which(y>0)))))
  binary_degree <- unlist(lapply(binary_degree, mean))


  ##############################################################################
  ## Set layout variables ######################################################
  ##############################################################################

  # define cut and min values
  vals <- overlaps$overlap[!is.na(overlaps$overlap)]
  vals <- vals[order(vals, decreasing=T)]
  if(!is.null(min.val)) min.val <- quantile(vals, 0.5)
  if(!is.null(cut.val)) cut.val <- quantile(vals, 0.9)

  # set color palette for nodes
  if(is.null(nodes.color)){
    c("#643BA6", "darkorange2", "grey")
  }


  # set network properties for each matrix
  for(i in 1:length(network_matrices)){

    # set node properties
    group_ids <- colnames(network_matrices[[i]])
    group_metadata <- data_tags_kitefins[data_tags_kitefins$ID %in% group_ids,]
    node_sizes <- as.numeric(plyr::mapvalues(group_ids, group_metadata$ID, group_metadata$length, warn_missing=F))
    node_sizes <- .rescale(node_sizes, from=range(data_tags_kitefins$length), c(1, 2.8))
    node_types <- as.numeric(plyr::mapvalues(group_ids, group_metadata$ID, as.numeric(group_metadata$sex), warn_missing=F))
    node_colors <- color_pal[node_types]
    edge_labels <- round(network_matrices[[i]],1)

    # set edge properties (overlap)
    edge_colors <- network_matrices[[i]]
    edge_colors[edge_colors>=0] <- NA
    colnames(edge_colors) <- plyr::mapvalues(colnames(edge_colors), group_metadata$ID, as.character(group_metadata$sex), warn_missing=F)
    rownames(edge_colors) <- plyr::mapvalues(rownames(edge_colors), group_metadata$ID, as.character(group_metadata$sex), warn_missing=F)
    for(e in 1:length(edge_colors)){
      index <- arrayInd(e, dim(edge_colors))
      group1 <- rownames(edge_colors)[index[,1]]
      group2 <- colnames(edge_colors)[index[,2]]
      if(group1==group2 & group1=="female"){edge_colors[e] <- color_pal[1]
      }else if(group1==group2 & group1=="male"){edge_colors[e] <- color_pal[2]
      }else{edge_colors[e] <- "black"}
    }

  }


  ##############################################################################
  ## Generate figure ###########################################################
  ##############################################################################

  # iterate over each matrix
  for(i in 1:length(network_matrices)){

    # overlap network(s) #######################################################
    qgraph::qgraph(network_matrices[[i]],
                   layout=network.layout,
                   cut = cut.val,
                   minimum = min.val,
                   edge.labels = edge_labels,
                   edge.color = edge.color,
                   edge.width = edge.width,
                   edge.label.cex = edge.label.cex,
                   edge.label.color = edge.label.color,
                   edge.label.font = edge.label.font,
                   edge.label.margin=0.005,
                   edge.arrow.size = edge.arrow.size,
                   edge.arrow.width = edge.arrow.width,
                   colFactor = 0.6,
                   label.color = nodes.label.color,
                   label.cex = nodes.label.cex,
                   node.width = node_sizes,
                   color = nodes.color,
                   curveAll = TRUE,
                   curveDefault = edge.curved,
                   curveShape = -1,
                   curveScale = TRUE,
                   curveScaleNodeCorrection = TRUE,
                   repulsion = 0.1,
                   trans = TRUE,
                   fade = TRUE,
                   usePCH = TRUE)

    # add title plus stats
    title(main=names(network_matrices)[i], cex.main=title.cex, line=4.6, font=2, xpd=T)
    network_metrics1 <- paste0("N\u00ba of individuals: ", length(group_ids))
    network_metrics2 <- paste0("N\u00ba of dyads: ", ndyads$ndyads[i])
    network_metrics3 <- paste0("Mean binary degree: ", sprintf("%.1f", binary_degree[i]))
    network_metrics4 <- paste0("Mean shared period: ", sprintf("%.0f", shared_period$days[i]), " days")
    legend("bottom", inset=c(0,-0.08), legend=c(network_metrics1, network_metrics2, network_metrics3, network_metrics4), bty="n", cex=1, xpd=T)
  }

  # plot random overlaps distribution #########################################
  if(!is.null(random.results)){
    par(mar=c(5,5,5,5))
    for(i in 1:nlevels(pairwise_table$type)){

      obs_overlap <- mean(random.results$pairwise_results, na.rm=T)
      type_indexes <- lapply(types, function(x) which(pairwise_stats$type==x))
      null_dist <- lapply(1:iterations, function(i) lapply(type_indexes, function(p) random_results[[i]][p]))
      null_dist <- unlist(lapply(null_dist, mean, na.rm=T))



      randomized_overlaps <-

      group_overlaps <- lapply(randomized_overlaps$randomized_overlaps, function(x) x[type_indices[[i]]])
      group_overlaps <- unlist(lapply(group_overlaps, mean, na.rm=T))
      group_overlaps <- group_overlaps[order(group_overlaps)]
      h <- hist(group_overlaps, breaks=60, plot=F)
      upper_freq <- max(h$counts)*1.1
      hist(group_overlaps, breaks=60, col=NA, border=F, main="", xlab="Overlap (%)", ylab="Frequency",
           cex.lab=1.2, yaxs="i", axes=F, ylim=c(0, upper_freq))
      rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col="grey96", border=NULL)
      hist(group_overlaps, breaks=60, col=adjustcolor(color_pal[i], alpha.f=0.7), lwd=0.05, main="", xlab="", ylab="", axes=F, las=1, yaxs="i", ylim=c(0, upper_freq), add=T)
      axis(1, at=grid::grid.pretty(range(h$breaks)), labels=sprintf("%.2f", grid::grid.pretty(range(h$breaks))), pos=0, las=1, cex.axis=1.2)
      axis(2, at=grid::grid.pretty(c(0, upper_freq)), labels=grid::grid.pretty(c(0, upper_freq)), las=1, cex.axis=1.2)
      xrange <- par('usr')[2] - par('usr')[1]
      if(observed_overlap<max(group_overlaps) & observed_overlap>min(group_overlaps)){
        line_val <- observed_overlap
      }else if(observed_overlap<max(group_overlaps)){
        line_val <- par('usr')[1]+xrange*0.01
      }else if(observed_overlap>min(group_overlaps)){
        line_val <- par('usr')[2]-xrange*0.01
      }
      abline(v=line_val, col="red2", lwd=3)
      if(pval>=0.05){pval_text <- paste("p =", sprintf("%.3f", pval))}
      if(pval<0.05){pval_text <- "p < 0.05"}
      if(pval<0.01){pval_text <- "p < 0.01"}
      if(pval<0.001){pval_text <- "p < 0.001"}
      # format legend info
      info1 <- paste0("Simulations = ", length(group_overlaps))
      info2 <- paste0("Obs. overlap = ", sprintf("%.2f", observed_overlap), "%")
      legend("topleft", legend=c(info1, info2, pval_text), bty="n", y.intersp=1.2, cex=1.2)
      box()
    }

  }






}




#######################################################################################################
#######################################################################################################
#######################################################################################################
