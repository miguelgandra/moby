#######################################################################################################
# Plot overlap network  ###############################################################################
#######################################################################################################

#' Network representation of pairwise overlaps.

#' @description This function plots a network of estimated pairwise overlaps between
#' individuals or a histogram showing the distribution of overlaps from null model permutations,
#' or both. If only `overlaps` are provided, it plots a network where each node represents
#' an individual, with edges indicating the extent of spatial overlap between pairs.
#' Thicker edges represent higher overlap. If only `random.results` are provided, it plots
#' a histogram showing the distribution of overlaps from the null model permutations.
#' This helps to identify whether the observed overlaps are significantly different from
#' what would be expected by chance. If both are provided, it plots both the network
#' and the histogram.

#' @param overlaps A data frame containing the pairwise overlap data. Should be the output from the `calculateOverlap` function.
#' @param random.results A list containing null model permutation results. Should be the output from the `randomizeOverlaps` function.
#' @param color.nodes.by Variable used to color nodes (individuals). It should be a categorical vector (factor or character).
#' Alternatively, it can be set to 'group' (default) to assign different colors to different ID groups.
#' @param scale.nodes.by Numeric vector to scale the nodes. Length should match the number of IDs in overlaps.
#' @param remove.missing Logical. If TRUE, excludes individuals without detections and null pairwise comparisons from the network.
#' Defaults to FALSE.
#' @param min.val Minimum value for overlap threshold. Defaults to the 50% quantile of the overlap values.
#' @param cut.val Cutoff value for overlap threshold. Defaults to the 90% quantile of the overlap values.
#' @param group.order Order of groups for plotting. Defaults to NULL.
#' @param plot.stats Logical. If TRUE, additional network metrics are plotted
#' below each network. Defaults to TRUE.
#' @param cex.title Size of the plot title. Defaults to 1.4.
#' @param cex.lab Size of the axis labels. Defaults to 1.2.
#' @param cex.axis Size of the axis text. Defaults to 1.0.
#' @param cex.legend Size of the legend text. Defaults to 1.0.
#' @param legend.inset Optional. Specifies how far the legend is inset from the
#' null model histogram margins. If a single value is given, it is used for both margins;
#' if two values are given, the first is used for x- distance, the second for y-distance.
#' Defaults to c(-0.02, 0).
#' @param network.layout Layout algorithm for the network (see \code{\link[qgraph]{qgraph}}).
#' Defaults to "spring".
#' @param nodes.color Color(s) for the nodes.
#' @param nodes.size Numeric vector of length 2, indicating the minimum and maximum
#' node sizes. Defaults to c(1, 2.2).
#' @param nodes.label.cex Font size for node labels. Defaults to 0.8.
#' @param nodes.label.scale Logical. If TRUE, scales node labels. Defaults to FALSE.
#' @param nodes.label.color Color of the node labels. Defaults to white.
#' @param edge.color Color(s) for the edges. Defaults to NULL, which averages node colors.
#' @param edge.curved Numeric value or logical indicating the curvature of the edges. Defaults to 0.5.
#' @param edge.width Numeric vector of length 2, indicating the minimum and maximum edge widths. Defaults to c(0.4, 3.5).
#' @param edge.label.cex Font size for edge labels. Defaults to 1.
#' @param edge.label.color Color of the edge labels. Defaults to black.
#' @param edge.label.font Font type for edge labels (1 is plain, 2 is bold, 3 is italic, 4 is bold and italic, 5 is symbol font). Defaults to 1.
#' @param background.col Background color of the plot. Defaults to "grey96".
#' @param overlap.line.color Color of the line representing observed overlap in the null model plot. Defaults to "red2".
#' @param overlap.line.lwd Line width of the observed overlap line. Defaults to 3.
#' @param overlap.line.lty Line type of the observed overlap line. Defaults to 1.
#' @param hist.side Position of the null model histogram relative to the network plot, either "bottom" or "right". Defaults to "bottom".
#' @param standardize.freqs Logical. If TRUE, standardizes frequencies in the null model plot. Defaults to FALSE.
#' @param cols Number of columns for the plot layout. Defaults to NULL.
#' @param ... Further arguments passed to \code{\link[qgraph]{qgraph}}.
#'
#' @seealso `calculateOverlap`, `randomizeOverlaps`
#' @export


plotOverlap <- function(overlaps = NULL,
                        random.results = NULL,
                        color.nodes.by = "group",
                        scale.nodes.by = NULL,
                        remove.missing = FALSE,
                        min.val = NULL,
                        cut.val = NULL,
                        group.order = NULL,
                        plot.stats = TRUE,
                        cex.title = 1.4,
                        cex.lab = 1.2,
                        cex.axis = 1.0,
                        cex.legend = 1.0,
                        legend.inset = c(-0.02, 0),
                        network.layout = "spring",
                        nodes.color = NULL,
                        nodes.size = c(1.2, 2.2),
                        nodes.label.cex = 1,
                        nodes.label.scale = FALSE,
                        nodes.label.color = "white",
                        edge.color = NULL,
                        edge.curved = 0.5,
                        edge.width = c(0.4, 3.5),
                        edge.label.cex = 1.4,
                        edge.label.color = "black",
                        edge.label.font = 1,
                        background.col = "grey96",
                        overlap.line.color = "red2",
                        overlap.line.lwd = 3,
                        overlap.line.lty = 1,
                        hist.side = "bottom",
                        standardize.freqs = FALSE,
                        cols=NULL,
                        ...){


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # validate parameters
  errors <- c()

  # check if at least one of 'overlaps' or 'random.results' is provided
  if (is.null(overlaps) && is.null(random.results)) errors <- c(errors, "At least one of 'overlaps' or 'random.results' must be provided.")

   # validate 'overlaps' if it is provided
  if (!is.null(overlaps)) {
    if (!requireNamespace("qgraph", quietly=TRUE)) errors <- c(errors, "The 'qgraph' package is required to plot the overlap network but is not installed. Please install 'qgraph' using install.packages('qgraph') and try again.")
    if (!is.data.frame(overlaps) || !"ids" %in% names(attributes(overlaps))) errors <- c(errors, "'overlaps' format not recognized. Please make sure to use the output from the 'calculateOverlap' function.")
    if (!is.null(scale.nodes.by) && length(scale.nodes.by) != length(attributes(overlaps)$ids)) errors <- c(errors, "scale.nodes.by variable must have the same length as the number of IDs in overlaps.")
    if (!is.null(color.nodes.by) && color.nodes.by != "group" && length(color.nodes.by) != length(attributes(overlaps)$ids)) errors <- c(errors, "color.nodes.by variable must have the same length as the number of IDs in overlaps.")
  }
  # validate 'random.results' if it is provided
  if (!is.null(random.results)) {
    if (!inherits(random.results, "list")) errors <- c(errors, "random.results must be a list, as expected from the output of the 'randomizeOverlaps' function.")
    if (!"metric" %in% names(attributes(random.results))) errors <- c(errors, "random.results format not recognized. Please make sure to use the output from the 'randomizeOverlaps' function.")
  }
  # validate 'hist.side'
  if(!hist.side %in% c("bottom", "right")) errors <- c(errors, "hist.side must be one of 'bottom' or 'right'.")
  # print errors if any
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # save the current par settings and ensure they are restored upon function exit
  original_par <- par(no.readonly=TRUE)
  on.exit(par(original_par))

  # set legend.inset
  if(length(legend.inset)==1)  legend.inset <- rep(legend.inset, 2)

  # get attributes
  if (!is.null(overlaps)){
    pairwise_overlaps <- overlaps
    id.groups <- attributes(overlaps)$id.groups
    complete_ids <- as.character(attributes(overlaps)$ids)
  }
  if(!is.null(random.results)) {
    pairwise_overlaps <- random.results$pairwise_results
    id.groups <- attributes(random.results)$id.groups
    complete_ids <- as.character(attributes(random.results)$ids)
    metric <- attributes(random.results)$metric
    iterations <- attributes(random.results)$iterations
  }

  # create checkup table for id.groups
  if(!is.null(id.groups)){
    id_lookup <- reshape2::melt(id.groups)
    colnames(id_lookup) <- c("id", "group")
    id_lookup$group <- as.factor(id_lookup$group)
  }

  # print to console
  .printConsole("Generating overlap network")


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # add missing ids if required
  if(!remove.missing){
    all_pairs <- as.data.frame(t(combn(complete_ids, 2)))
    colnames(all_pairs) <- c("id1", "id2")
    existing_pairs <- pairwise_overlaps[, c("id1", "id2")]
    reverse_existing_pairs <- pairwise_overlaps[, c("id2", "id1")]
    colnames(reverse_existing_pairs) <- c("id1", "id2")
    existing_pairs <- unique(rbind(existing_pairs, reverse_existing_pairs))
    missing_pairs <- all_pairs[!paste(all_pairs$id1, all_pairs$id2) %in% paste(existing_pairs$id1, existing_pairs$id2), ]
    if(nrow(missing_pairs)>0){
      if(!is.null(id.groups)){
        missing_pairs$group1 <- names(id.groups)[as.numeric(plyr::mapvalues(missing_pairs$id1, id_lookup$id, id_lookup$group, warn_missing=FALSE))]
        missing_pairs$group2 <- names(id.groups)[as.numeric(plyr::mapvalues(missing_pairs$id2, id_lookup$id, id_lookup$group, warn_missing=FALSE))]
        missing_pairs$type <- paste(missing_pairs$group1, "<->", missing_pairs$group2)
      }
      pairwise_overlaps <- plyr::rbind.fill(pairwise_overlaps, missing_pairs)
    }
  }

  # split overlaps by type (if required)
  if(!any(colnames(pairwise_overlaps)=="type")) pairwise_overlaps$type <- "All"
  if(!is.null(group.order)){
    pairwise_overlaps$type <- factor(pairwise_overlaps$type, levels=unique(pairwise_overlaps$type)[group.order])
  } else{
    pairwise_overlaps$type <- as.factor(pairwise_overlaps$type)
  }
  types <- levels(pairwise_overlaps$type)

  # assign name to color.nodes.by elements
  if (is.null(color.nodes.by)){
    color.nodes.by <- as.factor(rep(1, length(complete_ids)))
    names(color.nodes.by) <- complete_ids
  }else if(color.nodes.by=="group"){
    selected_ids <- reshape2::melt(id.groups)
    color.nodes.by <- factor(selected_ids$L1, levels=unique(selected_ids$L1))
    names(color.nodes.by) <- selected_ids$value
  }else{
    names(color.nodes.by) <- complete_ids
  }

  # convert color.nodes.by to factor if required
  if(!inherits(color.nodes.by, "factor")){
    color.nodes.by <- as.factor(color.nodes.by)
    warning("'color.nodes.by' vector converted to factor.", call.=FALSE)
  }

  # assign names to scale.nodes.by
  if (!is.null(scale.nodes.by) && is.null(names(scale.nodes.by))){
    names(scale.nodes.by) <- complete_ids
  }

  # split overlaps by type
  group_overlaps <- split(pairwise_overlaps, f=pairwise_overlaps$type)
  unique_ids <- lapply(group_overlaps, function(x) unique(c(x$id1, x$id2)))
  group_overlaps <- mapply(function(x,ids){x$id1<-factor(x$id1, levels=ids); return(x)}, x=group_overlaps, ids=unique_ids, SIMPLIFY=FALSE)
  group_overlaps <- mapply(function(x,ids){x$id2<-factor(x$id2, levels=ids); return(x)}, x=group_overlaps, ids=unique_ids, SIMPLIFY=FALSE)
  network_matrices <- lapply(group_overlaps, function(x) reshape2::dcast(x, formula="id1~id2", value.var="overlap", drop=FALSE))
  network_matrices <- lapply(network_matrices, function(x) {rownames(x)<-x$id1; return(x[,-1])})
  network_matrices <- lapply(network_matrices, as.matrix)

  # filter out discarded comparisons
  valid_overlaps <- pairwise_overlaps[!is.na(pairwise_overlaps$overlap),]

  # calculate average shared monitoring period
  shared_period <- aggregate(valid_overlaps$shared_monit_days, by=list(valid_overlaps$type), mean, na.rm=TRUE)
  colnames(shared_period) <- c("type", "days")

  # get nº of dyads
  n_dyads <- aggregate(valid_overlaps$overlap, by=list(valid_overlaps$type), length)
  colnames(n_dyads) <- c("type", "n_dyads")

  # get binary degree (number of edges each node has)
  binary_degree <- network_matrices
  binary_degree <- lapply(binary_degree, function(x) unlist(apply(x, 1, function(y) length(which(y>0)))))
  binary_degree <- unlist(lapply(binary_degree, mean))


  ##############################################################################
  ## Set network variables #####################################################
  ##############################################################################

  # set color palette for nodes
  if(is.null(nodes.color)){
    if(is.null(id.groups)){
      nodes.color <- "#008177"
    }else{
      if(length(id.groups)==2)  nodes.color <- c("#643BA6", "darkorange2")
      else if(length(id.groups)==3)  nodes.color <- c("#643BA6", "darkorange2", "#008177")
      else if(length(id.groups)==4)  nodes.color <- c("#643BA6", "darkorange2", "#008177", "#AF5C60")
      else nodes.color <- terrain.colors(length(id.groups))
    }
  }

  # initialize holding lists
  network.params <- vector("list", length(network_matrices))
  node.sizes.list <- vector("list", length(network_matrices))
  node.colors.list <- vector("list", length(network_matrices))
  edge.labels.list <- vector("list", length(network_matrices))
  edge.colors.list <- vector("list", length(network_matrices))

  # set network properties for each matrix
  for(i in 1:length(network_matrices)){

    # define cut and min values
    vals <- as.numeric(network_matrices[[i]])
    vals <- vals[!is.na(vals)]
    vals <- vals[order(vals, decreasing=TRUE)]
    if(is.null(min.val)) min.val <- quantile(vals, 0.5)
    if(is.null(cut.val)) cut.val <- quantile(vals, 0.9)
    network.params[[i]] <- c(min.val, cut.val)

    # set node properties
    network_ids <- colnames(network_matrices[[i]])
    if(!is.null(scale.nodes.by)){
      node_sizes <- as.numeric(plyr::mapvalues(network_ids, names(scale.nodes.by), scale.nodes.by, warn_missing=FALSE))
      node.sizes.list[[i]] <- .rescale(node_sizes, from=range(scale.nodes.by), nodes.size)
    }else{
      node.sizes.list[[i]] <- rep(mean(nodes.size), length(network_ids))
    }
    node_types <- as.numeric(plyr::mapvalues(network_ids, names(color.nodes.by), color.nodes.by, warn_missing=FALSE))
    node.colors.list[[i]] <- nodes.color[node_types]

    # set edge properties (overlap)
    edge.labels.list[[i]] <- round(network_matrices[[i]],1)
    edge_colors <- network_matrices[[i]]
    edge_colors[edge_colors>=0] <- NA
    # if edge.color is null, average node colors
    if(is.null(edge.color)){
      #colnames(edge_colors) <- node.colors.list[[i]]
      #rownames(edge_colors) <- node.colors.list[[i]]
      color_pairs <- expand.grid(node.colors.list[[i]], node.colors.list[[i]])
      edge_colors <- apply(color_pairs, 1, function(x) mix_colors(c(x[1], x[2])))
      #edge_colors <- unlist(mapply(function(x, y) mix_colors(c(x, y)), rownames(edge_colors), colnames(edge_colors), SIMPLIFY=FALSE, USE.NAMES=FALSE))
    }
    edge.colors.list[[i]] <- edge_colors
  }


  # convert NA values to 0 for the plots
  network_matrices <- lapply(network_matrices, function(x) {x[is.na(x)] <- 0; return(x)})


  ##############################################################################
  ## Calculate maximum frequency for standardization ###########################
  ##############################################################################

  freq_range <- NULL
  if (standardize.freqs && !is.null(random.results)) {
    max_freqs <- sapply(group_overlaps, function(x) {
      type_indices <- which(as.character(pairwise_overlaps$type) == as.character(x$type[1]))
      if (length(type_indices) == 0) stop(paste("No valid type indices found for type:", x$type[1]), call.=FALSE)
      null_dist <- lapply(1:iterations, function(i) random.results$randomized_overlaps[[i]][type_indices])
      null_dist <- unlist(lapply(null_dist, mean, na.rm=TRUE))
      if (length(null_dist) == 0) stop(paste("Null distribution is empty for type:", x$type[1]), call.=FALSE)
      hist(null_dist, breaks=60, plot=FALSE)$counts
    })
    freq_range <- range(max_freqs)
  }


  ##############################################################################
  ## Set layout params #########################################################
  ##############################################################################

  # count the number of plots to create
  n_groups <- length(types)
  plots_per_group <- sum(!is.null(overlaps), !is.null(random.results))
  cols_per_group <- ifelse(hist.side == "right" && plots_per_group == 2, 2, 1)
  total_plots <- n_groups * plots_per_group

  # determine number of rows and columns for layout
  if(is.null(cols)) {
    cols <- ifelse(hist.side=="bottom" || plots_per_group==1, n_groups, 2)
  }else{
    if(cols==1 && plots_per_group==2 && hist.side == "right") {cols<-2; warning(paste("Number of columns is set to 1, but hist.side is set to right side. An additional column will be added."), call.=FALSE)}
    if(cols>2 && !cols%%2==0 && plots_per_group==2 && hist.side == "right") stop(paste("Number of columns is set to", cols, "but hist.side is set to left side). This might lead to an unintended layout"), call.=FALSE)
  }

  # initialize a list to store the layout matrix rows and a vector for row heights
  layout_rows <- list()
  row_heights <- c()
  plot_index <- 1

  # iterate over each group
  for (i in 1:n_groups) {

    # add a blank plot at the top for the title of the group
    layout_rows[[length(layout_rows) + 1]] <- rep(plot_index, cols_per_group)
    row_heights <- c(row_heights, 0.8)
    plot_index <- plot_index + 1

    # add the network plot if overlaps were supplied
    if (!is.null(overlaps)) {
      if (hist.side == "right") {
        # network on left, blank space for histogram on right
        layout_rows[[length(layout_rows) + 1]] <- rep(plot_index, 2)
      } else {
        # network plot fills only one column
        layout_rows[[length(layout_rows) + 1]] <- plot_index
      }
      row_heights <- c(row_heights, 6)
      plot_index <- plot_index + 1
    }

    # add the histogram plot if random.results is supplied
    if (!is.null(random.results)) {
      if (hist.side == "right" && plots_per_group == 2) {
        # histogram on right, network on left
        layout_rows[[length(layout_rows)]] <- c(plot_index -1, plot_index)
      } else {
        # histogram fills the row
        layout_rows[[length(layout_rows) + 1]] <- plot_index
        row_heights <- c(row_heights, 6)
      }
      plot_index <- plot_index + 1
    }
  }

  # convert layout rows to a single vector and reshape based on 'cols' parameter
  layout_vector <- unlist(layout_rows)
  n_rows <- ceiling(length(layout_vector) / cols)
  #fill_by_row <- hist.side != "bottom" && plots_per_group != 1
  fill_by_row <- ifelse(hist.side == "bottom" || plots_per_group==1, FALSE, TRUE)
  layout_matrix <- matrix(c(layout_vector, rep(NA, n_rows *cols - length(layout_vector))), ncol=cols, byrow=fill_by_row)

  # fill any NA positions
  na_positions <- which(is.na(layout_matrix))
  if(length(na_positions)>0) layout_matrix[na_positions] <- max(layout_matrix, na.rm=TRUE) + 1

  # set layout
  graphics::layout(layout_matrix, heights=row_heights)

  # set margin settings
  if(plot.stats) margins_network <- c(10, 2.5, 2, 2.5)
  else margins_network <- c(4, 2.5, 2, 2.5)
  margins_hist <- c(5, 3.5, 1, 2)


  ##############################################################################
  ## Generate figure ###########################################################
  ##############################################################################

  par(oma=c(1,1,1,1), mgp=c(2.6, 1, 0))

  # iterate over each matrix
  for(i in 1:n_groups){

    ############################################################################
    # title ####################################################################

    par(mar=c(1, 0, 1, 0))
    plot.new()
    graphics::text(x=0.5, y=0.5, labels=names(network_matrices)[i], cex=cex.title, font=2)

    ############################################################################
    # overlap network(s) #######################################################

    if(!is.null(overlaps)){

      # plot network
      qgraph::qgraph(network_matrices[[i]],
                     mar = margins_network,
                     layout = network.layout,
                     directed = FALSE,
                     minimum = network.params[[i]][1],
                     cut = network.params[[i]][2],
                     color = node.colors.list[[i]],
                     node.width = node.sizes.list[[i]],
                     label.color = nodes.label.color,
                     label.cex = nodes.label.cex,
                     label.scale = nodes.label.scale,
                     edge.labels = edge.labels.list[[i]],
                     edge.color = edge.colors.list[[i]],
                     edge.width = 1,
                     edge.label.cex = edge.label.cex,
                     edge.label.color = edge.label.color,
                     edge.label.font = edge.label.font,
                     edge.label.margin=0.005,
                     colFactor = 0.6,
                     curveAll = TRUE,
                     curveDefault = edge.curved,
                     curveShape = -1,
                     curveScale = TRUE,
                     curveScaleNodeCorrection = TRUE,
                     repulsion = 0.1,
                     trans = TRUE,
                     fade = TRUE,
                     usePCH = TRUE,
                     ...)

      # add title plus stats
      if(plot.stats){
        network_metrics1 <- paste0("N\u00ba of individuals: ", ncol(network_matrices[[i]]))
        network_metrics2 <- paste0("N\u00ba of dyads: ", n_dyads$n_dyads[i])
        network_metrics3 <- paste0("Mean binary degree: ", sprintf("%.1f", binary_degree[i]))
        network_metrics4 <- paste0("Mean shared period: ", sprintf("%.0f", shared_period$days[i]), " days")
        legend("bottomleft", inset=c(0.06, 0.06), legend=c(network_metrics1, network_metrics2, network_metrics3, network_metrics4),
               bty="n", y.intersp=1, cex=cex.legend, xpd=NA)
      }

    }


    ############################################################################
    # null model distribution ##################################################
    if(!is.null(random.results)){

      # set margins
      par(mar=margins_hist)

      # estimate mean observed overlap
      obs_overlap <- mean(group_overlaps[[i]]$overlap, na.rm=TRUE)
      p_val <- as.numeric(random.results$summary$`P-value`[random.results$summary$Type==types[i]])

      # fetch null distribution
      type_indices <- which(pairwise_overlaps$type==types[[i]])
      null_dist <- lapply(1:iterations, function(i) random.results$randomized_overlaps[[i]][type_indices])
      null_dist <- unlist(lapply(null_dist, mean, na.rm=TRUE))

      # calculate barplot params
      h <- hist(null_dist, breaks=60, plot=FALSE)
      xrange <- range(h$breaks)
      xdelta <- xrange[2] - xrange[1]
      if(standardize.freqs==TRUE) ylim <- c(0, max(freq_range)*1.1)
      else  ylim <- c(0, max(h$counts)*1.1)
      bar.colors <- mix_colors(unique(node.colors.list[[i]]))

      # calculate axis and lim params based on observed overlap and null distribution ranges
      if(obs_overlap<xrange[1]){
        xlim <- c(xrange[1]-xdelta*(0.1), xrange[2])
        breakval <- xrange[1]-xdelta*(0.05)
        axis_range <- c(breakval, xrange[2])
        lineval <- xrange[1]-xdelta*(0.12)
      }else if(obs_overlap>xrange[2]){
        xlim <- c(xrange[1], xrange[2]+xdelta*(0.1))
        breakval <- xrange[2]+xdelta*(0.05)
        axis_range <- c(xrange[1], breakval)
        lineval <- xrange[2]+xdelta*(0.12)
      }else{
        xlim <- xrange
        breakval <- NULL
        axis_range <- xrange
        lineval <- obs_overlap
      }

      # generate barplot
      hist(null_dist, breaks=60, col=NA, border=FALSE, main="", xlab="", ylab="", yaxs="i", axes=FALSE, ylim=ylim, xlim=xlim)
      if(!is.null(background.col)) rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col="grey96", border=NULL)
      hist(null_dist, breaks=60, col=adjustcolor(bar.colors, alpha.f=0.7), lwd=0.05, main="", xlab="", ylab="",
           xlim=xlim, ylim=ylim, axes=FALSE, las=1, yaxs="i", add=TRUE)
      # add axes labels
      title(xlab="Overlap (%)", cex.lab=cex.lab, xpd=NA)
      title(ylab="Frequency", cex.lab=cex.lab, xpd=NA)
      # add axes
      axis(1, at=grid::grid.pretty(axis_range, n=6), labels=sprintf("%.2f", grid::grid.pretty(axis_range, n=6)), pos=0, las=1, cex.axis=cex.axis)
      axis(2, at=grid::grid.pretty(ylim), labels=grid::grid.pretty(ylim), las=1, cex.axis=cex.axis)
      # add axis break (if required)
      if(!is.null(breakval)) plotrix::axis.break(axis=1, breakpos=breakval, pos=NULL, bgcol="white", breakcol="black", style="slash", brw=0.02)
      # add observed overlap
      segments(x0=lineval, y0=par("usr")[3], y1=par("usr")[4], col=overlap.line.color, lwd=overlap.line.lwd, lty=overlap.line.lty, lend=3)
      # add legend
      p_val <- ifelse(p_val<0.001, "p < 0.001", paste("p =", sprintf("%.3f", p_val)))
      info1 <- paste0("Simulations = ", iterations)
      info2 <- paste0("Obs. overlap = ", sprintf("%.2f", obs_overlap), "%")
      legend("topleft", inset=legend.inset, legend=c(info1, info2, p_val), bty="n", y.intersp=1, cex=cex.legend)
      box()
    }
  }

}


################################################################################
# Define internal function to mix edge colors ##################################
################################################################################

#' Mix Colors
#'
#' @description The 'mix_colors' function blends any number of provided colors and returns the resulting color.
#' @return The resulting color from the average of the provided colors, in hexadecimal format.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

mix_colors <- function(colors) {
  # convert each color to RGB values
  rgb_values <- lapply(colors, col2rgb)
  # sum the RGB values
  sum_rgb <- Reduce(`+`, rgb_values)
  # calculate the average RGB values
  n_colors <- length(colors)
  avg_rgb <- sum_rgb / n_colors
  # create the mixed color
  rgb(avg_rgb[1], avg_rgb[2], avg_rgb[3], maxColorValue = 255)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
