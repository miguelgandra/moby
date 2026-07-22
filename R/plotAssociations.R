#######################################################################################################
# Plot association (co-occurrence) network ##########################################################
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

#' @param overlaps A data frame containing the pairwise overlap data. Should be the output from the `calculateAssociations` function (a mobyNetwork).
#' @param random.results The null-model output of \code{\link{randomizeAssociations}}.
#' @param color.by How to colour the nodes: "group" (default) assigns a colour per ID group, or
#' "single" uses a single colour for all nodes.
#' @param scale.nodes.by Numeric vector to scale the nodes. Length should match the number of IDs in overlaps.
#' @param discard.missing Logical. If TRUE, excludes individuals without detections and null pairwise comparisons from the network.
#' Defaults to FALSE.
#' @param min.val Minimum value for overlap threshold. Can be a single value (used for all networks)
#' or a vector with one value per network. Defaults to 0.
#' @param cut.val Cutoff value for overlap threshold. Can be a single value (used for all networks)
#' or a vector with one value per network. Defaults to the 90% quantile of the overlap values.
#' @param group.order A vector specifying the order in which each group comparison type should be plotted.
#' This parameter can also be used to select and plot only a subset of the available comparisons.
#' When set to `NULL` (the default), all comparisons are plotted in their original order.
#' @param plot.stats Logical. If TRUE, additional network metrics are plotted
#' below each network. Defaults to TRUE.
#' @param legend.inset Optional. Specifies how far the legend is inset from the
#' null model histogram margins. If a single value is given, it is used for both margins;
#' if two values are given, the first is used for x- distance, the second for y-distance.
#' Defaults to c(-0.02, 0).
#' @param network.layout Layout algorithm for the network (see \code{\link[qgraph]{qgraph}}).
#' Defaults to "spring".
#' @param network.repulsion A scalar controlling the repulsion radius in the spring layout
#' of the network visualization. This value is passed to the `repulsion` argument
#' of `qgraph::qgraph`. Lower values (e.g., 0.05) reduce node separation, while higher
#' values (e.g., 0.5) increase spacing to avoid clustering. Useful for generating
#' alternative spatial arrangements. Defaults to 0.1.
#' @param nodes.color Color(s) for the nodes.
#' @param nodes.size Numeric vector of length 2, indicating the minimum and maximum
#' node sizes. Defaults to c(1, 2.2).
#' @param nodes.label.scale Logical. If TRUE, scales node labels. Defaults to FALSE.
#' @param nodes.label.color Color of the node labels. Defaults to white.
#' @param edge.color Color(s) for the edges. Can be a single value (used for all networks)
#' or a vector with one value per network. Defaults to NULL, which averages node colors.
#' @param edge.curved Numeric value or logical indicating the curvature of the edges. Defaults to 0.5.
#' @param edge.width Numeric value scaling all edge widths uniformly (default = 1).
#' Larger values make all edges thicker, smaller values make them thinner.
#' @param edge.label.color Colour of the edge labels. If NULL (default), inherits `edge.color`.
#' @param edge.label.font Font type for edge labels (1 is plain, 2 is bold, 3 is italic, 4 is bold and italic, 5 is symbol font). Defaults to 1.
#' @param background.color Background color of the plot. Defaults to "grey96".
#' @param overlap.line.color Color of the line representing observed overlap in the null model plot. Defaults to "red2".
#' @param overlap.line.lwd Line width of the observed overlap line. Defaults to 3.
#' @param overlap.line.lty Line type of the observed overlap line. Defaults to 1.
#' @param hist.side Position of the null model histogram relative to the network plot, either "bottom" or "right". Defaults to "bottom".
#' @param standardize.edge.weights Logical. If TRUE, edge widths are standardized across all networks so they are directly comparable. Defaults to TRUE.
#' @param standardize.freqs Logical. If TRUE, standardizes frequencies in the null model plot. Defaults to FALSE.
#' @param main Optional overall title drawn above the panel grid.
#' @param cex Global expansion factor scaling all text (titles, node/edge labels, axes, legend). Defaults to 1.
#' @param ncol Number of columns for the plot layout. Defaults to NULL.
#' @template deviceArgs
#' @param ... Further arguments passed to \code{\link[qgraph]{qgraph}}.
#' @return Invisibly, a tidy per-comparison-type data frame of network statistics (number of dyads,
#' mean overlap, mean shared monitoring days, mean binary degree, and the edge display thresholds).
#' @seealso \code{\link{calculateAssociations}}, \code{\link{randomizeAssociations}}, \code{\link[qgraph]{qgraph}}
#' @examples
#' # Pairwise co-occurrence network from a wide time-bin x individual table
#' wide  <- createWideTable(rays, value.col = "station")
#' assoc <- calculateAssociations(wide)
#' if (requireNamespace("qgraph", quietly = TRUE)) {
#'   plotAssociations(assoc)
#' }
#'
#' # Overlay the null-model histogram from a permutation test
#' \donttest{
#' if (requireNamespace("qgraph", quietly = TRUE)) {
#'   rand <- randomizeAssociations(wide, assoc, iterations = 100, random.seed = 1)
#'   plotAssociations(assoc, rand)
#' }
#' }
#' @export


plotAssociations <- function(overlaps = NULL,
                        random.results = NULL,
                        color.by = c("group", "single"),
                        scale.nodes.by = NULL,
                        discard.missing = FALSE,
                        min.val = NULL,
                        cut.val = NULL,
                        group.order = NULL,
                        plot.stats = TRUE,
                        legend.inset = c(-0.02, 0),
                        network.layout = "spring",
                        network.repulsion = 0.1,
                        nodes.color = NULL,
                        nodes.size = c(1.2, 2.2),
                        nodes.label.scale = FALSE,
                        nodes.label.color = "white",
                        edge.color = NULL,
                        edge.curved = 0.5,
                        edge.width = 1,
                        edge.label.color = NULL,
                        edge.label.font = 1,
                        background.color = "grey96",
                        overlap.line.color = "red2",
                        overlap.line.lwd = 3,
                        overlap.line.lty = 1,
                        hist.side = c("bottom", "right"),
                        standardize.edge.weights = TRUE,
                        standardize.freqs = FALSE,
                        main = NULL,
                        cex = 1,
                        ncol=NULL,
                        file = NULL,
                        width = NULL,
                        height = NULL,
                        res = 300,
                        ...){

  color.by <- match.arg(color.by)
  hist.side <- match.arg(hist.side)
  # single global cex -> derived text sizes (replaces the old cex_title/lab/axis/legend/label sprawl)
  cex_title <- 1.4 * cex; cex_lab <- 1.2 * cex; cex_axis <- 1.0 * cex
  cex_legend <- 1.0 * cex; cex_nodelab <- 1.0 * cex; cex_edgelab <- 1.4 * cex


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  errors <- c()
  if(is.null(overlaps) && is.null(random.results))
    errors <- c(errors, "At least one of 'overlaps' or 'random.results' must be provided.")

  # resolve the data object and its stamped attributes up-front (random.results wins when both are
  # given) so downstream validation - notably edge.color - can see id.groups without erroring
  src <- if(!is.null(random.results)) random.results else overlaps
  pairwise_overlaps <- if(!is.null(random.results)) random.results$pairwise_results else overlaps
  id.groups    <- attributes(src)$id.groups
  complete_ids <- as.character(attributes(src)$ids)
  metric       <- attributes(random.results)$metric
  iterations   <- attributes(random.results)$iterations

  if(!is.null(overlaps)){
    if(!requireNamespace("qgraph", quietly=TRUE)) errors <- c(errors, "The 'qgraph' package is required to plot the overlap network (install.packages('qgraph')).")
    if(!is.data.frame(overlaps) || !"ids" %in% names(attributes(overlaps))) errors <- c(errors, "'overlaps' not recognised. Please supply the output of calculateAssociations().")
    if(!is.null(scale.nodes.by) && length(scale.nodes.by) != length(complete_ids)) errors <- c(errors, "'scale.nodes.by' must have the same length as the number of IDs in 'overlaps'.")
  }
  if(!is.null(random.results)) errors <- c(errors, .validateRandomResults(random.results))
  if(!is.null(group.order) && !is.numeric(group.order)) errors <- c(errors, "'group.order' must be a numeric vector giving the order of comparison groups.")
  if(!is.null(edge.color) && length(edge.color) > 1){
    ntypes_v <- length(unique(pairwise_overlaps$type))
    if(!is.null(id.groups)){ if(length(edge.color) != ntypes_v) errors <- c(errors, paste0("'edge.color' must be a single colour or one colour per network (n=", ntypes_v, ").")) }
    else warning("'edge.color' has multiple values but 'id.groups' is not provided; only the first is used.", call. = FALSE)
  }
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)

  # save the current par settings and ensure they are restored upon function exit
  original_par <- .savePar()
  on.exit(.restorePar(original_par))

  # set legend.inset
  if(length(legend.inset)==1)  legend.inset <- rep(legend.inset, 2)

  # define min.val, cut.val and edge.color (repeat if necessary)
  ntypes <- ifelse(!is.null(overlaps), nlevels(overlaps$type), nlevels(random.results$pairwise_results$type))
  if(ntypes==0 || is.na(ntypes)) ntypes <- 1
  min.val <- if(!is.null(min.val) && length(min.val) == 1) rep(min.val, ntypes) else min.val
  cut.val <- if(!is.null(cut.val) && length(cut.val) == 1) rep(cut.val, ntypes) else cut.val
  edge.color <- if(!is.null(edge.color) && length(edge.color) == 1) rep(edge.color, ntypes) else edge.color

  # id -> group lookup (the data object and its attributes were resolved in the checks block above)
  if(!is.null(id.groups)){
    id_lookup <- .meltList(id.groups)
    colnames(id_lookup) <- c("id", "group")
    id_lookup$group <- as.factor(id_lookup$group)
    ordered_types <- levels(pairwise_overlaps$type)
  }

  # issue warning
  if(!is.null(group.order) && is.null(id.groups)){
    group.order <- NULL
    warning("'group.order' was set but will not be used as 'id.groups' is missing.", call.=FALSE)
  }



  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # add missing ids
  all_pairs <- as.data.frame(t(combn(complete_ids, 2)))
  colnames(all_pairs) <- c("id1", "id2")
  existing_pairs <- pairwise_overlaps[, c("id1", "id2")]
  reverse_existing_pairs <- pairwise_overlaps[, c("id2", "id1")]
  colnames(reverse_existing_pairs) <- c("id1", "id2")
  existing_pairs <- unique(rbind(existing_pairs, reverse_existing_pairs))
  missing_pairs <- all_pairs[!paste(all_pairs$id1, all_pairs$id2) %in% paste(existing_pairs$id1, existing_pairs$id2), ]
  if(nrow(missing_pairs)>0){
    if(!is.null(id.groups)){
      missing_pairs$group1 <- .mapValues(missing_pairs$id1, from=id_lookup$id, to=as.character(id_lookup$group), warn_missing=FALSE)
      missing_pairs$group2 <- .mapValues(missing_pairs$id2, from=id_lookup$id, to=as.character(id_lookup$group), warn_missing=FALSE)
      missing_pairs$type <-  paste(missing_pairs$group1, "<->", missing_pairs$group2)
      missing_pairs$type <- ifelse(missing_pairs$type %in% ordered_types, missing_pairs$type,  paste(missing_pairs$group2, "<->", missing_pairs$group1))
      missing_pairs$type <- factor(missing_pairs$type, levels=ordered_types)
    }
    pairwise_overlaps <- .rbindFill(pairwise_overlaps, missing_pairs)
  }

  # split overlaps by type (if required)
  if(!any(colnames(pairwise_overlaps)=="type")) pairwise_overlaps$type <- factor("All")

  # node grouping factor (named by id): one colour per ID group ("group"), or a single group
  if(color.by == "group" && !is.null(id.groups)){
    sel <- .meltList(id.groups)
    node_group <- factor(sel$L1, levels=unique(sel$L1)); names(node_group) <- sel$value
  }else{
    node_group <- as.factor(rep(1, length(complete_ids))); names(node_group) <- complete_ids
  }


  # assign names to scale.nodes.by
  if (!is.null(scale.nodes.by) && is.null(names(scale.nodes.by))){
    names(scale.nodes.by) <- complete_ids
  }

  # split overlaps by comparison type; drop any group left empty (e.g. by discard.missing) so the
  # per-group dcast cannot hit a zero-row / zero-level crash
  group_overlaps <- split(pairwise_overlaps, f=pairwise_overlaps$type, drop=TRUE)
  if(discard.missing) group_overlaps <- lapply(group_overlaps, function(x) x[!is.na(x$association), , drop=FALSE])
  group_overlaps <- group_overlaps[vapply(group_overlaps, function(x) nrow(x) > 0, logical(1))]
  if(length(group_overlaps) == 0) stop("No pairwise overlaps left to plot.", call. = FALSE)

  # comparison types + plot order are defined by the (surviving) groups, so every downstream index
  # (network matrices, stats, histogram labels) refers to the same panel
  types <- names(group_overlaps)
  if(is.null(group.order)) group.order <- seq_along(group_overlaps)
  else group.order <- group.order[group.order >= 1 & group.order <= length(group_overlaps)]

  unique_ids <- lapply(group_overlaps, function(x) unique(c(as.character(x$id1), as.character(x$id2))))
  group_overlaps <- mapply(function(x,ids){x$id1<-factor(x$id1, levels=ids); return(x)}, x=group_overlaps, ids=unique_ids, SIMPLIFY=FALSE)
  group_overlaps <- mapply(function(x,ids){x$id2<-factor(x$id2, levels=ids); return(x)}, x=group_overlaps, ids=unique_ids, SIMPLIFY=FALSE)
  network_matrices <- lapply(group_overlaps, function(x) .castWide(x, "id1", "id2", "association", fun.aggregate=function(z) z[1L], fill=NA))
  network_matrices <- lapply(network_matrices, function(x) {rownames(x)<-x$id1; return(as.matrix(x[,-1, drop=FALSE]))})

  # per-comparison-type statistics, ALIGNED to group_overlaps (were aggregated over a table that
  # dropped all-NA types, so a panel's legend could show another panel's numbers)
  n_dyads       <- vapply(group_overlaps, function(x) sum(!is.na(x$association)), numeric(1))
  shared_period <- vapply(group_overlaps, function(x) mean(x$shared_monit_days[!is.na(x$association)], na.rm=TRUE), numeric(1))
  mean_overlap  <- vapply(group_overlaps, function(x) mean(x$association, na.rm=TRUE), numeric(1))
  binary_degree <- vapply(network_matrices, function(x) mean(apply(x, 1, function(y) sum(y > 0, na.rm=TRUE))), numeric(1))


  ##############################################################################
  ## Set network variables #####################################################
  ##############################################################################

  # colourblind-safe node colours (one per ID group); single teal when there is only one group
  if(is.null(nodes.color)){
    nodes.color <- if(nlevels(node_group) == 1) "#008177" else .okabe_ito_pal(nlevels(node_group))
  }

  # get global overlap min and max
  global_max <- lapply(network_matrices, max, na.rm=TRUE)
  global_max <- max(unlist(global_max))

  # initialize holding lists
  network.params <- vector("list", length(network_matrices))
  node.sizes.list <- vector("list", length(network_matrices))
  node.colors.list <- vector("list", length(network_matrices))
  edge.labels.list <- vector("list", length(network_matrices))
  edge.colors.list <- vector("list", length(network_matrices))
  edge.label.colors.list <- vector("list", length(network_matrices))


  # set network properties for each matrix
  for(i in seq_along(network_matrices)){

    # define cut and min values
    vals <- as.numeric(network_matrices[[i]])
    vals <- vals[!is.na(vals)]
    vals <- vals[order(vals, decreasing=TRUE)]
    min_val <- ifelse(!is.null(min.val[i]), min.val[i], 0)
    cut_val <- ifelse(!is.null(cut.val[i]), cut.val[i], quantile(vals, 0.9))
    max_val <- ifelse(standardize.edge.weights, global_max, max(vals))
    network.params[[i]] <- c(min_val, cut_val, max_val)

    # set node properties
    network_ids <- colnames(network_matrices[[i]])
    if(!is.null(scale.nodes.by)){
      node_sizes <- as.numeric(.mapValues(network_ids, names(scale.nodes.by), scale.nodes.by, warn_missing=FALSE))
      node.sizes.list[[i]] <- .rescale(node_sizes, from=range(scale.nodes.by), nodes.size)
    }else{
      node.sizes.list[[i]] <- rep(mean(nodes.size), length(network_ids))
    }
    # index nodes.color by the factor codes of node_group; going through the labels
    # would coerce them to NA (as.numeric("A")) and blank out every node colour
    node_types <- as.integer(node_group[match(network_ids, names(node_group))])
    node.colors.list[[i]] <- nodes.color[node_types]

    # set edge labels
    edge.labels.list[[i]] <- round(network_matrices[[i]],1)

    # edge (and edge-label) colours: average the two endpoint node colours, or use a supplied colour.
    # edge.label.color defaults to inheriting edge.color (both NULL -> the averaged node colours).
    n_edge_cells <- length(network_matrices[[i]])
    node_pairs <- expand.grid(node.colors.list[[i]], node.colors.list[[i]])
    mixed <- apply(node_pairs, 1, function(x) .mixColors(c(x[1], x[2])))
    elc <- if(!is.null(edge.label.color)) edge.label.color else edge.color
    edge.colors.list[[i]]       <- if(is.null(edge.color)) mixed else rep(edge.color[[i]], n_edge_cells)
    edge.label.colors.list[[i]] <- if(is.null(elc))        mixed else rep(elc[[i]], n_edge_cells)
  }

  # convert NA values to 0 for the plots
  network_matrices <- lapply(network_matrices, function(x) {x[is.na(x)] <- 0; return(x)})

  # tidy per-comparison-type statistics (also the invisible return) + console summary
  stats_df <- data.frame(type = types, n_dyads = n_dyads, mean_overlap = mean_overlap,
                         shared_monit_days = shared_period, mean_binary_degree = binary_degree,
                         edge_min = vapply(network.params, `[`, numeric(1), 1),
                         edge_cut = vapply(network.params, `[`, numeric(1), 2),
                         edge_max = vapply(network.params, `[`, numeric(1), 3),
                         row.names = NULL, stringsAsFactors = FALSE)
  .printAssociationsSummary(stats_df, n_ids = length(complete_ids), metric = metric,
                            has_network = !is.null(overlaps), has_null = !is.null(random.results))


  ##############################################################################
  ## Calculate maximum frequency for standardization ###########################
  ##############################################################################

  freq_range <- NULL
  if (standardize.freqs && !is.null(random.results)) {
    max_freqs <- sapply(group_overlaps, function(x) {
      # fetch null distribution
      type_pairs <- pairwise_overlaps$pair[which(pairwise_overlaps$type==x$type[1])]
      nm <- random.results$randomized_overlaps
      null_dist <- colMeans(nm[rownames(nm) %in% type_pairs, , drop = FALSE], na.rm = TRUE)
      hist(null_dist, breaks=60, plot=FALSE)$counts
    })
    freq_range <- range(max_freqs)
  }


  ##############################################################################
  ## Set layout params #########################################################
  ##############################################################################

  # count the number of plots to create
  n_groups <- length(network_matrices)
  plots_per_group <- sum(!is.null(overlaps), !is.null(random.results))
  cols_per_group <- ifelse(hist.side == "right" && plots_per_group == 2, 2, 1)
  total_plots <- n_groups * plots_per_group

  # determine number of rows and columns for layout
  if(is.null(ncol)) {
    ncol <- ifelse(hist.side=="bottom" || plots_per_group==1, n_groups, 2)
    ncol <- ifelse(ncol>5, 5, ncol)
  }else{
    if(ncol==1 && plots_per_group==2 && hist.side == "right") {ncol<-2; warning(paste("Number of columns is set to 1, but hist.side is set to right side. An additional column will be added."), call.=FALSE)}
    if(ncol>2 && !ncol%%2==0 && plots_per_group==2 && hist.side == "right") stop(paste("Number of columns is set to", ncol, "but hist.side is set to left side). This might lead to an unintended layout"), call.=FALSE)
  }

  # initialize a list to store the layout matrix rows and a vector for row heights
  layout_rows <- list()
  row_heights <- c()
  plot_index <- 1

  # iterate over each group
  for (i in seq_len(n_groups)) {

    # add a blank plot at the top for the title of the group
    layout_rows[[length(layout_rows) + 1]] <- rep(plot_index, cols_per_group)
    row_heights <- c(row_heights, 0.8)
    plot_index <- plot_index + 1

    # add the network plot if overlaps were supplied
    if(!is.null(overlaps)) {
      layout_rows[[length(layout_rows) + 1]] <- plot_index
      row_heights <- c(row_heights, 6)
      plot_index <- plot_index + 1
    }
    # add the histogram plot if random.results is supplied
    if(!is.null(random.results)) {
      if(hist.side == "right" && plots_per_group == 2) {
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

  # convert layout rows to a single vector and reshape based on 'ncol' parameter
  layout_vector <- unlist(layout_rows)
  n_rows <- ceiling(length(layout_vector) / ncol)
  # split the vector by multiples of 3
  split_indices <- which(layout_vector %% (plots_per_group+1) == 0) + 1
  group_rows <- split(layout_vector, cumsum(seq_along(layout_vector) %in% split_indices))
  group_rows <- lapply(group_rows, function(x) matrix(x, ncol=cols_per_group,  byrow=TRUE))
  #group_rows <- split(layout_vector, ceiling(seq_along(layout_vector)/(plots_per_group+1)))
  fill_by_row <- ifelse(hist.side == "bottom" || plots_per_group==1, FALSE, TRUE)
  layout_matrix <- do.call(cbind, group_rows)

  # split the matrix into groups of x columns
  n <- ceiling(ncol(layout_matrix) / ncol)
  column_groups <- lapply(seq_len(n), function(i) {
    start_col <- (i - 1) * ncol + 1
    end_col <- min(i * ncol, ncol(layout_matrix))
    group <- layout_matrix[, start_col:end_col, drop=FALSE]
    # pad with NA to ensure the group has exactly 'ncol' columns
    if (ncol(group) < ncol) group <- cbind(group, matrix(NA, nrow = nrow(group), ncol = ncol - ncol(group)))
    group
  })

  # stack each group below the previous group
  layout_matrix <- do.call(rbind, column_groups)

  # fill any NA positions
  na_positions <- which(is.na(layout_matrix))
  if(length(na_positions)>0) layout_matrix[na_positions] <- max(layout_matrix, na.rm=TRUE) + 1

  # optional file output: size the device to the layout grid (before graphics::layout below)
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule=list(base=2, slope=3.2, n=ncol(layout_matrix), lo=5, hi=30),
                     h.rule=list(base=1, slope=2.8, n=nrow(layout_matrix), lo=4, hi=30),
                     crowd.unit="panels")
    on.exit(grDevices::dev.off(), add=TRUE, after=FALSE)
  }

  # set layout
  graphics::layout(layout_matrix, heights=row_heights)

  # set margin settings
  if(plot.stats) {
    margins_network <- c(10, 2.5, 2, 2.5)
  }else {
    margins_network <- c(4, 2.5, 2, 2.5)}
  margins_hist <- c(5, 3.5, 1, 2)


  ##############################################################################
  ## Generate figure ###########################################################
  ##############################################################################

  par(oma=c(1,1,1,1), mgp=c(2.6, 1, 0))

  # iterate over each matrix
  for(i in group.order){

    ############################################################################
    # title ####################################################################

    par(mar=c(1, 0, 1, 0))
    plot.new()
    if(length(group.order) > 1) graphics::text(x=0.5, y=0.5, labels=names(network_matrices)[i], cex=cex_title, font=2)

    ############################################################################
    # overlap network(s) #######################################################

    if(!is.null(overlaps)){

      # plot network
      qgraph::qgraph(network_matrices[[i]],
                     mar = margins_network,
                     layout = network.layout,
                     repulsion = network.repulsion,
                     directed = FALSE,
                     minimum = network.params[[i]][1],
                     maximum = network.params[[i]][3],
                     cut = network.params[[i]][2],
                     color = node.colors.list[[i]],
                     node.width = node.sizes.list[[i]],
                     label.color = nodes.label.color,
                     label.cex = cex_nodelab,
                     label.scale = nodes.label.scale,
                     edge.labels = edge.labels.list[[i]],
                     edge.color = edge.colors.list[[i]],
                     edge.width = edge.width,
                     edge.label.cex = cex_edgelab,
                     edge.label.color =  edge.label.colors.list[[i]],
                     edge.label.font = edge.label.font,
                     edge.label.margin=0.005,
                     colFactor = 0.6,
                     curveAll = TRUE,
                     curveDefault = edge.curved,
                     curveShape = -1,
                     curveScale = TRUE,
                     curveScaleNodeCorrection = TRUE,
                     trans = TRUE,
                     fade = TRUE,
                     usePCH = TRUE,
                     ...)

      # add title plus stats
      if(plot.stats){
        network_metrics1 <- paste0("No. of dyads: ", n_dyads[i])
        network_metrics2 <- paste0("Mean shared period: ", sprintf("%.0f", shared_period[i]), " days")
        network_metrics3 <- paste0("Mean overlap: ", sprintf("%.2f %%", mean(group_overlaps[[i]]$association, na.rm=TRUE)))
        network_metrics4 <- paste0("Edge display minimum: ", sprintf("%.2f %%", network.params[[i]][1]))
        network_metrics5 <- paste0("Edge display cutoff: ", sprintf("%.2f %%", network.params[[i]][2]))
        network_metrics6 <- paste0("Mean binary degree: ", sprintf("%.1f", binary_degree[i]))
        network_metrics <- c(network_metrics1, network_metrics2, network_metrics3, network_metrics4, network_metrics5, network_metrics6)
        legend("bottomleft", inset=c(0.06, 0.06), legend=network_metrics, bty="n", y.intersp=1, cex=cex_legend, xpd=NA)
      }

    }


    ############################################################################
    # null model distribution ##################################################
    if(!is.null(random.results)){

      # set margins
      par(mar=margins_hist)

      # estimate mean observed overlap
      obs_overlap <- mean(group_overlaps[[i]]$association, na.rm=TRUE)
      p_val <- as.numeric(random.results$summary$`P-value`[random.results$summary$Type==types[i]])

      # fetch null distribution
      type_pairs <- pairwise_overlaps$pair[which(pairwise_overlaps$type==types[i])]
      nm <- random.results$randomized_overlaps
      null_dist <- colMeans(nm[rownames(nm) %in% type_pairs, , drop = FALSE], na.rm = TRUE)

      # calculate barplot params
      h <- hist(null_dist, breaks=60, plot=FALSE)
      xrange <- range(h$breaks)
      xdelta <- xrange[2] - xrange[1]
      if(standardize.freqs==TRUE) ylim <- c(0, max(freq_range)*1.1)
      else  ylim <- c(0, max(h$counts)*1.1)
      bar.colors <- .mixColors(unique(node.colors.list[[i]]))

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
      if(!is.null(background.color)) rect(par('usr')[1], par('usr')[3], par('usr')[2], par('usr')[4], col=background.color, border=NA)
      hist(null_dist, breaks=60, col=adjustcolor(bar.colors, alpha.f=0.7), lwd=0.05, main="", xlab="", ylab="",
           xlim=xlim, ylim=ylim, axes=FALSE, las=1, yaxs="i", add=TRUE)
      # add axes labels
      title(xlab="Overlap (%)", cex.lab=cex_lab, xpd=NA)
      title(ylab="Frequency", cex.lab=cex_lab, xpd=NA)
      # add axes
      axis(1, at=grid::grid.pretty(axis_range, n=6), labels=sprintf("%.2f", grid::grid.pretty(axis_range, n=6)), pos=0, las=1, cex.axis=cex_axis)
      axis(2, at=grid::grid.pretty(ylim), labels=grid::grid.pretty(ylim), las=1, cex.axis=cex_axis)
      # add axis break (if required)
      if(!is.null(breakval)) .axisBreak(breakval, brw=0.02, bgcol="white", breakcol="black")
      # add observed overlap
      segments(x0=lineval, y0=par("usr")[3], y1=par("usr")[4], col=overlap.line.color, lwd=overlap.line.lwd, lty=overlap.line.lty, lend=3)
      # add legend
      p_val <- ifelse(p_val<0.001, "p < 0.001", paste("p =", sprintf("%.3f", p_val)))
      info1 <- paste0("Simulations = ", iterations)
      info2 <- paste0("Obs. overlap = ", sprintf("%.2f", obs_overlap), "%")
      legend("topleft", inset=legend.inset, legend=c(info1, info2, p_val), bty="n", y.intersp=1, cex=cex_legend)
      box()
    }
  }

  if(!is.null(main)) mtext(main, side = 3, outer = TRUE, font = 2, cex = cex_title, line = -0.5)
  invisible(stats_df)
}


#######################################################################################################
## Internal helper ####################################################################################

#' @keywords internal
#' @noRd
.printAssociationsSummary <- function(stats_df, n_ids, metric, has_network, has_null){
  kv <- .kv
  .summaryOpen("Association network")
  kv("Individuals", n_ids)
  kv("Comparisons", paste(stats_df$type, collapse = ", "))
  kv("Metric", if(is.null(metric)) "association index" else metric)
  kv("Total dyads", sum(stats_df$n_dyads))
  kv("Shows", paste(c(if(has_network) "network", if(has_null) "null-model histogram"), collapse = " + "))
  .summaryClose()
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
