#######################################################################################################
## Network metrics (shared across association & movement networks) ###################################
#######################################################################################################

# compute node- and network-level metrics for a single igraph graph.
# 'w' are the (positive) edge weights interpreted as connection strengths; for shortest-path
# measures (betweenness) they are inverted, since stronger ties = shorter distances.
.graphMetrics <- function(g, directed, community = TRUE) {
  w <- igraph::E(g)$weight
  n_nodes <- igraph::vcount(g)
  n_edges <- igraph::ecount(g)
  inv_w <- if (n_edges > 0) 1 / w else NULL

  if (directed) {
    nm <- data.frame(
      node = igraph::V(g)$name,
      in_degree   = as.integer(igraph::degree(g, mode = "in")),
      out_degree  = as.integer(igraph::degree(g, mode = "out")),
      in_strength = igraph::strength(g, mode = "in", weights = w),
      out_strength = igraph::strength(g, mode = "out", weights = w),
      betweenness = igraph::betweenness(g, weights = inv_w, directed = TRUE),
      row.names = NULL, stringsAsFactors = FALSE)
  } else {
    nm <- data.frame(
      node = igraph::V(g)$name,
      degree = as.integer(igraph::degree(g)),
      strength = igraph::strength(g, weights = w),
      betweenness = igraph::betweenness(g, weights = inv_w, directed = FALSE),
      eigenvector = tryCatch(igraph::eigen_centrality(g, directed = FALSE, weights = w)$vector,
                             error = function(e) rep(NA_real_, n_nodes)),
      clustering = igraph::transitivity(g, type = "local", weights = w, isolates = "zero"),
      row.names = NULL, stringsAsFactors = FALSE)
  }

  # community detection (walktrap supports weights and both directions)
  modularity <- NA_real_; n_comm <- NA_integer_; membership <- rep(NA_integer_, n_nodes)
  if (community && n_edges > 0) {
    cl <- tryCatch(igraph::cluster_walktrap(g, weights = w), error = function(e) NULL)
    if (!is.null(cl)) {
      membership <- as.integer(igraph::membership(cl))
      n_comm <- length(unique(membership))
      modularity <- tryCatch(igraph::modularity(g, membership, weights = w), error = function(e) NA_real_)
    }
  }
  nm$community <- membership

  net <- data.frame(
    n_nodes = n_nodes, n_edges = n_edges,
    density = igraph::edge_density(g),
    mean_strength = if (n_nodes > 0) mean(igraph::strength(g, weights = w)) else NA_real_,
    modularity = modularity, n_communities = n_comm,
    row.names = NULL, stringsAsFactors = FALSE)
  if (directed) net$reciprocity <- igraph::reciprocity(g)
  else net$transitivity <- igraph::transitivity(g, type = "global", weights = w)

  list(nodes = nm, network = net)
}


#' Compute node- and network-level metrics for a moby network
#'
#' @description Computes a consistent set of graph-theoretic metrics for either an association
#' (individual co-occurrence) network or a movement (location-transition) network produced by
#' \code{\link{calculateAssociations}} / \code{\link{calculateTransitions}}. Metrics are computed
#' per group/subset when the network carries one.
#'
#' Edge weights are interpreted as connection strengths (association index for association
#' networks; number of movements for movement networks). For shortest-path metrics
#' (betweenness) weights are internally inverted, so that stronger ties act as shorter distances
#' and high-betweenness nodes identify brokers / movement corridors.
#'
#' @param network A \code{\link{mobyNetwork}} object (type `"association"` or `"movement"`).
#' @param weight Optional name of the edge column to use as weight. Defaults to `"association"`
#' (association networks) or `"n_movements"` (movement networks).
#' @param community Logical; run weighted community detection (walktrap) and report modularity
#' and community membership. Defaults to TRUE.
#'
#' @return An object of class `mobyNetworkMetrics`: a list with
#' \item{nodes}{A data frame of per-node metrics (one row per node and group/subset). For
#' association networks: `degree`, `strength`, `betweenness`, `eigenvector`, `clustering`,
#' `community`. For movement networks: `in_degree`, `out_degree`, `in_strength`, `out_strength`,
#' `betweenness`, `community`.}
#' \item{network}{A data frame of network-level metrics per group/subset (`n_nodes`, `n_edges`,
#' `density`, `mean_strength`, `modularity`, `n_communities`, plus `transitivity` for
#' association networks or `reciprocity` for movement networks).}
#'
#' @seealso \code{\link{calculateAssociations}}, \code{\link{calculateTransitions}}
#'
#' @examples
#' data(rays)
#' trans <- calculateTransitions(rays, spatial.col = "station")
#' # node- and network-level graph metrics
#' metrics <- networkMetrics(trans)
#' metrics
#' head(metrics$nodes)
#'
#' @export

networkMetrics <- function(network, weight = NULL, community = TRUE) {

  if (!inherits(network, "mobyNetwork")) {
    stop("'network' must be a 'mobyNetwork' object (see calculateAssociations() / calculateTransitions()).", call. = FALSE)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Computing network metrics requires the 'igraph' package. Install it with install.packages('igraph').", call. = FALSE)
  }
  type <- attr(network, "network.type")
  edges <- networkEdges(network)

  node_rows <- list(); net_rows <- list()

  if (type == "movement") {
    wcol <- if (is.null(weight)) "n_movements" else weight
    if (!wcol %in% colnames(edges)) stop(paste0("Weight column '", wcol, "' not found in the network edges."), call. = FALSE)
    nodes_attr <- networkNodes(network)
    for (g in unique(edges$group)) {
      e_g <- edges[edges$group == g & !is.na(edges[[wcol]]) & edges[[wcol]] > 0,
                   c("from", "to", wcol), drop = FALSE]
      colnames(e_g)[3] <- "weight"
      v_g <- unique(as.character(nodes_attr$site[nodes_attr$group == g]))
      gr <- igraph::graph_from_data_frame(e_g, vertices = data.frame(name = v_g), directed = TRUE)
      m <- .graphMetrics(gr, directed = TRUE, community = community)
      node_rows[[g]] <- cbind(group = g, m$nodes, stringsAsFactors = FALSE)
      net_rows[[g]] <- cbind(group = g, m$network, stringsAsFactors = FALSE)
    }

  } else if (type == "association") {
    wcol <- if (is.null(weight)) "association" else weight
    if (!wcol %in% colnames(edges)) stop(paste0("Weight column '", wcol, "' not found in the network edges."), call. = FALSE)
    ids <- as.character(attr(network, "ids"))
    has_subset <- "subset" %in% colnames(edges)
    groups <- if (has_subset) unique(edges$subset) else "all"
    for (g in groups) {
      e_all <- if (has_subset) edges[edges$subset == g, , drop = FALSE] else edges
      e_g <- e_all[!is.na(e_all[[wcol]]) & e_all[[wcol]] > 0, c("id1", "id2", wcol), drop = FALSE]
      colnames(e_g)[3] <- "weight"
      gr <- igraph::graph_from_data_frame(e_g, vertices = data.frame(name = ids), directed = FALSE)
      m <- .graphMetrics(gr, directed = FALSE, community = community)
      node_rows[[g]] <- cbind(group = g, m$nodes, stringsAsFactors = FALSE)
      net_rows[[g]] <- cbind(group = g, m$network, stringsAsFactors = FALSE)
    }

  } else {
    stop(paste0("Unknown network type: '", type, "'."), call. = FALSE)
  }

  result <- list(nodes = do.call(rbind, node_rows), network = do.call(.rbindFill, net_rows))
  rownames(result$nodes) <- NULL; rownames(result$network) <- NULL
  attr(result, "network.type") <- type
  class(result) <- "mobyNetworkMetrics"
  result
}


#' @export
print.mobyNetworkMetrics <- function(x, ...) {
  cat("<mobyNetworkMetrics>", attr(x, "network.type"), "network\n")
  cat("  network-level metrics:\n")
  print(x$network, row.names = FALSE)
  cat("  node-level metrics:", nrow(x$nodes), "rows x", ncol(x$nodes) - 1, "metrics (see $nodes)\n")
  invisible(x)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
