#######################################################################################################
## mobyNetwork object #################################################################################
#######################################################################################################

# moby represents the two telemetry network types it supports through a single lightweight S3
# object, `mobyNetwork`:
#   * "association" networks: nodes are individuals, edges are co-occurrence-based association
#     indices (see calculateAssociations()).
#   * "movement" networks: nodes are locations/receivers, edges are directed transitions
#     (see calculateTransitions(); added in a later release).
#
# As with `mobyData`, a `mobyNetwork` IS a data.frame (its rows are the edge list), so it can be
# used wherever an edge data frame is expected, while carrying the node table and network metadata
# as attributes. Edge and node tables are retrieved with networkEdges() / networkNodes().

#' Construct a mobyNetwork object (internal)
#'
#' @param edges A data frame of edges (one row per dyad / transition).
#' @param nodes A data frame describing the network nodes.
#' @param type Network type: `"association"` or `"movement"`.
#' @param ... Further attributes to attach (e.g. `ids`, `id.groups`, `metric`, `subset`).
#' @return A `mobyNetwork` object.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.mobyNetwork <- function(edges, nodes, type = c("association", "movement"), ...) {
  type <- match.arg(type)
  edges <- as.data.frame(edges)
  attr(edges, "nodes") <- nodes
  attr(edges, "network.type") <- type
  extra <- list(...)
  for (nm in names(extra)) attr(edges, nm) <- extra[[nm]]
  class(edges) <- unique(c("mobyNetwork", "data.frame"))
  edges
}


#' Inspect a moby network object
#'
#' @description Accessors and predicates for `mobyNetwork` objects (returned by
#' \code{\link{calculateAssociations}} and, in future, the movement-network functions).
#' A `mobyNetwork` is a `data.frame` of edges that additionally carries the node table and
#' network metadata as attributes.
#'
#' @param x An object.
#' @return `is_mobyNetwork()` returns a logical; `networkEdges()` returns the edge data frame;
#' `networkNodes()` returns the node data frame; `networkType()` returns `"association"` or
#' `"movement"`.
#' @name mobyNetwork
#' @seealso \code{\link{calculateAssociations}}
#' @examples
#' data(rays)
#' trans <- calculateTransitions(rays, spatial.col = "station")
#' is_mobyNetwork(trans)
#' networkType(trans)
#' head(networkEdges(trans))
#' head(networkNodes(trans))
NULL

#' @rdname mobyNetwork
#' @export
is_mobyNetwork <- function(x) inherits(x, "mobyNetwork")

#' @rdname mobyNetwork
#' @export
networkEdges <- function(x) {
  if (!is_mobyNetwork(x)) stop("'x' is not a mobyNetwork object.", call. = FALSE)
  out <- x
  attr(out, "nodes") <- NULL
  attr(out, "network.type") <- NULL
  class(out) <- "data.frame"
  out
}

#' @rdname mobyNetwork
#' @export
networkNodes <- function(x) {
  if (!is_mobyNetwork(x)) stop("'x' is not a mobyNetwork object.", call. = FALSE)
  attr(x, "nodes")
}

#' @rdname mobyNetwork
#' @export
networkType <- function(x) attr(x, "network.type")


#' @export
print.mobyNetwork <- function(x, ...) {
  type <- attr(x, "network.type")
  nodes <- attr(x, "nodes")
  cat("<mobyNetwork>", if (!is.null(type)) type else "", "network\n")
  cat("  nodes: ", if (!is.null(nodes)) nrow(nodes) else length(attr(x, "ids")),
      "  |  edges: ", nrow(x), "\n", sep = "")
  if (!is.null(attr(x, "metric"))) cat("  edge metric: ", attr(x, "metric"), "\n", sep = "")
  if (!is.null(attr(x, "subset"))) cat("  subset: ", paste(attr(x, "subset"), collapse = ", "), "\n", sep = "")
  if (!is.null(attr(x, "id.groups"))) cat("  groups: ", paste(names(attr(x, "id.groups")), collapse = ", "), "\n", sep = "")
  cat("  edge table (use networkEdges()/networkNodes() to extract):\n")
  print(utils::head(networkEdges(x)), ...)
  invisible(x)
}


#' Plot a moby network
#'
#' @description S3 `plot` method for \code{\link{mobyNetwork}} objects. For an `"association"`
#' network it draws the individual co-occurrence/association network (delegating to
#' \code{\link{plotAssociations}}). For a `"movement"` network it draws the directed transition
#' network between locations over a map (delegating to \code{\link{plotMovements}}), placing nodes
#' at their geographic coordinates when available.
#'
#' @param x A `mobyNetwork` object.
#' @param ... Further arguments passed to the underlying plotting routine:
#' \code{\link{plotAssociations}} for association networks (e.g. `color.by`, `nodes.size`,
#' `edge.color`) or \code{\link{plotMovements}} for movement networks (e.g. `land.shape`,
#' `epsg.code`, `background.layer`, `edge.metric`, `repel.nodes`).
#' @return Called for its side effect (a plot); invisibly returns `NULL`.
#' @seealso \code{\link{calculateAssociations}}, \code{\link{calculateTransitions}},
#' \code{\link{plotAssociations}}, \code{\link{plotMovements}}
#'
#' @examples
#' data(rays)
#' # association network (individual co-occurrences)
#' wide <- createWideTable(rays, value.col = "station")
#' assoc <- calculateAssociations(wide)
#' if (requireNamespace("qgraph", quietly = TRUE)) {
#'   plot(assoc)
#' }
#' # movement network (transitions between locations)
#' trans <- calculateTransitions(rays, spatial.col = "station")
#' plot(trans)
#'
#' @export

plot.mobyNetwork <- function(x, ...) {
  type <- attr(x, "network.type")
  if (is.null(type)) stop("This mobyNetwork has no 'network.type' attribute.", call. = FALSE)

  if (type == "association") {
    plotAssociations(overlaps = x, ...)
  } else if (type == "movement") {
    plotMovements(network = x, ...)
  } else {
    stop(paste0("Unknown network type: '", type, "'."), call. = FALSE)
  }
  invisible(NULL)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
