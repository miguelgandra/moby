#######################################################################################################
## Randomize movement transitions (null movement network) ############################################
#######################################################################################################

# (.empiricalPvalues has moved to helpers-stats.R, now shared by randomizeTransitions and
#  randomizeAssociations so both network randomisers use one correct, tested inference routine.)


#' Test movement-network transitions against a null model
#'
#' @description Tests whether the directed transitions in a movement network (from
#' \code{\link{calculateTransitions}}) occur more (or less) often than expected under random
#' movement, using a within-individual permutation null model. For each individual, the
#' time-ordered sequence of visited locations is randomly permuted and the transition network
#' recomputed; repeating this many times builds a null distribution for each edge. Empirical,
#' multiple-comparison-corrected p-values then identify transitions (and the overall network)
#' that are stronger or weaker than expected.
#'
#' The null permutes the order in which each individual visited its set of locations, preserving
#' that set (and hence each animal's site fidelity) while randomising connectivity. It therefore
#' tests the non-randomness of *movement order / connectivity*, not of site use itself.
#' Note that movement-network null models are an active area of research; this is one defensible,
#' transparent choice rather than a universally agreed standard.
#'
#' @param network A `mobyNetwork` object of type `"movement"`, from \code{\link{calculateTransitions}}.
#' @param iterations Number of permutations. Defaults to 1000.
#' @param alternative Alternative hypothesis: `"two.sided"` (default), `"greater"` or `"less"`.
#' @param p.adjust.method Multiple-comparison correction across edges, passed to
#' \code{\link[stats]{p.adjust}}. Defaults to `"fdr"`; use `"none"` for raw p-values.
#' @param conf.level Confidence level used for the significance label. Defaults to 0.95.
#' @param random.seed Optional integer for reproducibility.
#' @param verbose Logical; print progress (a permutation progress bar). Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A list with
#' \item{edges}{A data frame of per-transition results: `group`, `from`, `to`, `n_movements`
#' (observed), `mean_null`, `p_value`, `p_adjusted`, and `association` (`"more"`, `"less"` or
#' `"non-significant"`, based on the adjusted p-value).}
#' \item{network}{A data frame with a network-level test per group: the observed number of
#' realised transitions versus the null, with its empirical p-value.}
#'
#' @seealso \code{\link{calculateTransitions}}, \code{\link{randomizeAssociations}}
#'
#' @examples
#' \donttest{
#' data(rays)
#' trans <- calculateTransitions(rays, spatial.col = "station")
#' # test each transition against a within-individual permutation null
#' # (iterations kept low here for speed; use the default 1000 in practice)
#' rand <- randomizeTransitions(trans, iterations = 100, random.seed = 1)
#' head(rand$edges)
#' rand$network
#' }
#'
#' @export

randomizeTransitions <- function(network,
                                 iterations = 1000,
                                 alternative = c("two.sided", "greater", "less"),
                                 p.adjust.method = "fdr",
                                 conf.level = 0.95,
                                 random.seed = NULL,
                                 verbose = getOption("moby.verbose", TRUE)) {

  if (!inherits(network, "mobyNetwork") || !identical(attr(network, "network.type"), "movement")) {
    .mobyAbort("'network' must be a movement 'mobyNetwork' object (see calculateTransitions()).")
  }
  alternative <- match.arg(alternative)
  if (length(p.adjust.method) != 1 || !p.adjust.method %in% stats::p.adjust.methods) {
    .mobyAbort("Invalid 'p.adjust.method'. Must be one of: ", paste(stats::p.adjust.methods, collapse = ", "), ".")
  }
  if (!is.null(random.seed)) set.seed(random.seed)
  alpha <- 1 - conf.level

  records <- attr(network, "transition_records")
  edges <- networkEdges(network)
  groups <- unique(edges$group)

  edge_results <- list()
  net_results <- list()

  .mobyInform(sprintf("Running %d permutations across %d group(s)...", iterations, length(groups)),
              verbose = verbose)

  for (g in groups) {
    e_g <- edges[edges$group == g, , drop = FALSE]
    rec_g <- records[[g]]
    if (nrow(e_g) == 0 || is.null(rec_g) || nrow(rec_g) == 0) next
    pb <- .progressBar(iterations, verbose)

    # reconstruct each individual's time-ordered sequence of visited (run) sites:
    # for consecutive transitions f1->t1, f2->t2, ... the run sequence is c(f1, t1, t2, ...)
    seqs <- lapply(split(rec_g, rec_g$id), function(r) c(r$from[1], r$to))

    obs_keys <- paste(e_g$from, e_g$to, sep = "\r")
    obs_counts <- e_g$n_movements
    null_mat <- matrix(0L, nrow = length(obs_keys), ncol = iterations, dimnames = list(obs_keys, NULL))
    net_null <- integer(iterations)

    for (it in seq_len(iterations)) {
      cnt <- integer(0)
      for (s in seqs) {
        sp <- rle(sample(s))$values
        if (length(sp) < 2) next
        k <- paste(sp[-length(sp)], sp[-1], sep = "\r")
        tb <- table(k)
        nm <- names(tb)
        cnt[nm] <- (if (length(cnt)) ifelse(nm %in% names(cnt), cnt[nm], 0L) else 0L) + as.integer(tb)
      }
      hit <- intersect(obs_keys, names(cnt))
      if (length(hit) > 0) null_mat[hit, it] <- cnt[hit]
      net_null[it] <- sum(cnt)
      .progressSet(pb, it)
    }
    .progressEnd(pb)

    pv <- .empiricalPvalues(obs_counts, null_mat, alternative, p.adjust.method, alpha)
    edge_results[[g]] <- data.frame(group = g, from = e_g$from, to = e_g$to,
                                    n_movements = e_g$n_movements, pv,
                                    row.names = NULL, stringsAsFactors = FALSE)

    # network-level test: observed total transitions vs the null
    obs_total <- sum(e_g$n_movements)
    p_net <- if (alternative == "less") (sum(net_null <= obs_total) + 1) / (iterations + 1)
             else if (alternative == "greater") (sum(net_null >= obs_total) + 1) / (iterations + 1)
             else min(2 * min((sum(net_null >= obs_total) + 1) / (iterations + 1),
                              (sum(net_null <= obs_total) + 1) / (iterations + 1)), 1)
    net_results[[g]] <- data.frame(group = g, n_transitions = obs_total,
                                   mean_null = mean(net_null), p_value = p_net,
                                   row.names = NULL, stringsAsFactors = FALSE)
  }

  result <- list(edges = do.call(rbind, edge_results),
                 network = do.call(rbind, net_results))
  attr(result, "iterations") <- iterations
  attr(result, "alternative") <- alternative
  attr(result, "p.adjust.method") <- p.adjust.method
  result
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
