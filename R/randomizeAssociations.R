#######################################################################################################
# Randomize associations (null model permutation tests) ###########################################
#######################################################################################################

#' Test association-network co-occurrences against a null model

#' @description This function tests the null hypothesis of temporally independent space use
#' (i.e., each animal occurs independently of the other) using Monte Carlo permutation tests.
#' The algorithm permutes entries within each column, keeping the total number of detections
#' of each individual and the relative occurrence frequencies across receivers unchanged.
#' Permutations can be constrained by given factors (e.g., diel phase) and are performed
#' only within the monitoring period of each animal (interval between its release
#' and the date of its last detection). Empirical p-values are calculated by comparing
#' the observed spatiotemporal overlaps against their null distribution. The p-values can
#' be one-tailed or two-tailed based on the specified alternative hypothesis (a
#' continuity correction is applied to avoid p-values of exactly 0 or 1).
#'
#' @param table A data frame containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding to the receiver/station with the
#' highest number of detections), as returned by \code{\link{createWideTable}}.
#' @param overlaps Data frame containing paiwise overlaps, as returned by \code{\link{calculateAssociations}}.
#' @param constraint.by Optional. Variable(s) to constrain permutations, e.g. to account
#' for potential diel or seasonal trends in animal occurrences. If supplied, permutations
#' across time bins will be restricted to the same combination of levels of these variables.
#' @param alternative Character string specifying the alternative hypothesis. Must be one of "two.sided",
#' "less", or "greater". "two.sided" tests for deviation in either direction, "less" tests if the observed
#' value is significantly less than the null distribution, and "greater" tests if the observed value is
#' significantly greater than the null distribution. Defaults to "two.sided".
#' @param iterations Number of Monte Carlo iterations (simulated datasets). Defaults to 1000.
#' @param conf.level Confidence level for the test. Defaults to 0.95 (95%).
#' @param p.adjust.method Method used to correct the per-dyad p-values for multiple
#' comparisons, passed to \code{\link[stats]{p.adjust}}. Because every dyad in the network is
#' tested simultaneously, some correction is recommended to control false positives. Defaults
#' to `"fdr"` (Benjamini-Hochberg false discovery rate); use `"none"` to retain raw p-values, or
#' any other method accepted by \code{p.adjust} (e.g. `"bonferroni"`, `"holm"`). The dyad-level
#' `Association` label is based on the adjusted p-value, while both the raw (`P-value`) and
#' adjusted (`Adjusted p-value`) values are reported.
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations.
#' Run \code{parallel::detectCores()} to check the number of available cores.
#' @param random.seed Optional. Set the seed for a reproducible randomization.
#' See \code{\link[base]{set.seed}}.
#' @param verbose Logical; print progress (a permutation progress bar). Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{\strong{summary}}{A summary table with overall results.}
#'   \item{\strong{pairwise_results}}{A data frame similar to that returned by the \code{\link{calculateAssociations}}
#'  function, but with additional columns indicating the mean of the null distribution,
#'  empirical p-values, and association type for each animal dyad.}
#'   \item{\strong{randomized_overlaps}}{A numeric matrix of simulated pairwise overlaps, with one row
#'   per animal dyad (the row names give the dyad, e.g. \code{"A-B"}) and one column per iteration. If a
#'   subset variable was specified in the supplied overlaps, this is instead a named list of such
#'   matrices, one per level of the subset variable.}
#' }
#'
#' The summary table includes the following columns:
#' \itemize{
#'   \item {Type}: The type of comparison. This column is included only when \code{id.groups} are identified in the supplied overlaps.
#'   \item {Subset}: This column (e.g., a temporal subset) is included only when a subset variable was used to calculate the supplied overlaps.
#'   \item {N dyads}: The number of dyads (pairs of individuals).
#'   \item {Mean interval (d)}: The mean monitoring period (in days).
#'   \item {Mean overlap (%)}: The mean observed overlap percentage with standard deviation.
#'   \item {Mean null distr (%)}: The mean overlap percentage from the null distribution with standard deviation.
#'   \item {P-value}: The estimated p-value(s), indicating the significance of the observed overlaps.
#'   \item {Association}: The association type (e.g., positive, negative, non-significant) based on the p-value.
#'   \item {Pairs Non Sig}: The number of pairs with non-significant.
#'   \item {Pairs > Random}: The number of pairs with observed overlap greater than the null distribution.
#'   \item {Pairs < Random}: The number of pairs with observed overlap less than the null distribution.
#' }
#'
#' @seealso \code{\link{calculateAssociations}}, \code{\link{plotAssociations}}
#'
#' @references
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns. Ecology, 81(9), 2606-2621.
#'
#' Farine, D. R. (2017). A guide to null models for animal social network analysis. Methods Ecol Evol 8: 1309-1320.
#'
#' @examples
#' \donttest{
#' data(rays)
#' wide <- createWideTable(rays, value.col = "station")
#' assoc <- calculateAssociations(wide)
#' # test the observed co-occurrences against a permutation null model
#' # (iterations kept low here for speed; use the default 1000 in practice)
#' rand <- randomizeAssociations(wide, assoc, iterations = 100, random.seed = 1)
#' rand$summary
#' }
#'
#' @export


randomizeAssociations <- function(table,
                             overlaps,
                             constraint.by = NULL,
                             iterations = 1000,
                             alternative = c("two.sided"),
                             conf.level = 0.95,
                             p.adjust.method = "fdr",
                             cores = 1,
                             random.seed = NULL,
                             verbose = getOption("moby.verbose", TRUE)) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # validate parameters
  errors <- c()
  if (!is.data.frame(table)) errors <- c(errors, "The 'table' argument must be a data frame.")
  if(!c("ids") %in% names(attributes(table))) errors <- c(errors, "The supplied table does not seem to be in the wide format. Please use the output of the 'createWideTable' function.")
  if (!is.data.frame(overlaps) || !c("ids") %in% names(attributes(overlaps))) errors <- c(errors, "'overlaps' format not recognized. Please make sure to use the output from the 'calculateAssociations' function.")
  if (!is.null(constraint.by) && !all(constraint.by %in% colnames(table))) errors <- c(errors,  "Constraint variable(s) not found in the supplied table.")
  if (!is.null(random.seed) && !inherits(random.seed, "numeric")) errors <- c(errors, "The random seed must be an integer.")
  if (!alternative %in% c("two.sided", "less", "greater")) errors <- c(errors, "Invalid value for argument 'alternative'. Must be one of 'two.sided', 'less', or 'greater'.")
  if (length(p.adjust.method)!=1 || !p.adjust.method %in% stats::p.adjust.methods) errors <- c(errors, paste0("Invalid 'p.adjust.method'. Must be one of: ", paste(stats::p.adjust.methods, collapse=", "), "."))
  if (!is.numeric(cores) || cores < 1 || cores %% 1 != 0) errors <- c(errors, "Cores must be a positive integer.")
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # parallel computing degrades gracefully to serial when its (Suggested) packages are unavailable
  if (cores > 1) {
    pkgs <- c("foreach", "doSNOW", "doRNG", "parallel")
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      .mobyWarn("Parallel computing needs the package(s) ", paste(missing_pkgs, collapse = ", "),
                "; running on a single core instead. Install them to enable 'cores > 1'.")
      cores <- 1
    } else if (parallel::detectCores() < cores) {
      .mobyWarn("Only ", parallel::detectCores(), " core(s) available; reducing 'cores'.")
      cores <- parallel::detectCores()
    }
  }

  # get table attributes
  complete_ids <- as.character(attributes(table)$ids)
  timebin.col <- attributes(table)$timebin.col
  start_dates <- attributes(table)$start.dates
  end_dates <- attributes(table)$end.dates

  # get overlap results attributes
  id.groups <- attributes(overlaps)$id.groups
  group.comparisons <- attributes(overlaps)$group.comparisons
  metric <- attributes(overlaps)$metric
  subset <- attributes(overlaps)$subset


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # retrieve number of animals
  n_ids <- length(complete_ids)
  animal_cols <- which(colnames(table) %in% complete_ids)

  # create constraint variable
  if(!is.null(constraint.by)){
    if(length(constraint.by)==1){table$constraint <- table[,constraint.by]}
    if(length(constraint.by)>1){table$constraint <- apply(table[,constraint.by], 1 , paste, collapse="//")}
  }else{
    table$constraint <- 1
  }

  # create subset variable if not provided
  if(is.null(subset)){
    table$subset <- 1
    overlaps$subset <- 1
  } else {
    table$subset <- interaction(table[, subset], drop = TRUE)
  }

  # split table by subset
  subset_list <- list("table"= split(table, f=table$subset),
                      "overlaps"=split(overlaps, f=overlaps$subset))


  ##############################################################################
  ## Run permutations ##########################################################
  ##############################################################################

  # at least two individuals are required to form pairs
  if (n_ids < 2) stop("At least two individuals are required to randomize pairwise overlaps.", call.=FALSE)

  # generate all pairwise combinations
  pairwise_combinations <- combn(seq_len(n_ids), 2)
  unique_pairs <- apply(pairwise_combinations, 2, function(x) paste(complete_ids[x], collapse="-"))
  n_pairs <- length(unique_pairs)

  # initialize final results list
  final_results <- vector("list", length(subset_list$table))

  # set alpha value
  alpha <- 1 - conf.level

  # set random seed for reproducibility
  if(!is.null(random.seed)) set.seed(random.seed)

  ##############################################################################
  # iterate over each subset   #################################################
  for(s in seq_along(subset_list$table)){

    subset_table <- subset_list$table[[s]]
    subset_overlaps <- subset_list$overlaps[[s]]

    # assign pair to current data
    if(!any(colnames(subset_overlaps)=="pair")){
      subset_overlaps$pair <- paste0(subset_overlaps$id1, "-", subset_overlaps$id2)
    }

    # build the [pairs x iterations] null-overlap matrix (rows = unique_pairs). The serial engine is
    # bit-identical to the previous implementation; cores > 1 splits the iterations across workers.
    if (!is.null(subset)) .mobyInform("Running permutations - ", names(subset_list$table)[s], verbose = verbose)
    else if (cores > 1 && s == 1) .mobyInform("Starting parallel computation: ", cores, " cores", verbose = verbose)

    null_mat <- .nullOverlapMatrix(subset_table, complete_ids, timebin.col, start_dates, end_dates,
                                   pairwise_combinations, unique_pairs, id.groups, group.comparisons,
                                   metric, iterations, cores, random.seed, verbose)
    obs_overlap <- subset_overlaps$association[match(unique_pairs, subset_overlaps$pair)]
    pv <- .empiricalPvalues(obs_overlap, null_mat, alternative, p.adjust.method, alpha,
                            labels = c(more = "positive", less = "negative"))

    pairwise_stats <- data.frame(
      pair = unique_pairs, mean_null = pv$mean_null,
      p_value = ifelse(is.na(pv$p_value), NA_character_, sprintf("%.3f", pv$p_value)),
      significance = pv$association,
      p_adjusted = ifelse(is.na(pv$p_adjusted), NA_character_, sprintf("%.3f", pv$p_adjusted)),
      stringsAsFactors = FALSE)
    pairwise_stats <- .joinKeep(pairwise_stats, subset_overlaps, by="pair", type="left")


    ############################################################################
    # calculate population p-values ############################################

    # if no specific group IDs are provided, set all pairwise comparisons to a single "All" group.
    if(is.null(id.groups)){
      pairwise_stats$type <- factor("All")
    }

    # extract unique types to iterate through for p-value calculations.
    types <- levels(pairwise_stats$type)
    # get pairs in 'pairwise_stats' for each group comparison type
    type_pairs <- lapply(types, function(x) pairwise_stats$pair[which(pairwise_stats$type==x)])
    # extract randomized overlaps for each group
    null_dist <- lapply(seq_len(iterations), function(i) lapply(type_pairs, function(p) null_mat[unique_pairs %in% p, i]))
    # calculate the mean of the null distribution values for each group, handling NA values.
    null_dist <- lapply(null_dist, function(x) lapply(x, mean, na.rm=TRUE))
    null_dist <- lapply(seq_along(types), function(x) unlist(lapply(null_dist, function(y) y[[x]])))
    names(null_dist) <- types
    # reshape the list into a data frame
    null_dist <- .meltList(null_dist)
    colnames(null_dist) <- c("mean_overlap", "type")
    # ensure 'type' is an ordered factor
    null_dist$type <- factor(null_dist$type, levels=types)
    # calculate the mean of the null distribution for each group comparison type
    mean_null <- stats::aggregate(null_dist$mean_overlap, by=list(null_dist$type), mean, na.rm=TRUE)
    colnames(mean_null) <- c("type", "mean_null")
    mean_null <- mean_null[order(mean_null$type),]
    # start building the final population statistics table
    population_stats <- mean_null
    # add standard deviations of the null distribution means to 'population_stats'
    population_stats$sd <- stats::aggregate(null_dist$mean_overlap, by=list(null_dist$type), sd, na.rm=TRUE)$x
    population_stats$mean_null <-  paste(sprintf("%.2f", population_stats$mean_null), "\u00b1", sprintf("%.2f", population_stats$sd))
    population_stats <- population_stats[,-3]
    # set column names to make the output more interpretable
    colnames(population_stats) <- c("Type", "Mean null distr (%)")

    # population-level test: observed mean overlap vs the distribution of per-iteration mean nulls,
    # one raw-p test per comparison type (a deliberately coarser null than the per-dyad test, and
    # NOT cross-type adjusted). Uses the same shared inference helper for correct tails + capping.
    obs_overlap <- stats::aggregate(pairwise_stats$association, by=list(pairwise_stats$type), mean, na.rm=TRUE, drop=FALSE)$x
    null_by_type <- lapply(types, function(x) null_dist$mean_overlap[null_dist$type==x])
    null_mat_pop <- do.call(rbind, null_by_type)
    pv_pop <- .empiricalPvalues(obs_overlap, null_mat_pop, alternative, p.adjust.method = "none", alpha,
                                labels = c(more = "positive", less = "negative"))
    population_stats$"P-value" <- sprintf("%.3f", pv_pop$p_value)
    population_stats$Association <- pv_pop$association
    population_stats$Subset <- names(subset_list$overlaps)[s]


    ##############################################################################
    # Generate summary table #####################################################

    # add no dyads
    valid_overlaps <- pairwise_stats[!is.na(pairwise_stats$association),]
    summary_table <- stats::aggregate(valid_overlaps$pair, by=list(valid_overlaps$type), function(x) length(unique(x)), drop=FALSE)
    summary_table$x[is.na(summary_table$x)] <- 0
    colnames(summary_table) <- c("Type", "N dyads")
    # add mean shared period
    summary_table$period <- round(stats::aggregate(valid_overlaps$shared_monit_days, by=list(valid_overlaps$type), mean, na.rm=TRUE, drop=FALSE)$x)
    colnames(summary_table)[3] <- "Mean interval (d)"
    # add mean overlap ± standard deviation
    summary_table$overlap <- stats::aggregate(pairwise_stats$association, by=list(pairwise_stats$type), mean, na.rm=TRUE, drop=FALSE)$x
    summary_table$overlap <- sprintf("%.2f", summary_table$overlap)
    summary_table$sd <- stats::aggregate(pairwise_stats$association, by=list(pairwise_stats$type), function(x) sd(x, na.rm=TRUE), drop=FALSE)$x
    summary_table$sd[is.na(summary_table$sd)] <- 0
    summary_table$sd <- sprintf("%.2f", summary_table$sd)
    summary_table$overlap <- paste(summary_table$overlap, "\u00b1", summary_table$sd)
    colnames(summary_table)[4] <- "Mean overlap (%)"
    summary_table <- .dropCols(summary_table, "sd")
    # add mean null distribution (%) + p-values
    summary_table <- .joinKeep(summary_table, population_stats, by="Type", type="left")
    pairwise_signif <- as.data.frame.matrix(table(pairwise_stats$type, pairwise_stats$significance))
    labels <- c("Pairs Non Sig","Pairs > Random","Pairs < Random")
    colnames(pairwise_signif) <- gsub("non-significant", labels[1], colnames(pairwise_signif), fixed=TRUE)
    colnames(pairwise_signif) <- gsub("positive", labels[2], colnames(pairwise_signif), fixed=TRUE)
    colnames(pairwise_signif) <- gsub("negative", labels[3], colnames(pairwise_signif), fixed=TRUE)
    pairwise_signif[labels[!(labels %in% colnames(pairwise_signif))]] <- 0
    pairwise_signif$Type <- rownames(pairwise_signif)
    rownames(pairwise_signif) <- NULL
    summary_table <- .joinKeep(summary_table, pairwise_signif, by="Type", type="left")
    summary_table <- summary_table[summary_table$`N dyads`>0,]
    summary_table$`Mean overlap (%)`[summary_table$`Mean interval (d)`==0] <- "-"
    summary_table$`P-value`[summary_table$`Mean interval (d)`==0] <- "-"
    summary_table$`Mean null distr (%)`[summary_table$`Mean interval (d)`==0] <- "-"
    summary_table[is.na(summary_table)] <- "-"

    # clean up and format final results
    count_cols <- which(colnames(summary_table) %in% labels)
    pairwise_stats$mean_null <- round(pairwise_stats$mean_null, 2)
    original_cols <- match(colnames(overlaps), colnames(pairwise_stats))
    new_cols <- setdiff(seq_len(ncol(pairwise_stats)), original_cols)
    pairwise_stats <- pairwise_stats[,c(original_cols, new_cols)]
    pairwise_stats <- pairwise_stats[!is.na(pairwise_stats$id1) & !is.na(pairwise_stats$id2),]


    ##############################################################################
    # Save group results #########################################################

    final_results[[s]] <- list("summary"=summary_table, "pairwise_results"=pairwise_stats, "randomized_overlaps"=null_mat)
  }



    #######################################################################################################
    # Return results ######################################################################################
    #######################################################################################################

    if(is.null(subset)){
      final_results <- final_results[[1]]
      final_results$summary <- .dropCols(final_results$summary, "Subset")
    }else{
      snames <- names(subset_list$table)
      out <- list(summary          = do.call("rbind", lapply(final_results, `[[`, "summary")),
                  pairwise_results = do.call("rbind", lapply(final_results, `[[`, "pairwise_results")),
                  randomized_overlaps = setNames(lapply(final_results, `[[`, "randomized_overlaps"), snames))
      rownames(out$summary) <- NULL
      out$summary <- out$summary[,c(1,8,2:11)]
      out$summary <- out$summary[order(out$summary$Type, out$summary$Subset),]
      rownames(out$pairwise_results) <- NULL
      final_results <- out
    }

    # print run time
    .reportRuntime(start.time, verbose)

    # assemble final list and add attributes
    attr(final_results, 'ids') <- complete_ids
    attr(final_results, 'id.groups') <- id.groups
    attr(final_results, 'metric') <- metric
    attr(final_results, 'subset') <- subset
    attr(final_results, 'constraint.by') <- constraint.by
    attr(final_results, 'group.comparisons') <- group.comparisons
    attr(final_results, 'iterations') <- iterations
    attr(final_results, 'alternative') <- alternative
    attr(final_results, 'random.seed') <- random.seed

    # return results
    return(final_results)
  }


  ################################################################################
  # Internal permutation-null engine (matrix-based) ##############################
  ################################################################################

  #' Precompute the iteration-invariant structures for the permutation null model
  #' @description Builds an integer station-code matrix, the per-(constraint group, individual) row
  #' indices to permute, and each pair's shared-window row indices and validity. Everything here is
  #' constant across Monte Carlo iterations, so it is computed once rather than every iteration.
  #' @note This function is intended for internal use within the 'moby' package.
  #' @keywords internal
  #' @noRd

.nullPrecompute <- function(sub, ids, timebin.col, start_dates, end_dates,
                            pairwise_combinations, id.groups, group.comparisons, metric) {
  n_ids <- length(ids); n_pairs <- ncol(pairwise_combinations)
  # integer station codes via one shared factor (fast equality, NA preserved); columns aligned to 'ids'
  vals <- lapply(ids, function(id) as.character(sub[[id]]))
  lev <- sort(unique(unlist(vals, use.names = FALSE)))
  M <- matrix(NA_integer_, nrow = nrow(sub), ncol = n_ids, dimnames = list(NULL, ids))
  for (j in seq_len(n_ids)) M[, j] <- match(vals[[j]], lev)
  # numeric time + monitoring windows (avoids per-comparison POSIXct S4 dispatch)
  tnum <- as.numeric(sub[[timebin.col]])
  sdates <- setNames(as.numeric(start_dates), names(start_dates))[ids]
  edates <- setNames(as.numeric(end_dates), names(end_dates))[ids]
  # constraint groups, in split()'s order (factor levels, else sorted unique)
  cons <- sub$constraint
  groups <- if (is.factor(cons)) levels(factor(cons)) else sort(unique(as.character(cons)))
  cons <- as.character(cons)
  grp_rows <- lapply(groups, function(g) which(cons == g))
  # per-(group, individual) window row indices to permute; which() drops NA-date comparisons, exactly
  # as the previous .sampleCols() did, so the permutation is applied to the identical set of rows
  win <- lapply(grp_rows, function(gr) { tg <- tnum[gr]
    lapply(seq_len(n_ids), function(j) gr[which(tg >= sdates[j] & tg <= edates[j])]) })
  # per-pair shared-window rows + validity (id.group filter, monitoring-window overlap, detections)
  grp_of <- if (!is.null(id.groups))
    vapply(ids, function(id) { w <- which(vapply(id.groups, function(g) id %in% g, logical(1)))
      if (length(w)) w[1] else NA_integer_ }, integer(1)) else NULL
  pair_rows <- vector("list", n_pairs); valid <- logical(n_pairs)
  for (p in seq_len(n_pairs)) {
    a <- pairwise_combinations[1, p]; b <- pairwise_combinations[2, p]
    if (!is.null(id.groups)) {
      if (is.na(grp_of[a]) || is.na(grp_of[b])) next
      if (group.comparisons == "within" && grp_of[a] != grp_of[b]) next
      if (group.comparisons == "between" && grp_of[a] == grp_of[b]) next
    }
    if (is.na(edates[a]) || is.na(edates[b])) next
    st <- max(sdates[a], sdates[b]); en <- min(edates[a], edates[b])
    if (is.na(st) || is.na(en) || en <= st) next
    r <- which(tnum >= st & tnum <= en); if (!length(r)) next
    pair_rows[[p]] <- r; valid[p] <- TRUE
  }
  list(M = M, win = win, n_groups = length(groups), n_ids = n_ids, n_pairs = n_pairs,
       pc = pairwise_combinations, pair_rows = pair_rows, valid_pairs = which(valid),
       w = if (metric == "half-weight") 0.5 else 1)
}

  #' Run a chunk of permutation iterations on a precomputed structure
  #' @description Returns an `[n_pairs x iterations]` null-overlap matrix. Each iteration permutes each
  #' individual's detections within its monitoring-window intersection with each constraint group, then
  #' computes every valid pair's association index over its precomputed shared window. `sample()` is
  #' drawn in the same (group, individual) order as the previous engine, so the serial path is
  #' bit-identical to the earlier implementation.
  #' @note This function is intended for internal use within the 'moby' package.
  #' @keywords internal
  #' @noRd

.nullIterateChunk <- function(pc, iterations, pb = NULL) {
  out <- matrix(NA_real_, nrow = pc$n_pairs, ncol = iterations)
  for (it in seq_len(iterations)) {
    Mr <- pc$M
    for (gi in seq_len(pc$n_groups)) { wi <- pc$win[[gi]]
      for (j in seq_len(pc$n_ids)) { idx <- wi[[j]]
        if (length(idx) > 1L) Mr[idx, j] <- sample(Mr[idx, j]) } }
    for (p in pc$valid_pairs) {
      r <- pc$pair_rows[[p]]; va <- Mr[r, pc$pc[1, p]]; vb <- Mr[r, pc$pc[2, p]]
      pres <- !is.na(va) & !is.na(vb); both <- sum(pres)
      matching <- sum(pres & va == vb)
      incomplete <- sum(is.na(va) != is.na(vb))
      tot <- both + incomplete * pc$w                    # both + incomplete*w == matching+mismatch+incomplete*w
      out[p, it] <- if (tot > 0) matching / tot * 100 else NaN
    }
    if (!is.null(pb)) .progressSet(pb, it)
  }
  out
}

  #' Build the `pairs x iterations` null-overlap matrix for one subset (serial or parallel)
  #' @description Orchestrates the permutation null model. The serial path draws from the current RNG
  #' stream (seeded once by the caller) and is bit-identical to the previous implementation; the
  #' parallel path splits the iterations across `cores`, each worker precomputing locally and drawing a
  #' reproducible independent RNG stream (via doRNG) - statistically equivalent, not bit-identical to
  #' serial (as before).
  #' @note This function is intended for internal use within the 'moby' package.
  #' @keywords internal
  #' @noRd

.nullOverlapMatrix <- function(sub, ids, timebin.col, start_dates, end_dates, pairwise_combinations,
                               unique_pairs, id.groups, group.comparisons, metric, iterations, cores,
                               random.seed, verbose) {
  if (cores == 1L) {
    pc <- .nullPrecompute(sub, ids, timebin.col, start_dates, end_dates, pairwise_combinations,
                          id.groups, group.comparisons, metric)
    pb <- .progressBar(iterations, verbose)
    null_mat <- .nullIterateChunk(pc, iterations, pb)
    .progressEnd(pb)
  } else {
    # split the iterations into 'cores' chunks; each worker precomputes locally (cheap) and runs its
    # chunk, so only the subset table is shipped to the workers, not a matrix per iteration
    k <- NULL                                          # bound by foreach() at runtime (quiets R CMD check)
    precompute <- .nullPrecompute; iterate <- .nullIterateChunk  # export as locals (no ::: on own package)
    chunk <- as.integer(cut(seq_len(iterations), cores, labels = FALSE))
    sizes <- as.integer(table(factor(chunk, levels = seq_len(cores)))); sizes <- sizes[sizes > 0]
    `%dorng%` <- doRNG::`%dorng%`
    cl <- parallel::makeCluster(cores); doSNOW::registerDoSNOW(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    pb <- .progressBar(length(sizes), verbose)
    opts <- if (!is.null(pb)) list(progress = function(n) setTxtProgressBar(pb, n)) else list()
    parts <- foreach::foreach(k = seq_along(sizes), .options.snow = opts, .options.RNG = random.seed,
                              .export = c("precompute", "iterate"), .packages = "moby") %dorng% {
      pc <- precompute(sub, ids, timebin.col, start_dates, end_dates, pairwise_combinations,
                       id.groups, group.comparisons, metric)
      iterate(pc, sizes[k])
    }
    .progressEnd(pb)
    null_mat <- do.call(cbind, parts)
  }
  rownames(null_mat) <- unique_pairs
  null_mat
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
