#######################################################################################################
# Randomize overlaps  #################################################################################
#######################################################################################################

#' Test significance of estimated overlaps through null model randomization tests

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
#' @param overlaps Data frame containing paiwise overlaps, as returned by \code{\link{calculateOverlap}}.
#' @param constraint.by Optional. Variable(s) to constrain permutations, e.g. to account
#' for potential diel or seasonal trends in animal occurrences. If supplied, permutations
#' across time bins will be restricted to the same combination of levels of these variables.
#' @param id.groups Optional. A list containing ID groups for calculating stats independently
#' within each group and comparing relationships between IDs of different groups.
#' @param alternative Character string specifying the alternative hypothesis. Must be one of "two.sided",
#' "less", or "greater". "two.sided" tests for deviation in either direction, "less" tests if the observed
#' value is significantly less than the null distribution, and "greater" tests if the observed value is
#' significantly greater than the null distribution. Defaults to "two.sided".
#' @param iterations Number of Monte Carlo iterations (simulated datasets). Defaults to 1000.
#' @param conf.level Confidence level for the test. Defaults to 0.95 (95%).
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations.
#' Run \code{parallel::detectCores()} to check the number of available cores.
#' @param random.seed Optional. Set the seed for a reproducible randomization.
#' See \code{\link[base]{set.seed}}.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{\strong{summary}}{A summary table with overall results.}
#'   \item{\strong{pairwise_results}}{A data frame similar to that returned by the \code{\link{calculateOverlap}}
#'  function, but with additional columns indicating the mean of the null distribution,
#'  empirical p-values, and association type for each animal dyad.}
#'   \item{\strong{randomized_overlaps}}{A list of all randomized pairwise overlaps, containing an element for each iteration}
#' }\cr
#'
#' The summary table includes the following columns:
#' \itemize{
#'   \item {Type}: The type of comparison. This column is included only when \code{id.groups} are supplied (more than one comparison).
#'   \item {Nº dyads}: The number of dyads (pairs of individuals).
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
#' @references
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence patterns. Ecology, 81(9), 2606-2621.\cr\cr
#' Farine, D. R. (2017). A guide to null models for animal social network analysis. Methods Ecol Evol 8: 1309–1320.
#'
#' @importFrom foreach %dopar%
#' @export


randomizeOverlaps <- function(table, overlaps, constraint.by=NULL, id.groups=NULL, iterations=1000,
                              alternative=c("two.sided"), conf.level=0.95, cores=1, random.seed=NULL) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # validate parameters
  errors <- c()
  if (!is.data.frame(table)) errors <- c(errors, "The 'table' argument must be a data frame")
  if(!c("ids") %in% names(attributes(table))) errors <- c(errors, "The supplied table does not seem to be in the wide format. Please use the output of the 'createWideTable' function.")
  if (!is.data.frame(overlaps) || !c("ids") %in% names(attributes(overlaps))) errors <- c(errors, "'overlaps' format not recognized. Please make sure to use the output from the 'calculateOverlap' function.")
  if (!is.null(constraint.by) && !all(constraint.by %in% colnames(table))) errors <- c(errors,  "Constraint variable(s) not found in the supplied table.")
  if (is.null(id.groups) && any(colnames(overlaps) %in% c("group1", "group2"))) errors <- c(errors,  "ID groups not specified, but used in the supplied overlaps.")
  if (!is.null(random.seed) && !inherits(random.seed, "numeric")) errors <- c(errors, "The random seed must be an integer.")
  if (!alternative %in% c("two.sided", "less", "greater")) errors <- c(errors, "Invalid value for argument 'alternative'. Must be one of 'two.sided', 'less', or 'greater'.")
  if (!is.numeric(cores) || cores < 1 || cores %% 1 != 0) errors <- c(errors, "Cores must be a positive integer.")
  if(parallel::detectCores()<cores)  errors <- c(errors, paste("Please choose a different number of cores for parallel computing (only", parallel::detectCores(), "available)."))
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # get table attributes
  complete_ids <- as.character(attributes(table)$ids)
  timebin.col <- attributes(table)$timebin.col
  start_dates <- attributes(table)$start.dates
  end_dates <- attributes(table)$end.dates

  # get overlap results attributes
  metric <- attributes(overlaps)$metric
  subset <- attributes(overlaps)$subset
  group.comparisons <- attributes(overlaps)$group.comparisons

  # reorder ID levels if ID groups are defined
  if (!is.null(id.groups)) {
    if (any(duplicated(unlist(id.groups)))) stop("Repeated ID(s) in id.groups", call.=FALSE)
    if (any(!unlist(id.groups) %in% complete_ids)) stop("Some of the ID(s) in id.groups were not found in the supplied table", call.=FALSE)

    # remove individuals not present in ID groups
    selected_ids <- unique(unlist(id.groups))
    discard_indexes <- which(!complete_ids %in% selected_ids)
    id_cols <- which(colnames(table) %in% complete_ids)
    discard_cols <- id_cols[!colnames(table)[id_cols] %in% selected_ids]
    if (length(discard_cols) > 0) {
      table <- table[, -discard_cols]
      complete_ids <- complete_ids[-discard_indexes]
      start_dates <- start_dates[-discard_indexes]
      end_dates <- end_dates[-discard_indexes]
    }
  }

  # remove individuals with no detections
  missing_individuals <- names(which(is.na(end_dates)))
  if (length(missing_individuals) > 0) {
    table <- table[,-which(colnames(table) %in% missing_individuals)]
    complete_ids <- complete_ids[-which(complete_ids %in% missing_individuals)]
    start_dates <- start_dates[-which(names(start_dates) %in% missing_individuals)]
    end_dates <- end_dates[-which(names(end_dates) %in% missing_individuals)]
  }


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # set random seed for reproducibility
  if(!is.null(random.seed)) set.seed(random.seed)

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

  # split data by subset variable
  if(!is.null(subset)){

  }



  ##############################################################################
  ## Randomize data and calculate overlaps #####################################
  ##############################################################################

  # generate all pairwise combinations
  n_ids <- length(complete_ids)
  pairwise_combinations <- combn(1:n_ids, 2)
  n_pairs <- ncol(pairwise_combinations)

  # initiate results list
  random_results <- list()

  # print to console
  .printConsole("Running null model permutation tests")


  ###################################################################
  # randomize overlaps using the default method (single core)
    if(cores==1){

      # set progress bar
      pb <- txtProgressBar(max=iterations, initial=0, style=3)

      # start iterations
      for (i in 1:iterations) {
        # randomize data table
        randomized_table <- randomize(table, id.cols=animal_cols, start_dates, end_dates)
        # calculate overlap
        random_overlaps <- lapply(1:n_pairs, function(p){
          pairwiseOverlapLite(p, table=randomized_table, pairwise_combinations, complete_ids,
                              id.groups, group.comparisons, start_dates, end_dates, metric)})
        # save results
        random_results[[i]] <- unlist(random_overlaps)

        # show progress in console
        setTxtProgressBar(pb, i)
      }

      # Close progress bar
      close(pb)


  #############################################################
  # else use parallel computing to speed up calculations
  }else{

    # print to console
    cat(paste0("Starting parallel computation: ", cores, " cores\n"))

    # initialize cluster
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)

    # set progress bar
    pb <- txtProgressBar(max=iterations, initial=0, style=3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))

    # start iterations (using parallelization)
    random_results <- foreach::foreach(i=1:iterations, .options.snow=opts,
                                       .export=c("randomize", "sampleCols", "pairwiseOverlapLite")) %dopar% {
      randomized_table <- randomize(table, id.cols=animal_cols, start_dates, end_dates)
      unlist(lapply(1:n_pairs, function(p){pairwiseOverlapLite(p, table=randomized_table, pairwise_combinations, complete_ids,
                                                               id.groups, group.comparisons, start_dates, end_dates, metric)}))
    }

    # close progress bar and stop cluster
    close(pb)
    on.exit(parallel::stopCluster(cl))
  }


  ##############################################################################
  # Calculate pairwise p-values ################################################
  ##############################################################################

  # pairwise overlap significance
  unique_pairs <- apply(pairwise_combinations, 2, function(x) paste(complete_ids[x], collapse="-"))
  overlaps$pair <- paste0(overlaps$id1, "-", overlaps$id2)
  pairwise_stats <- list()
  alpha <- 1-conf.level

  # iterate over each dyad and calculate p-values
  for (i in 1:n_pairs) {
    null_dist <- unlist(lapply(random_results,  function(x) as.numeric(x)[i]))
    obs_overlap <- overlaps$overlap[overlaps$pair==unique_pairs[i]]
    if(is.na(obs_overlap)) {pairwise_stats[[i]]<-NA; next}
    if(all(is.na(null_dist))) {pairwise_stats[[i]]<-NA; next}
    # one-tailed test (greater than or equal to the observed value)
    if (alternative=="greater") {
      pval <- (sum(obs_overlap >= null_dist)+1)/(length(null_dist)+1)
      signif <- if (pval < alpha) "positive" else "non-significant"
    # one-tailed test (less than or equal to the observed value)
    } else if (alternative == "less") {
      pval <- (sum(obs_overlap <= null_dist)+1)/(length(null_dist)+1)
      signif <- if (pval < alpha) "negative" else "non-significant"
    # two-tailed test
    } else if (alternative == "two.sided") {
      pval_left <- (sum(obs_overlap >= null_dist)+1)/(length(null_dist)+1)
      pval_right <- (sum(obs_overlap <= null_dist)+1)/(length(null_dist)+1)
      pval <- 2 * min(pval_left, pval_right)
      if (pval < alpha) {
        signif <- if (obs_overlap > mean(null_dist)) "positive" else "negative"
      } else {
        signif <- "non-significant"
      }
    }
    pairwise_stats[[i]] <- data.frame("pair"=unique_pairs[i], "mean_null"=mean(null_dist),
                                      "p_value"=sprintf("%.3f", pval), "association"=signif)
  }

  # summary stats
  pairwise_stats <- pairwise_stats[unlist(lapply(pairwise_stats, length))>1]
  pairwise_stats <- do.call("rbind", pairwise_stats)
  pairwise_stats <- plyr::join(overlaps, pairwise_stats, by="pair", type="left")


  ##############################################################################
  # Calculate population p-values ##############################################
  ##############################################################################

  if(!is.null(id.groups)){
    pairwise_stats$type <- paste(pairwise_stats$group1, "<->", pairwise_stats$group2)
  }else{
    pairwise_stats$type <- "All"
  }

  pairwise_stats$type <- paste(pairwise_stats$group1, "<->", pairwise_stats$group2)
  types <- unique(pairwise_stats$type)
  type_indexes <- lapply(types, function(x) which(pairwise_stats$type==x))
  null_dist <- lapply(1:iterations, function(i) lapply(type_indexes, function(p) random_results[[i]][p]))
  null_dist <- lapply(null_dist, function(x) lapply(x, mean, na.rm=T))
  null_dist <- lapply(1:length(types), function(x) unlist(lapply(null_dist, function(y) y[[x]])))
  names(null_dist) <- types
  null_dist <- reshape2::melt(null_dist)
  colnames(null_dist) <- c("mean_overlap", "type")
  mean_null <- stats::aggregate(null_dist$mean_overlap, by=list(null_dist$type), mean, na.rm=T)
  colnames(mean_null) <- c("type", "mean_null")
  population_stats <- mean_null
  population_stats$sd <- stats::aggregate(null_dist$mean_overlap, by=list(null_dist$type), sd, na.rm=T)$x
  population_stats$mean_null <-  paste(sprintf("%.2f", population_stats$mean_null), "±", sprintf("%.2f", population_stats$sd))
  population_stats <- population_stats[,-3]
  colnames(population_stats) <- c("Type", "Mean null distr (%)")

  # calculate p-values
  obs_overlap <- stats::aggregate(pairwise_stats$overlap, by=list(pairwise_stats$type), mean, na.rm=T)$x
  null_dist <- lapply(types, function(x) null_dist$mean_overlap[null_dist$type==x])
  all_pvals <- numeric(length(types))
  all_signifs <- character(length(types))
  for (i in 1:length(types)) {
    null_dist_type <- null_dist[[i]]
    if (alternative == "greater") {
      pval <- (sum(obs_overlap[i]>=null_dist_type)+1)/(length(null_dist_type)+1)
      signif <- if (pval < alpha) "positive" else "non-significant"
    } else if (alternative == "less") {
      pval <- (sum(obs_overlap[i]<=null_dist_type)+1)/(length(null_dist_type)+1)
      signif <- if (pval < alpha) "negative" else "non-significant"
    } else { # "two.sided"
      pval_left <- (sum(obs_overlap[i] >= null_dist_type) + 1) / (length(null_dist_type) + 1)
      pval_right <- (sum(obs_overlap[i] <= null_dist_type) + 1) / (length(null_dist_type) + 1)
      pval <- 2 * min(pval_left, pval_right)
      signif <- if (pval < alpha) ifelse(obs_overlap[i] > mean(null_dist_type), "positive", "negative") else "non-significant"
    }
    all_pvals[i] <- pval
    all_signifs[i] <- signif
  }
  population_stats$"P-value" <- sprintf("%.3f", unlist(all_pvals))
  population_stats$Association <- unlist(all_signifs)


  ##############################################################################
  # Generate summary table #####################################################
  ##############################################################################

  # add nº dyads
  pairwise_stats <- pairwise_stats[!is.na(pairwise_stats$overlap),]
  summary_table <- stats::aggregate(pairwise_stats$pair, by=list(pairwise_stats$type), function(x) length(unique(x)))
  colnames(summary_table) <- c("Type", "Nº dyads")
  # add mean shared period
  summary_table$period <- round(stats::aggregate(pairwise_stats$shared_monit_days, by=list(pairwise_stats$type), mean, na.rm=T)$x)
  colnames(summary_table)[3] <- "Mean interval (d)"
  # add mean overlap ± standard deviation
  summary_table$overlap <- stats::aggregate(pairwise_stats$overlap, by=list(pairwise_stats$type), mean, na.rm=T)$x
  summary_table$overlap <- sprintf("%.2f", summary_table$overlap)
  summary_table$sd <- stats::aggregate(pairwise_stats$overlap, by=list(pairwise_stats$type), function(x) sd(x, na.rm=T))$x
  summary_table$sd <- sprintf("%.2f", summary_table$sd)
  summary_table$overlap <- paste(summary_table$overlap, "±", summary_table$sd)
  colnames(summary_table)[4] <- "Mean overlap (%)"
  summary_table <- summary_table[,-which(colnames(summary_table)=="sd")]
  # add mean null distribution (%) + p-values
  summary_table <- join(summary_table, population_stats, by="Type", type="left")
  pairwise_signif <- as.data.frame.matrix(table(pairwise_stats$type, pairwise_stats$association))
  labels <- c("Pairs Non Sig","Pairs > Random","Pairs < Random")
  colnames(pairwise_signif) <- gsub("non-significant", labels[1], colnames(pairwise_signif), fixed=T)
  colnames(pairwise_signif) <- gsub("positive", labels[2], colnames(pairwise_signif), fixed=T)
  colnames(pairwise_signif) <- gsub("negative", labels[3], colnames(pairwise_signif), fixed=T)
  pairwise_signif[labels[!(labels %in% colnames(pairwise_signif))]] <- 0
  pairwise_signif$Type <- rownames(pairwise_signif)
  rownames(pairwise_signif) <- NULL
  summary_table <- plyr::join(summary_table, pairwise_signif, by="Type", type="left")

  # clean up and format final results
  count_cols <- which(colnames(summary_table) %in% labels)
  summary_table <- summary_table[rowSums(summary_table[,count_cols])>0,]
  pairwise_stats <- pairwise_stats[,-which(colnames(pairwise_stats)=="pair")]
  if(is.null(id.groups)) {
    pairwise_stats <- pairwise_stats[,-which(colnames(pairwise_stats) %in% "type")]
  }
  pairwise_stats$mean_null <- round(pairwise_stats$mean_null, 2)


  #######################################################################################################
  # Return results ######################################################################################

  # print run time
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n"))

  # assemble final list and add attributes
  results_list <- list("summary"=summary_table, "pairwise_results"=pairwise_stats, "randomized_overlaps"=random_results)
  attr(results_list, 'metric') <- metric
  attr(results_list, 'subset') <- subset
  attr(results_list, 'constraint.by') <- constraint.by
  attr(results_list, 'group.comparisons') <- group.comparisons
  attr(results_list, 'iterations') <- iterations
  attr(results_list, 'alternative') <- alternative
  attr(results_list, 'random.seed') <- random.seed

  # return results
  return(results_list)
}


################################################################################
# Define internal function to randomize detections #############################
################################################################################

#' Randomize detections
#' @description This function randomizes the detection data while preserving the structure within
#' constraint groups.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
randomize <- function(table, id.cols, start_dates, end_dates) {
  # split the data by the constraint column
  data_subset <- split(table, table$constraint)
  # apply sampleCols to each subset and combine them into one data frame
  data_random <- do.call(rbind, lapply(data_subset, sampleCols, id.cols, start_dates, end_dates))
  # order the result by timebin
  data_random <- data_random[order(data_random$timebin), ]
  return(data_random)
}

#' Helper function - sample columns of a data frame within monitoring periods
#' @description This function samples each column (except the timebin column) of a data frame
#' within the monitoring periods of each animal, specified by start and end dates.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
sampleCols <- function(x, id.cols, start_dates, end_dates) {
  for (i in id.cols) {
    id <- colnames(x)[i]
    # filter rows within the monitoring period
    valid_rows <- which(x[,1]>=start_dates[id] &  x[,1]<=end_dates[id])
    if (length(valid_rows) > 1) {
      x[valid_rows, i] <- sample(x[valid_rows, i])
    }
  }
  return(x)
}

################################################################################
# Define core pairwise overlap function ########################################
################################################################################

#' Calculate pairwise spatio-temporal overlaps (lite version)
#' @description This internal function calculates the pairwise spatio-temporal overlaps
#' (association index) between two individuals. It is a streamlined adaptation of the internal
#' pairwiseOverlap function used within the calculateOverlap function. The function assumes
#' that joint space usage occurs whenever individuals overlap in space (receiver) and time (time bin).

#' @description This internal function calculates the pairwise spatio-temporal overlaps
#' (association index) between two individuals. It is a streamlined adaptation of the internal
#' `pairwiseOverlap` function used within the \code{\link{calculateOverlap}} function.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal

pairwiseOverlapLite <- function(p, table, pairwise_combinations, complete_ids, id.groups,
                                group.comparisons, start_dates, end_dates, metric) {

  # get animal IDs
  a <- pairwise_combinations[1, p]
  b <- pairwise_combinations[2, p]
  id1 <- complete_ids[a]
  id2 <- complete_ids[b]

  # discard comparisons between individuals belonging to the same group or
  # skip comparisons between different groups, if required
  if (!is.null(id.groups)) {
    group1 <- which(unlist(lapply(id.groups, function(x) id1 %in% x)))
    group2 <- which(unlist(lapply(id.groups, function(x) id2 %in% x)))
    if (group.comparisons == "within" & group1 != group2) return(NA)
    if (group.comparisons == "between" & group1 == group2) return(NA)
  }

  # if any of the individuals doesn't have detections, jump to next pair
  if (any(is.na(end_dates[c(id1, id2)]))) return(NA)

  # if monitoring periods don't overlap, jump to next pair
  start <- max(start_dates[names(start_dates) %in% c(id1, id2)])
  end  <- min(end_dates[names(end_dates) %in% c(id1, id2)])
  window_size <- as.numeric(difftime(end , start, units = "days"))
  if (window_size <= 0) return(NA)

  # get pair detections
  Ta <- table[, id1][table$timebin >= start & table$timebin <= end]
  Tb <- table[, id2][table$timebin >= start & table$timebin <= end]
  if (length(Ta) == 0 || length(Tb) == 0) return(NA)

  # Number of matching detections
  matching_rows <- length(which(Ta == Tb))
  # number of different detections
  mismatch_rows <- length(which(Ta != Tb))
  # number of incomplete rows
  incomplete_rows <- length(which(is.na(Ta) & !is.na(Tb) | !is.na(Ta) & is.na(Tb)))

  # calculate overlap %
  if (metric == "simple-ratio") total <- matching_rows + mismatch_rows + incomplete_rows
  else if (metric == "half-weight") total <- matching_rows + mismatch_rows + incomplete_rows * 0.5
  overlap <- (matching_rows / total) * 100

  # return results
  return(overlap)
}

##################################################################################################
##################################################################################################
##################################################################################################
