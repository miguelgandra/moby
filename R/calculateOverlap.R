#######################################################################################################
# Calculate pairwise overlaps  ########################################################################
#######################################################################################################

#' Calculate pairwise spatio-temporal overlaps (association index)
#'
#' @description Estimates the extent of association between each pair of individuals,
#' assuming that joint space usage occurs whenever individuals overlap in space (receiver) and time (bin).
#' The analysis only considers the total shared periods of detection for a given pair, i.e., the time bins from the
#' latest release date to the earliest last detection between each pair. This truncation step ensures that the absence
#' of an individual in the receiver array is not due to transmitter failure or premature death. The pairwise overlap
#' index can be calculated using one of two metrics: the simple ratio association index
#' or the half-weight index. These indices range from 0% (no overlap) to 100% (complete overlap).
#'
#' @param table A data frame containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding to the receiver/station with the
#' highest number of detections), as returned by \code{\link{createWideTable}}.
#' @param id.groups Optional. A list containing ID groups, used to calculate stats independently
#' within each group, as well as comparing relationships between ids of different groups.
#' @param subset If defined, overlaps are calculated independently for each level of
#' this variable. This can either be a single column name (variable) or a vector
#' of column names, corresponding to variables contained in the table. In the case of
#' multiple columns, their interaction is used for grouping. If left NULL,
#' single overlap indices are calculated for the whole monitoring period.
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' One of "within", "between" or "all". Useful to discard comparisons between individuals
#' belonging to the same group or skip comparisons between different groups, when these are
#' not required (less computing time). Defaults to "all".
#' @param metric One of "simple-ratio" or "half-weight".
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations.
#' Run \code{parallel::detectCores()} to check the number of available cores.
#'
#' @details The two metrics for calculating association indices are:
#
#' 1. Simple Ratio Association Index:
#' This index provides a measure of the proportion of time two individuals spend together.
#' It is calculated as the ratio of the number of time bins where both individuals are
#' detected together to the total number of time bins where at least one of the individuals
#' is detected. This straightforward metric gives a clear indication of the extent
#' to which two individuals' space usage overlaps.
#'
#' 2. Half-Weight Index:
#' This index modifies the simple ratio by giving half-weight to the occurrences
#' where only one of the individuals is detected. It accounts for the fact that sometimes
#' individuals might be in the same area but one might not be detected due to various reasons.
#' This index is useful when detections are not perfectly reliable and helps mitigate
#' the impact of missed detections.
#'
#' It is important to acknowledge that any fish pair detected simultaneously may
#' be anywhere from a few centimetres to hundreds of meters apart, depending on
#' the maximum detection ranges of the transmitters. Therefore, caution should
#' be taken when interpreting the results, especially when making inferences about
#' biotic or fish−habitat relationships.
#'
#' Useful references:
#' \itemize{
#'   \item Cairns, S. J., & Schwager, S. J. (1987). A comparison of association indices. Animal Behaviour, 35(5), 1454-1469.
#'   \item Ginsberg, J. R., & Young, T. P. (1992). Measuring association between individuals or groups in behavioural studies. Animal Behaviour, 44(2), 377-379.
#'   \item Farine, D. R., & Whitehead, H. (2015). Constructing, conducting and interpreting animal social network analysis. Journal of Animal Ecology, 84(5), 1144-1163.
#'   \item Hoppitt, W. J., & Farine, D. R. (2018). Association indices for quantifying social relationships: how to deal with missing observations of individuals or groups. Animal Behaviour, 136, 227-238.
#' }
#'
#' @return A data frame containing the following columns:
#' \describe{
#'   \item{id1}{ID of the first individual in the pair.}
#'   \item{id2}{ID of the second individual in the pair.}
#'   \item{overlap}{Pairwise overlap percentage between the two individuals.}
#'   \item{co_occurrences}{Total number of time bins where both individuals were detected together.}
#'   \item{shared_monit_days}{Number of days the two individuals were monitored simultaneously.}
#'   \item{start}{Start date of the dyad shared monitoring period.}
#'   \item{end}{End date of the dyad shared monitoring period.}
#'   \item{mean_consec_overlap}{Average number of consecutive time bins where both individuals were detected together.}
#'   \item{max_consec_overlap}{Maximum number of consecutive time bins where both individuals were detected together.}
#' }
#' If `subset` variable is defined, these results are nested within each of the grouping levels.
#' If `id.groups` is defined, the function will add the respective groups containing each one of the IDs.
#'
#' @return A list containing similarity matrices with pairwise overlaps in %,
#' time windows in days (period used to compute overlap for each pair),
#' the average number of time-bins spent by each pair consecutively in the same receiver,
#' the maximum number of time-bins spent by each pair consecutively in the same receiver and
#' the total number of co-occurrences, as well as the number of co-occurrences by location level.
#' If id.groups or a subset variable is defined, these results are nested within each of the grouping levels.
#' @export

calculateOverlap <- function(table, id.groups=NULL, subset=NULL, metric="simple-ratio",
                             group.comparisons="all", cores=1) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # validate parameters
  errors <- c()
  if(!c("ids") %in% names(attributes(table))) errors <- c(errors, "The supplied table does not seem to be in the wide format. Please use the output of the 'createWideTable' function.")
  if (!is.data.frame(table)) errors <- c(errors, "The 'table' argument must be a data frame.")
  if (!is.null(subset) && !all(subset %in% colnames(table))) errors <- c(errors,  "Subset variable(s) not found in the supplied table.")
  if (!metric %in% c("simple-ratio", "half-weight")) errors <- c(errors, "Metric must be one of 'simple-ratio' or 'half-weight'.")
  if (!group.comparisons %in% c("all", "within", "between")) errors <- c(errors, "Group comparisons must be one of 'all', 'within' or 'between'.")
  if (!is.numeric(cores) || cores < 1 || cores %% 1 != 0) errors <- c(errors, "The 'cores' parameter must be a positive integer.")
  if (cores>1 && requireNamespace("foreach", quietly=TRUE)) errors <- c(errors, "The 'foreach' package is required for parallel computing but is not installed. Please install 'foreach' using install.packages('foreach') and try again.")
  if (cores>1 && requireNamespace("doSNOW", quietly=TRUE)) errors <- c(errors, "The 'doSNOW' package is required for parallel computing but is not installed. Please install 'doSNOW' using install.packages('doSNOW') and try again.")
  if (cores>1 && requireNamespace("parallel", quietly=TRUE)){
    errors <- c(errors, "The 'parallel' package is required for parallel computing but is not installed. Please install 'parallel' using install.packages('parallel') and try again.")
  }else if(parallel::detectCores()<cores){
    errors <- c(errors, paste("Please choose a different number of cores for parallel computing (only", parallel::detectCores(), "available)."))
  }
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
  }

  # get table attributes
  complete_ids <- as.character(attributes(table)$ids)
  timebin.col <- attributes(table)$timebin.col
  start_dates <- attributes(table)$start.dates
  end_dates <- attributes(table)$end.dates

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

  # define local %dopar%
  if(cores>1) `%dopar%` <- foreach::`%dopar%`

  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # measure running time
  start.time <- Sys.time()

  # number of animals
  n_ids <- length(complete_ids)

  # select core columns
  keep_cols <- which(colnames(table) %in% complete_ids)
  missing_individuals <- names(which(is.na(end_dates)))
  keep_cols <- c(1, keep_cols)

  # initialize list
  data_list <- list()

  # if no subset variable is defined, calculate overlaps for the entire study duration
  if (is.null(subset)) {
    data_list[[1]] <- table[, keep_cols]
    names(data_list) <- "complete monitoring duration"
  # else, split data and calculate overlaps independently for each group
  } else {
    keep_cols <- c(keep_cols, which(colnames(table) %in% subset))
    data_list <- split(table[,keep_cols], f=table[,subset])
    data_list <- lapply(data_list, function(x) x[,-which(colnames(x) %in% subset)])
  }

  ##############################################################################
  # Calculate pairwise overlaps ################################################
  ##############################################################################

  # generate all pairwise combinations
  pairwise_combinations <- combn(1:n_ids, 2)

  # initiate results list
  results <- list()

  # iterate over each data subset
  for (i in 1:length(data_list)) {

    # set variables
    data <- data_list[[i]]
    .printConsole(paste0("Calculating overlap - ", names(data_list)[i]))

    # create list of variables to export to each worker
    args <- list(pairwise_combinations=pairwise_combinations, complete_ids=complete_ids,
                 start_dates=start_dates, end_dates=end_dates,
                 id.groups=id.groups, group.comparisons=group.comparisons,
                 timebin.col=timebin.col, data=data, metric=metric)

    ###################################################################
    # calculate pairwise results using the default method (single core)
    if (cores == 1) {

      # Set progress bar
      pb <- txtProgressBar(max=ncol(pairwise_combinations), initial=0, style=3)

      # Calculate pairwise overlaps
      results[[i]] <- lapply(1:ncol(pairwise_combinations), function(p) pairwiseOverlap(p, progressbar=pb, args=args))
      results[[i]] <- plyr::rbind.fill(results[[i]])

      # Close progress bar
      close(pb)

      #############################################################
      # else use parallel computing to speed up calculations
    } else {

      # print to console
      cat(paste0("Starting parallel computation: ", cores, " cores\n"))

      # initialize cluster
      cl <- parallel::makeCluster(cores)
      doSNOW::registerDoSNOW(cl)

      # set progress bar
      pb <- txtProgressBar(max=ncol(pairwise_combinations), initial=0, style=3)
      opts <- list(progress = function(n) setTxtProgressBar(pb, n))

      # calculate pairwise overlaps (multi core)
      results[[i]] <- foreach::foreach(
        p = 1:ncol(pairwise_combinations),
        .combine = rbind,
        .options.snow = opts,
        .export = "pairwiseOverlap"
        ) %dopar% {
        pairwiseOverlap(p, progressbar = NULL, args =args)
        }


      # close progress bar and stop cluster
      close(pb)
      on.exit(parallel::stopCluster(cl))
    }
  }

  # format final results
  if (!is.null(subset)){
    for(i in 1:length(data_list)) results[[i]]$subset <- names(data_list)[i]
    results <- do.call("rbind", results)
  }else{
    results <- results[[1]]
  }

  # remove skipped pairwise combinations
  results <- results[!is.na(results$id1),]

  # remove dyads containing individuals without detections
  if (length(missing_individuals) > 0){
    results <- results[-which(results$id1 %in% missing_individuals),]
    results <- results[-which(results$id2 %in% missing_individuals),]
    cat(paste0(length(missing_individuals), " individual(s) with no detections\n"))
  }

  # create new attributes to save relevant variables
  attr(results, 'ids') <- complete_ids
  attr(results, 'id.groups') <- id.groups
  attr(results, 'metric') <- metric
  attr(results, 'subset') <- subset
  attr(results, 'group.comparisons') <- group.comparisons

  # print run time
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste("Total execution time:", sprintf("%.02f", as.numeric(time.taken)), base::units(time.taken), "\n"))

  # return results
  return(results)
}


################################################################################
# Define core pairwise overlap function ########################################
################################################################################

#' Calculate pairwise spatio-temporal overlaps
#' @description This internal function calculates the spatio-temporal overlap
#' (association index) between two individuals. The function assumes
#' that joint space usage occurs whenever individuals overlap in space
#' (receiver) and time (time bin).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

pairwiseOverlap <- function(p, progressbar=NULL, args=NULL) {

  # update progress bar
  if (!is.null(progressbar)) setTxtProgressBar(progressbar, p)

  # assign variables (needed for parallelization)
  if (!is.null(args)) list2env(args, envir = environment())

  # initialize output data frame
  result <- data.frame(matrix(ncol=9, nrow=1))
  colnames(result) <- c("id1", "id2", "overlap", "co_occurrences", "shared_monit_days", "start", "end",
                        "mean_consec_overlap", "max_consec_overlap")

  # get animal IDs
  a <- pairwise_combinations[1, p]
  b <- pairwise_combinations[2, p]
  id1 <- complete_ids[a]
  id2 <- complete_ids[b]
  result$id1 <- id1
  result$id2 <- id2

  # discard comparisons between individuals belonging to the same group or
  # skip comparisons between different groups, if required
  if (!is.null(id.groups)) {
    group1 <- which(unlist(lapply(id.groups, function(x) id1 %in% x)))
    group2 <- which(unlist(lapply(id.groups, function(x) id2 %in% x)))
    result$group1 <- names(id.groups)[group1]
    result$group2 <- names(id.groups)[group2]
    result$type <- paste(result$group1, "<->", result$group2)
    if (group.comparisons == "within" & group1 != group2) return()
    if (group.comparisons == "between" & group1 == group2) return()
  }

  # set default window size
  result$shared_monit_days <- 0

  # if any of the individuals doesn't have detections, jump to next pair
  if (any(is.na(end_dates[c(id1, id2)]))) return(result)

  # if monitoring periods don't overlap, jump to next pair
  start <- max(start_dates[names(start_dates) %in% c(id1, id2)])
  end  <- min(end_dates[names(end_dates) %in% c(id1, id2)])
  window_size <- as.numeric(difftime(end , start, units = "days"))
  if (window_size <= 0) return(result)
  else result$shared_monit_days <- round(window_size, 1)
  result$start <- lubridate::floor_date(start, "day")
  result$end <- lubridate::floor_date(end, "day")

  # get pair detections
  Ta <- data[, id1][data[, timebin.col] >= start & data[, timebin.col] <= end]
  Tb <- data[, id2][data[, timebin.col] >= start & data[, timebin.col] <= end]
  if (length(Ta) == 0 || length(Tb) == 0) return(result)

  # Number of matching detections
  matching_rows <- length(which(Ta == Tb))
  result$co_occurrences <- matching_rows
  # number of different detections
  mismatch_rows <- length(which(Ta != Tb))
  # number of incomplete rows
  incomplete_rows <- length(which(is.na(Ta) & !is.na(Tb) | !is.na(Ta) & is.na(Tb)))

  # calculate overlap %
  if (metric == "simple-ratio") total <- matching_rows + mismatch_rows + incomplete_rows
  else if (metric == "half-weight") total <- matching_rows + mismatch_rows + incomplete_rows * 0.5
  result$overlap <- (matching_rows / total) * 100
  result$overlap <- round(result$overlap, 2)

  # calculate lengths and values of consecutive co-occurrences
  if (result$co_occurrences > 0) {
    runs <- rle(as.character(Ta) == as.character(Tb))
    result$mean_consec_overlap <- mean(runs$lengths[runs$values==TRUE], na.rm=TRUE)
    result$max_consec_overlap <- suppressWarnings(max(runs$lengths[runs$values==TRUE], na.rm=TRUE))
  }

  # return results
  return(result)
}

##################################################################################################
##################################################################################################
##################################################################################################
