#######################################################################################################
# Randomize overlaps  #################################################################################
#######################################################################################################

#' Test significance of estimated overlaps through null model randomization tests

#' @description Tests the null hypothesis of temporally independent space use
#' (i.e. each animal occurs independently of the other) using Monte Carlo permutation tests.
#' The algorithm permutes entries within each column, so that the total number of detections
#' of each individual and the relative occurrence frequencies across receivers are kept unchanged.
#'
#' @param table A data frame object containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding to the receiver with the
#' highest number of detections), as returned by \code{\link{createWideTable}}.
#' @param overlaps Similarity matrix containing pairwise overlaps, as returned by \code{\link{calculateOverlap}}.
#' @param id.metadata A data frame containing individual's metadata. It should contain
#' at least  animal IDs  and a tagging/release date column in POSIXct format.
#' @param contraint.by Variable or variables used to constraint permutations, e.g.
#' to account for potential diel or seasonal trends in animal occurrences. If supplied,
#' permutations across time bins will be restricted to the same combination of
#' levels of these variables.
#' @param id.groups Optional. A list containing ID groups, used to calculate stats independently
#' within each group, as well as comparing relationships between ids of different groups.
#' @param id.col Name of the column in id.metadata containing animal IDs. Defaults to 'ID'.
#' @param tagdate.col Name of the column in id.metadata containing animal tagging/release dates in
#' POSIXct format. Defaults to 'tagging_date'.
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' One of "within", "between" or "all". Useful to discard comparisons between individuals
#' belonging to the same group or skip comparisons between different groups, when these are
#' not required (less computing time). Defaults to "all".
#' @param iterations Number of Monte Carlo iterations (simulated datasets).
#' @param conf.level Confidence level. Defaults to 95%.
#' @param cores Number of CPU cores for parallel computing (increase processing power). Defaults to 1.
#' @param random.seed Set the seed for a reproducible randomization.
#' @importFrom foreach %dopar%
#' @export


randomizeOverlaps <- function(table, overlaps, id.metadata, contraint.by=NULL, id.groups=NULL,
                              id.col="ID", tagdate.col="tagging_date", group.comparisons="all",
                              iterations=100, conf.level=0.95, cores=1, random.seed=NULL) {

  ######################################################################
  ## Initial checks ####################################################

  if(!id.col %in% colnames(id.metadata)){
    stop("ID column not found in id.metadata. Please specify the correct column using 'id.col'")
  }

  if(!tagdate.col %in% colnames(id.metadata)){
    stop("Tagging date column not found in id.metadata. Please specify the correct column using 'tagdate.col'")
  }

  if(class(id.metadata[,id.col])!="factor"){
    cat("Converting ids to factor\n")
    id.metadata[,id.col] <- as.factor(id.metadata[,id.col])
  }

  if(!is.null(contraint.by) & !all(contraint.by %in% colnames(table))){
    stop("contraint.by variable(s) not found in the supplied table")
  }

  if(!group.comparisons %in% c("all", "within", "between")){
    stop("Wrong group.comparisons argument, please select one of: 'within', 'between' or 'all'")
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {
      stop("Repeated ID(s) in id.groups")
    }
    if(any(!unlist(id.groups) %in% levels(id.metadata[,id.col]))){
      stop("Some of the ID(s) in id.groups were not found in the supplied metadata")
    }
  }


  # check parallel computing
  if(cores>1){
    if(parallel::detectCores()<cores){
      stop(paste("Please choose a different number of cores for parallel computing (only", parallel::detectCores(), "available)"))
    }
  }



  ####################################################################
  ## Prepare data ####################################################

  # measure running time
  start.time <- Sys.time()

  # set random seed for reproducibility
  if(!is.null(random.seed)){
    set.seed(random.seed)
  }

  # extract overlap matrix (if not directly supplied)
  if(class(overlaps)=="list"){overlaps<-overlaps$overlap}

  # retrieve number of animals
  n_individuals <- unique(dim(overlaps))
  unique_ids <- names(overlaps)
  animal_cols <- which(colnames(table) %in% unique_ids)

  # retrieve last detections
  last.detections <- suppressWarnings(apply(table[,animal_cols], 2, function(x) table$timebin[max(which(!is.na(x)))]))
  last.detections <- as.POSIXct(last.detections, origin='1970-01-01', tz="UTC")

  # create constraint variable
  if(!is.null(contraint.by)){
    if(length(contraint.by)==1){table$subset <- table[,contraint.by]}
    if(length(contraint.by)>1){table$subset <- apply(table[,contraint.by], 1 , paste, collapse="//")}
  }else{
    table$subset <- 1
  }


  ####################################################################
  ## Calculate overlaps ##############################################

  # initiate results list
  random_results <- list()

  #############################################################
  # calculate results using the default method (single core
  if(cores==1){

    # print message to console
    cat("Randomizing overlaps...\n")

    # set progress bar
    pb <- txtProgressBar(min=0, max=iterations, initial=0, style=3)

    # start iterations
    for (i in 1:iterations) {
      # randomize data table
      randomized_table <- randomize(table, id.cols=animal_cols)
      # calculate overlap
      random_results[[i]]  <- apply(combn(n_individuals, 2), 2, calculateOverlap2, data=randomized_table,
                                    id.metadata, id.col, tagdate.col, last.detections, id.groups, group.comparisons)
      # show progress in console
      setTxtProgressBar(pb, i)
    }

    # close progress bar
    close(pb)

  #############################################################
  # else use parallel computing to speed up calculations
  }else{
    cat(paste0("Starting parallel computation: ", cores, " cores\n"))
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    cat("Randomizing overlaps...\n")

    random_results <- foreach::foreach(i=1:iterations, .export=c("randomize", "calculateOverlap2")) %dopar% {
      # randomize data table
      randomized_table <- randomize(table, id.cols=animal_cols)

      # calculate overlap
      apply(combn(n_individuals, 2), 2, calculateOverlap2, data=randomized_table,
            id.metadata, id.col, tagdate.col, last.detections, id.groups, group.comparisons)
    }
    #stop cluster
    parallel::stopCluster(cl)
  }



  #######################################################################################################
  # Test pairwise overlap significance ##################################################################

  # pairwise overlap significance
  pair_combn <- combn(n_individuals, 2)
  n_pairs <- ncol(pair_combn)
  unique_pairs <- apply(pair_combn, 2, function(x) paste(unique_ids[x], collapse="-"))
  random_pair_overlap <- list()
  pairwise_stats <- list()
  sign_tails <- (1-conf.level)/2

  for (i in 1:n_pairs) {
    random_pair_overlap[[i]] <- unlist(lapply(random_results,  function(x) as.numeric(x)[i]))
    if(is.na(overlaps[i])) {pairwise_stats[[i]]<-NA; next}
    if(all(is.na(random_pair_overlap[[i]]))) {pairwise_stats[[i]]<-NA; next}
    lower_bound <- quantile(random_pair_overlap[[i]], probs=sign_tails)
    upper_bound <- quantile(random_pair_overlap[[i]], probs=1-sign_tails)
    pairwise_stats[[i]] <- ifelse(overlaps[i]<lower_bound, "-", ifelse(overlaps[i]>upper_bound, "+", "ns"))
    if(is.nan(pairwise_stats[[i]])){pairwise_stats[[i]]<-NA}
  }

  # summary stats
  pairwise_stats <- sapply(pairwise_stats, function(x) ifelse(is.null(x), NA, x), simplify=F)
  pairwise_stats <- data.frame("pair"=1:n_pairs, "ids"=unique_pairs, "overlap"=as.numeric(overlaps), "p.val"=unlist(pairwise_stats))


  #######################################################################################################
  # Generate summary table ##############################################################################

  pairwise_summary <- pairwise_stats
  pairwise_summary$id1 <- unique_ids[combn(n_individuals, 2)[1,]]
  pairwise_summary$id2 <- unique_ids[combn(n_individuals, 2)[2,]]

  if(!is.null(id.groups)){
    pairwise_summary$group1 <- sapply(pairwise_summary$id1, function(x) which(unlist(lapply(id.groups, function(y) x %in% y))))
    pairwise_summary$group1 <- names(id.groups)[pairwise_summary$group1]
    pairwise_summary$group2 <- sapply(pairwise_summary$id2, function(x) which(unlist(lapply(id.groups, function(y) x %in% y))))
    pairwise_summary$group2 <- names(id.groups)[pairwise_summary$group2]
    pairwise_summary$type <- paste(pairwise_summary$group1, "<->", pairwise_summary$group2)
  }else{
    pairwise_summary$type <- "All"
  }

  summary_table <- aggregate(pairwise_summary$overlap, by=list(pairwise_summary$type), mean, na.rm=T)
  summary_table$x <- sprintf("%.2f", summary_table$x)
  colnames(summary_table) <- c("Type", "Mean overlap (%)")
  summary_table$error <- aggregate(pairwise_summary$overlap, by=list(pairwise_summary$type), function(x) plotrix::std.error(x, na.rm=T))$x
  summary_table$error <- sprintf("%.2f", summary_table$error)
  summary_table$`Mean overlap (%)` <- paste(summary_table$`Mean overlap (%)`, "±", summary_table$error)
  summary_table <- summary_table[,-which(colnames(summary_table)=="error")]

  pairwise_summary <- as.data.frame.matrix(table(pairwise_summary$type, pairwise_summary$p.val))
  pairwise_summary$Type <- rownames(pairwise_summary)
  rownames(pairwise_summary) <- NULL
  summary_table <- plyr::join(summary_table, pairwise_summary, by="Type", type="left")

  count_cols <- which(colnames(summary_table) %in% c("+", "-", "ns"))
  summary_table <- summary_table[rowSums(summary_table[,count_cols])>0,]


  #######################################################################################################
  # Return results ######################################################################################

  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(list("summary"=summary_table, "randomized_overlaps"=random_results, "pairwise_significance"=pairwise_stats))

}


################################################################################
# Helper function I - randomize detections #####################################

randomize <- function(data, id.cols) {

  sampleCols <- function(x){
    if(nrow(x)<2){return(x)}
    else{return(data.frame("timebin"=x[,1], apply(x[,-1], 2, sample), check.names=F))}
  }

  data_subset <- split(data, f=data$subset)
  data_subset <- lapply(data_subset, function(x) x[,c(1,id.cols)])
  data_random <- lapply(data_subset, sampleCols)
  data_random <- dplyr::bind_rows(data_random)
  data_random <- data_random[order(data_random$timebin),]
  return(as.data.frame(data_random))
}


################################################################################
# Helper function II - simplified overlap function (for better performance) ####

calculateOverlap2 <- function(x, data, id.metadata, id.col, tagdate.col, last.detections,
                              id.groups, group.comparisons) {

  a <- x[1]
  b <- x[2]
  unique_ids <- levels(id.metadata[,id.col])
  id1 <- unique_ids[a]
  id2 <- unique_ids[b]

  # if any of the individuals doesn't have detections, jump to next pair
  if(any(is.na(last.detections[c(id1,id2)]))){return(NA)}

  # discard comparisons between individuals belonging to the same group or
  # skip comparisons between different groups, if required
  if(!is.null(id.groups)){
    group1 <- which(unlist(lapply(id.groups, function(x) id1 %in% x)))
    group2 <- which(unlist(lapply(id.groups, function(x) id2 %in% x)))
    if(group.comparisons=="within" & group1!=group2){return(NA)}
    if(group.comparisons=="between" & group1==group2){return(NA)}
  }

  # if transmitters don't overlap jump to next pair
  id_indexes <- which(id.metadata[,id.col] %in% c(id1, id2))
  start <- max(id.metadata[id_indexes, tagdate.col])
  end <- min(last.detections[names(last.detections) %in% c(id1, id2)])
  if(as.numeric(difftime(end, start, units = "days"))<=0) {return(NA)}

  # get data
  Ta <- data[,id1][data$timebin>=start & data$timebin<=end]
  Tb <- data[,id2][data$timebin>=start & data$timebin<=end]
  if(length(Ta)==0 || length(Tb)==0) {return(NA)}

  # number of matching detections
  matching_rows <- length(which(Ta==Tb))
  # number of different detections
  mismatch_rows <- length(which(Ta!=Tb))
  # number of incomplete rows
  incomplete_rows <- length(which(is.na(Ta) & !is.na(Tb) | !is.na(Ta) & is.na(Tb)))

  # calculate overlap %
  overlap <- matching_rows / (mismatch_rows + incomplete_rows + matching_rows)*100
  if(is.nan(overlap)){overlap <- 0}
  return(overlap)
}


##################################################################################################
##################################################################################################
##################################################################################################
