#######################################################################################################
# Calculate pairwise overlaps  ########################################################################
#######################################################################################################

#' Calculate pairwise spatio-temporal overlaps (association index)

#' @description Estimates the extent of association between each two individuals,
#' assuming that joint space usage occurs whenever individuals overlap in
#' space (receiver) and time (bin). Only the total shared periods of detection
#' for a given pair are considered for this analysis, i.e. the time bins from the
#' latest release date to the earliest last detection between each pair. This
#' truncation step ensures that the absence of an individual in the receiver array
#' is not due to transmitter failure or premature death. This pairwise overlap
#' index is analogous to the simple ratio association index (Cairns & Schwager 1987,
#' Ginsberg & Young 1992), ranging from 0% (no overlap) to 100% (complete overlap).
#'
#' @param table A data frame object containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding to the receiver with the
#' highest number of detections), as returned by \code{\link{createWideTable}}.
#' @param id.metadata A data frame containing individual's metadata. It should contain
#' at least  animal IDs  and a tagging/release date column in POSIXct format.
#' @param id.groups Optional. A list containing ID groups, used to calculate stats independently
#' within each group, as well as comparing relationships between ids of different groups.
#' @param subset If defined, overlaps are calculated independently for each level
#' of this variable (usually a temporal class, e.g., diel phase or annual season).
#' If left NULL, single overlap indices are calculated for the whole monitoring period.
#' @param id.col Name of the column in id.metadata containing animal IDs. Defaults to 'ID'.
#' @param tagdate.col Name of the column in id.metadata containing animal tagging/release dates in
#' POSIXct format. Defaults to 'tagging_date'.
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' One of "within", "between" or "all". Useful to discard comparisons between individuals
#' belonging to the same group or skip comparisons between different groups, when these are
#' not required (less computing time). Defaults to "all".
#'
#' @return A list containing similarity matrices with pairwise overlaps in %,
#' time windows in days (period used to compute overlap for each pair),
#' the average number of time-bins spent by each pair consecutively in the same receiver,
#' the maximum number of time-bins spent by each pair consecutively in the same receiver and
#' the total number of co-occurrences, as well as the number co-occurrences by location level.
#' If id.groups or a subset variable is defined, these results are nested within each of the grouping levels.
#' @export


calculateOverlap <- function(table, id.metadata, id.groups=NULL, subset=NULL,
                             id.col="ID", tagdate.col="tagging_date", group.comparisons="all") {


  ######################################################################
  ## Initial checks ####################################################

  if(!is.null(subset) & !all(subset %in% colnames(table))){
    stop("subset variable(s) not found in the supplied table")
  }

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

  if(class(id.metadata[,tagdate.col])[1]!="POSIXct"){
    stop("Please convert the tagging dates to POSIXct format")
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {
      stop("Repeated ID(s) in id.groups")
    }
    if(any(!unlist(id.groups) %in% levels(id.metadata[,id.col]))){
      stop("Some of the ID(s) in id.groups were not found in the supplied metadata")
    }

    # remove individuals not present in ID groups
    complete_ids <- unique(id.metadata[,id.col])
    selected_ids <- unique(unlist(id.groups))
    id_cols <- which(colnames(table) %in% complete_ids)
    discard_cols <- id_cols[!colnames(table)[id_cols] %in% selected_ids]
    if(length(discard_cols)>0){
      table <- table[,-discard_cols]
      id.metadata <- id.metadata[id.metadata[,id.col] %in% selected_ids,]
      id.metadata[,id.col] <- droplevels(id.metadata[,id.col])
    }
  }



  ######################################################################
  ## Prepare data ######################################################

  # measure running time
  start.time <- Sys.time()

  #
  unique_ids <- colnames(data_table)[colnames(data_table) %in% levels(id.metadata[,id.col])]


  # retrieve last detection of each individual
  keep_cols <- which(colnames(table) %in% levels(id.metadata[,id.col]))
  tagging_dates <- id.metadata[,tagdate.col]
  last_detections <- suppressWarnings(apply(table[,keep_cols], 2, function(x) table$timebin[max(which(!is.na(x)))]))
  last_detections <- as.POSIXct(last_detections, origin='1970-01-01', tz="UTC")
  missing_individuals <- names(which(is.na(last_detections)))
  keep_cols <- c(1, keep_cols)


  # retrieve animal IDs
  unique_ids <- levels(id.metadata[,id.col])
  n_ids <- length(unique_ids)

  # check if overlaps need to be calculated for different groups
  multiple <- !is.null(subset)
  data_list <- list()


  ##########################################################################
  # if not, calculate overlaps for the entire study duration ###############
  if(multiple==F){
    data_list[[1]] <- table[,keep_cols]
    names(data_list) <- "entire duration"

  # else, split data and calculate overlaps independently for each group ####
  } else {
    keep_cols <- c(keep_cols, which(colnames(table)==subset))
    data_list <- split(table[,keep_cols], f=table[,subset])
    data_list <- lapply(data_list, function(x) x[,-which(colnames(x)==subset)])
  }


  ###########################################################################
  # calculate pairwise overlaps #############################################

  # generate all pairwise combinations
  pairs <- combn(1:n_ids, 2)
  # initiate template matrix and results list
  template_matrix <- matrix(NA, n_ids, n_ids)
  colnames(template_matrix) <- unique_ids
  rownames(template_matrix) <- unique_ids
  results <- list()

  for(i in 1:length(data_list)){

    # set variables
    data <- data_list[[i]]
    cat(paste0("Calculating overlap - ", names(data_list)[i], "\n"))

    # initiate data objects to hold the results
    overlap_matrix <- proxy::as.simil(template_matrix)
    window_size <- proxy::as.simil(template_matrix)
    mean_run <- proxy::as.simil(template_matrix)
    max_run <- proxy::as.simil(template_matrix)
    matching_detections <- proxy::as.simil(template_matrix)
    station_counts <- list()
    no_data <- c()

    # set progress bar
    pb <- txtProgressBar(min=0, max=ncol(pairs), initial=0, style=3)

    # calculate pairwise overlaps
    for (p in 1:ncol(pairs)) {

      # show progress in console
      setTxtProgressBar(pb,p)

      # select transmitters
      a <- pairs[1,p]
      b <- pairs[2,p]
      id1 <- unique_ids[a]
      id2 <- unique_ids[b]

      # if any of the individuals doesn't have detections, jump to next pair
      if(any(is.na(last_detections[c(a,b)]))){window_size[p]<-0; next}

      # discard comparisons between individuals belonging to the same group or
      # skip comparisons between different groups, if required
      if(!is.null(id.groups)){
        group1 <- which(unlist(lapply(id.groups, function(x) id1 %in% x)))
        group2 <- which(unlist(lapply(id.groups, function(x) id2 %in% x)))
        if(group.comparisons=="within" & group1!=group2){window_size[p]<-0; next}
        if(group.comparisons=="between" & group1==group2){window_size[p]<-0; next}
      }

      # calculate overlap period (time window)
      id_indexes <- which(id.metadata[,id.col] %in% c(id1, id2))
      start <- max(id.metadata[id_indexes,tagdate.col])
      end <- min(last_detections[names(last_detections) %in% c(id1, id2)])
      window_size[p] <- as.numeric(difftime(end, start, units = "days"))
      # if transmitters don't overlap jump to next pair
      if(window_size[p]<=0) next

      # extract data from the overlapping time window
      Ta <- data[,id1][data$timebin>=start & data$timebin<=end]
      Tb <- data[,id2][data$timebin>=start & data$timebin<=end]
      if(length(Ta)==0 || length(Tb)==0) {no_data <-c(no_data, p); next}

      # number of matching values
      matching_rows <- length(which(Ta==Tb))
      # number of different values
      mismatch_rows <- length(which(Ta!=Tb))
      # number of incomplete rows
      incomplete_a_rows <- length(which(is.na(Ta) & !is.na(Tb)))
      incomplete_b_rows <- length(which(!is.na(Ta) & is.na(Tb)))

      # calculate overlap (%)
      total <- mismatch_rows + incomplete_a_rows + incomplete_b_rows + matching_rows
      overlap_matrix[p] <- (matching_rows / total)*100

      # calculate lengths and values of runs of equal stations
      runs <- rle(as.character(Ta)==as.character(Tb))
      mean_run[p] <- mean(runs$lengths[runs$values==T], na.rm=T)
      max_run[p] <- suppressWarnings(max(runs$lengths[runs$values==T], na.rm=T))

      # save number of matching values
      matching_detections[p] <- matching_rows

      # matching detection counts per station
      station_counts[[p]] <- table(Ta[which(Ta==Tb)])
      names(station_counts)[p] <- paste0(id1,"<->",id2)
    }

    # close progress bar
    close(pb)

    # replace meaningless values
    window_size[window_size<0] <- 0
    overlap_matrix[c(which(window_size==0), no_data)] <- NA
    mean_run[c(which(window_size==0), no_data)] <- NA
    mean_run[is.nan(mean_run)] <- 0
    max_run[c(which(window_size==0), no_data)] <- NA
    max_run[is.nan(max_run)] <- 0
    max_run[matching_detections==0] <- 0


    # add all results to a list
    results[[i]] <- list("overlap"=overlap_matrix, "window_size"=window_size, "mean_run"=mean_run,
                         "max_run"=max_run, "matching_detections"=matching_detections,
                         "station_counts"=station_counts)
  }

  # return final results
  if(multiple==T) {names(results) <- names(data_list)}
  if(multiple==F) {results <- unlist(results, recursive=F)}
  if(length(missing_individuals)>0) {cat(paste0(length(missing_individuals), " individual(s) with no detections\n"))}
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(results)
}

##################################################################################################
##################################################################################################
##################################################################################################
