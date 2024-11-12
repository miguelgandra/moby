#######################################################################################################
# Calculate detection frequencies/presence probabilities  #############################################
#######################################################################################################

#' Calculate hourly detection/presence frequencies

#' @description Estimates hourly detection frequencies / presence probabilities of each
#' individual. If a subset variable is defined, frequencies are calculated
#' independently by group.
#'
#' @param table A data frame object containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding to the number of detections),
#' as returned by \code{\link{createWideTable}} with value.col="detections".
#' @param n.individuals Number of tagged individuals.
#' @param tagging.dates A POSIXct vector containing the tag/release date of each animal.
#' @param type Either "detections" or "presences".
#' @param subset If defined, hourly frequencies are calculated independently for each level
#' of this variable. If left NULL, single detection/presence frequencies are calculated for the
#' whole monitoring period.
#' @return Hourly frequences for each individual.
#' @export


calculateFrequencies <- function(table,
                                 n.individuals,
                                 tagging.dates,
                                 type = "detections",
                                 subset = NULL) {

  # initial checks
  if(!class(table[,2]) %in% c("integer","numeric")) {stop("Supplied table needs to contain number of detections")}
  if(length(type)!=1 | !type %in% c("detections", "presences")) {stop("Wrong type argument")}

  # retrieve last detection of each individual
  out_cols <- setdiff(2:ncol(table), 1:n.individuals+1)
  last.detections <- suppressWarnings(apply(table[,-c(1,out_cols)], 2, function(x) table$timebin[max(which(!is.na(x)))]))
  last.detections <- as.POSIXct(last.detections, origin='1970-01-01', tz="UTC")
  missing_individuals <- names(which(is.na(last.detections)))

  # get time bins interval (in hours)
  interval <- unique(as.numeric(difftime(dplyr::lead(table$timebin), table$timebin, units="mins")))
  interval <- interval[!is.na(interval)]

  # retrieve animal IDs
  unique_ids <- colnames(table[,-c(1,out_cols)])

  # if type equal to "presence", convert detections to presences
  if(type=="presences") {table[,-c(1,out_cols)] <- apply(table[,-c(1,out_cols)], 2, function(x){x[x>1]<-1; return(x)})}

  # check if frequencies need to be calculated for different groups
  multiple <- !is.null(subset)
  data_list <- list()

  ##########################################################################
  # if not, calculate frequencies for the entire study duration ############
  if(!multiple){
    if(length(out_cols)>0) {
      data_list[[1]] <- table[,-out_cols]
    }else{
      data_list[[1]] <- table
    }
    names(data_list) <- "entire duration"
  }

  ###########################################################################
  # else, split data and calculate frequencies independently for each group #
  if(multiple){

    # check if subset groups are available in the supplied data
    if(any(!subset %in% colnames(table))) {stop("supplied grouping variables could not be found", call.=FALSE)}

    # split data
    table_subset <- sapply(subset, function(x) split(table, f=(table[,x])), simplify=FALSE, USE.NAMES=TRUE)
    group_names <- unlist(lapply(table_subset, names))
    table_subset <- unlist(table_subset, recursive=FALSE, use.names=TRUE)
    names(table_subset) <- group_names
    data_subset <- lapply(table_subset, function(x) x[,-out_cols])
    data_list <- data_subset
  }

  ###########################################################################
  # calculate frequences ####################################################

  results <- list()
  for(i in 1:length(data_list)){

    # set variables
    data <- data_list[[i]]
    string <- ifelse(type=="detections", "detection frequencies", "presence probabilities")
    cat(paste0("Calculating ", string, " - ", names(data_list)[i], "\n"))
    hourly_frequences <- c()

    # set progress bar
    pb <- txtProgressBar(min=0, max=n.individuals, initial=0, style=3)

    for(f in 1:n.individuals) {
      setTxtProgressBar(pb, f)
      individual_data <- data[,f+1][data$timebin>=tagging.dates[f] & data$timebin<=last.detections[f]]
      if(length(individual_data)==0) {hourly_frequences[f]<-NA; next}
      detections <- sum(individual_data, na.rm=TRUE)
      hourly_frequences[f] <- detections/(length(individual_data)*interval/60)
    }

    # close progress bar and save results
    close(pb)
    names(hourly_frequences) <- unique_ids
    results[[i]] <- hourly_frequences
  }

  # return final results
  if (multiple) {
    names(results) <- names(data_list)
  } else {
    results <- unlist(results, recursive=FALSE)
  }

  if(length(missing_individuals)>0) {cat(paste0(length(missing_individuals), " individual(s) with no detections\n"))}
  return(results)
}

##################################################################################################
##################################################################################################
##################################################################################################
