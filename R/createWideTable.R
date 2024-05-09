#######################################################################################################
## Plot receiver statistics ###########################################################################
#######################################################################################################

#' Create detections table in wide format.

#' @description Returns a data frame object containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding either to the receiver > detections
#' or the number of detections).
#'
#' @param data A data frame containing animal detections with corresponding time-bins
#' (must contain a 'timebin' column).
#' @param value.col The values to be assigned to each entry, should match a column
#' in the supplied data (except if set to "detections"). If set to "detections" and no detections
#' column is found in the dataset, the function automatically assumes one detection per row.
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param start.dates Optional. A POSIXct vector containing the start date of the monitoring period of each
#' animal (e.g.tag/release date of each individual). If left null, defaults to each individual's first detection.
#' If a single value is provided, it will be used to all IDs.
#' @param end.dates Optional. A POSIXct vector containing the tag/release date
#' of each animal.If left null, defaults to the individual's last detections.
#' If a single value is provided, it will be used to all IDs.
#' @param agg.fun Function to assign the final value across each time bin.
#' By default (if left null), it automatically adjusts based on the class type of the value.col.
#' If this column is of type character or factor (e.g., 'station' or 'habitat'),
#' the value is assigned based on the most frequently occurring observation within each timebin (i.e., mode).
#' If it is numeric (e.g., number of detections), the assigned value will correspond to
#' the sum of all observations within each timebin.
#' @param round.dates Boolean to indicate if start and end dates should be rounded
#' (floor of earliest tagging and ceiling of last detection).
#' @export


createWideTable <- function(data, value.col, id.col="ID", start.dates=NULL, end.dates=NULL,
                            agg.fun=NULL, round.dates=F) {


  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  if(!id.col %in% colnames(data)) {
    stop("'id.col' variable not found in the supplied data")
  }

  if(!value.col %in% colnames(data)) {
    if(value.col=="detections"){
      cat("Warning: No 'detections' column found, assuming one detection per row\n")
      data$detections <- 1
    }else{
      stop("'value.col' variable not found in the supplied data")
    }
  }

  if(class(data[,id.col])!="factor") {
    data[,id.col] <- as.factor(data[,id.col])
    cat("Warning: 'id.col' converted to factor\n")
  }

  if(!c("timebin") %in% colnames(data)) {
    stop("'timebin' column not found in the supplied data")
  }

  if(!is.null(start.dates)){
    if(!grepl("POSIXct", paste(class(start.dates), collapse=" "))){
      stop("'start.dates' must be provided in POSIXct format")
    }
    if(length(start.dates)>1 & length(start.dates)!=nlevels(data[,id.col])){
      stop("Incorrect number of start.dates. Must be either a single value or
           a vector containing a start date for each individual")
    }
    if(length(start.dates)==1){
      start.dates <- rep(start.dates, nlevels(data[,id.col]))
    }
  }

  if(!is.null(end.dates)){
    if(!grepl("POSIXct", paste(class(end.dates), collapse=" "))){
      stop("'end.dates' must be provided in POSIXct format")
    }
    if(length(end.dates)>1 & length(end.dates)!=nlevels(data[,id.col])){
      stop("Incorrect number of end.dates. Must be either a single value or
             a vector containing a end date for each individual")
    }
    if(length(end.dates)==1){
      end.dates <- rep(end.dates, nlevels(data[,id.col]))
    }
  }

  ##############################################################################
  # Declare function to check ties #############################################
  ##############################################################################

  checkTies <- function(x){
    mode <- max(table(x))
    if(length(which(table(x)==mode))>1){
      counts <- reshape2::melt(table(x)[table(x)==mode])
      counts$value <- paste0(counts[,1], " (", counts[,2], ")")
      return(paste(counts$value, collapse=" | "))
    }else{NA_character_}
  }

  ##############################################################################
  # Aggregate values per time bin ##############################################
  ##############################################################################

  # check class type
  value_type <- class(data[,value.col])

  # if no agg.fun has been defined, proceed with the default behaviour
  if(is.null(agg.fun)){

    # aggregate values (assign most frequent one)
    if(value_type %in% c("character", "factor")){
      assignVal <- function(x){names(which.max(table(x)))}
      data_table <- reshape2::dcast(data, formula=paste0("timebin~", id.col), value.var=value.col,
                                    fill=NA_character_, fun.aggregate=assignVal, drop=F)
      ties <- reshape2::dcast(data, formula=paste0("timebin~", id.col), value.var=value.col,
                              fill=NA_character_, fun.aggregate=checkTies, drop=F)
      ties$count <- apply(ties[,-1], 1, function(x) length(which(!is.na(x))))
      ties <- ties[ties$count>0,]
      if(nrow(ties)>0){
        ties <- reshape2::melt(ties[,-ncol(ties)], id.vars="timebin")
        colnames(ties) <- c("timebin", "ID", "ties")
        ties <- ties[!is.na(ties$ties),]
        cat(paste0("Warning: ", nrow(ties), " instances with value ties\n"))
        cat("First value assigned:\n")
        if(nrow(ties)<=10){print(ties)}else{print(head(ties))}
      }
    }

    # aggregate values (assign sum)
    if(value_type %in% c("numeric", "integer")){
      data_table <- reshape2::dcast(data, formula=paste0("timebin~", id.col), value.var=value.col,
                                    fill=0, fun.aggregate=sum, drop=F)
    }

  # else use supplied function to aggregate values
  }else{

    data_table <- reshape2::dcast(data, formula=paste0("timebin~", id.col), value.var=value.col,
                                  fill=0, fun.aggregate=agg.fun, drop=F)
}


  ##############################################################################
  # Create final table #########################################################
  ##############################################################################

  # if start.dates were not provided, retrieve first detections dates
  if(is.null(start.dates)){
    start.dates <- tapply(X=data$timebin, INDEX=data[,id.col], FUN=min)
    start.dates <- as.POSIXct(start.dates, origin='1970-01-01', tz="UTC")
  }

  # if end.dates were not provided, retrieve last detections dates
  if(is.null(end.dates)){
    end.dates <- tapply(X=data$timebin, INDEX=data[,id.col], FUN=max)
    end.dates <- as.POSIXct(end.dates, origin='1970-01-01', tz="UTC")
  }

  # get time bins interval (in minutes)
  interval <- difftime(data$timebin, data.table::shift(data$timebin), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=T))


  # create complete sequence of timebins for the study duration
  # round dates if needed
  if(round.dates==T){
    start <-  lubridate::floor_date(min(start.dates, na.rm=T), unit="day")
    end <- lubridate::ceiling_date(max(end.dates, na.rm=T), unit="day")-60*60*(interval/60)
  } else {
    start <- min(start.dates, na.rm=T)
    end <- max(end.dates, na.rm=T)
  }
  time_seq <- seq.POSIXt(from=start, to=end, by=interval*60)
  time_seq <- as.POSIXct(setdiff(time_seq, data_table$timebin), origin='1970-01-01', tz="UTC")
  empty_bins <- data.frame(time_seq, matrix(NA, ncol=nlevels(data[,id.col])))
  colnames(empty_bins) <- c("timebin", levels(data[,id.col]))
  # create complete matrix
  data_table <- rbind(data_table, empty_bins)
  data_table <- data_table[with(data_table, order(timebin)), ]

  # set table entries before start date and after end date to NA
  if(value_type %in% c("numeric", "integer")){
    data_table[is.na(data_table)] <- 0
    detected_ids <- names(table(data[,id.col]))[table(data[,id.col])>0]
    missing_ids <- names(table(data[,id.col]))[table(data[,id.col])==0]

    if(length(missing_ids)>=1){data_table[,missing_ids] <- NA_integer_}
    for (i in detected_ids) {
      index <- which(levels(data[,id.col])==i)
      data_table[data_table$timebin<start.dates[index],i] <- NA_integer_
      data_table[data_table$timebin>end.dates[index],i] <- NA_integer_
    }
  }

  # discard row names
  rownames(data_table) <- NULL

  # return formatted table
  return(data_table)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################