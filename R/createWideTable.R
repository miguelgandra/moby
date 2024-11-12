#######################################################################################################
## Plot receiver statistics ###########################################################################
#######################################################################################################

#' Create detections table in wide format

#' @description This function generates a data frame containing binned detections in a wide format
#' (time bin x individual matrix). The values in the matrix are determined by the specified column
#' in the input data frame, which can represent, for instance, the receiver where most detections were
#' registered within that time frame or the number of detections.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections with corresponding time-bins.
#' @param value.col The column name in the data frame whose values will be assigned
#'  to each entry in the resulting wide-format table. If set to "detections" and no
#'  such column exists in the dataset, the function will assume one detection per row
#'  and create a corresponding column.
#' @param start.dates Optional. A POSIXct vector containing the start date of the monitoring period of each
#' animal (e.g., tag/release date of each individual). If left null, defaults to each individual's first detection.
#' If a single value is provided, it will be used for all IDs.
#' @param end.dates Optional. A POSIXct vector containing a cut-off date for each animal.
#' If left null, defaults to the individual's last detection.
#' If a single value is provided, it will be used for all IDs.
#' @param agg.fun A function to determine how values within each time bin are aggregated.
#' If left NULL, the aggregation method is automatically chosen based on the type of `value.col`:
#' for character or factor columns (e.g., 'station'), the most frequent value (mode) is used;
#' for numeric columns (e.g., 'detections'), the sum of values is used.
#' @param round.dates Logical. If TRUE, the start and end dates are rounded to the nearest day:
#' the earliest start date is floored to the beginning of the day, and the latest end date is
#' rounded up to the end of the day. This ensures that the time bins cover whole days. Defaults to FALSE.
#' @param verbose Output ties info to console? Defaults to TRUE.
#' @return A data frame in wide format with time bins as rows and individuals as columns. Each cell
#' in the matrix represents the value from `value.col` for that time bin and individual, aggregated
#' according to the specified aggregation function.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   timebin = as.POSIXct(c('2023-01-01 00:00:00', '2023-01-01 00:05:00', '2023-01-01 00:10:00',
#'                          '2023-01-01 00:00:00', '2023-01-01 00:10:00', '2023-01-01 00:15:00',
#'                          '2023-01-01 00:05:00', '2023-01-01 00:10:00', '2023-01-01 00:20:00')),
#'   id = c('A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C'),
#'   detections = c(2, 3, 1, 5, 1, 4, 6, 2, 3)
#' )
#'
#' start.dates <- as.POSIXct(c('2023-01-01 00:00:00', '2023-01-01 00:00:00', '2023-01-01 00:00:00'))
#' end.dates <- as.POSIXct(c('2023-01-01 00:10:00', '2023-01-01 00:15:00', '2023-01-01 00:20:00'))
#'
#' createWideTable(data, value.col="detections", start.dates=start.dates, end.dates=end.dates)
#' }

createWideTable <- function(data,
                            id.col = getDefaults("id"),
                            timebin.col = getDefaults("timebin"),
                            value.col,
                            start.dates = NULL,
                            end.dates = NULL,
                            agg.fun = NULL,
                            round.dates = FALSE,
                            verbose = TRUE) {


  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # check if data contains value.col
  if(!value.col %in% colnames(data)) {
    if(value.col=="detections"){
      warning("No 'detections' column found, assuming one detection per row", call.=FALSE)
      data$detections <- 1
    }else{
      stop("Value column not found. Please specify the correct column using 'value.col'", call.=FALSE)
    }
  }

  # check and replicate start.dates if it is a single value
  if(!is.null(start.dates) && length(start.dates)==1){
      start.dates <- rep(start.dates, nlevels(data[,id.col]))
  }

  # check and replicate end.dates if it is a single value
  if(!is.null(end.dates) && length(end.dates)==1){
      end.dates <- rep(end.dates, nlevels(data[,id.col]))
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
      data_table <- reshape2::dcast(data, formula=paste0(timebin.col, "~", id.col), value.var=value.col,
                                    fill=NA_character_, fun.aggregate=assignVal, drop=FALSE)
      ties <- reshape2::dcast(data, formula=paste0(timebin.col, "~", id.col), value.var=value.col,
                              fill=NA_character_, fun.aggregate=checkTies, drop=FALSE)
      ties$count <- apply(ties[,-1], 1, function(x) length(which(!is.na(x))))
      ties <- ties[ties$count>0,]

      if(nrow(ties)>0){
        ties <- reshape2::melt(ties[,-ncol(ties)], id.vars=timebin.col)
        colnames(ties) <- c("timebin", "ID", "ties")
        ties <- ties[!is.na(ties$ties),]
        warning(paste0(nrow(ties), " instances with value ties\n"), call.=FALSE)
        if(verbose){
          cat(paste0("Warning: ", nrow(ties), " instances with value ties\n"))
          cat("First value assigned:\n")
          if(nrow(ties)<=10){print(ties)}else{print(head(ties))}
        }
      }
    }

    # aggregate values (assign sum)
    if(value_type %in% c("numeric", "integer")){
      data_table <- reshape2::dcast(data, formula=paste0(timebin.col, "~", id.col), value.var=value.col,
                                    fill=0, fun.aggregate=sum, drop=FALSE)
    }

  # else use supplied function to aggregate values
  }else{
    data_table <- reshape2::dcast(data, formula=paste0(timebin.col, "~", id.col), value.var=value.col,
                                  fill=0, fun.aggregate=agg.fun, drop=FALSE)
}


  ##############################################################################
  # Create final table #########################################################
  ##############################################################################

  # if start.dates were not provided, retrieve first detections dates
  if(is.null(start.dates)){
    start.dates <- tapply(X=data[,timebin.col], INDEX=data[,id.col], FUN=min)
    start.dates <- as.POSIXct(start.dates, origin='1970-01-01', tz="UTC")
  }

  # if end.dates were not provided, retrieve last detections dates
  if(is.null(end.dates)){
    end.dates <- tapply(X=data[,timebin.col], INDEX=data[,id.col], FUN=max)
    end.dates <- as.POSIXct(end.dates, origin='1970-01-01', tz="UTC")
  }

  # get time bins interval (in minutes)
  interval <- difftime(data[,timebin.col], dplyr::lag(data[,timebin.col]), units="min")
  interval <- as.numeric(min(interval[interval>0], na.rm=TRUE))


  # create complete sequence of timebins for the study duration
  # round dates if needed
  if(round.dates){
    start <-  lubridate::floor_date(min(start.dates, na.rm=TRUE), unit="day")
    end <- lubridate::ceiling_date(max(end.dates, na.rm=TRUE), unit="day")-60*60*(interval/60)
  } else {
    start <- min(start.dates, na.rm=TRUE)
    end <- max(end.dates, na.rm=TRUE)
  }
  time_seq <- seq.POSIXt(from=start, to=end, by=interval*60)
  time_seq <- as.POSIXct(setdiff(time_seq, data_table[,timebin.col]), origin='1970-01-01', tz="UTC")
  empty_bins <- data.frame(time_seq, matrix(NA, ncol=nlevels(data[,id.col])))
  colnames(empty_bins) <- c(timebin.col, levels(data[,id.col]))
  # create complete matrix
  data_table <- rbind(data_table, empty_bins)
  data_table <- data_table[order(data_table[,timebin.col]), ]

  # set table entries before start date and after end date to NA
  if(value_type %in% c("numeric", "integer")){
    data_table[is.na(data_table)] <- 0
    detected_ids <- names(table(data[,id.col]))[table(data[,id.col])>0]
    missing_ids <- names(table(data[,id.col]))[table(data[,id.col])==0]

    if(length(missing_ids)>=1){data_table[,missing_ids] <- NA_integer_}
    for (i in detected_ids) {
      index <- which(levels(data[,id.col])==i)
      data_table[data_table[,timebin.col]<start.dates[index],i] <- NA_integer_
      data_table[data_table[,timebin.col]>end.dates[index],i] <- NA_integer_
    }
  }

  # discard row names
  rownames(data_table) <- NULL

  # create new attributes to save relevant variables
  attr(data_table, 'ids') <- levels(data[,id.col])
  attr(data_table, 'timebin.col') <- as.character(timebin.col)
  names(start.dates) <- levels(data[,id.col])
  names(end.dates) <- levels(data[,id.col])
  attr(data_table, 'start.dates') <- start.dates
  attr(data_table, 'end.dates') <- end.dates

  # return formatted table
  return(data_table)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
