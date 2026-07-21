#######################################################################################################
## Plot receiver statistics ###########################################################################
#######################################################################################################

#' Create detections table in wide format

#' @description This function generates a data frame containing binned detections in a wide format
#' (time bin x individual matrix). The values in the matrix are determined by the specified column
#' in the input data frame, which can represent, for instance, the receiver where most detections were
#' registered within that time frame or the number of detections.
#'
#' @inheritParams as_moby
#' @param data A data frame containing animal detections and including a time bin column
#' (as specified by the `timebin.col` argument). Time bins can be created using the
#' \code{\link{getTimeBins}} function.
#' @param value.col The column name in the data frame whose values will be assigned
#' to each entry in the resulting wide-format table. If set to "detections" and no
#' such column exists in the dataset, the function will assume one detection per row
#' and create a corresponding column.
#' @param start.dates Optional. A POSIXct vector containing the start date of the monitoring period of each
#' animal (e.g., tag/release date of each individual). This parameter must be either:
#' - A single POSIXct value, which will be applied to all unique animal IDs; or
#' - A named POSIXct vector, where the names correspond to the animal IDs in the `id.col` column.
#' If multiple tag durations are provided, the vector must include all IDs and will be reordered to align with the levels of `id.col`.
#' If left null, defaults to each individual's first detection.
#' @param end.dates Optional. A POSIXct vector containing a cut-off date for each animal.
#' This parameter must be either:
#' - A single POSIXct value, which will be applied to all unique animal IDs; or
#' - A named POSIXct vector, where the names correspond to the animal IDs in the `id.col` column.
#' If multiple tag durations are provided, the vector must include all IDs and will be reordered to align with the levels of `id.col`.
#' If left null, defaults to the individual's last detection.
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
#' data <- data.frame(
#'   timebin = as.POSIXct(c('2023-01-01 00:00:00', '2023-01-01 00:05:00', '2023-01-01 00:10:00',
#'                          '2023-01-01 00:00:00', '2023-01-01 00:10:00', '2023-01-01 00:15:00',
#'                          '2023-01-01 00:05:00', '2023-01-01 00:10:00', '2023-01-01 00:20:00'),
#'                        tz = "UTC"),
#'   id = c('A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C'),
#'   detections = c(2, 3, 1, 5, 1, 4, 6, 2, 3)
#' )
#'
#' start.dates <- as.POSIXct(rep('2023-01-01 00:00:00', 3), tz = "UTC")
#' end.dates <- as.POSIXct(c('2023-01-01 00:10:00', '2023-01-01 00:15:00', '2023-01-01 00:20:00'),
#'                         tz = "UTC")
#'
#' createWideTable(data, id.col = "id", timebin.col = "timebin", value.col = "detections",
#'                 start.dates = start.dates, end.dates = end.dates)

createWideTable <- function(data,
                            id.col = NULL,
                            timebin.col = NULL,
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
  start.dates <- reviewed_params$start.dates
  end.dates <- reviewed_params$end.dates


  # check if data contains value.col
  if(!value.col %in% colnames(data)) {
    if(value.col=="detections"){
      warning("- No 'detections' column found, assuming one detection per row", call.=FALSE)
      data$detections <- 1
    }else{
      stop("Value column not found. Please specify the correct column using 'value.col'", call.=FALSE)
    }
  }


  ##############################################################################
  # Declare function to check ties #############################################
  ##############################################################################

  checkTies <- function(x){
    mode <- max(table(x))
    if(length(which(table(x)==mode))>1){
      tb <- table(x)[table(x)==mode]
      counts <- data.frame(category=names(tb), value=as.integer(tb), stringsAsFactors=FALSE)
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
      data_table <- .castWide(data, timebin.col, id.col, value.col, fun.aggregate=assignVal, fill=NA_character_)
      ties <- .castWide(data, timebin.col, id.col, value.col, fun.aggregate=checkTies, fill=NA_character_)
      ties$count <- rowSums(!is.na(as.matrix(ties[, -1])))
      if(nrow(ties)>0){
        wide <- ties[, -ncol(ties)]
        idcols <- setdiff(names(wide), timebin.col)
        ties <- data.frame(wide[[timebin.col]][rep(seq_len(nrow(wide)), length(idcols))],
                           rep(idcols, each=nrow(wide)),
                           unlist(wide[idcols], use.names=FALSE),
                           stringsAsFactors=FALSE)
        colnames(ties) <- c("timebin", "ID", "ties")
        ties <- ties[!is.na(ties$ties),]
        if(nrow(ties)>0){
          .mobyWarn(nrow(ties), " (ID, time-bin) combination(s) had multiple differing values; ",
                    "the first was kept. Aggregate upstream (e.g. calculateCOAs) to control this.")
          if(verbose){
            .mobyInform("Tied (ID, time-bin) instances (first value kept):", verbose = verbose)
            if(nrow(ties)<=10){print(ties)}else{print(head(ties))}
          }
        }
      }
    }

    # aggregate values (assign sum)
    if(value_type %in% c("numeric", "integer")){
      data_table <- .castWide(data, timebin.col, id.col, value.col, fun.aggregate=sum, fill=0)
    }

  # else use supplied function to aggregate values
  }else{
    data_table <- .castWide(data, timebin.col, id.col, value.col, fun.aggregate=agg.fun, fill=0)
}


  ##############################################################################
  # Create final table #########################################################
  ##############################################################################

  # timezone of the supplied time bins, preserved across numeric round-trips below
  tz <- .dataTZ(data[,timebin.col])

  # if start.dates were not provided, retrieve first detections dates
  if(is.null(start.dates)){
    start.dates <- tapply(X=data[,timebin.col], INDEX=data[,id.col], FUN=min)
    start.dates <- as.POSIXct(start.dates, origin='1970-01-01', tz=tz)
  }

  # if end.dates were not provided, retrieve last detections dates
  if(is.null(end.dates)){
    end.dates <- tapply(X=data[,timebin.col], INDEX=data[,id.col], FUN=max)
    end.dates <- as.POSIXct(end.dates, origin='1970-01-01', tz=tz)
  }

  # get time bins interval (in minutes)
  interval <- difftime(data[,timebin.col], .lag(data[,timebin.col]), units="min")
  interval <- suppressWarnings(as.numeric(min(interval[interval>0], na.rm=TRUE)))
  if(!is.finite(interval)) {
    stop(paste0("Could not determine the time-bin interval: at least two distinct time bins are ",
                "required in '", timebin.col, "'. Please check the input data."), call.=FALSE)
  }


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
  time_seq <- as.POSIXct(setdiff(time_seq, data_table[,timebin.col]), origin='1970-01-01', tz=tz)
  # append empty (all-NA) rows only for time bins not already present; skip when the
  # observed bins already cover the whole sequence (otherwise the row counts mismatch)
  if(length(time_seq) > 0){
    empty_bins <- data.frame(time_seq, matrix(NA, nrow=length(time_seq), ncol=nlevels(data[,id.col])))
    colnames(empty_bins) <- c(timebin.col, levels(data[,id.col]))
    data_table <- rbind(data_table, empty_bins)
  }
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
