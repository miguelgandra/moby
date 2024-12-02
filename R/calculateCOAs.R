#######################################################################################################
# Function to calculate COAs (Centers of Activity Positions) ##########################################
#######################################################################################################

#' Calculate Centers of Activity (COAs)
#'
#' @description This function calculates position estimates (Centers of Activity, COAs)
#' from presence data collected at multiple receivers. It uses weighted means of longitude
#' and latitude based on detection counts during specified time bins. Additionally, it aggregates
#' all remaining columns dynamically. For numeric columns, the mean is calculated, while for
#' character or factor columns, unique values are concatenated and separated by "|".
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections and including a time bin column
#' (as specified by the `timebin.col` argument). Time bins can be created using the
#' \code{\link{getTimeBins}} function.
#'
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item The unique identifier (`id.col`) for each individual.
#'   \item The time bin (`timebin.col`) for which the COA was calculated.
#'   \item Mean longitude and latitude (`lon.col` and `lat.col`) for each ID and time bin.
#'   \item The number of detections (`detections`) for each ID and time bin.
#'   \item The number of unique stations visited (`stations`) for each ID and time bin.
#'   \item For numeric columns: Mean values for each ID and time bin.
#'   \item For character or factor columns: Concatenated unique values (separated by "|") for each ID and time bin.
#' }
#'
#' @export

calculateCOAs <- function(data,
                          id.col = getDefaults("ID"),
                          timebin.col = getDefaults("timebin"),
                          lon.col = getDefaults("lon"),
                          lat.col = getDefaults("lat"),
                          station.col = NULL){

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # add a column to count individual detections
  data$detections <- 1


  ##############################################################################
  ## Aggregate data by ID and time bin  ########################################
  ##############################################################################

  # calculate mean longitude and latitude for each ID and time bin
  formula1 <- as.formula(paste(sprintf("cbind(%s, %s)", lon.col, lat.col), sprintf("~ %s + %s", id.col, timebin.col)))
  coas <- stats::aggregate(formula1, data=data, mean, na.rm=TRUE)
  # count the number of detections for each ID and time bin
  formula2 <- as.formula(paste("detections", sprintf("~ %s + %s", id.col, timebin.col)))
  detections <- stats::aggregate(formula2, data=data, length)

  # list data frames to merge
  data_list <- list(coas, detections)

  # count the number of unique stations visited per time bin
  if(!is.null(station.col)){
    formula3 <- as.formula(paste(station.col, sprintf("~ %s + %s", id.col, timebin.col)))
    stations <- stats::aggregate(formula3, data=data, function(x) length(unique(x)))
    colnames(stations)[3] <- "stations"
    data_list <- c(data_list, list(stations))
  }

  # dynamically aggregate all remaining columns (excepting those of class POSIXct)
  non_posix_cols <- setdiff(colnames(data), colnames(data)[sapply(data, inherits, what="POSIXct")])
  remaining_cols <- setdiff(non_posix_cols, c(id.col, timebin.col, lon.col, lat.col, "detections"))
  if (length(remaining_cols) > 0) {
    for (col in remaining_cols) {
      if (is.numeric(data[[col]])) {
        # aggregate numeric columns by mean
        formula_numeric <- as.formula(paste(col, sprintf("~ %s + %s", id.col, timebin.col)))
        numeric_mean <- stats::aggregate(formula_numeric, data = data, mean, na.rm = TRUE)
        data_list <- c(data_list, list(numeric_mean))
      } else {
        # aggregate character or factor columns by concatenating unique values
        formula_char <- as.formula(paste(col, sprintf("~ %s + %s", id.col, timebin.col)))
        char_concat <- stats::aggregate(formula_char, data = data, function(x) paste(unique(x), collapse = "|"))
        data_list <- c(data_list, list(char_concat))
      }
    }
  }

  # merge all data frames using Reduce and plyr::join
  result <- Reduce(function(x, y) plyr::join(x, y, by=c(id.col, timebin.col), type="left"), data_list)

  # return results
  return(result)

}


#######################################################################################################
#######################################################################################################
#######################################################################################################
