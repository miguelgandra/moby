#######################################################################################################
# Function to calculate COAs (Centers of Activity Positions) ##########################################
#######################################################################################################

#' Calculate Centers of Activity (COAs)
#'
#' @description This function calculates position estimates (Centers of Activity, COAs)
#' from presence data collected at multiple receivers. It uses weighted means of longitude
#' and latitude based on detection counts during specified time bins.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing the detections data.
#' @param sensor.cols Character vector (optional). Names of additional sensor columns
#' to include in the COA calculation. Default is `NULL`.
#'
#' @return A data frame containing the following columns:
#' \itemize{
#'   \item The unique identifier (`id.col`) for each individual.
#'   \item The time bin (`timebin.col`) for which the COA was calculated.
#'   \item Mean longitude and latitude (`lon.col` and `lat.col`) for each ID and time bin.
#'   \item The number of detections (`detections`) for each ID and time bin.
#'   \item The number of unique stations visited (`station.col`) for each ID and time bin.
#'   \item Mean values for any additional sensor columns (if `sensor.cols` is provided).
#' }
#'
#' @export

calculateCOAs <- function(data,
                          id.col = getDefaults("ID"),
                          timebin.col = getDefaults("timebin"),
                          lon.col = getDefaults("lon"),
                          lat.col = getDefaults("lat"),
                          station.col = NULL,
                          sensor.cols = NULL){

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # check if sensor columns exist in the data (if provided)
  if(!is.null(sensor.cols)) {
    missing_sensors <- setdiff(sensor.cols, colnames(data))
    if(length(missing_sensors)>0) {
      stop("The following sensor columns are missing from the data: ", paste(missing_sensors, collapse = ", "), call.=FALSE)
    }
  }

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
    data_list <- append(data_list, stations)
  }

  # aggregate sensor data (if provided)
  if(!is.null(sensor.cols)){
    for(s in sensor.cols){
      sensor_formula <- as.formula(paste(s, sprintf("~ %s + %s", id.col, timebin.col)))
      sensor_mean <- stats::aggregate(sensor_formula, data=data, mean, na.rm=TRUE)
      data_list <- append(data_list, sensor_mean)
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
