#######################################################################################################
# Function to assign time-bins ########################################################################
#######################################################################################################

#' Assign time bins
#'
#' @description This function assigns time bins to a dataset based on a specified time interval,
#' enabling the grouping of detection data into fixed time periods. Time binning is particularly
#' valuable in telemetry studies, where detection events are often aggregated over regular intervals
#' (e.g., hourly or every 30 minutes).
#'
#' The function determines time bins by rounding datetime values to the nearest specified interval
#' using one of three methods: `"floor"` (round down), `"ceiling"` (round up), or `"round"` (round to the nearest).
#' It essentially serves as a wrapper for the `round_date()`, `floor_date()`, and `ceiling_date()`
#' functions from the `lubridate` package (Grolemund & Wickham, 2011).
#'
#' @param datetimes A POSIXct vector containing datetimes.
#' @param interval A character string specifying the time unit or multiple of a unit to round to. Valid base units are
#' `"second"`, `"minute"`, `"hour"`, `"day"`, `"week"`, `"month"`, `"bimonth"`, `"quarter"`, `"season"`, `"halfyear"`, and `"year"`.
#' Arbitrary unique English abbreviations as in the `lubridate::period()` constructor are also allowed. Defaults to `"30 mins"`.
#' @param rounding.method The method for assigning time bins. Options are `"floor"` (default), `"ceiling"`, or `"round"`.
#'
#'
#' @return A vector of datetime values in POSIXct format representing the assigned time bins,
#' rounded according to the specified `interval` and `rounding.method`.
#'
#' @references Grolemund, G., & Wickham, H. (2011). Dates and times made easy with lubridate. Journal of statistical software, 40, 1-25.
#'
#' @seealso \code{\link[lubridate]{round_date}}
#'
#' @examples
#' # Sample dataset with datetime column
#' data <- data.frame(
#'   id = 1:6,
#'   datetime = as.POSIXct(c(
#'     "2024-11-15 08:23:45",
#'     "2024-11-15 08:45:00",
#'     "2024-11-15 09:05:30",
#'     "2024-11-15 09:20:00",
#'     "2024-11-15 09:59:59",
#'     "2024-11-15 10:10:00"
#'   ))
#' )
#'
#' # Assign time bins using a 30-minute interval
#' \dontrun{
#' data$timebin <- getTimeBins(datetimes = data$datetime, interval = "30 mins")
#' }
#'
#' @export

getTimeBins <- function(datetimes,
                        interval = "30 mins",
                        rounding.method = "floor"){


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # validate parameters
  errors <- c()
  if(!inherits(datetimes, "POSIXct")) errors <- c(errors, "Datetimes must be provided in POSIXct format.")
  if(!rounding.method %in% c("floor", "ceiling", "round")) errors <- c(errors, "Invalid 'rounding.method' argument. Must be one of 'floor', 'ceiling', or 'round'.")
  if(!is.character(interval)) errors <- c(errors, "'interval' must be a character string specifying a valid time unit or period.")
  # print errors if any
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }


  ##############################################################################
  ## Calculate time bins #######################################################
  ##############################################################################

  # assign time bins based on the specified rounding method
  timebins <- switch(rounding.method,
    "floor" = lubridate::floor_date(datetimes, unit=interval),
    "ceiling" = lubridate::ceiling_date(datetimes, unit=interval),
    "round" = lubridate::round_date(datetimes, unit=interval)
  )

  # return results
  return(timebins)

}

#######################################################################################################
#######################################################################################################
#######################################################################################################
