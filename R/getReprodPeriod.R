#######################################################################################################
# Function to determine spawning/resting season ######################################################
#######################################################################################################

#' Assign reproductive status
#'
#' @description Assign spawning vs resting periods based on given start and end interval.
#' @param datetimes A POSIXct object containing the respective datetimes or time-bins.
#' @param spawning.start Start of the spawning season. Can be supplied as a
#' POSIXct object or as month (e.g., "08" for August).
#' @param spawning.end End of the spawning season. Can be supplied as a
#' POSIXct object or as month (e.g., "09" for September).
#' @param format Character string giving the datetimes-time format of the supplied
#' spawning.start and spawning.end. See \code{\link[base]{strptime}} details for
#' further info on available formats. Defaults to "%m" (month).
#' @param tz A character string specifying the time zone to be used for the conversion.
#' Defaults to "UTC".
#' @return A factor indicating the reproductive state (resting vs spawning).
#'
#' @examples
#' # Using integer month format (spawning period between June and September)
#' datetimes <- as.POSIXct("2024-05-30")
#' getReprodPeriod(datetimes, spawning.start="05", spawning.end="09", format="%m")
#'
#' # Using abbreviated month format (spawning period between June and September)
#' datetimes <- as.POSIXct("2024-05-30")
#' getReprodPeriod(datetimes, spawning.start="Jun", spawning.end="Sep", format="%b")
#'
#' # Using day/month format (spawning period between 15th August and 31th August)
#' datetimes <- as.POSIXct("2024-09-01")
#' getReprodPeriod(datetimes, "15/08", "31/08", format="%d/%m")
#'
#' # Using POSIXct objects
#' datetimes <- as.POSIXct("2024-05-30")
#' spawning.start <- as.POSIXct("2024-04-01")
#' spawning.end <- as.POSIXct("2024-09-30")
#' getReprodPeriod(datetimes, spawning.start, spawning.end, format="%Y-%m-%d")
#'
#'@export


getReprodPeriod <- function(datetimes, spawning.start, spawning.end, format="%m", tz="UTC") {

  # extract year from the input datetimes
  year <- strftime(datetimes, "%Y", tz=tz)

  # add dummy year if missing in the format
  if(!grepl("%y|%Y", format, fixed=T)){
    format <- paste0(format, "/%Y")
    spawning.start <- paste0(spawning.start, "/", year)
    spawning.end <- paste0(spawning.end, "/", year)
  }

  # add day of month if missing in the format
  if(!grepl("%d", format, fixed=T)){
    format <- paste0("%d/", format)
    spawning.start <- paste0("01/", spawning.start)
    end.month <- as.POSIXct(paste0("01/", spawning.end), format, tz=tz)
    month.days <-  as.integer(lubridate::days_in_month(end.month))
    spawning.end <- paste0(month.days, "/", spawning.end)
  }

  spawning.start <- as.POSIXct(spawning.start, format, tz=tz)
  spawning.end <- as.POSIXct(spawning.end, format, tz=tz)
  reprod_period <- ifelse(datetimes >= spawning.start & datetimes <= spawning.end, "spawning", "resting")
  reprod_period <- factor(reprod_period, levels=c("resting", "spawning"))
  return(reprod_period)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
