#######################################################################################################
# Function to determine spawning/resting season ######################################################
#######################################################################################################

#' Assign spawning periods
#'
#' @description Assign spawning vs resting periods based on given start and end interval.
#' @param date A POSIXct object containing the respective datetimes or time-bins.
#' @param spawning.start Start of the spawning season. Can be supplied as
#' a POSIXct object or as month.
#' @param spawning.end End of the spawning season. Can be supplied as
#' a POSIXct object or as month.
#' @param format Character string giving the date-time format of the supplied
#' spawning.start and spawning.end.
#' @return Factor with reproductive season (resting vs spawning).
#' @export


getReprodPeriod <- function(date, spawning.start, spawning.end, format="%m") {

  year <- strftime(date, "%Y", tz="UTC")

  # add day of month if missing
  if(!grepl("%d", format, fixed=T)){
    spawning.start <- paste0("01/", spawning.start)
    spawning.end <- paste0("31/", spawning.end)
    format <- paste0("%d/", format)
  }

  # add dummy year
  if(!grepl("%y|%Y", format, fixed=T)){
    spawning.start <- paste0(spawning.start, "/", year)
    spawning.end <- paste0(spawning.end, "/", year)
    format <- paste0(format, "/%Y")
  }

  spawning.start <- as.POSIXct(spawning.start, format, tz="UTC")
  spawning.end <- as.POSIXct(spawning.end, format, tz="UTC")
  reprod_period <- ifelse(date >= spawning.start & date <= spawning.end, "spawning", "resting")
  reprod_period <- factor(reprod_period, levels=c("resting", "spawning"))
  return(reprod_period)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
