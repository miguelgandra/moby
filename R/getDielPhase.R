#######################################################################################################
# Estimate diel phase at a given datetime and location #################################################
#######################################################################################################

#' Estimate diel phase
#'
#' @description Function to retrieve diel phase (e.g. day/night) for given datetimes and coordinates.\cr
#'
#' The number of retrieved levels can be set by the 'phases' argument:\cr
#'  • phases=2:  day | night\cr
#'  • phases=3:  day | crepuscule | night\cr
#'  • phases=4:  dawn | day | dusk | night\cr
#'
#' Crepuscular periods can be defined based on different solar elevation angles:\cr
#'  • solar.depth=6: civil twilight\cr
#'  • solar.depth=12: nautical twilight\cr
#'  • solar.depth=18: astronomical twilight\cr
#'
#' @param datetime A POSIXct object containing the respective datetimes or time-bins.
#' @param coordinates A SpatialPoints or matrix object containing longitude and
#'  latitude coordinates (in that order) at which to estimate sunrise and sunset times.
#' @param phases Integer indicating the number of diel phases to return (2, 3, or 4).
#' @param solar.depth Numeric value indicating the angle of the sun below the horizon (in degrees).
#'  Passed to the \code{\link[suntools]{crepuscule}} function.
#' @return A factor indicating the diel phase.
#'
#' @examples
#' datetime <- as.POSIXct("2024-05-30 12:00:00")
#' coordinates <- c(-7.997, 37.008)
#' getDielPhase(datetime, coordinates, phases=4, solar.depth=12)
#'
#' @export


getDielPhase <- function (datetime, coordinates, phases=2, solar.depth=18) {

  # validate number of phases
  if(!phases %in% c(2,3,4)) stop("Number of phases should be between 2, 3 and 4")


  # calculate sunrise and sunset times for the given coordinates
  coordinates <- matrix(coordinates, ncol=2)
  sunrise <- suntools::sunriset(coordinates, datetime, POSIXct.out=T, direction="sunrise")$time
  sunset <- suntools::sunriset(coordinates, datetime, POSIXct.out=T, direction="sunset")$time

  # directly return day/night without further calculations
  if(phases==2) {
    timeofday <- ifelse(datetime >= sunrise & datetime <= sunset, "day", "night")
    return(factor(timeofday, levels=c("day", "night")))
  }

  # otherwise calculate dawn/dusk times
  dusk <- suntools::crepuscule(coordinates, datetime, POSIXct.out=T, solarDep=solar.depth, direction="dusk")$time
  prev_dusk <- suntools::crepuscule(coordinates, datetime-60*60*24, POSIXct.out=T, solarDep=solar.depth, direction="dusk")$time
  dawn <- suntools::crepuscule(coordinates, datetime, POSIXct.out=T, solarDep=solar.depth, direction="dawn")$time

  # determine diel phase using a streamlined conditional structure
  timeofday <- ifelse(datetime > prev_dusk & datetime < dawn, "night",
                      ifelse(datetime >= dawn & datetime <= sunrise, "dawn",
                             ifelse(datetime > sunrise & datetime < sunset, "day",
                                    ifelse(datetime >= sunset & datetime <= dusk, "dusk", "night"))))

  # if 3 phases rename dusk and dawn to crepuscule
  if(phases == 3) {
    timeofday[timeofday %in% c("dawn", "dusk")] <- "crepuscule"
    return(factor(timeofday, levels=c("day", "crepuscule", "night")))
  }

  # else return the 4 phases
  if(phases == 4) {
    return(factor(timeofday, levels=c("dawn", "day", "dusk", "night")))
  }

}

#######################################################################################################
#######################################################################################################
#######################################################################################################
