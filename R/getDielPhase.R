#######################################################################################################
# Estimate diel phase at a given datetime and location #################################################
#######################################################################################################

#' Estimate diel phase

#' @description Function to retrieve diel phase (e.g. day/night) in which detections occurred.
#' Number of retrieved levels can be set by the 'phases' argument:
#'  • phases=2:  day | night
#'  • phases=3:  day| crepuscule | night
#'  • phases=4:  dawn | day | dusk | night
#' Crepuscular periods can be defined based on the different solar elevation angles:
#'  • solar.depth=6: civil twilight
#'  • solar.depth=12: nautical twilight
#'  • solar.depth=18: astronomical twilight
#'
#' @param datetime A POSIXct object containing the respective datetimes or time-bins.
#' @param coordinates A SpatialPoints or matrix object containing longitude and
#'  latitude coordinates (in that order) at which to estimate sunrise and sunset times.
#' @param phases Number of diel phases to return.
#' @param solar.depth Angle of the sun below the horizon (in degrees). Passed the
#' solarDep argument in \code{\link[suntools]{crepuscule}} function.
#'
#' @return Factor with diel phase
#' @export


getDielPhase <- function (datetime, coordinates, phases=2, solar.depth=18) {

  # check if number of phases is valid
  if(!phases %in% c(2,3,4)){stop("Number of phases should be between 2, 3 and 4")}

  # calculate sunrise and sunset times for the given coordinates
  coordinates <- matrix(coordinates, ncol=2)
  sunrise <- suntools::sunriset(coordinates, datetime, POSIXct.out=T, direction="sunrise")$time
  sunset <- suntools::sunriset(coordinates, datetime, POSIXct.out=T, direction="sunset")$time

  # return day/night results if only 2 phases are needed
  if(phases==2) {
    timeofday <- ifelse(datetime >= sunrise & datetime <= sunset, "day", "night")
    timeofday <- factor(timeofday, levels=c("day", "night"))
    return(timeofday)
  }

  # otherwise calculate dawn/dusk times
  dusk <- suntools::crepuscule(coordinates, datetime, POSIXct.out=T, solarDep=solar.depth, direction="dusk")$time
  prev_dusk <- suntools::crepuscule(coordinates, datetime-60*60*24, POSIXct.out=T, solarDep=solar.depth, direction="dusk")$time
  dawn <- suntools::crepuscule(coordinates, datetime, POSIXct.out=T, solarDep=solar.depth, direction="dawn")$time

  timeofday <- ifelse(datetime > prev_dusk & datetime < dawn, "night",
                      ifelse(datetime >= dawn & datetime <= sunrise, "dawn",
                             ifelse(datetime > sunrise & datetime < sunset, "day",
                                    ifelse(datetime >= sunset & datetime <=dusk, "dusk",
                                           ifelse(datetime > dusk, "night", NA)))))

  # if 3 phases rename dusk and dawn as crepescule
  if(phases==3) {
    timeofday[timeofday %in% c("dawn", "dusk")] <- "crepuscule"
    timeofday <- factor(timeofday, levels=c("day", "crepuscule", "night"))
  }
  if(phases==4) {
    timeofday <- factor(timeofday, levels=c("dawn", "day", "dusk", "night"))
  }

  #return result
  return(timeofday)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
