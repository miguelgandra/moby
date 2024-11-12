#######################################################################################################
# Function to get table of diel phase times (in hours) ################################################
#######################################################################################################

#' Retrieve diel phase boundary times (hours)

#' @description Returns a table containing average sunrise, sunset and twilight times
#' for a given period (study duration).
#' '
#' @param sunriset.coords A SpatialPoints, matrix or numeric object containing longitude and
#' latitude coordinates (in that order) at which to estimate sunrise and sunset times.
#' @param start.time A POSIXct object containing the earliest date of the monitoring period.
#' @param end.time A POSIXct object containing the latest date of the monitoring period.
#' @param by Date-time format (as defined by \code{\link[base]{strptime}}) containing the
#' time frame used to average sunrise, sunset and twilight times. Defaults to month ("%m").
#' @param solar.depth Angle of the sun below the horizon (in degrees). Passed the
#' solarDep argument in \code{\link[suntools]{crepuscule}} function.
#' Defaults to 18 (astronomical twilight).
#' @return Data frame with diel phase' boundary times (in hours)
#' @export


getSunTimes <- function(sunriset.coords, start.time, end.time, by="%m", solar.depth=18) {

  if(class(sunriset.coords)[1]=="numeric" & length(sunriset.coords)==2){
    sunriset.coords <- matrix(sunriset.coords, ncol=2)
  }

  if(!class(sunriset.coords)[1] %in% c("matrix", "SpatialPoints")){
    stop("sunriset.coords should be supplied as SpatialPoints, matrix or numeric format\n")
  }

  timebins <- seq.POSIXt(from=start.time, to=end.time, by=60*60*24)
  sunrises <- suntools::sunriset(sunriset.coords, timebins, POSIXct.out=TRUE, direction="sunrise")$time
  sunsets <- suntools::sunriset(sunriset.coords, timebins, POSIXct.out=TRUE, direction="sunset")$time
  dusks <- suntools::crepuscule(sunriset.coords, timebins, POSIXct.out=TRUE, solarDep=solar.depth, direction="dusk")$time
  dawns <- suntools::crepuscule(sunriset.coords, timebins, POSIXct.out=TRUE, solarDep=solar.depth, direction="dawn")$time

  intervals <- strftime(timebins, by)
  sunrises <- as.POSIXlt(sunrises)
  sunrises <- sunrises$hour + sunrises$min/60 + sunrises$sec/3600
  sunsets <- as.POSIXlt(sunsets)
  sunsets <- sunsets$hour + sunsets$min/60 + sunsets$sec/3600
  dusks <- as.POSIXlt(dusks)
  dusks <- dusks$hour + dusks$min/60 + dusks$sec/3600
  dawns <- as.POSIXlt(dawns)
  dawns <- dawns$hour + dawns$min/60 + dawns$sec/3600

  daytimes_table <- data.frame(intervals, dawns, sunrises, sunsets, dusks)
  daytimes_table$order <- 1:nrow(daytimes_table)
  row_order <- stats::aggregate(daytimes_table$order, by=list(daytimes_table$intervals), FUN=min)$x
  daytimes_table <- stats::aggregate(daytimes_table[,-1], by=list(daytimes_table$intervals), FUN=mean)
  daytimes_table <- daytimes_table[order(row_order),]
  daytimes_table <- daytimes_table[,-ncol(daytimes_table)]
  colnames(daytimes_table)[1] <- "interval"
  return(daytimes_table)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
