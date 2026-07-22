#######################################################################################################
# Function to get table of diel phase times (in hours) ################################################
#######################################################################################################

#' Retrieve diel phase boundary times (hours)

#' @description Returns a table containing average sunrise, sunset and twilight times
#' for a given period (study duration).
#' '
#' @param coords A SpatialPoints, matrix or numeric object containing longitude and
#' latitude coordinates (in that order) at which to estimate sunrise and sunset times.
#' @param start.time A POSIXct object containing the earliest date of the monitoring period.
#' @param end.time A POSIXct object containing the latest date of the monitoring period.
#' @param by Date-time format (as defined by \code{\link[base]{strptime}}) containing the
#' time frame used to average sunrise, sunset and twilight times. Defaults to month ("%m").
#' @param solar.depth Angle of the sun below the horizon (in degrees). Passed the
#' solarDep argument in \code{\link[suntools]{crepuscule}} function.
#' Defaults to 18 (astronomical twilight).
#' @return Data frame with diel phase' boundary times (in hours)
#' @examples
#' # average diel-phase boundary times (hours) off SW Portugal, early May
#' getSunTimes(c(-9, 38.4),
#'             as.POSIXct("2023-05-01", tz = "UTC"),
#'             as.POSIXct("2023-05-07", tz = "UTC"))
#' @export


getSunTimes <- function(coords, start.time, end.time, by="%m", solar.depth=18) {

  if(is.numeric(coords) && length(coords)==2){
    coords <- matrix(coords, ncol=2)
  } else if(is.data.frame(coords) && ncol(coords)==2){
    # accept a 2-column data.frame (longitude, latitude), consistent with getDielPhase()
    coords <- as.matrix(coords)
  }

  if(!inherits(coords, c("matrix", "SpatialPoints"))){
    .mobyAbort("'coords' must be a SpatialPoints, matrix, numeric, or 2-column data.frame.")
  }

  # diel times are reported as decimal hours in the timezone of the supplied
  # period (start.time/end.time), so that they align with the local-time axes
  # used by the plotting functions (suntools returns instants in UTC)
  tz <- .dataTZ(start.time)

  timebins <- seq.POSIXt(from=start.time, to=end.time, by=60*60*24)
  sunrises <- suntools::sunriset(coords, timebins, POSIXct.out=TRUE, direction="sunrise")$time
  sunsets <- suntools::sunriset(coords, timebins, POSIXct.out=TRUE, direction="sunset")$time
  dusks <- suntools::crepuscule(coords, timebins, POSIXct.out=TRUE, solarDep=solar.depth, direction="dusk")$time
  dawns <- suntools::crepuscule(coords, timebins, POSIXct.out=TRUE, solarDep=solar.depth, direction="dawn")$time

  intervals <- strftime(timebins, by, tz=tz)
  sunrises <- as.POSIXlt(sunrises, tz=tz)
  sunrises <- sunrises$hour + sunrises$min/60 + sunrises$sec/3600
  sunsets <- as.POSIXlt(sunsets, tz=tz)
  sunsets <- sunsets$hour + sunsets$min/60 + sunsets$sec/3600
  dusks <- as.POSIXlt(dusks, tz=tz)
  dusks <- dusks$hour + dusks$min/60 + dusks$sec/3600
  dawns <- as.POSIXlt(dawns, tz=tz)
  dawns <- dawns$hour + dawns$min/60 + dawns$sec/3600

  daytimes_table <- data.frame(intervals, dawns, sunrises, sunsets, dusks)
  daytimes_table$order <- seq_len(nrow(daytimes_table))
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
