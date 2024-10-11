#######################################################################################################
## Set Defaults ###############################################################################
#######################################################################################################

#' Set Default Configuration Parameters
#'
#' @description This function updates the default parameters for the 'moby' package.
#' It accesses a dedicated environment to store user-defined settings, such as column names
#' for animal IDs, date-times, and time-bins This allows for streamlined usage across
#' multiple functions within the package, as these parameters only need to be defined once.
#'
#' @param id.col Name of the column containing animal IDs. Defaults to 'ID'.
#' @param datetime.col Name of the column containing datetimes in POSIXct format. Defaults to 'datetime'.
#' @param timebin.col Name of the column containing time bins (in POSIXct format). Defaults to 'timebin'.
#' @param station.col Name of the column containing station/receiver IDs. Defaults to 'station'.
#' @param lon.col Name of the column containing longitude values. Defaults to 'lon'.
#' @param lat.col Name of the column containing latitude values. Defaults to 'lat'.
#' @param epsg.code The EPSG code (integer) representing the coordinate reference system (CRS) to be used
#' for projecting positions/layers.
#' @param tagging.dates A POSIXct vector containing the tag/release date of each animal.
#' The length of this vector should match the number of unique animal IDs, and the dates must be in the same order as the ID levels.
#' Alternatively, if a single value is provided, it will be applied to all IDs.
#'
#' @examples
#' # Setting default parameters using setDefaults
#' setDefaults(
#'   id.col = "animal_id",
#'   datetime.col = "timestamp",
#'   timebin.col = "time_bin",
#'   tagging.dates = c("2020-01-01", "2020-06-01"),
#' )
#'
#' # Retrieving default parameters using getDefaults
#' # Get the default column name for animal IDs
#' id_column <- getDefaults("id")
#' print(id_column)  # Should print "animal_id"
#'
#' @export


setDefaults <- function(id.col=NULL, datetime.col=NULL, timebin.col=NULL,
                        station.col=NULL, lon.col=NULL, lat.col=NULL,
                        epsg.code=NULL, tagging.dates=NULL) {

  # store new default settings
  if(!is.null(id.col)) mobyEnv$defaults$id.col <- id.col
  if(!is.null(datetime.col)) mobyEnv$defaults$datetime.col <- datetime.col
  if(!is.null(timebin.col)) mobyEnv$defaults$timebin.col <- timebin.col
  if(!is.null(station.col)) mobyEnv$defaults$station.col <- station.col
  if(!is.null(lon.col)) mobyEnv$defaults$lon.col <- lon.col
  if(!is.null(lat.col)) mobyEnv$defaults$lat.col <- lat.col
  if(!is.null(epsg.code)) mobyEnv$defaults$epsg.code <- epsg.code
  if(!is.null(tagging.dates)) mobyEnv$defaults$tagging.dates <- tagging.dates
}


