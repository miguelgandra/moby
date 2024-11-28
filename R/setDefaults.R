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
#' @param lon.col Name of the column containing longitude values. Defaults to 'lon'. Coordinates can be supplied
#' either in geographic (unprojected) or projected format. In both cases, the appropriate projected CRS
#' should be provided via the `epsg.code` argument.
#' @param lat.col Name of the column containing latitude values. Defaults to 'lat'. Coordinates can be supplied
#' either in geographic (unprojected) or projected format. In both cases, the appropriate projected CRS
#' should be provided via the `epsg.code` argument.
#' @param epsg.code An integer representing the EPSG code of the coordinate reference system (CRS). This parameter
#' serves two purposes: (1) to specify the projection system used when coordinates are already projected, or
#' (2) to define the target projection system when coordinates are supplied in geographic format. Note that
#' this must always refer to a **projected CRS**
#' @param tagging.dates A POSIXct vector specifying the tagging or release dates for each animal.
#' This parameter must be either:
#' - A single POSIXct value, which will be applied to all unique animal IDs; or
#' - A named POSIXct vector, where the names correspond to the animal IDs in the `id.col` column.
#' If multiple tagging dates are provided, the vector must include all IDs and will be reordered to align with the levels of `id.col`.
#'
#' @examples
#' # Setting default parameters using setDefaults
#' setDefaults(
#'   id.col = "animal_id",
#'   datetime.col = "timestamp",
#'   timebin.col = "time_bin",
#'   tagging.dates = c("2020-01-01"),
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


