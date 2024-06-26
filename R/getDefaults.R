#######################################################################################################
## Get Defaults #######################################################################################
#######################################################################################################

#' Get Default Configuration Parameters
#'
#' @description This function retrieves the current default parameters for the 'moby' package.
#' It accesses a dedicated environment where user-defined settings are stored, such as column
#' names for animal IDs, date-times, and time-bins.
#'
#' @param var A character string specifying the variable name for which to retrieve the default setting.
#' Accepted values are "id.col", "datetime.col", "timebin.col", "epsg.code" and "tagging.dates".
#' Alternatively, these values can also be accessed using their abbreviated forms: "id", "datetime", "timebin", "epsg" and "tagdates".
#'
#' @return The current default setting for the specified variable.
#'
#' @examples
#' # Retrieve the default column name for animal IDs
#' getDefaults("id")
#'
#' # Retrieve the default column name for date-times
#' getDefaults("datetime")
#'
#' @export

getDefaults <- function(var) {

  if(tolower(var)=="id") var <- "id.col"
  if(tolower(var)=="datetime") var <- "datetime.col"
  if(tolower(var)=="timebin") var <- "timebin.col"
  if(tolower(var)=="station") var <- "station.col"
  if(tolower(var)=="lon") var <- "lon.col"
  if(tolower(var)=="lat") var <- "lat.col"
  if(tolower(var)=="epsg") var <- "epsg.code"
  if(tolower(var)=="tagdates") var <- "tagging.dates"

  if(!var %in% c("id.col", "datetime.col", "timebin.col", "station.col",
                 "lon.col", "lat.col", "epsg.code", "tagging.dates")){
    stop("Invalid variable specified.", call. = FALSE)
  }

  # retrieve default setting
  return(unlist(mobyEnv$defaults[var]))
}


