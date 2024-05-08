#######################################################################################################
# Function to determine season (meteorological definition) ############################################
#######################################################################################################

#' Estimate annual season.
#'
#' @description In temperate and sub-polar regions, 4 seasons based on the Gregorian calendar
#' (meteorological definition) are generally recognized:

#' Northern Hemisphere:
#'  • Spring: March - May
#'  • Summer: June - August
#'  • Autumn: September - November
#'  • Winter: December - February
#'
#' Southern Hemisphere
#'  • Spring: September - November
#'  • Summer: December - February
#'  • March - May
#'  • Winter: June - August
#'
#' @param date A POSIXct object containing the respective datetimes or time-bins.
#' @param hemisphere Earth hemisphere for which to calculate seasons.
#' @return Factor with season
#' @export


getSeason <- function(date, hemisphere="Northern") {

  # convert dates to year/month and add one month (1/12)
  yq <- zoo::as.yearqtr(zoo::as.yearmon(date, "%m/%d/%Y") + 1/12)
  if(hemisphere=="Northern") {
    season <- factor(format(yq, "%q"), levels=1:4, labels=c("winter", "spring", "summer", "autumn"))}
  if(hemisphere=="Southern") {
    season <- factor(format(yq, "%q"), levels=1:4, labels=c("summer", "autumn", "winter", "spring"))}
  return(season)
}




#######################################################################################################
#######################################################################################################
#######################################################################################################
