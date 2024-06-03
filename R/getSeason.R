#######################################################################################################
# Function to determine season (meteorological definition) ############################################
#######################################################################################################

#' Estimate annual season
#'
#' @description This function determines the meteorological season for a given date based on the
#' meteorological definition, which divides the year into four seasons of three months each:\cr
#'
#' Northern Hemisphere:\cr
#'  • Spring: March - May\cr
#'  • Summer: June - August\cr
#'  • Autumn: September - November\cr
#'  • Winter: December - February\cr
#'
#' Southern Hemisphere:\cr
#'  • Spring: September - November\cr
#'  • Summer: December - February\cr
#'  • Autumn: March - May\cr
#'  • Winter: June - August\cr
#'
#' @param date A POSIXct object containing the respective datetimes or time-bins.
#' @param hemisphere A character string specifying the hemisphere ("Northern" or "Southern").
#' @return A factor indicating the season.
#' @examples
#' date <- as.POSIXct("2024-05-30")
#' getSeason(date, hemisphere="Northern")
#' @export

getSeason <- function(date, hemisphere="Northern") {

  # validate hemisphere input
  if (!hemisphere %in% c("Northern", "Southern")) stop("Hemisphere must be 'Northern' or 'Southern'")

  # extract month
  month <- as.integer(format(date, "%m"))

  # determine season based on hemisphere
  if (hemisphere == "Northern") {
    season <- c("winter", "winter", "spring", "spring", "spring", "summer",
                "summer", "summer", "autumn", "autumn", "autumn", "winter")[month]
  } else {
    season <- c("summer", "summer",  "autumn", "autumn", "autumn", "winter",
                "winter", "winter", "spring", "spring", "spring", "summer")[month]
  }

  return(factor(season, levels=c("spring", "summer", "autumn", "winter")))

}



#######################################################################################################
#######################################################################################################
#######################################################################################################
