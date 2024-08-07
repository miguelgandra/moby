% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSeason.R
\name{getSeason}
\alias{getSeason}
\title{Estimate annual season}
\usage{
getSeason(datetimes, hemisphere = "Northern")
}
\arguments{
\item{datetimes}{A POSIXct object containing the respective datetimes or time-bins.}

\item{hemisphere}{A character string specifying the hemisphere ("Northern" or "Southern").}
}
\value{
A factor indicating the season.
}
\description{
This function determines the meteorological season for a given date based on the
meteorological definition, which divides the year into four seasons of three months each:\cr

Northern Hemisphere:\cr
• Spring: March - May\cr
• Summer: June - August\cr
• Autumn: September - November\cr
• Winter: December - February\cr

Southern Hemisphere:\cr
• Spring: September - November\cr
• Summer: December - February\cr
• Autumn: March - May\cr
• Winter: June - August\cr
}
\examples{
datetimes <- as.POSIXct("2024-05-30")
getSeason(datetimes, hemisphere="Northern")
}
