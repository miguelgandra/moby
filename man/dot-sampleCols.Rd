% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomizeOverlap.R
\name{.sampleCols}
\alias{.sampleCols}
\title{Helper function - sample columns of a data frame within monitoring periods}
\usage{
.sampleCols(x, id.cols, start_dates, end_dates)
}
\description{
This function samples each column (except the timebin column) of a data frame
within the monitoring periods of each animal, specified by start and end dates.
}
\note{
This function is intended for internal use within the 'moby' package.
}
\keyword{internal}
