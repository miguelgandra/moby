% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validateArguments.R
\name{.validateArguments}
\alias{.validateArguments}
\title{Validate Arguments}
\usage{
.validateArguments()
}
\value{
A list containing potentially modified versions of \code{data}, \code{tagging.dates}, and \code{tag.durations}.
If validation fails, the function will stop and return an error message.
}
\description{
This internal function performs a comprehensive set of validation checks on the arguments
provided to functions within the \code{moby} package. It ensures that the data and various parameters
are in the correct format and meet the expected criteria, issuing errors and warnings as necessary.
The function tries to issue as many errors and warnings at the same time as possible, to speed up the
completion of required changes and reduce the time it takes to get the function call right.
}
\details{
The \code{.validateArguments} function performs rigorous checks to ensure that all input arguments
are in the correct format, of the correct type, and consistent with the requirements of the \code{moby} package.
This includes checks for data structure, column specifications, and additional parameters. The validation
aims to provide detailed error and warning messages to guide users in correcting their input.

#' Specifically, the following checks are performed:
\describe{
\item{\code{data}}{Ensures the first argument is a data frame.}
\item{\code{id.col}}{Validates that the ID column exists in the data frame, is unique for each individual, and is converted to a factor if necessary.}
\item{\code{id.groups}}{Ensures ID groups are unique and correctly referenced in the data.}
\item{\code{datetime.col}}{Checks that the datetime column exists and is in POSIXct format.}
\item{\code{timebin.col}}{Validates the existence of the time-bin column and checks that it is in POSIXct format.}
\item{\code{lon.col} and \code{lat.col}}{Ensures the longitude and latitude columns exist and contain numeric values.}
\item{\code{station.col}}{Checks for the existence of the station column.}
\item{\code{spatial.col}}{Validates spatial column(s), ensuring they exist and match the required format for spatial analyses.}
\item{\code{dist.col}}{Checks the distance column, ensuring it contains numeric values.}
\item{\code{split.by}}{Validates the grouping variable(s) used to split the data.}
\item{\code{subset}}{Checks that subset criteria are valid and appropriately referenced in the data.}
\item{\code{variable}}{Ensures a single variable is specified and exists in the data.}
\item{\code{variables}}{Validates that all specified variables exist in the data.}
\item{\code{color.by}}{Checks that the color-by column exists, contains no missing values, and is converted to a factor if needed.}
\item{\code{tagging.dates}}{Ensures tagging dates are in POSIXct format, match the unique IDs, and are appropriately named if provided as a vector.}
\item{\code{tag.durations}}{Validates tag durations, ensuring they are numeric and consistent with the unique IDs.}
\item{\code{start.dates} and \code{end.dates}}{Checks start and end dates, ensuring they are in POSIXct format and consistent with the data.}
\item{\code{cutoff.dates}}{Validates cutoff dates, ensuring they are in POSIXct format and logically set.}
\item{\code{last.monitoring.date}}{Ensures the last monitoring date is in POSIXct format and appropriately set relative to the data.}
\item{\code{land.shape}}{Checks that the provided spatial object is either an \code{sf} object or convertible to one.}
\item{\code{color.pal}}{Validates the color palette, ensuring it is a character vector or a function.}
\item{\code{diel.lines}}{Checks that diel line values are valid (e.g., 0, 2, or 4).}
\item{\code{polygons}}{Ensures polygons are validly set to \code{'diel'}, \code{'season'}, or \code{FALSE}.}
\item{\code{background.color}}{Checks that the background color is a valid named color or hexadecimal code.}
\item{\code{style}}{Validates that the style parameter is either \code{'raster'} or \code{'points'}.}
\item{\code{cores}}{Ensures that the number of cores for parallel processing is specified as a positive integer. If this parameter is missing or invalid, a default value of 1 is used.}
\item{\code{cols}}{Checks that the \code{cols} parameter is correctly formatted as a single positive integer.}
}
}
\note{
This function is intended for internal use within the \code{moby} package. It helps streamline
argument validation across multiple functions, ensuring consistency and reducing redundancy.
}
\keyword{internal}
