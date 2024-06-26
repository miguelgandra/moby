% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setDefaults.R
\name{setDefaults}
\alias{setDefaults}
\title{Set Default Configuration Parameters}
\usage{
setDefaults(
  id.col = NULL,
  datetime.col = NULL,
  timebin.col = NULL,
  station.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  epsg.code = NULL,
  tagging.dates = NULL
)
}
\arguments{
\item{id.col}{Name of the column containing animal IDs. Defaults to 'ID'.}

\item{datetime.col}{Name of the column containing datetimes in POSIXct format. Defaults to 'datetime'.}

\item{timebin.col}{Name of the column containing time bins (in POSIXct format). Defaults to 'timebin'.}

\item{station.col}{Name of the column containing station/receiver IDs. Defaults to 'station'.}

\item{lon.col}{Name of the column containing longitude values (unprojected). Defaults to 'lon'.}

\item{lat.col}{Name of the column containing latitude values (unprojected). Defaults to 'lat'.}

\item{epsg.code}{Coordinate reference system used to project positions/layers (class 'CRS').}

\item{tagging.dates}{A POSIXct vector containing the tag/release date of each animal.
The length of this vector should match the number of unique animal IDs.
Alternatively, if a single value is provided, it will be applied to all IDs.}
}
\description{
This function updates the default parameters for the 'moby' package.
It accesses a dedicated environment to store user-defined settings, such as column names
for animal IDs, date-times, and time-bins This allows for streamlined usage across
multiple functions within the package, as these parameters only need to be defined once.
}
\examples{
# Setting default parameters using setDefaults
setDefaults(
  id.col = "animal_id",
  datetime.col = "timestamp",
  timebin.col = "time_bin",
  tagging.dates = c("2020-01-01", "2020-06-01"),
)

# Retrieving default parameters using getDefaults
# Get the default column name for animal IDs
id_column <- getDefaults("id")
print(id_column)  # Should print "animal_id"

}
