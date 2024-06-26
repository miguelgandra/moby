% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFFTs.R
\name{plotFFTs}
\alias{plotFFTs}
\title{Analyze fine-scale rhythmic patterns in detections}
\usage{
plotFFTs(
  data,
  tagging.dates,
  type = "detections",
  axis.periods = c(48, 24, 12, 8, 6),
  id.col = "ID",
  min.days = 10,
  detrend = T,
  cols = 2,
  same.scale = T,
  id.groups = NULL
)
}
\arguments{
\item{data}{A data frame containing binned animal detections (must contain a "timebin" column).}

\item{tagging.dates}{A POSIXct vector containing the tag/release date of each animal.}

\item{type}{Type of response used to calculate FFTs. One of  "detections" or "presences".
Presences are estimated by time bin (0 - 1), independently of the number of detections.}

\item{axis.periods}{Periods (in hours) to display below the plots.}

\item{id.col}{Name of the column containing animal IDs. Defaults to 'ID'.}

\item{min.days}{Discard individuals that were detected in less than x days.}

\item{detrend}{Detrend time series using differences (\code{\link[base]{diff}})
rather than the actual detections.}

\item{cols}{Number of columns in the final panel (passed to the mfrow argument).}

\item{same.scale}{Forces same spectral scale (ylims) across all plots,
allowing for density comparison between individuals. If set to false, all y-axis are displayed.
Otherwise they are only displayed in the left-most plots to save space.}

\item{id.groups}{Optional. A list containing ID groups, used to
visually aggregate animals belonging to the same class (e.g. different species).}
}
\description{
Function that analyzes and plots periodic patterns in detections using a
FFT (Fast Fourier Transform) framework. By decomposing data series into the
frequency domain, FFT allows the identification of dominant spectral peaks
that may reflect, for example, tidal (6 h–12 h) or diel (24 h) cyclic patterns in habitat use.
}
