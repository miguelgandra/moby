% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCWTs.R
\name{plotCWTs}
\alias{plotCWTs}
\title{Continuous wavelet transform (CWT)}
\usage{
plotCWTs(
  data,
  variable,
  plot.title = paste("Wavelet Power Spectrum -", tools::toTitleCase(variable)),
  id.col = getDefaults("ID"),
  timebin.col = getDefaults("timebin"),
  id.groups = NULL,
  wavelet.type = "MORLET",
  same.scale = FALSE,
  period.range = c(3, 48),
  axis.periods = c(6, 12, 16, 24, 48),
  time.unit = "hours",
  color.pal = NULL,
  min.days = NULL,
  detrend = FALSE,
  date.format = "\%d/\%b",
  date.interval = 4,
  date.start = 1,
  cex.main = 1.2,
  cex.lab = 1.1,
  cex.axis = 1,
  cex.legend = 0.9,
  legend.xpos = c(0.895, 0.91),
  legend.ypos = c(0.15, 0.85),
  cols = 1,
  cores = 1,
  ...
)
}
\arguments{
\item{data}{A data frame containing binned animal detections or other time-based measurements.
The time series does not need to be regular; gaps between time bins are allowed.
If there are gaps between measurements (i.e., missing time bins), the function will
automatically assume a value of zero for those missing time steps.}

\item{variable}{Name of the column containing the numeric variable to be analyzed.}

\item{plot.title}{A string specifying the title of the plot. By default, this title is
automatically generated as "Wavelet Power Spectrum" followed by the \code{variable} name.}

\item{id.col}{Name of the column containing animal IDs. Defaults to 'ID'.}

\item{timebin.col}{Name of the column containing time bins (in POSIXct format). Defaults to 'timebin'.}

\item{id.groups}{Optional. A list containing ID groups, used to
visually aggregate animals belonging to the same class (e.g. different species).}

\item{wavelet.type}{A character string specifying the wavelet type to be used in the continuous wavelet transform.
This is passed to the \code{wname} argument in the \code{\link[wavScalogram]{cwt_wst}} function.
Possible values are: "MORLET", "DOG", "PAUL", "HAAR", or "HAAR2". The default is "MORLET".}

\item{same.scale}{Forces same spectral scale (zlims) across all plots,
allowing for density comparison between individuals.}

\item{period.range}{The range of period scales (y-axis limits) to be considered,
specified in the units defined by \code{time.unit}. Defaults to c(3, 48) in hours.}

\item{axis.periods}{Periods to include/highlight on the y-axis, specified in the units
defined by \code{time.unit}. Defaults to c(6,12,16,24,48) in hours.}

\item{time.unit}{Time unit for \code{period.range} and \code{axis.periods}.
Options are "mins", "hours" and "days". Defaults to "hours".}

\item{color.pal}{Color palette. Defaults to the Jet colormap.}

\item{min.days}{Discard individuals that were detected in less than x days.}

\item{detrend}{Detrend time series using differences (\code{\link[base]{diff}})
rather than the actual values. Defaults to False.}

\item{date.format}{Date-time format (as used in \code{\link[base]{strptime}}),
defining the x-axis labels. Defaults to month ("\%d/\%b").}

\item{date.interval}{Number defining the interval between each
displayed date (x-axis label). Defaults to 4.}

\item{date.start}{Integer defining the first displayed date (can be used in combination
with 'date.interval" to better control the x-axis labels). Defaults to 1.}

\item{cex.main}{Determines the size of the title(s). Defaults to 1.2.}

\item{cex.lab}{Determines the size of the y-axis and y-axis labels. Defaults to 1.1.}

\item{cex.axis}{Determines the size of the text labels on the axes. Defaults to 1.}

\item{cex.legend}{Determines the size of the text labels on the color legend. Defaults to 0.9.}

\item{legend.xpos}{Relative position of left and right edge of color bar on first axis (0-1).
Defaults to c(0.90, 0.915).}

\item{legend.ypos}{Relative position of left and right edge of color bar on second axis (0-1).
Defaults to c(0.15, 0.85).}

\item{cols}{Number of columns in the final panel (passed to the mfrow argument).}

\item{cores}{Number of CPU cores to use for the computations. Defaults to 1, which
means no parallel computing (single core).  If set to a value greater than 1,
the function will use parallel computing to speed up calculations.
Run \code{parallel::detectCores()} to check the number of available cores.}

\item{...}{Further arguments passed to the \code{\link[wavScalogram]{cwt_wst}} function.
(e.g. border_effects, waverad).}
}
\description{
Function that analyzes and plots periodic patterns in a time series using a
Continuous Wavelet Transform (CWT) framework. CWT analysis provides an alternative method
to the Fast Fourier Transform (FFT) or other time-frequency decomposition techniques,
enabling the examination of periodicities over time. This function serves mostly as a
wrapper for the \code{\link[wavScalogram]{cwt_wst}} function from the \code{wavScalogram} package
(Bolos & Benitez, 2022).
}
\references{
Bolos, V. J., & Benitez, R. (2022). wavScalogram: an R package with scalogram wavelet tools for time series analysis. The R Journal, 14(2), 164-185.
}
