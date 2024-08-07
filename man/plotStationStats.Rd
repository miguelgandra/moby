% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotStationStats.R
\name{plotStationStats}
\alias{plotStationStats}
\title{Calculate and plot receivers' statistics}
\usage{
plotStationStats(
  data,
  overlaps = NULL,
  type = c("detections", "co-occurrences"),
  id.col = getDefaults("id"),
  spatial.col = getDefaults("station"),
  number.stations = F,
  rotate.labels = F,
  show.percentage = T,
  lab.threshold = NULL,
  legend.pos = "top",
  legend.style = "horizontal",
  id.groups = NULL,
  color.pal = NULL,
  group.comparisons = "all",
  cols = 1
)
}
\arguments{
\item{data}{A data frame containing animal detections with corresponding time-bins.}

\item{overlaps}{List containing counts of stations where animals co-occurred, as returned by
\code{\link{calculateOverlap}}. Only required if type contains "co-occurrences".}

\item{type}{Type of statistic to calculate/plot. If more than one type is supplied a grouped barplot
is generated. Possible values are "detections", "average detections", "individuals" and "co-occurrences".}

\item{id.col}{Name of the column containing animal IDs. Defaults to 'ID'.}

\item{spatial.col}{Should match the name of a variable in the supplied data containing the
spatial categories to plot (e.g. receivers or broader areas/locations). Defaults to "station".}

\item{number.stations}{Replace receiver names with numbers? Useful if receiver names are
too long to be displayed correctly.}

\item{rotate.labels}{Another alternative to address long receiver names.}

\item{show.percentage}{If true, percentages are displayed on top of each bar.}

\item{lab.threshold}{Optional. If defined, text labels are only only displayed to
frequencies above this threshold (e.g., if lab.threshold=0.05 count values are only displayed
in bars with frequencies >= 5\%).}

\item{legend.pos}{Legend position.}

\item{legend.style}{Legend style (horizontal vs vertical).}

\item{id.groups}{Optional. A list containing ID groups, used to generate different
plots to each class (e.g. different species).}

\item{color.pal}{Color palette for the bars. Should match the number of types defined (1 or 2).}

\item{group.comparisons}{Controls the type of comparisons to be run, when id.groups are defined.
One of "within", "between" or "all". Defaults to "all".}

\item{cols}{Number of columns in the final panel (passed to the mfrow argument).}
}
\description{
Calculates several receivers' statistics and generates corresponding barplots.
Implemented statistics include total detections, average detection frequency (per individual),
number of animals detected and number of co-occurrences.
}
