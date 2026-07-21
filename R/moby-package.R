#' moby: Streamlined Acoustic Telemetry Analyses and Data Visualization
#'
#' @description
#' moby provides a comprehensive suite of functions and methods tailored for processing,
#' visualizing, and interpreting passive acoustic telemetry data, with a primary focus on
#' marine environments. It offers tools spanning data cleaning, temporal and spatial
#' classification, movement and home-range analysis, residency metrics, spatiotemporal
#' overlap and social-network analysis, and a wide range of publication-ready plots.
#'
#' Most functions accept user-defined column names (e.g. `id.col`, `datetime.col`,
#' `station.col`). To avoid repeating these on every call, wrap your detection data once
#' with \code{\link{as_moby}}: the resulting `mobyData` object carries the column mapping,
#' coordinate reference system, tagging dates and ID groups as metadata, which downstream
#' functions read automatically. Plain data frames remain fully supported.
#'
#' @section Getting started:
#' See \code{help(package = "moby")} for the full list of functions, and the package
#' vignette (\code{browseVignettes("moby")}) for an introductory workflow.
#'
#' @keywords internal
"_PACKAGE"
