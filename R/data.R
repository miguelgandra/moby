#######################################################################################################
## Toy teaching dataset ##############################################################################
#######################################################################################################

#' Toy acoustic-telemetry dataset (two rays in an MPA)
#'
#' @description A small, **fully synthetic** acoustic-telemetry dataset used throughout the moby
#' tutorials. It mimics two demersal elasmobranchs - thornback ray (*Raja clavata*) and common
#' stingray (*Dasyatis pastinaca*) - monitored by a six-receiver array inside a small coastal
#' marine protected area (MPA). Coordinates, dates and detection patterns are simulated (no real
#' protected-species locations are involved), but are designed to be ecologically plausible enough
#' to exercise the whole moby workflow, including per-species analyses via `id.groups`.
#'
#' The objects come as a small family so that each tutorial can start from the appropriate
#' checkpoint without re-running earlier ones:
#' \itemize{
#'   \item `rays_detections` - raw, harmonised detections (a plain data frame).
#'   \item `rays_tags` - tag/animal metadata (one row per animal).
#'   \item `rays_deployments` - receiver-deployment log (one row per station).
#'   \item `rays` - the cleaned, analysis-ready \code{\link{mobyData}} checkpoint (detections with a
#'   `timebin` column and metadata attached, including the two-species `id.groups`).
#' }
#' Raw CSV versions for the import tutorial are also shipped: the detections in a generic
#' hand-mapped layout at \code{system.file("extdata", "rays_detections.csv", package = "moby")}, and
#' the receiver-deployment log in raw ETN (European Tracking Network) export format at
#' \code{system.file("extdata", "rays_deployments.csv", package = "moby")} (harmonise it with
#' \code{importDeployments(..., source = "etn")}).
#'
#' @format
#' `rays_detections`: a data frame with one row per detection and columns `ID`, `datetime`
#' (POSIXct, UTC), `station`, `lon`, `lat`, `species`, `receiver`, `transmitter`.
#'
#' `rays_tags`: a data frame with columns `ID`, `transmitter`, `species`, `tagging_date` (POSIXct),
#' `tagging_station`.
#'
#' `rays_deployments`: a data frame with columns `receiver`, `station`, `lon`, `lat`, `deploy`,
#' `recover` (POSIXct) and `depth` - the canonical deployment-log schema (as produced by
#' \code{\link{importDeployments}}) used by \code{\link{checkDeployments}},
#' \code{\link{matchDeployments}} and \code{\link{plotDeployments}}.
#'
#' `rays`: a \code{\link{mobyData}} (data-frame subclass) with the columns of `rays_detections` plus
#' a `timebin` column, and a `"moby"` metadata attribute (column mapping, `epsg.code` (32629, UTM 29N),
#' tagging dates and the per-species `id.groups`).
#'
#' @source Simulated with `data-raw/make-toy-data.R` (seeded, reproducible). Not derived from any
#' real tracking study.
#'
#' @name rays
#' @aliases rays rays_detections rays_tags rays_deployments
#' @docType data
#' @keywords datasets
#' @examples
#' # the analysis-ready mobyData checkpoint (metadata attached)
#' data(rays)
#' is_moby(rays)
#' mobyMeta(rays)$id.groups        # the two-species grouping
#' head(rays)
#'
#' # the accompanying raw tables
#' head(rays_detections)           # harmonised detections (plain data frame)
#' head(rays_tags)                 # one row per tagged animal
#' head(rays_deployments)          # receiver-deployment log
"rays"

#' @rdname rays
"rays_detections"

#' @rdname rays
"rays_tags"

#' @rdname rays
"rays_deployments"
