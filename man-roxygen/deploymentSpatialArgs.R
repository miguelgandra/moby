#' @param deployment.station.col,deployment.lon.col,deployment.lat.col Names of the station,
#' longitude and latitude columns in the receiver-deployment log (`deployments`). Default to the
#' canonical `"station"`/`"lon"`/`"lat"` produced by \code{\link{importDeployments}}; set them when a
#' hand-made log uses other names (e.g. `deployment.lon.col = "Longitude"`). The `deployment.` prefix
#' marks these as deployment-log columns, keeping them distinct from the bare `*.col` arguments,
#' which always refer to the detection dataset. The `receiver` column is the canonical join key and
#' is always taken as-is.
