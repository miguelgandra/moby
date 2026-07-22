#' @param deployment.deploy.col,deployment.recover.col Names of the deployment and recovery
#' date-time columns in the receiver-deployment log (`deployments`). Default to the canonical
#' `"deploy"`/`"recover"` produced by \code{\link{importDeployments}}; set them when a hand-made log
#' uses other names (e.g. `deployment.deploy.col = "deploy_date"`). The `deployment.` prefix marks
#' these as deployment-log columns, keeping them distinct from the bare `*.col` arguments, which
#' always refer to the detection dataset.
