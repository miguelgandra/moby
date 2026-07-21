## -----------------------------------------------------------------------------
#| message: false
library(moby)
data(rays_detections)
data(rays_tags)


## -----------------------------------------------------------------------------
#| eval: false
# dep_file <- system.file("extdata", "rays_deployments.csv", package = "moby")   # raw ETN export
# deployments <- importDeployments(dep_file, source = "etn")
# head(deployments)


## -----------------------------------------------------------------------------
#| eval: false
# qc <- checkDeployments(deployments)
# qc   # flags overlaps, impossible dates/coords, duplicates, coverage gaps


## -----------------------------------------------------------------------------
#| eval: false
# checkDeployments(deployments, detections = rays_detections, datetime.col = "datetime")


## -----------------------------------------------------------------------------
#| eval: false
# checkDeployments(deployments, checks = c("dates", "overlaps", "coordinates"))


## -----------------------------------------------------------------------------
#| eval: false
# plotDeployments(deployments)


## -----------------------------------------------------------------------------
#| eval: false
# detections <- matchDeployments(rays_detections, deployments,
#                                station.col = "station")


## -----------------------------------------------------------------------------
#| eval: false
# filtered <- filterDetections(
#   rays_detections,
#   id.col        = "ID",
#   datetime.col  = "datetime",
#   lon.col       = "lon", lat.col = "lat",
#   tagging.dates = setNames(rays_tags$tagging_date, rays_tags$ID),
#   nominal.delay = 120,         # s; enables the short-interval false-detection filter (~30x delay)
#   max.speed     = 2,           # m/s; great-circle distances (no land.shape needed)
#   speed.unit    = "m/s"
# )
# filtered                        # a printable <mobyFilter> summary
# clean <- filtered$data          # the cleaned detections (with a qc_flag column)


## -----------------------------------------------------------------------------
#| eval: false
# clean$timebin <- getTimeBins(clean$datetime, interval = "1 hour")

