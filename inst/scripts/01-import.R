## -----------------------------------------------------------------------------
#| message: false
library(moby)
detections_csv <- system.file("extdata", "rays_detections.csv", package = "moby")


## -----------------------------------------------------------------------------
#| eval: false
# raw <- read.csv(detections_csv, stringsAsFactors = FALSE)
# raw$timestamp <- as.POSIXct(raw$timestamp, tz = "UTC")
# 
# detections <- importDetections(
#   raw,
#   source = "generic",
#   col.map = list(
#     ID          = "animal_id",
#     datetime    = "timestamp",
#     station     = "station_name",
#     lon         = "deploy_longitude",
#     lat         = "deploy_latitude",
#     transmitter = "transmitter"      # kept so assignAnimalIDs() can join tag metadata
#   )
# )


## -----------------------------------------------------------------------------
#| eval: false
# tags <- importTags(rays_tags, source = "generic",
#                    col.map = list(ID = "ID", tagging_date = "tagging_date"))
# 
# detections <- assignAnimalIDs(detections, tags)


## -----------------------------------------------------------------------------
#| eval: false
# # one named list element per species -> drives per-species outputs everywhere
# id_groups <- split(rays_tags$ID, rays_tags$species)
# 
# dataset <- as_moby(
#   detections,
#   id.col        = "ID",
#   datetime.col  = "datetime",
#   station.col   = "station",
#   lon.col       = "lon",
#   lat.col       = "lat",
#   tagging.dates = setNames(rays_tags$tagging_date, rays_tags$ID),
#   id.groups     = id_groups
# )
# 
# mobyMeta(dataset)   # inspect the attached metadata

