## -----------------------------------------------------------------------------
#| message: false
library(moby)
detections_csv <- system.file("extdata", "rays_detections.csv", package = "moby")


## -----------------------------------------------------------------------------
#| eval: false
# detections <- importDetections(
#   detections_csv,
#   source = "generic",
#   col.map = list(
#     ID          = "animal_id",
#     datetime    = "timestamp",
#     station     = "station_name",
#     lon         = "deploy_longitude",
#     lat         = "deploy_latitude",
#     transmitter = "transmitter"      # kept so matchTags() can join tag metadata
#   )
# )


## -----------------------------------------------------------------------------
#| eval: false
# tags <- importTags(rays_tags, source = "generic",
#                    col.map = list(ID = "ID", tagging_date = "tagging_date"))
# 
# detections <- matchTags(detections, tags)


## -----------------------------------------------------------------------------
#| eval: false
# # one named list element per species -> drives per-species outputs everywhere
# id_groups <- split(rays_tags$ID, rays_tags$species)
# 
# dataset <- as_moby(detections, tags = tags, id.groups = id_groups)
# 
# mobyMeta(dataset)   # tagging dates, derived from 'tags' at this call
# dataset             # the print method summarises what is attached

