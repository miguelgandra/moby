## -----------------------------------------------------------------------------
#| message: false
library(moby)
data(rays)
id_groups <- mobyMeta(rays)$id.groups


## -----------------------------------------------------------------------------
#| eval: false
# coas <- calculateCOAs(rays)   # aggregates the detections in each time bin (rays already carries 'timebin')


## -----------------------------------------------------------------------------
#| eval: false
# # great-circle step distances (supply a land.shape for shortest in-water paths)
# tracks <- calculateStepDistances(coas)
# getTrajectories(tracks)   # the reconstructed path geometries (an attribute)


## -----------------------------------------------------------------------------
#| eval: false
# rom <- calculateROM(tracks)                 # total distance + mean/max rate of movement
# li  <- calculateLinearityIndex(tracks)      # movement directness (0 = convoluted, 1 = straight)


## -----------------------------------------------------------------------------
#| eval: false
# kud <- calculateUDs(coas, id.groups = id_groups)   # 'coas' carries the CRS forward from 'rays'
# 
# plotMaps(coas, uds = kud, animal.tracks = tracks,
#          id.groups = id_groups)


## -----------------------------------------------------------------------------
#| eval: false
# movementTable(tracks, ud.results = kud, id.groups = id_groups)

