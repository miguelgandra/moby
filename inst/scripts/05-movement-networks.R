## -----------------------------------------------------------------------------
#| message: false
library(moby)
data(rays)
id_groups <- mobyMeta(rays)$id.groups


## -----------------------------------------------------------------------------
#| eval: false
# mov_net <- calculateTransitions(rays, spatial.col = "station", id.groups = id_groups)
# mov_net                       # a mobyNetwork object (edge list + node table)


## -----------------------------------------------------------------------------
#| eval: false
# networkMetrics(mov_net)       # node- and network-level graph metrics


## -----------------------------------------------------------------------------
#| eval: false
# mov_rand <- randomizeTransitions(mov_net, iterations = 1000)


## -----------------------------------------------------------------------------
#| eval: false
# plotMovements(mov_net, id.groups = id_groups)   # also: plot(mov_net)
# transitionsTable(mov_net)

