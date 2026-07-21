## -----------------------------------------------------------------------------
#| message: false
library(moby)
data(rays)
id_groups <- mobyMeta(rays)$id.groups


## -----------------------------------------------------------------------------
#| eval: false
# rays_wide <- createWideTable(rays, value.col = "station")   # time-bin x individual co-occurrence table
# assoc_net <- calculateAssociations(rays_wide, id.groups = id_groups)
# assoc_net                       # a mobyNetwork object


## -----------------------------------------------------------------------------
#| eval: false
# networkMetrics(assoc_net)


## -----------------------------------------------------------------------------
#| eval: false
# assoc_rand <- randomizeAssociations(rays_wide, assoc_net, iterations = 1000,
#                                     p.adjust.method = "fdr")


## -----------------------------------------------------------------------------
#| eval: false
# plotAssociations(assoc_net)        # also: plot(assoc_net)
# plotAssociationMatrix(assoc_rand)  # the matrix needs the null-model (randomised) result
# plotGroupSizeDistribution(rays, id.groups = id_groups)

