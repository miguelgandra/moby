# 06 - Association (social) networks

> **In this module**
>
> Build an **association network** where nodes are individuals and edges
> are co-occurrence-based association indices (how often two animals are
> detected together in space and time). The workflow is deliberately
> identical to module 05 — **calculate → metrics → randomise → plot** —
> over the *same* `mobyNetwork` object and
> [`networkMetrics()`](https://miguelgandra.github.io/moby/reference/networkMetrics.md)
> function.
>
> **Prerequisites:** none. Loads the `rays` checkpoint.

> **Interpretation**
>
> Co-occurrence indices (SRI/HWI) measure *shared space-and-time use*,
> **not** confirmed social interaction. Always pair them with a
> randomisation test, and report associations with multiple-comparison
> control.

## Learning objectives

1.  Build an association network with
    [`calculateAssociations()`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md).
2.  Quantify it with
    [`networkMetrics()`](https://miguelgandra.github.io/moby/reference/networkMetrics.md)
    (same function as module 05).
3.  Test it against a datastream-permutation null with
    [`randomizeAssociations()`](https://miguelgandra.github.io/moby/reference/randomizeAssociations.md).
4.  Visualise it with
    [`plotAssociations()`](https://miguelgandra.github.io/moby/reference/plotAssociations.md),
    [`plotAssociationMatrix()`](https://miguelgandra.github.io/moby/reference/plotAssociationMatrix.md)
    and
    [`plotGroupSizeDistribution()`](https://miguelgandra.github.io/moby/reference/plotGroupSizeDistribution.md).

## Setup

``` r

library(moby)
data(rays)
id_groups <- mobyMeta(rays)$id.groups
```

## 1. Build the network

Associations are usually computed *within* species; passing `id.groups`
builds one network per group.

``` r

rays_wide <- createWideTable(rays, value.col = "station")   # time-bin x individual co-occurrence table
assoc_net <- calculateAssociations(rays_wide, id.groups = id_groups)
assoc_net                       # a mobyNetwork object
```

## 2. Network metrics

``` r

networkMetrics(assoc_net)
```

## 3. Null model (with multiple-comparison control)

``` r

assoc_rand <- randomizeAssociations(rays_wide, assoc_net, iterations = 1000,
                                    p.adjust.method = "fdr")
```

## 4. Visualise

``` r

plotAssociations(assoc_net)        # also: plot(assoc_net)
plotAssociationMatrix(assoc_rand)  # the matrix needs the null-model (randomised) result
plotGroupSizeDistribution(rays, id.groups = id_groups)
```

## Recap

You have now run both moby network pipelines and seen how the shared
`mobyNetwork` object and
[`networkMetrics()`](https://miguelgandra.github.io/moby/reference/networkMetrics.md)
make movement and association analyses mirror each other.

🎉 That completes the core tutorial series. For pulling your own data
from a live database, see the optional [ETN
appendix](https://miguelgandra.github.io/moby/articles/00-etn-import.md).
