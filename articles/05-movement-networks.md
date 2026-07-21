# 05 - Movement (spatial) networks

> **In this module**
>
> Build a **movement network** where nodes are locations
> (receivers/sites) and directed edges are transitions between them.
> moby treats the two network types — movement and association —
> symmetrically, so the workflow here is the same four-verb rhythm you
> will reuse in module 06: **calculate → metrics → randomise → plot**.
>
> **Prerequisites:** none. Loads the `rays` checkpoint.

## Learning objectives

1.  Build a movement network with
    [`calculateTransitions()`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md).
2.  Quantify it with
    [`networkMetrics()`](https://miguelgandra.github.io/moby/reference/networkMetrics.md).
3.  Test it against a null model with
    [`randomizeTransitions()`](https://miguelgandra.github.io/moby/reference/randomizeTransitions.md).
4.  Map it with
    [`plotMovements()`](https://miguelgandra.github.io/moby/reference/plotMovements.md)
    and tabulate it with
    [`transitionsTable()`](https://miguelgandra.github.io/moby/reference/transitionsTable.md).

## Setup

``` r

library(moby)
data(rays)
id_groups <- mobyMeta(rays)$id.groups
```

## 1. Build the network

A transition is a movement between two *different* consecutive
**residence events** (visits) in each animal’s time-ordered sequence.
Visits are segmented by `max.gap` (default 48 h): consecutive detections
at the same station are one visit unless separated by a longer absence,
in which case the later return is treated as a new visit rather than one
continuous stay. Supplying `id.groups` builds one network per species.

``` r

mov_net <- calculateTransitions(rays, spatial.col = "station", id.groups = id_groups)
mov_net                       # a mobyNetwork object (edge list + node table)
```

The node table reports `n_residence` (the number of visits to each
station) and `mean_residence_h` (mean visit duration); edge
`mean_duration_h` is the mean transit time over continuously-observed
moves only (transits that spanned a `> max.gap` absence are kept for
connectivity but flagged `crossed_gap` and excluded from that mean).

> **Tip**
>
> If you want the discrete visits themselves — one row per
> arrival/departure — rather than the network summary, call
> `calculateVisits(rays, spatial.col = "station")`. It shares the same
> `max.gap` segmentation, so its counts and durations line up with the
> network’s node table.

## 2. Network metrics

``` r

networkMetrics(mov_net)       # node- and network-level graph metrics
```

## 3. Null model

``` r

mov_rand <- randomizeTransitions(mov_net, iterations = 1000)
```

## 4. Visualise & tabulate

``` r

plotMovements(mov_net, id.groups = id_groups)   # also: plot(mov_net)
transitionsTable(mov_net)
```

> **Tip**
>
> [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on any
> `mobyNetwork` dispatches to the right routine automatically —
> [`plotMovements()`](https://miguelgandra.github.io/moby/reference/plotMovements.md)
> for movement networks,
> [`plotAssociations()`](https://miguelgandra.github.io/moby/reference/plotAssociations.md)
> for association networks.

## Recap & what’s next

You built, measured, tested and mapped a movement network per species.

➡️ **Next:** [06 — Association (social)
networks](https://miguelgandra.github.io/moby/articles/06-association-networks.md)
— the same four steps, applied to co-occurrence between individuals.
