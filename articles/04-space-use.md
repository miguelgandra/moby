# 04 - Movement & space use

> **In this module**
>
> Turn detections into positions, then quantify how much animals move
> and the areas they use: centres of activity, step distances, rates of
> movement, and home ranges (kernel utilisation distributions), all per
> species.
>
> **Prerequisites:** none. Loads the `rays` checkpoint.

## Learning objectives

1.  Estimate positions (**centres of activity**) with
    [`calculateCOAs()`](https://miguelgandra.github.io/moby/reference/calculateCOAs.md).
2.  Reconstruct movement paths and step distances with
    [`calculateStepDistances()`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).
3.  Summarise movement with
    [`calculateROM()`](https://miguelgandra.github.io/moby/reference/calculateROM.md)
    and
    [`calculateLinearityIndex()`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md).
4.  Estimate home ranges with
    [`calculateUDs()`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
    (autocorrelated KDE) and map them with
    [`plotMaps()`](https://miguelgandra.github.io/moby/reference/plotMaps.md).

## Setup

``` r

library(moby)
data(rays)
id_groups <- mobyMeta(rays)$id.groups
```

## 1. Centres of activity (positions)

Acoustic detections are presence-at-a-receiver events;
[`calculateCOAs()`](https://miguelgandra.github.io/moby/reference/calculateCOAs.md)
converts them into short-interval position estimates suitable for
movement and home-range analysis.

``` r

coas <- calculateCOAs(rays)   # aggregates the detections in each time bin (rays already carries 'timebin')
```

## 2. Step distances and movement tracks

``` r

# great-circle step distances (supply a land.shape for shortest in-water paths)
tracks <- calculateStepDistances(coas)
getTrajectories(tracks)   # the reconstructed path geometries (an attribute)
```

## 3. Rates of movement & directness

``` r

rom <- calculateROM(tracks)                 # total distance + mean/max rate of movement
li  <- calculateLinearityIndex(tracks)      # movement directness (0 = convoluted, 1 = straight)
```

## 4. Home ranges (UDs) and maps

[`calculateUDs()`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
defaults to **autocorrelated kernel density estimation (AKDE)** via
`ctmm`, which accounts for the autocorrelation inherent in tracking
data.

``` r

kud <- calculateUDs(coas, id.groups = id_groups)   # 'coas' carries the CRS forward from 'rays'

plotMaps(coas, uds = kud, animal.tracks = tracks,
         id.groups = id_groups)
```

## 5. A movement summary table

``` r

movementTable(tracks, ud.results = kud, id.groups = id_groups)
```

## Recap & what’s next

You estimated positions, movement metrics and home ranges per species.

➡️ **Next:** [05 — Movement (spatial)
networks](https://miguelgandra.github.io/moby/articles/05-movement-networks.md).
