# Plot a moby network

S3 `plot` method for
[`mobyNetwork`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
objects. For an `"association"` network it draws the individual
co-occurrence/association network (delegating to
[`plotAssociations`](https://miguelgandra.github.io/moby/reference/plotAssociations.md)).
For a `"movement"` network it draws the directed transition network
between locations over a map (delegating to
[`plotMovements`](https://miguelgandra.github.io/moby/reference/plotMovements.md)),
placing nodes at their geographic coordinates when available.

## Usage

``` r
# S3 method for class 'mobyNetwork'
plot(x, ...)
```

## Arguments

- x:

  A `mobyNetwork` object.

- ...:

  Further arguments passed to the underlying plotting routine:
  [`plotAssociations`](https://miguelgandra.github.io/moby/reference/plotAssociations.md)
  for association networks (e.g. `color.by`, `nodes.size`, `edge.color`)
  or
  [`plotMovements`](https://miguelgandra.github.io/moby/reference/plotMovements.md)
  for movement networks (e.g. `land.shape`, `epsg.code`,
  `background.layer`, `edge.metric`, `repel.nodes`).

## Value

Called for its side effect (a plot); invisibly returns `NULL`.

## See also

[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md),
[`plotAssociations`](https://miguelgandra.github.io/moby/reference/plotAssociations.md),
[`plotMovements`](https://miguelgandra.github.io/moby/reference/plotMovements.md)

## Examples

``` r
data(rays)
# association network (individual co-occurrences)
wide <- createWideTable(rays, value.col = "station")
#> Warning: - 'id.col' converted to factor.
#> ── createWideTable() ─────────────────────────────────────────────────── moby ──
#> 
#> ℹ Reshaping detections into a time-bin × individual matrix
#> • Input: 1,643 records · 8 individuals
#> 
#> → Method
#>   • values  station
#> Warning: 3 (ID, time-bin) combination(s) had multiple differing values; the first was kept. Aggregate upstream (e.g. calculateCOAs) to control this.
#> Tied (ID, time-bin) instances (first value kept):
#>                  timebin  ID                ties
#> 1898 2023-06-14 11:00:00 D03 ST01 (1) | ST06 (1)
#> 2163 2023-04-24 05:00:00 D04 ST01 (4) | ST05 (4)
#> 5302 2023-06-25 16:00:00 R04 ST03 (1) | ST05 (1)
#> 
#> ✔ 2,130 time bins × 8 individuals
assoc <- calculateAssociations(wide)
#> ── calculateAssociations() ───────────────────────────────────────────── moby ──
#> 
#> ℹ Building a co-occurrence association network
#> • Input: 8 individuals · 2,130 time bins
#> 
#> → Method
#>   • metric  simple-ratio index (SRI)
if (requireNamespace("qgraph", quietly = TRUE)) {
  plot(assoc)
}
#> ── plotAssociations() ────────────────────────────────────────────────── moby ──
#> 
#> ℹ Mapping pairwise overlaps between individuals
#> • Individuals: 8
#> • Comparisons: All
#> • Total dyads: 28
#> 
#> → Method
#>   • Metric  association index

# movement network (transitions between locations)
trans <- calculateTransitions(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> ── calculateTransitions() ────────────────────────────────────────────── moby ──
#> 
#> ℹ Building a directed movement network between locations
#> • Input: 1,643 detections · 8 individuals · 6 nodes (station)
#> 
#> → Method
#>   • max.gap  48 hours (a longer absence starts a new visit; tune per system)
plot(trans)
#> ── plotMovements() ───────────────────────────────────────────────────── moby ──
#> 
#> ℹ Drawing the network of movements between sites
#> • Nodes/edges: 6 sites, 58 edges
#> • Groups:      2
#> 
#> → Method
#>   • Edge metric  movements
#>   • Projection   force-directed (no map)

```
