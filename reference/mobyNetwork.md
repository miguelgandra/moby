# Inspect a moby network object

Accessors and predicates for `mobyNetwork` objects (returned by
[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md)
and, in future, the movement-network functions). A `mobyNetwork` is a
`data.frame` of edges that additionally carries the node table and
network metadata as attributes.

## Usage

``` r
is_mobyNetwork(x)

networkEdges(x)

networkNodes(x)

networkType(x)
```

## Arguments

- x:

  An object.

## Value

`is_mobyNetwork()` returns a logical; `networkEdges()` returns the edge
data frame; `networkNodes()` returns the node data frame;
`networkType()` returns `"association"` or `"movement"`.

## See also

[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md)

## Examples

``` r
data(rays)
trans <- calculateTransitions(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> - Segmenting residence events with max.gap = 48 hours; tune/justify per system (Inf = split on location change only).
is_mobyNetwork(trans)
#> [1] TRUE
networkType(trans)
#> [1] "movement"
head(networkEdges(trans))
#>          group from   to n_movements n_individuals mean_duration_h
#> 1 Raja clavata ST01 ST02           1             1              NA
#> 2 Raja clavata ST01 ST03           5             3        20.67031
#> 3 Raja clavata ST01 ST04           2             2              NA
#> 4 Raja clavata ST01 ST05           1             1        32.49407
#> 5 Raja clavata ST01 ST06           1             1        32.77300
#> 6 Raja clavata ST02 ST01           2             1              NA
head(networkNodes(trans))
#>          group site n_detections n_individuals n_residence mean_residence_h
#> 1 Raja clavata ST01           60             4          10         3.030878
#> 2 Raja clavata ST02          192             3          29         4.721206
#> 3 Raja clavata ST03          244             4          33         5.355790
#> 4 Raja clavata ST04          165             4          30         3.642327
#> 5 Raja clavata ST05          130             4          15         6.976132
#> 6 Raja clavata ST06          120             4          18         3.488669
#>      lon    lat
#> 1 -9.020 38.454
#> 2 -9.008 38.464
#> 3 -8.996 38.456
#> 4 -8.985 38.466
#> 5 -8.972 38.452
#> 6 -8.990 38.442
```
