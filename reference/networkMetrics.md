# Compute node- and network-level metrics for a moby network

Computes a consistent set of graph-theoretic metrics for either an
association (individual co-occurrence) network or a movement
(location-transition) network produced by
[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md)
/
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md).
Metrics are computed per group/subset when the network carries one.

Edge weights are interpreted as connection strengths (association index
for association networks; number of movements for movement networks).
For shortest-path metrics (betweenness) weights are internally inverted,
so that stronger ties act as shorter distances and high-betweenness
nodes identify brokers / movement corridors.

## Usage

``` r
networkMetrics(network, weight = NULL, community = TRUE)
```

## Arguments

- network:

  A
  [`mobyNetwork`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
  object (type `"association"` or `"movement"`).

- weight:

  Optional name of the edge column to use as weight. Defaults to
  `"association"` (association networks) or `"n_movements"` (movement
  networks).

- community:

  Logical; run weighted community detection (walktrap) and report
  modularity and community membership. Defaults to TRUE.

## Value

An object of class `mobyNetworkMetrics`: a list with

- nodes:

  A data frame of per-node metrics (one row per node and group/subset).
  For association networks: `degree`, `strength`, `betweenness`,
  `eigenvector`, `clustering`, `community`. For movement networks:
  `in_degree`, `out_degree`, `in_strength`, `out_strength`,
  `betweenness`, `community`.

- network:

  A data frame of network-level metrics per group/subset (`n_nodes`,
  `n_edges`, `density`, `mean_strength`, `modularity`, `n_communities`,
  plus `transitivity` for association networks or `reciprocity` for
  movement networks).

## See also

[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md)

## Examples

``` r
data(rays)
trans <- calculateTransitions(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> - Segmenting residence events with max.gap = 48 hours; tune/justify per system (Inf = split on location change only).
# node- and network-level graph metrics
metrics <- networkMetrics(trans)
metrics
#> <mobyNetworkMetrics> movement network
#>   network-level metrics:
#>               group n_nodes n_edges   density mean_strength modularity
#>        Raja clavata       6      30 1.0000000      38.00000          0
#>  Dasyatis pastinaca       6      28 0.9333333      33.33333          0
#>  n_communities reciprocity
#>              1   1.0000000
#>              1   0.9285714
#>   node-level metrics: 12 rows x 7 metrics (see $nodes)
head(metrics$nodes)
#>          group node in_degree out_degree in_strength out_strength betweenness
#> 1 Raja clavata ST01         5          5           9           10           0
#> 2 Raja clavata ST02         5          5          23           22           0
#> 3 Raja clavata ST03         5          5          26           26          10
#> 4 Raja clavata ST04         5          5          28           28           7
#> 5 Raja clavata ST05         5          5          12           12           0
#> 6 Raja clavata ST06         5          5          16           16           0
#>   community
#> 1         1
#> 2         1
#> 3         1
#> 4         1
#> 5         1
#> 6         1
```
