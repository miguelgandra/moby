# Extract movement trajectories from a calculateStepDistances() result

Retrieves the per-individual movement trajectories (spatial-line
geometries) attached by
[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
to its distance-enriched output.

## Usage

``` r
getTrajectories(x)
```

## Arguments

- x:

  A distance-enriched data frame returned by
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).

## Value

A named list of trajectory geometries (one element per individual), or
`NULL` if the object carries no trajectories.

## See also

[`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md),
[`plotMaps`](https://miguelgandra.github.io/moby/reference/plotMaps.md)

## Examples

``` r
data(rays)
rays_dist <- calculateStepDistances(rays)
#> Warning: - 'id.col' converted to factor.
#>   |                                                                              |                                                                      |   0%Calculating linear paths between consecutive positions...
#>   |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> Total execution time: 0.05 secs 
traj <- getTrajectories(rays_dist)
length(traj)
#> [1] 8
```
