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
#> ── calculateStepDistances() ──────────────────────────────────────────── moby ──
#> 
#> ℹ Measuring distance between consecutive positions
#> • Input: 1,643 positions · 8 individuals
#> 
#> → Method
#>   • paths  straight-line (great-circle)
#> 
#> ✔ 1,635 steps measured
#> ℹ Median step: 0 m (0-4,196 m)
traj <- getTrajectories(rays_dist)
length(traj)
#> [1] 8
```
