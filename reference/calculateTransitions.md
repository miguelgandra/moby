# Build a movement (transition) network between locations

Builds a **movement network**, in which nodes are spatial units
(receivers, stations, habitats or any user-defined spatial class) and
directed edges represent transitions

- i.e. individuals moving from one location to another. Edge weights
  summarise how many movements occurred and how many distinct
  individuals performed them. This is the spatial counterpart to
  [`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
  and the foundation for movement-network metrics, visualisation and (in
  a forthcoming release) randomisation.

Transitions are defined between *consecutive distinct* residence events
(visits) in each individual's time-ordered sequence, where visits are
segmented by `max.gap` (see
[`calculateVisits`](https://miguelgandra.github.io/moby/reference/calculateVisits.md)):
consecutive detections at the same location are collapsed into one visit
unless separated by an absence longer than `max.gap`, so a transition is
genuine movement between two different nodes and a long absence is not
mistaken for one continuous stay.

## Usage

``` r
calculateTransitions(
  data,
  id.col = NULL,
  datetime.col = NULL,
  spatial.col = NULL,
  id.groups = NULL,
  max.gap = 48,
  max.gap.unit = c("hours", "days", "mins", "secs")
)
```

## Arguments

- data:

  A data frame (or
  [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md))
  of detections.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- datetime.col:

  Name of the column containing date-times in POSIXct format. Defaults
  to `"datetime"`.

- spatial.col:

  Name of the column defining the network nodes (e.g. receiver, station,
  habitat, region). If `NULL` (default), it is taken from the `mobyData`
  station column (or the canonical `"station"`); set it to build
  transitions between any other spatial unit.

- id.groups:

  Optional named list of ID groups; when supplied, an independent
  network is built for each group (carried in the `group` column of the
  node/edge tables).

- max.gap:

  Maximum tolerated gap between successive detections within a single
  visit (passed to
  [`calculateVisits`](https://miguelgandra.github.io/moby/reference/calculateVisits.md));
  a longer absence ends a stay, so a later return counts as a new visit.
  Defaults to 48 (hours). Use `Inf` to segment on location changes only.

- max.gap.unit:

  Units of `max.gap`: one of `"hours"` (default), `"days"`, `"mins"`,
  `"secs"`.

## Value

A
[`mobyNetwork`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
object of type `"movement"`. Its rows are the directed edges with
columns: `group`, `from`, `to`, `n_movements` (number of transitions),
`n_individuals` (distinct individuals performing them), and
`mean_duration_h` (mean transit time in hours, computed over
continuously-observed movements only — transits that spanned a `max.gap`
absence are excluded as their timing is unobserved). The node table
([`networkNodes`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md))
has one row per location with `group`, `site`, `n_detections`,
`n_individuals`, `n_residence` (number of residence events / visits),
`mean_residence_h` (mean visit duration in hours), and `lon`/`lat` (mean
coordinates, when available). The full per-transition records (departure
/ arrival times, hour, month, `duration_h`, and `crossed_gap`) are
stored in the `"transition_records"` attribute.

## See also

[`calculateVisits`](https://miguelgandra.github.io/moby/reference/calculateVisits.md),
[`calculateAssociations`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
[`mobyNetwork`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md),
[`transitionsTable`](https://miguelgandra.github.io/moby/reference/transitionsTable.md)

## Examples

``` r
data(rays)
# build a movement network with the receiver stations as nodes
trans <- calculateTransitions(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> - Segmenting residence events with max.gap = 48 hours; tune/justify per system (Inf = split on location change only).
trans
#> <mobyNetwork> movement network
#>   nodes: 12  |  edges: 58
#>   groups: Raja clavata, Dasyatis pastinaca
#>   edge table (use networkEdges()/networkNodes() to extract):
#>          group from   to n_movements n_individuals mean_duration_h
#> 1 Raja clavata ST01 ST02           1             1              NA
#> 2 Raja clavata ST01 ST03           5             3        20.67031
#> 3 Raja clavata ST01 ST04           2             2              NA
#> 4 Raja clavata ST01 ST05           1             1        32.49407
#> 5 Raja clavata ST01 ST06           1             1        32.77300
#> 6 Raja clavata ST02 ST01           2             1              NA
head(networkEdges(trans))
#>          group from   to n_movements n_individuals mean_duration_h
#> 1 Raja clavata ST01 ST02           1             1              NA
#> 2 Raja clavata ST01 ST03           5             3        20.67031
#> 3 Raja clavata ST01 ST04           2             2              NA
#> 4 Raja clavata ST01 ST05           1             1        32.49407
#> 5 Raja clavata ST01 ST06           1             1        32.77300
#> 6 Raja clavata ST02 ST01           2             1              NA
```
