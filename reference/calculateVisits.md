# Extract residence events (visits)

Segments each individual's time-ordered detections into discrete
**residence events** (also called *visits*): episodes of continuous
presence at a location, bounded by arrival and departure. A visit ends
when the animal moves to a different location *or* when the gap since
its previous detection exceeds `max.gap` — i.e. a long absence is
treated as the animal having left and (possibly) returned, rather than
as one uninterrupted stay.

This is the event-resolution counterpart to
[`calculateResidency`](https://miguelgandra.github.io/moby/reference/calculateResidency.md)
(which measures day-scale residency as a presence/absence ratio) and the
residence layer underlying
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md)
(which summarises the same visits as movement-network nodes). Both
functions derive their residence segmentation from the same internal
engine, so counts and durations are consistent across them.

## Usage

``` r
calculateVisits(
  data,
  spatial.col = NULL,
  id.col = NULL,
  datetime.col = NULL,
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

- spatial.col:

  Name of the column defining the locations visited (e.g. receiver,
  station, habitat, region). If `NULL` (default), it is taken from the
  `mobyData` station column (or the canonical `"station"`); set it to
  track visits to any other spatial unit.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- datetime.col:

  Name of the column containing date-times in POSIXct format. Defaults
  to `"datetime"`.

- id.groups:

  Optional named list of ID groups; when supplied, visits are computed
  within each group and the group is carried in the `group` column.

- max.gap:

  Maximum tolerated gap between successive detections within a single
  visit. Absences longer than this split the sequence into separate
  residence events. Defaults to 48 (hours). Use `Inf` to segment on
  location changes only.

- max.gap.unit:

  Units of `max.gap`: one of `"hours"` (default), `"days"`, `"mins"`,
  `"secs"`.

## Value

A data frame with one row per residence event (visit), ordered by group,
individual and arrival time:

- group:

  the ID group (`"all"` when `id.groups` is not supplied).

- id:

  the individual.

- site:

  the location visited (`spatial.col` value).

- arrival, departure:

  first and last detection times of the visit (POSIXct).

- residence_h:

  visit duration in hours (`departure - arrival`; 0 for a
  single-detection visit).

- n_detections:

  number of detections recorded during the visit.

## Details

The `max.gap` threshold is the one biological choice here: it encodes
how long an absence you treat as "left and came back" versus "still
present but temporarily undetected". There is no universal value — it
depends on the species, the detection range and the array layout — so it
is exposed as a first-class, reportable parameter. A principled way to
pick it is to inspect the distribution of intervals between successive
detections, which is often bimodal, with a natural trough separating
within-visit gaps (brief range dropouts) from between-visit gaps (true
excursions). Set `max.gap = Inf` to disable gap-splitting entirely (a
visit then ends only on a change of location).

## References

Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., &
Abecasis, D. (2023). Residency and space use estimation methods based on
passive acoustic telemetry data. Movement Ecology, 11(1), 12.

See also the residence-event extractors `glatos::residence_events` and
`VTrack::RunResidenceExtraction` for the same concept in neighbouring
toolkits.

## See also

[`calculateResidency`](https://miguelgandra.github.io/moby/reference/calculateResidency.md),
[`calculateTransitions`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md)

## Examples

``` r
data(rays)
# discrete visits to each receiver station (48 h gap threshold)
visits <- calculateVisits(rays, spatial.col = "station")
#> Warning: - 'id.col' converted to factor.
#> - Defining residence events with max.gap = 48 hours (a longer absence starts a new visit); tune/justify per system.
head(visits)
#>                group  id site             arrival           departure
#> 1 Dasyatis pastinaca D01 ST03 2023-04-08 15:05:10 2023-04-08 16:51:56
#> 2 Dasyatis pastinaca D01 ST06 2023-04-09 12:54:09 2023-04-09 13:26:11
#> 3 Dasyatis pastinaca D01 ST03 2023-04-13 15:21:31 2023-04-13 17:21:57
#> 4 Dasyatis pastinaca D01 ST02 2023-04-13 18:53:08 2023-04-13 21:48:08
#> 5 Dasyatis pastinaca D01 ST04 2023-04-15 20:14:38 2023-04-15 22:13:38
#> 6 Dasyatis pastinaca D01 ST06 2023-04-17 02:11:27 2023-04-17 06:33:47
#>   n_detections residence_h
#> 1            5   1.7793973
#> 2            3   0.5339714
#> 3            4   2.0072882
#> 4           10   2.9167023
#> 5            7   1.9833155
#> 6           13   4.3721619
```
