# Toy acoustic-telemetry dataset (two rays in an MPA)

A small, **fully synthetic** acoustic-telemetry dataset used throughout
the moby tutorials. It mimics two demersal elasmobranchs - thornback ray
(*Raja clavata*) and common stingray (*Dasyatis pastinaca*) - monitored
by a six-receiver array inside a small coastal marine protected area
(MPA). Coordinates, dates and detection patterns are simulated (no real
protected-species locations are involved), but are designed to be
ecologically plausible enough to exercise the whole moby workflow,
including per-species analyses via `id.groups`.

The objects come as a small family so that each tutorial can start from
the appropriate checkpoint without re-running earlier ones:

- `rays_detections` - raw, harmonised detections (a plain data frame).

- `rays_tags` - tag/animal metadata (one row per animal).

- `rays_deployments` - receiver-deployment log (one row per station).

- `rays` - the cleaned, analysis-ready
  [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
  checkpoint (detections with a `timebin` column and metadata attached,
  including the two-species `id.groups`).

Raw CSV versions for the import tutorial are also shipped: the
detections in a generic hand-mapped layout at
`system.file("extdata", "rays_detections.csv", package = "moby")`, and
the receiver-deployment log in raw ETN (European Tracking Network)
export format at
`system.file("extdata", "rays_deployments.csv", package = "moby")`
(harmonise it with `importDeployments(..., source = "etn")`).

## Usage

``` r
rays

rays_detections

rays_tags

rays_deployments
```

## Format

`rays_detections`: a data frame with one row per detection and columns
`ID`, `datetime` (POSIXct, UTC), `station`, `lon`, `lat`, `species`,
`receiver`, `transmitter`.

`rays_tags`: a data frame with columns `ID`, `transmitter`, `species`,
`tagging_date` (POSIXct), `tagging_station`.

`rays_deployments`: a data frame with columns `receiver`, `station`,
`lon`, `lat`, `deploy`, `recover` (POSIXct) and `depth` - the canonical
deployment-log schema (as produced by
[`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md))
used by
[`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md),
[`matchDeployments`](https://miguelgandra.github.io/moby/reference/matchDeployments.md)
and
[`plotDeployments`](https://miguelgandra.github.io/moby/reference/plotDeployments.md).

`rays`: a
[`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
(data-frame subclass) with the columns of `rays_detections` plus a
`timebin` column, and a `"moby"` metadata attribute (column mapping,
`epsg.code` (32629, UTM 29N), tagging dates and the per-species
`id.groups`).

An object of class `data.frame` with 1643 rows and 8 columns.

An object of class `data.frame` with 8 rows and 5 columns.

An object of class `data.frame` with 6 rows and 7 columns.

## Source

Simulated with `data-raw/make-toy-data.R` (seeded, reproducible). Not
derived from any real tracking study.

## Examples

``` r
# the analysis-ready mobyData checkpoint (metadata attached)
data(rays)
is_moby(rays)
#> [1] TRUE
mobyMeta(rays)$id.groups        # the two-species grouping
#> $`Raja clavata`
#> [1] "R01" "R02" "R03" "R04"
#> 
#> $`Dasyatis pastinaca`
#> [1] "D01" "D02" "D03" "D04"
#> 
head(rays)
#> <mobyData> 6 records x 9 columns
#>   individuals: 1  (id.col = 'ID')
#>   period: 2023-04-08 15:05:10 to 2023-04-09 12:54:09 (tz = UTC)
#>   columns: datetime.col='datetime', timebin.col='timebin', station.col='station', lon.col='lon', lat.col='lat'
#>   metadata: epsg=32629, tagging.dates (8), id.groups (2)

# the accompanying raw tables
head(rays_detections)           # harmonised detections (plain data frame)
#>    ID            datetime station    lon    lat            species  receiver
#> 1 D01 2023-04-08 15:05:10    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#> 2 D01 2023-04-08 15:14:59    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#> 3 D01 2023-04-08 15:52:10    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#> 4 D01 2023-04-08 16:51:29    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#> 5 D01 2023-04-08 16:51:56    ST03 -8.996 38.456 Dasyatis pastinaca VR2W-1003
#> 6 D01 2023-04-09 12:54:09    ST06 -8.990 38.442 Dasyatis pastinaca VR2W-1006
#>      transmitter
#> 1 A69-1602-30005
#> 2 A69-1602-30005
#> 3 A69-1602-30005
#> 4 A69-1602-30005
#> 5 A69-1602-30005
#> 6 A69-1602-30005
head(rays_tags)                 # one row per tagged animal
#>    ID    transmitter            species tagging_date tagging_station
#> 1 R01 A69-1602-30001       Raja clavata   2023-04-03            ST02
#> 2 R02 A69-1602-30002       Raja clavata   2023-04-08            ST03
#> 3 R03 A69-1602-30003       Raja clavata   2023-04-10            ST06
#> 4 R04 A69-1602-30004       Raja clavata   2023-04-01            ST03
#> 5 D01 A69-1602-30005 Dasyatis pastinaca   2023-04-07            ST03
#> 6 D02 A69-1602-30006 Dasyatis pastinaca   2023-04-02            ST02
head(rays_deployments)          # receiver-deployment log
#>    receiver station    lon    lat              deploy             recover depth
#> 1 VR2W-1001    ST01 -9.020 38.454 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 2 VR2W-1002    ST02 -9.008 38.464 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 3 VR2W-1003    ST03 -8.996 38.456 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 4 VR2W-1004    ST04 -8.985 38.466 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 5 VR2W-1005    ST05 -8.972 38.452 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 6 VR2W-1006    ST06 -8.990 38.442 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
```
