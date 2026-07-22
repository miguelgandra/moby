# Check receiver deployment metadata (quality control)

Performs quality control on receiver deployment metadata by identifying
inconsistencies that can compromise acoustic telemetry analyses. The
function checks the internal consistency of deployment records and,
optionally, cross-validates deployment metadata against a detection
dataset.

The function is intentionally non-destructive: it only reports potential
issues and never modifies input data. Reported issues can then be
reviewed and, where appropriate, resolved using
[`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md).

## Usage

``` r
checkDeployments(
  deployments,
  detections = NULL,
  id.col = NULL,
  datetime.col = NULL,
  station.col = NULL,
  deployment.station.col = "station",
  deployment.lon.col = "lon",
  deployment.lat.col = "lat",
  deployment.deploy.col = "deploy",
  deployment.recover.col = "recover",
  checks = "all",
  scope = c("all", "detected"),
  coord.tolerance = 500,
  gap.tolerance = 1,
  land.shape = NULL,
  epsg.code = NULL,
  land.tolerance = 500,
  verbose = TRUE
)
```

## Arguments

- deployments:

  A receiver-deployment data frame, e.g. from
  [`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md),
  with a `receiver` column plus station, deployment / recovery
  date-times and (optionally) longitude / latitude. Those columns may
  carry non-canonical names, resolved via the `deployment.*` arguments
  below.

- detections:

  Optional. A detection dataset (`mobyData` or data frame) with
  `receiver`, `station` and a date-time column, used for the
  detection-vs-metadata checks.

- id.col, datetime.col, station.col:

  Column names in the **detection** dataset (`detections`), resolved
  from its `mobyData` metadata or canonical defaults when `NULL`. Bare
  `*.col` arguments always refer to the detections; the deployment log's
  columns use the `deployment.*` arguments.

- deployment.station.col, deployment.lon.col, deployment.lat.col:

  Names of the station, longitude and latitude columns in the
  receiver-deployment log (`deployments`). Default to the canonical
  `"station"`/`"lon"`/`"lat"` produced by
  [`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md);
  set them when a hand-made log uses other names (e.g.
  `deployment.lon.col = "Longitude"`). The `deployment.` prefix marks
  these as deployment-log columns, keeping them distinct from the bare
  `*.col` arguments, which always refer to the detection dataset. The
  `receiver` column is the canonical join key and is always taken as-is.

- deployment.deploy.col, deployment.recover.col:

  Names of the deployment and recovery date-time columns in the
  receiver-deployment log (`deployments`). Default to the canonical
  `"deploy"`/`"recover"` produced by
  [`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md);
  set them when a hand-made log uses other names (e.g.
  `deployment.deploy.col = "deploy_date"`). The `deployment.` prefix
  marks these as deployment-log columns, keeping them distinct from the
  bare `*.col` arguments, which always refer to the detection dataset.

- checks:

  Character vector selecting which check groups to run (any of, or
  `"all"`, the default): `"dates"` (deploy/recover date integrity),
  `"overlaps"` (duplicate and overlapping deployment records), `"gaps"`
  (coverage gaps between consecutive deployments - often benign
  servicing, so the easiest to drop), `"coordinates"` (implausible
  coordinates, station-vs-coordinate naming consistency, and - when
  `land.shape` is supplied - receiver positions on land; skipped if
  `lon`/`lat` are absent) and `"detections"` (cross-check detections
  against deployment windows; requires the `detections` argument).

- scope:

  Character; how much of the deployment log the *metadata-internal*
  checks cover. `"all"` (default) audits every receiver in
  `deployments`. `"detected"` restricts those checks to receivers that
  appear in `detections` (i.e. that recorded at least one detection),
  which removes warnings about receivers that are irrelevant to the data
  being analysed - typically the bulk of the coverage-gap and
  duplicate-station noise. Each retained receiver keeps its full
  deployment timeline, so gap and overlap logic stays correct. `scope`
  only affects the metadata-internal checks; the `"detections"`
  cross-checks are inherently detection-scoped and always run against
  the complete metadata (so a receiver present in the detections but
  missing from the log is still flagged). `"detected"` needs
  `detections`; without it, the function warns and audits everything.

- coord.tolerance:

  Numeric. Distance (in the coordinate units of the data; metres for
  geographic coordinates) beyond which coordinates sharing a station
  name are flagged as inconsistent, and below which differently-named
  stations are flagged as possible duplicates. Defaults to 500.

- gap.tolerance:

  Numeric. Minimum gap (in days) between consecutive deployments of a
  receiver to report as a coverage gap. Defaults to 1.

- land.shape:

  Optional `sf` (or `SpatialPolygons*`) polygon layer of landmasses.
  When supplied, the `"coordinates"` group additionally flags receiver
  positions that fall on land. Off by default; if `NULL`, it is taken
  from the `detections` `mobyData` metadata when present (i.e. from
  `as_moby(detections, land.shape = ...)`). The layer must carry a
  coordinate reference system; if it lacks one, or cannot be reconciled
  with the deployment coordinates, the on-land check is skipped with a
  message (the audit never aborts).

- epsg.code:

  Optional EPSG code for the deployment coordinates, used only by the
  on-land check. Needed only when the coordinates are projected (not
  longitude/latitude); geographic coordinates are assumed to be WGS84.

- land.tolerance:

  Numeric. How far inside the coastline (in metres) a receiver position
  must lie to be reported as on land. This spares genuine near-shore
  receivers that a coarse coastline overlaps by a small margin; gross
  coordinate-entry errors (swapped lon/lat, a sign flip) fall kilometres
  inland and are still flagged. Defaults to 500.

- verbose:

  Logical; print a summary to the console. Defaults to TRUE.

## Value

An object of class `mobyQC`: a list with

- report:

  A tidy data frame of flagged issues (`type`, `receiver`, `station`,
  `first`, `last`, `n_detections`, `n_individuals`, `details`). Columns
  that do not apply to a given issue are `NA` (kept typed, so the frame
  stays usable programmatically); for a more readable manual export,
  write with an explicit placeholder, e.g.
  `write.csv(x$report, "report.csv", row.names = FALSE, na = "-")`.

- deployments:

  The input deployments with a cleaned `recover_clean` column.

- counts:

  Named integer vector of issue counts by type.

## Details

Internal metadata checks include missing or invalid deployment dates,
overlapping deployments, duplicate records, coordinate inconsistencies,
station naming issues, and gaps in receiver coverage. When a land layer
is supplied (`land.shape`), the function additionally flags receiver
positions that fall on land - a common consequence of a coordinate-entry
error (swapped longitude/latitude, a wrong sign, a decimal typo).

When detections are provided, the function additionally checks whether
detections are consistent with the deployment history, including unknown
receivers, detections outside deployment periods, and station
mismatches.

## See also

[`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md),
[`importDetections`](https://miguelgandra.github.io/moby/reference/importDetections.md),
[`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)

## Examples

``` r
# quality-control a receiver-deployment log
checkDeployments(rays_deployments)
#> - 'detections' checks requested but no 'detections' supplied; skipping the detection-vs-metadata checks.
#> <mobyQC> deployment metadata quality-control report
#>   6 deployment records | 6 receivers | 6 stations
#>   No issues flagged.

# also cross-check the detections against the deployment windows
data(rays)
checkDeployments(rays_deployments, detections = rays)
#> <mobyQC> deployment metadata quality-control report
#>   6 deployment records | 6 receivers | 6 stations
#>   No issues flagged.

# restrict the metadata checks to receivers that actually recorded detections
checkDeployments(rays_deployments, detections = rays, scope = "detected")
#> <mobyQC> deployment metadata quality-control report
#>   6 deployment records | 6 receivers | 6 stations
#>   No issues flagged.
```
