# Match detections to receiver deployments and back-fill metadata

Joins a receiver-deployment log (see
[`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md))
to a detection dataset, matching each detection to the deployment window
(`receiver` + time within `[deploy, recover]`) it falls in, and
back-filling missing coordinates and station names from the receiver
log. When a receiver was moved between stations and a detection matches
more than one window, the window whose station name agrees with the
detection is preferred. Coordinate disagreements between the detection
export and the metadata (a common VUE artefact) are flagged. This
applies the corrections for the issues reported by
[`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md),
and is a natural follow-up to it.

## Usage

``` r
matchDeployments(
  detections,
  deployments,
  datetime.col = NULL,
  station.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  deploy.col = "deploy",
  recover.col = "recover",
  coord.tolerance = 500,
  fill.coords = TRUE,
  fill.station = TRUE,
  drop.unmatched = FALSE,
  verbose = TRUE
)
```

## Arguments

- detections:

  A detection dataset (`mobyData` or data frame) with a `receiver`
  column.

- deployments:

  A receiver-deployment data frame from
  [`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md)
  (`receiver`, `station`, `lon`, `lat`, `deploy`, and optionally
  `recover`).

- datetime.col, station.col, lon.col, lat.col:

  Column names in `detections` (resolved from the `mobyData` metadata or
  canonical defaults when `NULL`).

- deploy.col, recover.col:

  Names of the deployment and recovery date-time columns in
  `deployments`. Default to the canonical `"deploy"`/`"recover"`; set
  them when a log uses other names (e.g. `"deploy_date"`).

- coord.tolerance:

  Numeric. Distance (metres for geographic coordinates) above which a
  detection's own coordinates are flagged as disagreeing with the
  deployment metadata. The metadata coordinates are treated as
  authoritative for back-filling. Defaults to 500.

- fill.coords, fill.station:

  Logical; back-fill missing detection coordinates / station names from
  the matched deployment. Default TRUE.

- drop.unmatched:

  Logical; drop detections that fall outside every deployment window.
  Defaults to FALSE (they are retained and flagged).

- verbose:

  Logical; print a short summary. Defaults to TRUE.

## Value

A [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
object (the detections) with back-filled `lon`/`lat`/station where
missing, plus two logical flag columns: `deployment_matched` (detection
fell within a valid window) and `coord_mismatch` (own vs metadata
coordinates differ beyond `coord.tolerance`).

## See also

[`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md),
[`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md),
[`importDetections`](https://miguelgandra.github.io/moby/reference/importDetections.md)

## Examples

``` r
# match detections to deployment windows, back-filling station/coordinates
matched <- matchDeployments(rays_detections, rays_deployments, station.col = "station")
#> matchDeployments: 1643/1643 detections matched a deployment window (0 unmatched).
#>    0 coordinate(s) back-filled from metadata; 0 coordinate mismatch(es) > 500.
#> Note: the following mapped column(s) are not present in the data and will only matter for functions that use them: timebin.
table(matched$deployment_matched)
#> 
#> TRUE 
#> 1643 
```
