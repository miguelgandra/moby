# Check receiver deployment metadata (quality control)

Quality-controls a receiver-deployment / station log and (optionally)
cross-checks it against a detection dataset, flagging the kinds of
issues that commonly corrupt acoustic telemetry analyses. This
complements
[`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
(which cleans the detections themselves) by auditing the *metadata*. It
generalises the receiver-log curation step that is otherwise done by
hand. The function only *reports* problems (it never edits data or
aborts on data issues), so it is a safe, explicit preprocessing step;
the companion
[`matchDeployments`](https://miguelgandra.github.io/moby/reference/matchDeployments.md)
applies the corrections.

Metadata-internal checks: missing deploy/recover dates; invalid date
ranges (recover before deploy); overlapping deployments on the same
receiver; duplicate deployment records; missing/implausible coordinates;
inconsistent station naming (one station name with divergent
coordinates, or near-identical coordinates under different names); and
gaps in a receiver's temporal coverage.

Detection-vs-metadata checks (when `detections` is supplied): receivers
present in the detections but absent from the metadata; detections
occurring before a receiver's first deployment, after its last recovery,
or within a gap between deployments; and station-name mismatches (a
known receiver associated with a station it was never deployed at).

## Usage

``` r
checkDeployments(
  deployments,
  detections = NULL,
  id.col = NULL,
  datetime.col = NULL,
  station.col = NULL,
  deploy.col = "deploy",
  recover.col = "recover",
  checks = "all",
  coord.tolerance = 500,
  gap.tolerance = 1,
  verbose = TRUE
)
```

## Arguments

- deployments:

  A receiver-deployment data frame, e.g. from
  [`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md),
  with columns `receiver`, `station`, `lon`, `lat`, `deploy` and
  (optionally) `recover`.

- detections:

  Optional. A detection dataset (`mobyData` or data frame) with
  `receiver`, `station` and a date-time column, used for the
  detection-vs-metadata checks.

- id.col, datetime.col, station.col:

  Column names in `detections` (resolved from the `mobyData` metadata or
  canonical defaults when `NULL`).

- deploy.col, recover.col:

  Names of the deployment and recovery date-time columns in
  `deployments`. Default to the canonical `"deploy"`/`"recover"` (as
  produced by
  [`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md));
  set them when a log uses other names (e.g. `"deploy_date"`).

- checks:

  Character vector selecting which check groups to run (any of, or
  `"all"`, the default): `"dates"` (deploy/recover date integrity),
  `"overlaps"` (duplicate and overlapping deployment records), `"gaps"`
  (coverage gaps between consecutive deployments - often benign
  servicing, so the easiest to drop), `"coordinates"` (implausible
  coordinates and station-vs-coordinate naming consistency; skipped if
  `lon`/`lat` are absent) and `"detections"` (cross-check detections
  against deployment windows; requires the `detections` argument).

- coord.tolerance:

  Numeric. Distance (in the coordinate units of the data; metres for
  geographic coordinates) beyond which coordinates sharing a station
  name are flagged as inconsistent, and below which differently-named
  stations are flagged as possible duplicates. Defaults to 500.

- gap.tolerance:

  Numeric. Minimum gap (in days) between consecutive deployments of a
  receiver to report as a coverage gap. Defaults to 1.

- verbose:

  Logical; print a summary to the console. Defaults to TRUE.

## Value

An object of class `mobyQC`: a list with

- report:

  A tidy data frame of flagged issues (`type`, `receiver`, `station`,
  `first`, `last`, `n_detections`, `n_individuals`, `details`).

- deployments:

  The input deployments with a cleaned `recover_clean` column.

- counts:

  Named integer vector of issue counts by type.

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
```
