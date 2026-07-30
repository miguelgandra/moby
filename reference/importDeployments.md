# Import and harmonise receiver deployment metadata

Reads a receiver-deployment / station log from a range of common sources
and harmonises it into a consistent schema for use with
[`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md)
and the rest of the `moby` workflow.

## Usage

``` r
importDeployments(
  x,
  source = c("vue", "glatos", "otn", "etn", "generic"),
  tz = "UTC",
  col.map = NULL,
  datetime.format = NULL,
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- x:

  A path to a `.csv`/`.xlsx` deployment log, or a data frame (e.g. the
  output of `etn::get_acoustic_deployments()`).

- source:

  One of `"vue"`, `"glatos"`, `"otn"`, `"etn"` or `"generic"`.

- tz:

  Time zone used to parse deploy/recover date-times. Defaults to
  `"UTC"`.

- col.map:

  Optional named list mapping canonical deployment fields to source
  column name(s), merged over the `source` preset. See
  [moby_import_schema](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
  for the full list of deployment fields and which are required.

- datetime.format:

  Optional explicit `strptime` format for the deploy/recover columns.
  When `NULL` (default) the layout is inferred, but only when the answer
  is unambiguous: a column such as `"07/06/2013"`, where day-first and
  month-first both fit and disagree, raises an error asking for this
  argument rather than silently choosing one. A column no layout can
  read is an error too, so an unreadable column can never travel on as
  silent `NA`s.

  Timezone handling depends on the input type. TEXT carries no zone, so
  `tz` is applied as the zone while parsing. A column that is ALREADY
  `POSIXct` (as when a data frame is passed in, e.g. from an API) is an
  absolute instant chosen by the caller: it is never reinterpreted, and
  `tz` changes only how it is displayed.

- verbose:

  Logical; print a summary of the operation. Defaults to
  `getOption("moby.verbose", TRUE)`.

## Value

A data frame with columns `receiver`, `station`, `lon`, `lat`, `deploy`
(POSIXct), `recover` (POSIXct) and, where available, `depth`; sorted by
receiver and deployment date.

## See also

[`moby_import_schema`](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
for the canonical field list;
[`importDetections`](https://miguelgandra.github.io/moby/reference/importDetections.md),
[`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md)

## Examples

``` r
# read a raw ETN deployment export and harmonise it
csv <- system.file("extdata", "rays_deployments.csv", package = "moby")
deployments <- importDeployments(csv, source = "etn")
#> ── importDeployments() ───────────────────────────────────────────────── moby ──
#> 
#> ℹ Reading and harmonising the receiver deployment log
#> • Input: rays_deployments.csv
#> 
#> → Method
#>   • source    etn preset
#>   • timezone  UTC
#> 
#> ✔ 6 deployment records imported across 6 receivers
#> ℹ Fields resolved: receiver, station, lat, lon, deploy, recover, depth
head(deployments)
#>    receiver station    lon    lat              deploy             recover depth
#> 1 VR2W-1001    ST01 -9.020 38.454 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 2 VR2W-1002    ST02 -9.008 38.464 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 3 VR2W-1003    ST03 -8.996 38.456 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 4 VR2W-1004    ST04 -8.985 38.466 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 5 VR2W-1005    ST05 -8.972 38.452 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 6 VR2W-1006    ST06 -8.990 38.442 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
```
