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
  datetime.format = NULL
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

  Optional named list mapping canonical fields (`receiver`, `station`,
  `lon`, `lat`, `deploy`, `recover`, `depth`) to source column name(s);
  merged over the `source` preset.

- datetime.format:

  Optional explicit `strptime` format for the deploy/recover columns.

## Value

A data frame with columns `receiver`, `station`, `lon`, `lat`, `deploy`
(POSIXct), `recover` (POSIXct) and, where available, `depth`; sorted by
receiver and deployment date.

## See also

[`importDetections`](https://miguelgandra.github.io/moby/reference/importDetections.md),
[`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md)

## Examples

``` r
# read a raw ETN deployment export and harmonise it
csv <- system.file("extdata", "rays_deployments.csv", package = "moby")
deployments <- importDeployments(csv, source = "etn")
head(deployments)
#>    receiver station    lon    lat              deploy             recover depth
#> 1 VR2W-1001    ST01 -9.020 38.454 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 2 VR2W-1002    ST02 -9.008 38.464 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 3 VR2W-1003    ST03 -8.996 38.456 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 4 VR2W-1004    ST04 -8.985 38.466 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 5 VR2W-1005    ST05 -8.972 38.452 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
#> 6 VR2W-1006    ST06 -8.990 38.442 2023-03-25 09:00:00 2023-07-05 15:00:00    NA
```
