# Validate Arguments

This internal function performs a comprehensive set of validation checks
on the arguments provided to functions within the `moby` package. It
ensures that the data and various parameters are in the correct format
and meet the expected criteria, issuing errors and warnings as
necessary. The function tries to issue as many errors and warnings at
the same time as possible, to speed up the completion of required
changes and reduce the time it takes to get the function call right.

## Usage

``` r
.validateArguments()
```

## Value

A list containing potentially modified versions of `data`,
`tagging.dates`, and `tag.durations`. If validation fails, the function
will stop and return an error message.

## Details

The `.validateArguments` function performs rigorous checks to ensure
that all input arguments are in the correct format, of the correct type,
and consistent with the requirements of the `moby` package. This
includes checks for data structure, column specifications, and
additional parameters. The validation aims to provide detailed error and
warning messages to guide users in correcting their input.

\#' Specifically, the following checks are performed:

- `data`:

  Ensures the first argument is a data frame.

- `id.col`:

  Validates that the ID column exists in the data frame, is unique for
  each individual, and is converted to a factor if necessary.

- `id.groups`:

  Ensures ID groups are unique and correctly referenced in the data.

- `datetime.col`:

  Checks that the datetime column exists and is in POSIXct format.

- `timebin.col`:

  Validates the existence of the time-bin column and checks that it is
  in POSIXct format.

- `lon.col` and `lat.col`:

  Ensures the longitude and latitude columns exist and contain numeric
  values.

- `station.col`:

  Checks for the existence of the station column.

- `spatial.col`:

  Validates spatial column(s), ensuring they exist and match the
  required format for spatial analyses.

- `dist.col`:

  Checks the distance column, ensuring it contains numeric values.

- `split.by`:

  Validates the grouping variable(s) used to split the data.

- `subset`:

  Checks that subset criteria are valid and appropriately referenced in
  the data.

- `variable`:

  Ensures a single variable is specified and exists in the data.

- `variables`:

  Validates that all specified variables exist in the data.

- `color.by`:

  Checks that the color-by column exists, contains no missing values,
  and is converted to a factor if needed.

- `tagging.dates`:

  Ensures tagging dates are in POSIXct format, match the unique IDs, and
  are appropriately named if provided as a vector.

- `tag.durations`:

  Validates tag durations, ensuring they are numeric and consistent with
  the unique IDs.

- `start.dates` and `end.dates`:

  Checks start and end dates, ensuring they are in POSIXct format and
  consistent with the data.

- `cutoff.dates`:

  Validates cutoff dates, ensuring they are in POSIXct format and
  logically set.

- `last.monitoring.date`:

  Ensures the last monitoring date is in POSIXct format and
  appropriately set relative to the data.

- `land.shape`:

  Checks that the provided spatial object is either an `sf` object or
  convertible to one.

- `color.pal`:

  Validates the color palette, ensuring it is a character vector or a
  function.

- `diel.lines`:

  Checks that diel line values are valid (e.g., 0, 2, or 4).

- `polygons`:

  Ensures polygons are validly set to `'diel'`, `'season'`, or `FALSE`.

- `background.color`:

  Checks that the background color is a valid named color or hexadecimal
  code.

- `style`:

  Validates that the style parameter is either `'raster'` or `'points'`.

- `cores`:

  Ensures that the number of cores for parallel processing is specified
  as a positive integer. If this parameter is missing or invalid, a
  default value of 1 is used.

- `ncol`:

  Checks that the `ncol` parameter is correctly formatted as a single
  positive integer.

## Note

This function is intended for internal use within the `moby` package. It
helps streamline argument validation across multiple functions, ensuring
consistency and reducing redundancy.
