# Validate and Process Animal-Specific Parameters

This function validates and processes an argument intended for
animal-specific parameters, ensuring it conforms to the required
structure and class. It handles single values, unnamed vectors, and
named vectors while issuing appropriate warnings or errors as needed.

## Usage

``` r
.checkAnimalParams(
  argument,
  arg_label,
  expected_class = c("POSIXct", "numeric"),
  data,
  id_col,
  allow.missing = FALSE
)
```

## Value

A list with three elements:

- `vector`: The processed argument (replicated, reordered, or unchanged
  as appropriate).

- `errors`: A character vector of error messages, or `NULL` if no errors
  are found.

- `warnings`: A character vector of warning messages, or `NULL` if no
  warnings are found.

## Note

This function is intended for internal use within the `moby` package.
