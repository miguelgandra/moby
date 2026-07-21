# Validate Column Presence and Format

This function checks whether a specified column name is valid, exists in
the provided dataset, and is correctly formatted as a single character
string. It returns informative error messages if the validation fails.

## Usage

``` r
.checkColumn(argument, col_label, data)
```

## Value

A character string containing an error message if validation fails, or
`NULL` if the column name is valid.

## Note

This function is intended for internal use within the `moby` package.
