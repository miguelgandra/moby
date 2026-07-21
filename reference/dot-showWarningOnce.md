# Display a warning only once per session

This function checks if a warning has already been shown in the current
session. If not, it displays the warning message and ensures that the
warning will not be shown again during the same session.

## Usage

``` r
.showWarningOnce(arg_name, message)
```

## Arguments

- message:

  A character string containing the warning message to be displayed.

## Note

This function is intended for internal use within the `moby` package.
