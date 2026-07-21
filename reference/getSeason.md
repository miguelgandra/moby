# Estimate annual season

This function determines the meteorological season for a given date
based on the meteorological definition, which divides the year into four
seasons of three months each:  

Northern Hemisphere:  

- Spring: March - May  

- Summer: June - August  

- Autumn: September - November  

- Winter: December - February  

Southern Hemisphere:  

- Spring: September - November  

- Summer: December - February  

- Autumn: March - May  

- Winter: June - August  

## Usage

``` r
getSeason(datetimes, hemisphere = "Northern")
```

## Arguments

- datetimes:

  A POSIXct object containing the respective datetimes or time-bins.

- hemisphere:

  A character string specifying the hemisphere ("Northern" or
  "Southern").

## Value

A factor indicating the season.

## Examples

``` r
datetimes <- as.POSIXct("2024-05-30")
getSeason(datetimes, hemisphere="Northern")
#> [1] spring
#> Levels: spring summer autumn winter
```
