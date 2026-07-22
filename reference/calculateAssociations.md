# Calculate pairwise spatiotemporal associations (co-occurrence network)

Builds an individual **association network**, in which nodes are tagged
individuals and edges quantify the strength of their spatiotemporal
**co-occurrence**. Two individuals are considered to co-occur whenever
they are detected at the same receiver within the same time bin. The
edge weight is a co-occurrence-based association index — either the
simple-ratio index (SRI) or the half-weight index (HWI) — ranging from
0% (never co-occur) to 100% (always co-occur).

Note that associations here are defined operationally from
spatiotemporal co-occurrence and represent **shared space-and-time use,
not confirmed interaction**: a co-occurring pair may be anywhere from a
few centimetres to the receiver's detection range apart. Results should
be interpreted as a proxy for potential association, and statistically
tested against a null model with
[`randomizeAssociations`](https://miguelgandra.github.io/moby/reference/randomizeAssociations.md).

Each pair is compared only over its shared monitoring period (time bins
from the later of the two release dates to the earlier of the two
last-detection dates), so that an individual's absence is not mistaken
for non-association when it was due to transmitter failure or death.

## Usage

``` r
calculateAssociations(
  data,
  id.groups = NULL,
  subset = NULL,
  metric = "simple-ratio",
  group.comparisons = "all",
  cores = 1,
  verbose = getOption("moby.verbose", TRUE)
)
```

## Arguments

- data:

  A data frame containing binned detections in the wide format (time bin
  x individual matrix, with values corresponding to the receiver/station
  with the highest number of detections), as returned by
  [`createWideTable`](https://miguelgandra.github.io/moby/reference/createWideTable.md).

- id.groups:

  Optional. A list containing ID groups, used to calculate stats
  independently within each group, as well as comparing relationships
  between ids of different groups.

- subset:

  If defined, overlaps are calculated independently for each level of
  this variable. This can either be a single column name (variable) or a
  vector of column names, corresponding to variables contained in the
  data. In the case of multiple columns, their interaction is used for
  grouping. If left NULL, single overlap indices are calculated for the
  whole monitoring period.

- metric:

  One of "simple-ratio" or "half-weight".

- group.comparisons:

  Controls the type of comparisons to be run, when id.groups are
  defined. One of "within", "between" or "all". Useful to discard
  comparisons between individuals belonging to the same group or skip
  comparisons between different groups, when these are not required
  (less computing time). Defaults to "all".

- cores:

  Number of CPU cores to use for the computations. Defaults to 1, which
  means no parallel computing (single core). If set to a value greater
  than 1, the function will use parallel computing to speed up
  calculations.

- verbose:

  Logical; print progress and a summary to the console. Defaults to
  `getOption("moby.verbose", TRUE)`. Run
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to check the number of available cores.

## Value

A
[`mobyNetwork`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
object of type `"association"`: a data frame of network edges (one row
per dyad) carrying the node table and metadata as attributes. The edge
columns are:

- id1, id2:

  The two individuals forming the dyad (network nodes).

- association:

  The co-occurrence-based association index (SRI or HWI) for the pair,
  in percent.

- co_occurrences:

  Total number of time bins in which both individuals were detected
  together.

- shared_monit_days:

  Number of days the pair was monitored simultaneously.

- start, end:

  Start and end of the dyad's shared monitoring period.

- mean_consec_overlap, max_consec_overlap:

  Mean and maximum run of consecutive co-occurring time bins.

When `subset` is defined, results are stacked with a `subset` column;
when `id.groups` is defined, a comparison `type` column (and
`group1`/`group2`) is added. Use
[`networkEdges`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
/
[`networkNodes`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)
to extract the components, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) to draw the
network.

## Details

The two metrics for calculating association indices are:

1.  Simple Ratio Association Index: This index provides a measure of the
    proportion of time two individuals spend together. It is calculated
    as the ratio of the number of time bins where both individuals are
    detected together to the total number of time bins where at least
    one of the individuals is detected. This straightforward metric
    gives a clear indication of the extent to which two individuals'
    space usage overlaps.

2.  Half-Weight Index: This index modifies the simple ratio by giving
    half-weight to the occurrences where only one of the individuals is
    detected. It accounts for the fact that sometimes individuals might
    be in the same area but one might not be detected due to various
    reasons. This index is useful when detections are not perfectly
    reliable and helps mitigate the impact of missed detections.

It is important to acknowledge that any fish pair detected
simultaneously may be anywhere from a few centimetres to hundreds of
meters apart, depending on the maximum detection ranges of the
transmitters. Therefore, caution should be taken when interpreting the
results, especially when making inferences about biotic or fish-habitat
relationships.

Useful references:

- Cairns, S. J., & Schwager, S. J. (1987). A comparison of association
  indices. Animal Behaviour, 35(5), 1454-1469.

- Ginsberg, J. R., & Young, T. P. (1992). Measuring association between
  individuals or groups in behavioural studies. Animal Behaviour, 44(2),
  377-379.

- Farine, D. R., & Whitehead, H. (2015). Constructing, conducting and
  interpreting animal social network analysis. Journal of Animal
  Ecology, 84(5), 1144-1163.

- Hoppitt, W. J., & Farine, D. R. (2018). Association indices for
  quantifying social relationships: how to deal with missing
  observations of individuals or groups. Animal Behaviour, 136, 227-238.

## See also

[`randomizeAssociations`](https://miguelgandra.github.io/moby/reference/randomizeAssociations.md),
[`plotAssociations`](https://miguelgandra.github.io/moby/reference/plotAssociations.md),
[`mobyNetwork`](https://miguelgandra.github.io/moby/reference/mobyNetwork.md)

## Examples

``` r
data(rays)
# associations require the wide (time bin x individual) table as input
wide <- createWideTable(rays, value.col = "station")
#> Warning: - 'id.col' converted to factor.
#> Warning: 3 (ID, time-bin) combination(s) had multiple differing values; the first was kept. Aggregate upstream (e.g. calculateCOAs) to control this.
#> Tied (ID, time-bin) instances (first value kept):
#>                  timebin  ID                ties
#> 1898 2023-06-14 11:00:00 D03 ST01 (1) | ST06 (1)
#> 2163 2023-04-24 05:00:00 D04 ST01 (4) | ST05 (4)
#> 5302 2023-06-25 16:00:00 R04 ST03 (1) | ST05 (1)
assoc <- calculateAssociations(wide)
#> Calculating overlap - complete monitoring duration
#> Total execution time: 0.15 secs
assoc
#> <mobyNetwork> association network
#>   nodes: 8  |  edges: 28
#>   edge metric: simple-ratio
#>   edge table (use networkEdges()/networkNodes() to extract):
#>   id1 id2 association co_occurrences shared_monit_days      start        end
#> 1 D01 D02        0.00              0              80.6 2023-04-08 2023-06-28
#> 2 D01 D03        0.00              0              80.3 2023-04-09 2023-06-28
#> 3 D01 D04        0.57              1              76.9 2023-04-08 2023-06-24
#> 4 D01 R01        0.00              0              81.1 2023-04-08 2023-06-28
#> 5 D01 R02        1.62              3              81.0 2023-04-08 2023-06-28
#> 6 D01 R03        0.00              0              68.2 2023-04-21 2023-06-28
#>   mean_consec_overlap max_consec_overlap
#> 1                  NA                 NA
#> 2                  NA                 NA
#> 3                   1                  1
#> 4                  NA                 NA
#> 5                   3                  3
#> 6                  NA                 NA
head(networkEdges(assoc))
#>   id1 id2 association co_occurrences shared_monit_days      start        end
#> 1 D01 D02        0.00              0              80.6 2023-04-08 2023-06-28
#> 2 D01 D03        0.00              0              80.3 2023-04-09 2023-06-28
#> 3 D01 D04        0.57              1              76.9 2023-04-08 2023-06-24
#> 4 D01 R01        0.00              0              81.1 2023-04-08 2023-06-28
#> 5 D01 R02        1.62              3              81.0 2023-04-08 2023-06-28
#> 6 D01 R03        0.00              0              68.2 2023-04-21 2023-06-28
#>   mean_consec_overlap max_consec_overlap
#> 1                  NA                 NA
#> 2                  NA                 NA
#> 3                   1                  1
#> 4                  NA                 NA
#> 5                   3                  3
#> 6                  NA                 NA
```
