# Changelog

## moby 1.0.0

First stable release of moby, a toolkit for processing, analysing and
visualising passive acoustic telemetry data. It supersedes the v0.1.0
beta (2024-12-03), which was distributed for internal testing only.
v1.0.0 is a ground-up rewrite: the spatial stack moved from
`raster`/`sp`/`gdistance` to `terra`/`igraph`, the dependency closure
was cut from 42 packages to 29, and the API was reorganised around a
`mobyData` object.

### Breaking changes

**v1.0.0 is not backward compatible with v0.1.0.** No deprecation shims
or [`.Defunct()`](https://rdrr.io/r/base/Defunct.html) aliases are
provided, so scripts written against the beta will fail with
`could not find function`. To reproduce an analysis run under v0.1.0,
install that tag explicitly:

``` r

remotes::install_github("miguelgandra/moby@v0.1.0")
```

Renamed functions:

| v0.1.0 | v1.0.0 |
|----|----|
| `calculateKUDs()` | [`calculateUDs()`](https://miguelgandra.github.io/moby/reference/calculateUDs.md) |
| `calculateOverlap()` | [`calculateAssociations()`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md) |
| `randomizeOverlap()` | [`randomizeAssociations()`](https://miguelgandra.github.io/moby/reference/randomizeAssociations.md) |
| `calculateTracks()` | [`calculateStepDistances()`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md) (trajectories via [`getTrajectories()`](https://miguelgandra.github.io/moby/reference/getTrajectories.md)) |
| `migrationsTable()` | [`transitionsTable()`](https://miguelgandra.github.io/moby/reference/transitionsTable.md) |
| `plotMigrations()` | [`plotMovements()`](https://miguelgandra.github.io/moby/reference/plotMovements.md) |
| `plotOverlap()` | [`plotAssociations()`](https://miguelgandra.github.io/moby/reference/plotAssociations.md) |
| `plotOverlapTable()` | [`plotAssociationMatrix()`](https://miguelgandra.github.io/moby/reference/plotAssociationMatrix.md) |
| `plotGroupSizes()` | [`plotGroupSizeDistribution()`](https://miguelgandra.github.io/moby/reference/plotGroupSizeDistribution.md) |
| `plotBoxplots()` | [`plotMetricComparison()`](https://miguelgandra.github.io/moby/reference/plotMetricComparison.md) |
| `plotCWTs()` | [`plotScalogram()`](https://miguelgandra.github.io/moby/reference/plotScalogram.md) |
| `plotFFTs()` | [`plotPeriodogram()`](https://miguelgandra.github.io/moby/reference/plotPeriodogram.md) |
| `setDefaults()`, `getDefaults()` | [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md), [`mobyMeta()`](https://miguelgandra.github.io/moby/reference/mobyMeta.md) |

Removed without a direct replacement: `calculateFrequencies()` and
`plotDetections()`. The nearest equivalents are
[`createWideTable()`](https://miguelgandra.github.io/moby/reference/createWideTable.md)
and
[`plotChronogram()`](https://miguelgandra.github.io/moby/reference/plotChronogram.md)/[`plotAbacus()`](https://miguelgandra.github.io/moby/reference/plotAbacus.md)
respectively.

### Correction affecting v0.1.0 results

In v0.1.0 the land layer was rasterised in a way that left all but the
first polygon passable, so least-cost in-water distances were computed
as if only one landmass blocked movement. Analyses run under v0.1.0 with
a **multi-polygon** land layer are affected — including in-water step
distances and anything derived from them (rate of movement, linearity
index). Single-polygon land layers are unaffected. v1.0.0 treats every
polygon as impassable and pins the behaviour with a regression test.

### Data model

- [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
  creates a `mobyData` object — a `data.frame` subclass that carries the
  dataset’s metadata (column mapping, coordinate reference system,
  tagging dates, ID groups, land layer) as attributes. Downstream
  functions read this metadata automatically, so column names, EPSG
  codes and tagging dates do not need to be repeated on every call.
  Plain data frames are also fully supported, falling back to
  conventional column names (`ID`, `datetime`, `timebin`, `station`,
  `lon`, `lat`). Helpers:
  [`mobyMeta()`](https://miguelgandra.github.io/moby/reference/mobyMeta.md),
  [`is_moby()`](https://miguelgandra.github.io/moby/reference/is_moby.md).

### Import, quality control and curation

- [`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md),
  [`importDeployments()`](https://miguelgandra.github.io/moby/reference/importDeployments.md)
  and
  [`importTags()`](https://miguelgandra.github.io/moby/reference/importTags.md)
  read and harmonise detections, receiver-deployment logs and tag
  metadata from common sources (Innovasea/VEMCO VUE, Innovasea
  VDAT/Fathom, GLATOS, OTN and ETN), plus a `generic` mode with
  user-supplied column maps.
- [`checkDeployments()`](https://miguelgandra.github.io/moby/reference/checkDeployments.md)
  audits receiver-deployment metadata (overlapping deployments,
  missing/invalid dates, implausible coordinates, duplicate records,
  inconsistent station naming, coverage gaps) and cross-checks it
  against detections, returning a structured report; a `checks` argument
  selects which check groups to run.
- [`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md)
  matches detections to deployment windows and back-fills coordinates
  and station names;
  [`assignAnimalIDs()`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
  joins tag metadata to assign animal IDs and tagging dates.
- [`plotArray()`](https://miguelgandra.github.io/moby/reference/plotArray.md)
  maps a receiver array from a deployment log — station positions,
  nominal detection ranges, deployment status and per-station detection
  effort.
- [`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
  removes spurious detections through a configurable pipeline:
  exact-duplicate removal, tagging/cut-off dates, a short-interval
  (`min_lag`) false-detection filter (Pincock 2012), a
  corroboration-shielded swim-speed filter (land-aware, in-water
  distances), and minimum-detection/day criteria, returning a classed
  `mobyFilter` object.
  [`plotMinLag()`](https://miguelgandra.github.io/moby/reference/plotMinLag.md)
  helps choose the `min_lag` threshold from a dataset’s empirical
  distribution.

### Temporal and spatial classification

- [`getTimeBins()`](https://miguelgandra.github.io/moby/reference/getTimeBins.md),
  [`getDielPhase()`](https://miguelgandra.github.io/moby/reference/getDielPhase.md),
  [`getSeason()`](https://miguelgandra.github.io/moby/reference/getSeason.md),
  [`getReprodPeriod()`](https://miguelgandra.github.io/moby/reference/getReprodPeriod.md)
  assign temporal classes;
  [`getSunTimes()`](https://miguelgandra.github.io/moby/reference/getSunTimes.md)
  returns diel boundary times.
- [`calculateCOAs()`](https://miguelgandra.github.io/moby/reference/calculateCOAs.md)
  (centres of activity);
  [`calculateStepDistances()`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
  reconstructs movement paths (shortest in-water or great-circle),
  returning the data with a `dist_m` column added and the trajectories
  attached as an attribute (retrieved with
  [`getTrajectories()`](https://miguelgandra.github.io/moby/reference/getTrajectories.md)).
- [`interpolateDistances()`](https://miguelgandra.github.io/moby/reference/interpolateDistances.md),
  [`calculateLandDists()`](https://miguelgandra.github.io/moby/reference/calculateLandDists.md),
  [`correctPositions()`](https://miguelgandra.github.io/moby/reference/correctPositions.md),
  [`createWideTable()`](https://miguelgandra.github.io/moby/reference/createWideTable.md).

### Summaries, residency and home range

- [`calculateResidency()`](https://miguelgandra.github.io/moby/reference/calculateResidency.md)
  returns residency indices (IR1, IR2, IWR) and their building blocks as
  a tidy numeric table;
  [`summaryTable()`](https://miguelgandra.github.io/moby/reference/summaryTable.md)
  formats these (and metadata, sensors) for publication.
- [`calculateVisits()`](https://miguelgandra.github.io/moby/reference/calculateVisits.md)
  extracts discrete residence events (arrival, departure, duration) per
  animal and site, with a configurable maximum within-visit gap.
- [`calculateROM()`](https://miguelgandra.github.io/moby/reference/calculateROM.md)
  (total distance, mean/max rate of movement) and
  [`calculateLinearityIndex()`](https://miguelgandra.github.io/moby/reference/calculateLinearityIndex.md)
  (movement directness) return tidy numeric per-animal tables;
  [`movementTable()`](https://miguelgandra.github.io/moby/reference/movementTable.md)
  formats these (with home-range areas) for publication.
- [`calculateUDs()`](https://miguelgandra.github.io/moby/reference/calculateUDs.md)
  estimates utilization distributions, defaulting to **autocorrelated
  kernel density estimation (AKDE)** with confidence intervals (via
  `ctmm`), with classic fixed-bandwidth KDE available;
  [`calculateUDOverlap()`](https://miguelgandra.github.io/moby/reference/calculateUDOverlap.md)
  quantifies pairwise home-range overlap;
  [`plotMaps()`](https://miguelgandra.github.io/moby/reference/plotMaps.md).

### Network analysis

- Two symmetric pipelines over a shared `mobyNetwork` object:
  - **Association networks** (individuals):
    [`calculateAssociations()`](https://miguelgandra.github.io/moby/reference/calculateAssociations.md),
    [`randomizeAssociations()`](https://miguelgandra.github.io/moby/reference/randomizeAssociations.md),
    [`plotAssociations()`](https://miguelgandra.github.io/moby/reference/plotAssociations.md),
    [`plotAssociationMatrix()`](https://miguelgandra.github.io/moby/reference/plotAssociationMatrix.md).
  - **Movement networks** (locations):
    [`calculateTransitions()`](https://miguelgandra.github.io/moby/reference/calculateTransitions.md),
    [`randomizeTransitions()`](https://miguelgandra.github.io/moby/reference/randomizeTransitions.md),
    [`plotMovements()`](https://miguelgandra.github.io/moby/reference/plotMovements.md),
    [`transitionsTable()`](https://miguelgandra.github.io/moby/reference/transitionsTable.md).
  - [`networkMetrics()`](https://miguelgandra.github.io/moby/reference/networkMetrics.md)
    computes node- and network-level graph metrics (degree/strength,
    betweenness, eigenvector, community/modularity, reciprocity) for
    both, and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
    draws either network.
- Co-occurrence-based association indices (SRI/HWI) are interpreted as
  shared space-and-time use, not confirmed interaction; randomisation
  tests use multiple-comparison correction.

### Visualization

A consistent family of base-graphics plotting functions sharing a
device-stable design: a single global `cex`, colourblind-safe palettes
(Okabe-Ito categorical, viridis continuous), an optional
`file`/`width`/`height`/`res` export, and an invisible tidy return.

- **Detection & effort over time:**
  [`plotAbacus()`](https://miguelgandra.github.io/moby/reference/plotAbacus.md)
  (detection abacus),
  [`plotDeployments()`](https://miguelgandra.github.io/moby/reference/plotDeployments.md)
  (receiver operating-period timeline),
  [`plotActograms()`](https://miguelgandra.github.io/moby/reference/plotActograms.md),
  [`plotChronogram()`](https://miguelgandra.github.io/moby/reference/plotChronogram.md),
  [`plotContours()`](https://miguelgandra.github.io/moby/reference/plotContours.md).
- **Per-station & distributions:**
  [`plotStationStats()`](https://miguelgandra.github.io/moby/reference/plotStationStats.md),
  [`plotGroupSizeDistribution()`](https://miguelgandra.github.io/moby/reference/plotGroupSizeDistribution.md),
  [`plotMetricComparison()`](https://miguelgandra.github.io/moby/reference/plotMetricComparison.md)
  (per-individual metrics across grouping levels, with
  design-appropriate repeated-measures tests and effect sizes).
- **Rhythms & periodicity:**
  [`plotPeriodogram()`](https://miguelgandra.github.io/moby/reference/plotPeriodogram.md)
  (FFT periodogram / opt-in Lomb-Scargle) and
  [`plotScalogram()`](https://miguelgandra.github.io/moby/reference/plotScalogram.md)
  (continuous wavelet scalogram).
- **Spatial & network:**
  [`plotArray()`](https://miguelgandra.github.io/moby/reference/plotArray.md)
  (receiver-array map),
  [`plotMaps()`](https://miguelgandra.github.io/moby/reference/plotMaps.md)
  (home-range utilisation isopleths),
  [`plotMovements()`](https://miguelgandra.github.io/moby/reference/plotMovements.md),
  [`plotAssociations()`](https://miguelgandra.github.io/moby/reference/plotAssociations.md),
  [`plotAssociationMatrix()`](https://miguelgandra.github.io/moby/reference/plotAssociationMatrix.md);
  [`shadeSeasons()`](https://miguelgandra.github.io/moby/reference/shadeSeasons.md)
  adds seasonal background shading.
