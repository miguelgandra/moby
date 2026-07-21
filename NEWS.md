# moby 1.0.0

First stable release of moby, a toolkit for processing, analysing and visualising passive acoustic
telemetry data. It supersedes the v0.1.0 beta (2024-12-03), which was distributed for internal
testing only. v1.0.0 is a ground-up rewrite: the spatial stack moved from `raster`/`sp`/`gdistance`
to `terra`/`igraph`, the dependency closure was cut from 42 packages to 29, and the API was
reorganised around a `mobyData` object.

## Breaking changes

**v1.0.0 is not backward compatible with v0.1.0.** No deprecation shims or `.Defunct()` aliases are
provided, so scripts written against the beta will fail with `could not find function`. To reproduce
an analysis run under v0.1.0, install that tag explicitly:

``` r
remotes::install_github("miguelgandra/moby@v0.1.0")
```

Renamed functions:

| v0.1.0 | v1.0.0 |
|---|---|
| `calculateKUDs()` | `calculateUDs()` |
| `calculateOverlap()` | `calculateAssociations()` |
| `randomizeOverlap()` | `randomizeAssociations()` |
| `calculateTracks()` | `calculateStepDistances()` (trajectories via `getTrajectories()`) |
| `migrationsTable()` | `transitionsTable()` |
| `plotMigrations()` | `plotMovements()` |
| `plotOverlap()` | `plotAssociations()` |
| `plotOverlapTable()` | `plotAssociationMatrix()` |
| `plotGroupSizes()` | `plotGroupSizeDistribution()` |
| `plotBoxplots()` | `plotMetricComparison()` |
| `plotCWTs()` | `plotScalogram()` |
| `plotFFTs()` | `plotPeriodogram()` |
| `setDefaults()`, `getDefaults()` | `as_moby()`, `mobyMeta()` |

Removed without a direct replacement: `calculateFrequencies()` and `plotDetections()`. The nearest
equivalents are `createWideTable()` and `plotChronogram()`/`plotAbacus()` respectively.

## Correction affecting v0.1.0 results

In v0.1.0 the land layer was rasterised in a way that left all but the first polygon passable, so
least-cost in-water distances were computed as if only one landmass blocked movement. Analyses run
under v0.1.0 with a **multi-polygon** land layer are affected — including in-water step distances and
anything derived from them (rate of movement, linearity index). Single-polygon land layers are
unaffected. v1.0.0 treats every polygon as impassable and pins the behaviour with a regression test.

## Data model

* `as_moby()` creates a `mobyData` object — a `data.frame` subclass that carries the dataset's
  metadata (column mapping, coordinate reference system, tagging dates, ID groups, land layer)
  as attributes. Downstream functions read this metadata automatically, so column names, EPSG
  codes and tagging dates do not need to be repeated on every call. Plain data frames are also
  fully supported, falling back to conventional column names (`ID`, `datetime`, `timebin`,
  `station`, `lon`, `lat`). Helpers: `mobyMeta()`, `is_moby()`.

## Import, quality control and curation

* `importDetections()`, `importDeployments()` and `importTags()` read and harmonise detections,
  receiver-deployment logs and tag metadata from common sources (Innovasea/VEMCO VUE, Innovasea
  VDAT/Fathom, GLATOS, OTN and ETN), plus a `generic` mode with user-supplied column maps.
* `checkDeployments()` audits receiver-deployment metadata (overlapping deployments, missing/invalid
  dates, implausible coordinates, duplicate records, inconsistent station naming, coverage gaps)
  and cross-checks it against detections, returning a structured report; a `checks` argument selects
  which check groups to run.
* `matchDeployments()` matches detections to deployment windows and back-fills coordinates and
  station names; `assignAnimalIDs()` joins tag metadata to assign animal IDs and tagging dates.
* `plotArray()` maps a receiver array from a deployment log — station positions, nominal detection
  ranges, deployment status and per-station detection effort.
* `filterDetections()` removes spurious detections through a configurable pipeline: exact-duplicate
  removal, tagging/cut-off dates, a short-interval (`min_lag`) false-detection filter (Pincock 2012),
  a corroboration-shielded swim-speed filter (land-aware, in-water distances), and
  minimum-detection/day criteria, returning a classed `mobyFilter` object. `plotMinLag()` helps
  choose the `min_lag` threshold from a dataset's empirical distribution.

## Temporal and spatial classification

* `getTimeBins()`, `getDielPhase()`, `getSeason()`, `getReprodPeriod()` assign temporal classes;
  `getSunTimes()` returns diel boundary times.
* `calculateCOAs()` (centres of activity); `calculateStepDistances()` reconstructs movement paths
  (shortest in-water or great-circle), returning the data with a `dist_m` column added and the
  trajectories attached as an attribute (retrieved with `getTrajectories()`).
* `interpolateDistances()`, `calculateLandDists()`, `correctPositions()`, `createWideTable()`.

## Summaries, residency and home range

* `calculateResidency()` returns residency indices (IR1, IR2, IWR) and their building blocks as a
  tidy numeric table; `summaryTable()` formats these (and metadata, sensors) for publication.
* `calculateVisits()` extracts discrete residence events (arrival, departure, duration) per animal
  and site, with a configurable maximum within-visit gap.
* `calculateROM()` (total distance, mean/max rate of movement) and
  `calculateLinearityIndex()` (movement directness) return tidy numeric per-animal tables;
  `movementTable()` formats these (with home-range areas) for publication.
* `calculateUDs()` estimates utilization distributions, defaulting to **autocorrelated kernel
  density estimation (AKDE)** with confidence intervals (via `ctmm`), with classic fixed-bandwidth
  KDE available; `calculateUDOverlap()` quantifies pairwise home-range overlap; `plotMaps()`.

## Network analysis

* Two symmetric pipelines over a shared `mobyNetwork` object:
  * **Association networks** (individuals): `calculateAssociations()`, `randomizeAssociations()`,
    `plotAssociations()`, `plotAssociationMatrix()`.
  * **Movement networks** (locations): `calculateTransitions()`, `randomizeTransitions()`,
    `plotMovements()`, `transitionsTable()`.
  * `networkMetrics()` computes node- and network-level graph metrics (degree/strength,
    betweenness, eigenvector, community/modularity, reciprocity) for both, and `plot()` draws
    either network.
* Co-occurrence-based association indices (SRI/HWI) are interpreted as shared space-and-time use,
  not confirmed interaction; randomisation tests use multiple-comparison correction.

## Visualization

A consistent family of base-graphics plotting functions sharing a device-stable design: a single global
`cex`, colourblind-safe palettes (Okabe-Ito categorical, viridis continuous), an optional
`file`/`width`/`height`/`res` export, and an invisible tidy return.

* **Detection & effort over time:** `plotAbacus()` (detection abacus), `plotDeployments()` (receiver
  operating-period timeline), `plotActograms()`, `plotChronogram()`, `plotContours()`.
* **Per-station & distributions:** `plotStationStats()`, `plotGroupSizeDistribution()`,
  `plotMetricComparison()` (per-individual metrics across grouping levels, with design-appropriate
  repeated-measures tests and effect sizes).
* **Rhythms & periodicity:** `plotPeriodogram()` (FFT periodogram / opt-in Lomb-Scargle) and
  `plotScalogram()` (continuous wavelet scalogram).
* **Spatial & network:** `plotArray()` (receiver-array map), `plotMaps()` (home-range utilisation
  isopleths), `plotMovements()`, `plotAssociations()`, `plotAssociationMatrix()`; `shadeSeasons()`
  adds seasonal background shading.
