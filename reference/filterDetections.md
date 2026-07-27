# Filter and quality-control acoustic detections

Cleans a detection dataset by applying a sequence of quality-control
filters, each of which is opt-in (nothing but the temporal bounds is
applied unless you ask for it). Removed detections are never discarded
silently: they are returned in `data_discarded` with the exact `reason`
for removal, and the filtered detections carry a `qc_flag` column.

## Usage

``` r
filterDetections(
  data,
  id.col = NULL,
  datetime.col = NULL,
  lon.col = NULL,
  lat.col = NULL,
  land.shape = NULL,
  epsg.code = NULL,
  tagging.dates = NULL,
  cutoff.dates = NULL,
  remove.duplicates = TRUE,
  nominal.delay = NULL,
  min.lag.factor = 30,
  min.lag.threshold = NULL,
  isolation.window = NULL,
  max.speed = NULL,
  speed.unit = "m/s",
  acoustic.range = 600,
  min.corroboration = 2L,
  max.iterations = 20L,
  min.detections = 0,
  min.days = 0,
  verbose = getOption("moby.verbose", TRUE),
  ...
)
```

## Arguments

- data:

  A data frame (or `mobyData`) of raw animal detections.

- id.col:

  Name of the column containing animal IDs. Defaults to `"ID"`.

- datetime.col:

  Name of the column containing date-times in POSIXct format. Defaults
  to `"datetime"`.

- lon.col:

  Name of the column containing longitude (or projected x) values.
  Defaults to `"lon"`.

- lat.col:

  Name of the column containing latitude (or projected y) values.
  Defaults to `"lat"`.

- land.shape:

  Optional. Coastline/landmass polygons (`sf`, or `SpatialPolygons*`),
  enabling shortest in-water distances in the speed filter. Passed to
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md).

- epsg.code:

  Optional integer EPSG code of a projected (metre-based) CRS used to
  project positions for the speed filter.

- tagging.dates:

  Optional POSIXct vector of tagging/release dates. Either a single
  value (applied to all individuals) or a named vector whose names match
  the animal IDs.

- cutoff.dates:

  Optional. Cut-off date(s) beyond which detections are discarded (tag
  expiry, last download, etc.). A single POSIXct applied to all IDs, or
  a named POSIXct vector keyed by `id.col`.

- remove.duplicates:

  Logical; drop exact-duplicate records (same animal ID, timestamp and
  station) before any other filter. Defaults to TRUE.

- nominal.delay:

  Transmitter nominal (mean) delay, in seconds, used to scale the
  min_lag false-detection window. A single value applied to all
  individuals, or a named numeric vector keyed by `id.col` (for mixed
  tag families). Read from the `mobyData` metadata (`nominal.delay`)
  when not supplied. When `NULL` (and absent from metadata) the adaptive
  min_lag filter is off, but a fixed `min.lag.threshold` still enables
  it.

- min.lag.factor:

  Multiplier defining the min_lag threshold as
  `min.lag.factor * nominal.delay` (in seconds). Defaults to 30 (Pincock
  2012 rule of thumb). Use
  [`plotMinLag()`](https://miguelgandra.github.io/moby/reference/plotMinLag.md)
  to check whether this default suits a given dataset.

- min.lag.threshold:

  Optional. Set the min_lag threshold directly, in seconds, selecting
  the per-receiver *fixed*-window mode (see the ‘Isolation-based
  filtering’ section). Enables the filter on its own (no `nominal.delay`
  needed) and applies to every tag, including foreign IDs absent from
  `nominal.delay`; when both are supplied it overrides
  `min.lag.factor * nominal.delay`.

- isolation.window:

  Optional residency / sporadic-visitor filter (in hours): a detection
  is removed when the gap to BOTH temporal neighbours exceeds this
  window. `NULL` (default) or `FALSE` disables it. This is the
  whole-array counterpart of the per-receiver `min_lag` filter, and
  complementary to it (see the ‘Isolation-based filtering’ section); it
  is a residency / coverage check, NOT a false-detection filter (use
  `nominal.delay` for that).

- max.speed:

  Maximum plausible swim speed. A single value applied to all
  individuals, or a named numeric vector keyed by `id.col` (for
  multi-species / size-specific limits). `NULL` disables the speed
  filter.

- speed.unit:

  Units of `max.speed`: either "m/s" or "km/h".

- acoustic.range:

  Assumed detection range of the receivers (metres), used in the speed
  filter's minimum-distance correction. Note this makes the speed filter
  blind to consecutive detections closer than `2 * acoustic.range`.
  Defaults to 600.

- min.corroboration:

  Integer; in the speed filter, a flagged spatial outlier is
  auto-removed only if its cluster has fewer than this many detections.
  A corroborated group of `min.corroboration` or more is flagged for
  review instead of deleted (prevents cascading deletion of real
  relocations). Defaults to 2 (remove isolated singletons only).

- max.iterations:

  Maximum passes of the speed-filter convergence loop. Defaults to 20;
  use `Inf` for unlimited.

- min.detections:

  Discard individuals with fewer than this many surviving detections. 0
  = off.

- min.days:

  Discard individuals detected on fewer than this many distinct days. 0
  = off.

- verbose:

  Logical; print progress and the filtering summary. Defaults to
  `getOption("moby.verbose", TRUE)`.

- ...:

  Further arguments passed to
  [`calculateStepDistances`](https://miguelgandra.github.io/moby/reference/calculateStepDistances.md)
  (e.g. `grid.resolution`, `mov.directions`, `cores`).

## Value

An object of class `mobyFilter` (a list, with a `print` method)
containing:

- `data`:

  the filtered detections as a `mobyData`, with an added `qc_flag`
  column (`"valid"`, or `"overspeed_review"` for speed-flagged
  detections retained for review);

- `data_discarded`:

  the removed detections, each with a `reason` column;

- `summary`:

  a per-individual data frame of raw counts and removals by filter.

The filtering parameters used are stored as the `"parameters"`
attribute.

## Details

Filters are applied in this order (a filter is skipped when its
controlling argument is left at its "off" value):

1.  **Duplicate removal** (`remove.duplicates`): exact-duplicate records
    (same animal, timestamp and station) are dropped up front so they
    cannot slip past the other filters.

2.  **Pre-tagging**: detections before an animal's `tagging.date`.

3.  **Cut-off**: detections after an optional `cutoff.date`.

4.  **False-detection (min_lag)**: the standard short-interval
    false-detection filter (Pincock 2012; Simpfendorfer et al. 2015).
    For each detection it computes the shortest interval to the nearest
    other detection of the *same tag on the same receiver*; a detection
    whose min_lag exceeds `min.lag.factor * nominal.delay` (default 30x
    the transmitter nominal delay), or that has no same-receiver
    companion at all, is treated as a spurious (code-collision) decode.
    Enabled by `nominal.delay` (adaptive window) or `min.lag.threshold`
    (fixed window); off otherwise.

5.  **Isolation** (`isolation.window`): an optional residency /
    sporadic-visitor tool that removes a detection when the time gap to
    BOTH temporal neighbours (across the whole array) exceeds the
    window. This is NOT a false-detection filter and is off by default.

6.  **Speed**: flags steps exceeding `max.speed` using a minimum
    inter-detection distance
    (`straight/least-cost distance - 2 * acoustic.range`, a conservative
    lower bound). For each flagged step the local +/-3-detection window
    is hierarchically clustered; an isolated spatial outlier (a
    singleton cluster) is removed, whereas a corroborated group of
    `min.corroboration` or more detections is FLAGGED
    (`qc_flag = "overspeed_review"`) and retained, so a genuine
    relocation or a systematic error (bad coordinates, clock drift,
    under-estimated range) is surfaced for review rather than cascaded
    away. Requires `max.speed`; off otherwise.

7.  **Minimum detections / minimum days**: drop individuals with too few
    surviving detections or too few days with detections.

When a `land.shape` is supplied, distances (for the speed filter) are
shortest in-water (least-cost) distances rather than straight lines.
`nominal.delay` and `max.speed` may be a single value applied to all
individuals or a named numeric vector keyed by `id.col` (for mixed tag
families or multi-species datasets); `nominal.delay` is additionally
read from the `mobyData` metadata when not supplied.

## Isolation-based filtering

Two of the filters above, `min_lag` and `isolation.window`, are the same
underlying operation: a both-sides temporal-isolation test that flags a
detection whose nearest companion, on each side, lies further away than
a window. They differ on two axes.

**Spatial scale** (the primary distinction). `min_lag` works *per
receiver* – "is this detection corroborated on *this* receiver?" – the
scale at which spurious code-collision decodes actually occur.
`isolation.window` works *across the whole array* – "does the individual
have *any* other detection, anywhere, within the window?" – a residency
/ coverage question.

**Window definition**. `min_lag` uses a window that adapts to tag
transmission timing (`min.lag.factor * nominal.delay`, about 1 h for a
120 s tag), or a fixed value via `min.lag.threshold`. `isolation.window`
uses a fixed window in hours. A transmission-scaled window is only
meaningful per receiver; array-wide gaps are governed by animal
movement, so hours is the natural unit there.

|  |  |  |
|----|----|----|
|  | **per receiver** (min_lag) | **whole array** (isolation.window) |
| **adaptive window** | default (min.lag.factor) | – |
| **fixed window** | min.lag.threshold | isolation.window |

(The per-receiver filter is enabled by supplying either `nominal.delay`
– for the adaptive window – or `min.lag.threshold` – for the fixed
window, in seconds. If both are given, the fixed `min.lag.threshold`
takes precedence and applies to every tag.)

**They are complementary, not alternatives.** A typical workflow applies
`min_lag` first, as routine cleaning (removing spurious decodes), and
*then* – optionally – `isolation.window` as a residency check (flagging
animals seen only in sporadic one-off blips). Neither replaces the
other: use `min_lag` whenever a `nominal.delay` is available, and add
`isolation.window` when you also want to screen for biologically
unsupported presence.

*Coming from ATfiltR?* `min_lag` corresponds to
`findSolo(per.receiver = TRUE)` and `isolation.window` to
`findSolo(per.receiver = FALSE)`. moby deliberately uses a different
natural window at each scale (transmission-scaled per receiver, fixed
hours array-wide), so do not carry a single hours value across the two.

## References

Pincock, D.G. (2012) False detections: what they are and how to remove
them from detection data. VEMCO Application Note DOC-004691.

Simpfendorfer, C.A., Huveneers, C., Steckenreuter, A., Tattersall, K.,
Hoenner, X., Harcourt, R. & Heupel, M.R. (2015) Ghosts in the data:
false detections in VEMCO pulse position modulation acoustic telemetry
monitoring equipment. Animal Biotelemetry, 3:55.

Kessel, S.T., Cooke, S.J., Heupel, M.R., Hussey, N.E., Simpfendorfer,
C.A., Vagle, S. & Fisk, A.T. (2014) A review of detection range testing
in aquatic passive acoustic telemetry studies. Reviews in Fish Biology
and Fisheries, 24:199-218.

## See also

[`plotMinLag()`](https://miguelgandra.github.io/moby/reference/plotMinLag.md)
to check the min_lag threshold empirically.

## Examples

``` r
data(rays)
# default: only the temporal bounds are applied (tagging dates read from the mobyData metadata)
filtered <- filterDetections(rays)
#> Warning: - 'id.col' converted to factor.
#> - No 'nominal.delay' or 'min.lag.threshold' supplied or found in metadata: the min_lag false-detection filter is OFF. Supply the transmitter nominal delay (s) for an adaptive window, or a fixed 'min.lag.threshold' (s), to enable it.
#> ── filterDetections() ────────────────────────────────────────────────── moby ──
#> 
#> ℹ Removing spurious detections
#> • Input: 1,643 detections · 8 individuals
#> 
#> ✔ 1,643 detections retained (0 removed, 0%)
#> → Breakdown by filter: print() the result
filtered                     # prints the summary
#> <mobyFilter> filtered acoustic detections
#>   1643 retained | 0 removed
#>   Per-individual summary:
#>   ID Raw detections Before tagging Total removed
#>  D01            249              0             -
#>  D02            154              0             -
#>  D03            160              0             -
#>  D04            169              0             -
#>  R01            283              0             -
#>  R02            160              0             -
#>  R03            207              0             -
#>  R04            261              0             -
#>   Inspect $data, $data_discarded and attr(, "parameters").
head(filtered$data_discarded[, c("ID", "datetime", "reason")])
#> [1] ID       datetime reason  
#> <0 rows> (or 0-length row.names)

# enable the short-interval false-detection filter (needs the transmitter nominal delay)
filtered2 <- filterDetections(rays, nominal.delay = 120)   # 120 s tags
#> Warning: - 'id.col' converted to factor.
#> ── filterDetections() ────────────────────────────────────────────────── moby ──
#> 
#> ℹ Removing spurious detections
#> • Input: 1,643 detections · 8 individuals
#> 
#> → Filtering criteria
#>   • min_lag  3,600 s (30 × nominal delay)
#> 
#> ✔ 1,556 detections retained (87 removed, 5%)
#> → Breakdown by filter: print() the result

# or, without a known nominal delay, a fixed per-receiver isolation window (1 h = 3600 s)
filtered2b <- filterDetections(rays, min.lag.threshold = 3600)
#> Warning: - 'id.col' converted to factor.
#> ── filterDetections() ────────────────────────────────────────────────── moby ──
#> 
#> ℹ Removing spurious detections
#> • Input: 1,643 detections · 8 individuals
#> 
#> → Filtering criteria
#>   • min_lag  3,600 s (fixed window)
#> 
#> ✔ 1,556 detections retained (87 removed, 5%)
#> → Breakdown by filter: print() the result

# \donttest{
# add a movement-speed filter (slower: computes step distances); on a 3-animal subset
sub <- rays[rays$ID %in% head(levels(factor(rays$ID)), 3), ]
filtered3 <- filterDetections(sub, max.speed = 5, speed.unit = "km/h")
#> Warning: - 'id.col' converted to factor.
#> - No 'nominal.delay' or 'min.lag.threshold' supplied or found in metadata: the min_lag false-detection filter is OFF. Supply the transmitter nominal delay (s) for an adaptive window, or a fixed 'min.lag.threshold' (s), to enable it.
#> ── filterDetections() ────────────────────────────────────────────────── moby ──
#> 
#> ℹ Removing spurious detections
#> • Input: 563 detections · 3 individuals
#> 
#> → Filtering criteria
#>   • speed  5 km/h (great-circle)
#> 
#> ✔ 563 detections retained (0 removed, 0%)
#> ! 9 flagged for review (over-speed, retained)
#> → Breakdown by filter: print() the result
#> ⏱ runtime: 1.6s
# }
```
