# The moby pipeline: from a raw export to an analysis-ready dataset

`moby`'s preparation functions look numerous, but they fill only **five
roles**. Learn the roles and the running order stops being something to
memorise. Two things are true throughout, and both simplify the picture:

- Every function also accepts an **ordinary data frame**. A
  [`mobyData`](https://miguelgandra.github.io/moby/reference/as_moby.md)
  is a convenience that saves you repeating arguments - never a gate you
  must pass through.

- Only the *enrich* and *clean* roles change your data. The audit role
  cannot.

## The five roles

- **1. Read** -
  [`importDetections`](https://miguelgandra.github.io/moby/reference/importDetections.md),
  [`importTags`](https://miguelgandra.github.io/moby/reference/importTags.md),
  [`importDeployments`](https://miguelgandra.github.io/moby/reference/importDeployments.md):

  Turn a vendor export into moby's canonical columns: renames columns,
  parses date-times, coerces coordinates. Each returns a plain data
  frame. Optional - skip them if your table already uses
  `ID`/`datetime`/`station`/`lon`/`lat` with a POSIXct datetime. See
  [`moby_import_schema`](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
  for every field.

- **2. Audit** -
  [`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md):

  Validates the receiver log (coverage gaps, overlapping deployments,
  duplicate stations, invalid dates, implausible or on-land coordinates,
  short monitoring effort) and cross-checks detections against it.
  Returns a report and **modifies nothing**. Recommended.

- **3. Enrich** -
  [`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md),
  [`matchDeployments`](https://miguelgandra.github.io/moby/reference/matchDeployments.md):

  Fold a side table into the detections.
  [`assignAnimalIDs()`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
  resolves transmitter codes to animals and attaches per-animal tagging
  dates and nominal delays;
  [`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md)
  matches detections to deployment windows and back-fills station names
  and coordinates. They key on different columns (transmitter vs
  receiver) and are **order-independent**.

- **4. Clean & bin** -
  [`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md),
  [`correctPositions`](https://miguelgandra.github.io/moby/reference/correctPositions.md),
  [`getTimeBins`](https://miguelgandra.github.io/moby/reference/getTimeBins.md):

  Where the data actually becomes analysis-ready: removing spurious
  detections, relocating positions that fall on land, and assigning
  regular time bins.

- **5. Declare** -
  [`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md):

  Records which column plays which role and attaches study metadata
  (`id.groups`, `epsg.code`, `land.shape`). Not really a stage - it is
  idempotent, so call it whenever you want to add metadata without
  losing what is already attached.

## A typical run

    # 1 - read
    det  <- importDetections(det_file, source = "vue")
    tags <- importTags(tag_file, source = "vue")

    # 2 - enrich: animal identity, plus tagging dates and nominal delays
    det <- assignAnimalIDs(det, tags)
    det <- det[!is.na(det$ID), ]          # drop transmitters absent from 'tags'
    det$ID <- droplevels(det$ID)

    # 3 - clean: note this returns a REPORT, not a data frame
    qc  <- filterDetections(det, max.speed = 2, speed.unit = "m/s")
    qc                                    # what was removed, and why
    det <- qc$data

    # 4 - bin, and declare anything only you know
    det$timebin <- getTimeBins(det$datetime, interval = "30 mins")
    dataset <- as_moby(det, id.groups = groups, epsg.code = 32629)

Add the receiver-log branch when you need it - audit, then match:

    deployments <- importDeployments(log_file, source = "vue")
    checkDeployments(deployments, detections = det)     # report only; fix what it flags
    det <- matchDeployments(det, deployments)

Already hold a tidy data frame? No moby object is required at any point:

    qc <- filterDetections(df, tagging.dates = tagging_dates, max.speed = 2)

## What happens if you skip a step

|  |  |  |
|----|----|----|
| **function** | **status** | **skip it and...** |
| [`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md) | optional | nothing breaks; supply canonical columns yourself |
| [`importTags()`](https://miguelgandra.github.io/moby/reference/importTags.md) | optional | nothing breaks; you must supply tagging dates another way |
| [`importDeployments()`](https://miguelgandra.github.io/moby/reference/importDeployments.md) | optional | nothing breaks; only needed if you use the receiver log |
| [`checkDeployments()`](https://miguelgandra.github.io/moby/reference/checkDeployments.md) | recommended | nothing breaks; you may match against a log with real errors |
| [`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md) | optional | coordinates stay as exported, unvalidated and ungapfilled |
| [`assignAnimalIDs()`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md) | **essential in practice** | [`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md) **errors** without tagging dates; the min_lag filter is unavailable |
| [`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md) | **essential** | false detections, pre-tagging records and impossible speeds stay in |
| [`correctPositions()`](https://miguelgandra.github.io/moby/reference/correctPositions.md) | optional | nothing breaks unless coordinates sit on land |
| [`getTimeBins()`](https://miguelgandra.github.io/moby/reference/getTimeBins.md) | recommended | no per-bin analyses (COAs, wide tables, chronograms) |
| [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md) | recommended | no `id.groups`, projected CRS or land-aware distances |

## Ordering rules (there are only two)

1.  **Audit before you delete.**
    [`checkDeployments`](https://miguelgandra.github.io/moby/reference/checkDeployments.md)
    is read-only and idempotent, so it has no single correct slot - run
    it early to find problems, and again after corrections to confirm
    they landed. The one position that matters: run it *before*
    `matchDeployments(drop.unmatched = TRUE)`, which removes the records
    that four of its five detection-side checks rely on. Under the
    default (`FALSE`) nothing is lost and order is free.

2.  **Enrich before you filter.**
    [`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
    needs tagging dates, which
    [`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
    supplies. The two enrich steps themselves may run in either order.

## Three easy traps

- [`filterDetections`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
  and
  [`correctPositions`](https://miguelgandra.github.io/moby/reference/correctPositions.md)
  return **lists**, not data frames - take `$data`. Writing
  `det <- filterDetections(det)` leaves you holding a `mobyFilter`
  object and the next call fails.

- [`getTimeBins`](https://miguelgandra.github.io/moby/reference/getTimeBins.md)
  returns a **vector**; assign it to a column yourself.

- [`assignAnimalIDs`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
  leaves `NA` IDs (with a warning) for transmitters that are not in your
  tag table - typically other projects' tags. Remove them, or they
  travel into the analysis.

## See also

[`moby_import_schema`](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
for the canonical column names;
[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md)
for the data model; the package vignettes for worked examples.
