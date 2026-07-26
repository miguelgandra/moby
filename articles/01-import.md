# 01 - From raw detections to a mobyData object

> **In this module**
>
> Read raw detections from a file, harmonise them into moby’s expected
> layout, attach tag metadata, and wrap everything into a **`mobyData`**
> object that downstream functions understand. You will also define the
> **`id.groups`** object that drives per-species analyses throughout the
> tutorials.
>
> **Prerequisites:** none. Uses the bundled CSV at
> `system.file("extdata", "rays_detections.csv", package = "moby")`.
> (For pulling data from a live database, see the optional [ETN
> appendix](https://miguelgandra.github.io/moby/articles/00-etn-import.md).)

## Learning objectives

1.  Read a generic detection export and map its columns with
    [`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md).
2.  Attach tag metadata with
    [`importTags()`](https://miguelgandra.github.io/moby/reference/importTags.md)
    /
    [`assignAnimalIDs()`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md).
3.  Understand the `mobyData` data model and create one with
    [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md).
4.  Define `id.groups` (two species) once, for reuse everywhere.

## Setup

``` r

library(moby)
detections_csv <- system.file("extdata", "rays_detections.csv", package = "moby")
```

## 1. Read and harmonise the detections

The raw file uses source-specific column names (`animal_id`,
`timestamp`, `station_name`, …). The generic importer maps them onto
moby’s canonical names — and parses the date-times for you, so there is
no need to convert them beforehand.

``` r

detections <- importDetections(
  detections_csv,
  source = "generic",
  col.map = list(
    ID          = "animal_id",
    datetime    = "timestamp",
    station     = "station_name",
    lon         = "deploy_longitude",
    lat         = "deploy_latitude",
    transmitter = "transmitter"      # kept so assignAnimalIDs() can join tag metadata
  )
)
```

> **Tip**
>
> moby also ships dedicated importers for common sources —
> `importDetections(source = "vue")`, `"vdat"`, `"glatos"`, `"otn"`,
> `"etn"` — so you usually do **not** need to hand-map columns. See
> [`?moby_import_schema`](https://miguelgandra.github.io/moby/reference/moby_import_schema.md)
> for every canonical field and which ones are required.

> **Importing reshapes;
> [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
> declares**
>
> [`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md)
> returns a **plain data frame** in moby’s canonical layout — it renames
> columns, parses dates and derives `ID`, but attaches no study
> metadata. Steps 2 and 3 below turn it into a `mobyData`. If your table
> is *already* tidy (your own column names, dates already parsed), skip
> the importer and go straight to
> [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md).

## 2. Attach tag / animal metadata

``` r

tags <- importTags(rays_tags, source = "generic",
                   col.map = list(ID = "ID", tagging_date = "tagging_date"))

detections <- assignAnimalIDs(detections, tags)
```

[`assignAnimalIDs()`](https://miguelgandra.github.io/moby/reference/assignAnimalIDs.md)
resolves transmitter codes to animals — several codes can belong to one
animal, typically because a sensor tag transmits on more than one code
space — and it **derives per-animal metadata** from the tag table:
tagging dates, and transmitter nominal delays where present. Those are
attached to the returned object, so you do not pass them again later:
[`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
picks them up on its own.

> **Tip**
>
> This step is independent of
> [`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md)
> (module 02): one resolves transmitters to animals, the other resolves
> receivers to deployment windows. Either order gives the same result.

## 3. Build the `mobyData` object

A `mobyData` *is* a data frame, but it carries metadata (column mapping,
tagging dates, `id.groups`) as an attribute, so you don’t repeat those
arguments on every call.

The object returned above is already a `mobyData` with its column roles
and tagging dates attached. Call
[`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
to add what only you can supply — here the species grouping (and, when
you have them, a projected CRS and a coastline). Because the columns
already use moby’s canonical names, no column arguments are needed.

``` r

# one named list element per species -> drives per-species outputs everywhere
id_groups <- split(rays_tags$ID, rays_tags$species)

dataset <- as_moby(detections, id.groups = id_groups)

mobyMeta(dataset)   # tagging dates are already there, from assignAnimalIDs()
dataset             # the print method summarises what is attached
```

## Recap & what’s next

You produced a `mobyData` object and defined `id.groups`. The bundled
`rays` object is exactly this, already cleaned — later modules load it
directly.

➡️ **Next:** [02 — Quality control, deployment matching &
filtering](https://miguelgandra.github.io/moby/articles/02-quality-control.md).
