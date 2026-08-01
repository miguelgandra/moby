# Canonical column names for the moby import functions

The import functions
[`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md),
[`importTags()`](https://miguelgandra.github.io/moby/reference/importTags.md)
and
[`importDeployments()`](https://miguelgandra.github.io/moby/reference/importDeployments.md)
read a raw export and rename its columns onto a small, fixed set of
*canonical* field names, so the rest of moby always sees the same layout
whatever the source system was. The `col.map` argument is how you
declare that mapping for a non-standard file: a named list of
`canonical_field = "your column name"`. For the built-in `source`
presets (`"vue"`, `"glatos"`, `"otn"`, `"etn"`, ...) the mapping is
filled in for you, so `col.map` only needs the fields a preset missed.

Column matching is *tolerant*: names are compared ignoring case, spaces
and punctuation, so `"Station Name"`, `"station.name"` and
`"STATION_NAME"` all match the canonical `station`. A field may list
several candidate names, and the first one present is used:
`list(lon = c("Longitude", "deploy_longitude"))`.

Only a handful of fields are **required**; every other field is optional
and is simply absent downstream if you do not supply it (functions that
need a missing field say so when you call them).

## Detection fields — [`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md)

|  |  |  |
|----|----|----|
| Field | Description | Required |
| `datetime` | Detection timestamp, parsed to POSIXct (see `datetime.format`) | **yes** |
| `transmitter` | Full tag/transmitter code — the key used to join tag metadata | **yes** (or `ID`) |
| `ID` | Animal identifier | no — set from `transmitter` if absent |
| `transmitter_codespace` | Code space; joined as `codespace-transmitter` when both are present | no |
| `transmitter_name` | Human-readable tag name | no |
| `receiver` | Receiver serial / ID | no |
| `station` | Receiver station name | no |
| `lon`, `lat` | Coordinates (coerced to numeric) | no |
| `sensor_value`, `sensor_unit` | Sensor reading and its unit | no |

If `ID` is absent it is initialised from `transmitter` (with a message);
assign true animal IDs afterwards with
[`matchTags()`](https://miguelgandra.github.io/moby/reference/matchTags.md).

## Tag fields — [`importTags()`](https://miguelgandra.github.io/moby/reference/importTags.md)

|  |  |  |
|----|----|----|
| Field | Description | Required |
| `transmitter` | Full tag code — the key joined to detections | **yes** |
| `ID` | Animal identifier | no — derived from `serial`/`transmitter` if absent |
| `serial` | Tag serial number | no |
| `tagging_date` | Release / tagging date (POSIXct) | no |
| `nominal_delay` | Transmitter mean delay, seconds (scales the min_lag filter) | no |
| `min_delay`, `max_delay` | Delay range; `nominal_delay` is set to their midpoint when it is absent | no |
| `species`, `sex`, `length`, `weight`, `tagging_location` | Biometrics / metadata | no |

## Deployment fields — [`importDeployments()`](https://miguelgandra.github.io/moby/reference/importDeployments.md)

|  |  |  |
|----|----|----|
| Field | Description | Required |
| `receiver` | Receiver serial / ID | **yes** |
| `station` | Station name | **yes** |
| `deploy` | Deployment date-time (POSIXct) | **yes** |
| `recover` | Recovery date-time (POSIXct); open (ongoing) deployments are handled downstream | no |
| `lon`, `lat` | Receiver coordinates (numeric) | no |
| `depth` | Station / bottom depth | no |

## What happens when an optional field is omitted

- No `ID` — initialised from `transmitter`; attach real animal IDs with
  [`matchTags()`](https://miguelgandra.github.io/moby/reference/matchTags.md).

- No `station` / `lon` / `lat`, but a `receiver` is present —
  [`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md)
  can back-fill them from the deployment log.

- No `receiver` — deployment matching
  ([`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md),
  [`checkDeployments()`](https://miguelgandra.github.io/moby/reference/checkDeployments.md))
  is unavailable.

- No coordinates at all — the import still succeeds, but spatial
  functions (COAs, maps, step distances) error when you later call them.

- No `nominal_delay` (and no `min_delay`/`max_delay`) — the min_lag
  false-detection filter in
  [`filterDetections()`](https://miguelgandra.github.io/moby/reference/filterDetections.md)
  stays off unless you set `min.lag.threshold` (see that function's
  ‘Isolation-based filtering’ section).

## Which function do I use?

The import functions and
[`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
do different jobs and are normally used **in sequence**, not as
alternatives:

|  |  |  |
|----|----|----|
|  | [`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md) etc. | [`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md) |
| **Job** | *Reshape* a raw export into moby's canonical layout | *Declare* a tidy table as a `mobyData` |
| **Does** | renames columns, parses date-times, derives `ID` | records which column plays which role; attaches metadata |
| **Column names** | replaced by canonical names | kept as they are (`id.col = "animal_id"` is fine) |
| **Carries** | nothing — returns a plain data frame | `id.groups`, `tagging.dates`, `nominal.delay`, CRS, land layer |
| **Use when** | data comes from a receiver system, or its columns/dates need fixing | your table is already tidy, or you are attaching study metadata |

Roughly: importing is to
[`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
what [`read.csv()`](https://rdrr.io/r/utils/read.table.html) is to
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) — a reader,
then a constructor. The usual pipeline runs both:

    det <- importDetections(file, source = "vue")      # 1. reshape -> data frame
    det <- matchTags(det, tags)                        # 2. join    -> data frame
    det <- as_moby(det, tags = tags,                   # 3. declare -> mobyData
                   id.groups = groups, epsg.code = 3395, land.shape = coast)

If you **already have a tidy data frame in R** (your own column names,
date-times already parsed), skip the importers and call
[`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
directly — it keeps your names and can carry the full metadata. Reach
for `source = "generic"` only when you also want the reshaping:
canonical renaming, date parsing and `ID` derivation.

## Worked `col.map` examples

    # 1. A non-standard VUE export where the transmitter is split across two columns
    importDetections(file, source = "vue",
                     col.map = list(transmitter_codespace = "Code.Space",
                                    transmitter = "ID", receiver = "Receiver.Name"))

    # 2. A GLATOS data frame already loaded in R — the preset handles it, no col.map needed
    importDetections(glatos_df, source = "glatos")

    # 3. A plain generic CSV with your own column names
    importDetections(file, source = "generic",
                     col.map = list(ID = "animal", datetime = "time_utc", station = "site",
                                    lon = "x", lat = "y", transmitter = "tag"))

    # 4. A tag database that stores a delay RANGE rather than a nominal delay
    importTags(file, source = "generic",
               col.map = list(transmitter = "Codespace", tagging_date = "Tagdeployed",
                              min_delay = "Minoff_sec", max_delay = "Maxoff_sec"))

## See also

[`importDetections()`](https://miguelgandra.github.io/moby/reference/importDetections.md),
[`importTags()`](https://miguelgandra.github.io/moby/reference/importTags.md),
[`importDeployments()`](https://miguelgandra.github.io/moby/reference/importDeployments.md),
[`matchTags()`](https://miguelgandra.github.io/moby/reference/matchTags.md),
[`matchDeployments()`](https://miguelgandra.github.io/moby/reference/matchDeployments.md),
[`as_moby()`](https://miguelgandra.github.io/moby/reference/as_moby.md)
