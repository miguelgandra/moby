#######################################################################################################
## Standardized importers / harmonisers ##############################################################
#######################################################################################################

# Source presets: for each canonical field, a vector of candidate source column names
# (first match present in the data wins). Matching is case-insensitive and ignores
# spaces/dots/underscores, so header variants are tolerated.

.detectionPresets <- function() {
  list(
    vue = list(
      datetime         = c("Date and Time (UTC)", "Date.and.Time..UTC.", "Date Time", "Date.Time", "datetime"),
      transmitter      = c("Transmitter"),
      transmitter_name = c("Transmitter Name"),
      receiver         = c("Receiver"),
      station          = c("Station Name", "Station"),
      lat              = c("Latitude", "Station Latitude"),
      lon              = c("Longitude", "Station Longitude"),
      sensor_value     = c("Sensor Value", "Sensor.1"),
      sensor_unit      = c("Sensor Unit", "Units.1")
    ),
    vdat = list(
      datetime         = c("Time", "Device Time (UTC)", "time", "devicetimeutc"),
      transmitter      = c("Full ID", "fullid", "ID", "id"),
      receiver         = c("Serial Number", "serialnumber", "Serial"),
      station          = c("Station Name", "stationname", "Station"),
      lat              = c("Latitude", "latitude"),
      lon              = c("Longitude", "longitude"),
      sensor_value     = c("Sensor Value", "Data", "rawdata"),
      sensor_unit      = c("Sensor Unit", "Units")
    ),
    glatos = list(
      ID               = c("animal_id"),
      datetime         = c("detection_timestamp_utc"),
      transmitter      = c("transmitter_id"),
      transmitter_codespace = c("transmitter_codespace"),
      transmitter_name = c("tag_serial_number"),
      receiver         = c("receiver_sn", "ins_serial_no"),
      station          = c("station", "glatos_array"),
      lat              = c("deploy_lat", "deploy_latitude"),
      lon              = c("deploy_long", "deploy_longitude"),
      sensor_value     = c("sensor_value"),
      sensor_unit      = c("sensor_unit")
    ),
    otn = list(
      ID               = c("catalognumber", "animal_id"),
      datetime         = c("datecollected", "detection_timestamp_utc"),
      transmitter      = c("tagname", "fieldnumber", "transmitter_id"),
      receiver         = c("receiver", "receiver_sn", "collectornumber"),
      station          = c("station", "station_name"),
      lat              = c("latitude", "deploy_lat"),
      lon              = c("longitude", "deploy_long"),
      sensor_value     = c("sensorvalue", "sensor_value"),
      sensor_unit      = c("sensorunit", "sensor_unit")
    ),
    etn = list(
      ID               = c("animal_id"),
      datetime         = c("date_time"),
      transmitter      = c("acoustic_tag_id"),
      transmitter_name = c("tag_serial_number"),
      receiver         = c("receiver_id"),
      station          = c("station_name"),
      lat              = c("deploy_latitude"),
      lon              = c("deploy_longitude"),
      sensor_value     = c("sensor_value"),
      sensor_unit      = c("sensor_unit")
    )
  )
}

.deploymentPresets <- function() {
  list(
    vue = list(
      receiver = c("Receiver"), station = c("Station", "Station Name"),
      lat = c("Latitude"), lon = c("Longitude"),
      deploy = c("Deploymentdate", "Deployment Date", "Deploy Date", "deploy_date_time"),
      recover = c("Dateout", "Date Out", "Recover Date", "recover_date_time"),
      depth = c("Stationdepth", "Bottomdepth", "Station Depth", "Bottom Depth")
    ),
    glatos = list(
      receiver = c("ins_serial_no", "receiver_sn"), station = c("station", "glatos_array"),
      lat = c("deploy_lat", "deploy_latitude"), lon = c("deploy_long", "deploy_longitude"),
      deploy = c("deploy_date_time"), recover = c("recover_date_time"),
      depth = c("bottom_depth")
    ),
    otn = list(
      receiver = c("receiver_sn", "receiver"), station = c("station_name", "station"),
      lat = c("deploy_lat", "latitude"), lon = c("deploy_long", "longitude"),
      deploy = c("deploy_date_time"), recover = c("recover_date_time"),
      depth = c("bottom_depth")
    ),
    etn = list(
      receiver = c("receiver_id"), station = c("station_name"),
      lat = c("deploy_latitude"), lon = c("deploy_longitude"),
      deploy = c("deploy_date_time"), recover = c("recover_date_time"),
      depth = c("bottom_depth", "deploy_depth")
    )
  )
}

# normalise a header for tolerant matching: lowercase, strip non-alphanumerics
.normHeader <- function(x) tolower(gsub("[^a-z0-9]", "", tolower(x)))

# find the first candidate column present in the data (tolerant match); returns name or NA
.matchColumn <- function(candidates, data_cols) {
  norm_data <- .normHeader(data_cols)
  for (cand in candidates) {
    hit <- which(norm_data == .normHeader(cand))
    if (length(hit) > 0) return(data_cols[hit[1]])
  }
  NA_character_
}

# robust POSIXct parsing (handles common acoustic-telemetry datetime layouts).
# Tries a sequence of explicit formats and keeps the one parsing the most values,
# avoiding lubridate's multi-order regex (which can fail on some PCRE2 builds).
.parseDatetime <- function(x, tz, format = NULL) {
  if (inherits(x, "POSIXct")) {
    attr(x, "tzone") <- tz
    return(x)
  }
  x <- as.character(x)
  if (!is.null(format)) return(as.POSIXct(x, format = format, tz = tz))
  formats <- c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M",
               "%Y/%m/%d %H:%M:%S", "%m/%d/%Y %H:%M:%S", "%d/%m/%Y %H:%M:%S", "%Y-%m-%d")
  non_na <- !is.na(x) & nzchar(x)
  n_target <- sum(non_na)
  best <- as.POSIXct(rep(NA_real_, length(x)), tz = tz)
  best_ok <- -1L
  for (fmt in formats) {
    parsed <- suppressWarnings(as.POSIXct(x, format = fmt, tz = tz))
    n_ok <- sum(!is.na(parsed) & non_na)
    if (n_ok > best_ok) { best_ok <- n_ok; best <- parsed }
    if (best_ok == n_target) break
  }
  best
}

# core harmoniser shared by importDetections / importDeployments
.harmonise <- function(data, mapping, datetime_fields, tz, datetime.format, keep.extra) {
  data_cols <- colnames(data)
  out <- data.frame(row.names = seq_len(nrow(data)))
  used <- character(0)
  mapped <- character(0)
  for (field in names(mapping)) {
    src <- .matchColumn(mapping[[field]], data_cols)
    if (!is.na(src)) {
      out[[field]] <- data[[src]]
      used <- c(used, src)
      mapped <- c(mapped, field)
    }
  }
  # parse datetimes
  for (field in intersect(datetime_fields, names(out))) {
    out[[field]] <- .parseDatetime(out[[field]], tz = tz, format = datetime.format)
  }
  # coerce coordinates / depth to numeric where present
  for (field in intersect(c("lon", "lat", "depth"), names(out))) {
    out[[field]] <- suppressWarnings(as.numeric(out[[field]]))
  }
  # retain unmapped columns if requested
  if (keep.extra) {
    extra <- setdiff(data_cols, used)
    for (e in extra) if (!e %in% names(out)) out[[e]] <- data[[e]]
  }
  # Record which canonical fields were genuinely resolved to a source column. The importers report
  # this to the user, and names(out) cannot answer it: under keep.extra an UNMAPPED source column is
  # copied through under its own name, so a canonical field whose col.map entry matched nothing would
  # still appear in names(out) whenever the source happens to carry a literal column of that name -
  # masking the very col.map typo the report exists to expose.
  attr(out, "mapped_fields") <- mapped
  out
}

# read a CSV or xlsx file path, or pass through a data.frame
.readSource <- function(x) {
  if (is.data.frame(x)) return(as.data.frame(x))
  if (!is.character(x) || length(x) != 1 || !file.exists(x)) {
    stop("'x' must be a data frame or a path to an existing .csv/.xlsx file.", call. = FALSE)
  }
  ext <- tolower(tools::file_ext(x))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Reading Excel files requires the 'readxl' package. Install it with install.packages('readxl'), or export to CSV.", call. = FALSE)
    }
    return(as.data.frame(readxl::read_excel(x)))
  }
  utils::read.csv(x, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}


# ---- console helpers shared by the three importers ------------------------------------------------
# One-line description of what was read and at what scale: the file's basename (never the full path,
# which would wrap the header) or "data frame", plus the incoming rows x columns.
.importSource <- function(x) {
  # basename() so a long path cannot wrap the header line
  if (is.data.frame(x)) "data frame" else basename(as.character(x)[1])
}

# The methodological choices of an import: which mapping was applied, and how timestamps were read.
.importCriteria <- function(source, col.map, tz, datetime.format) {
  crit <- c(source = if (source == "generic") "generic (user col.map)"
                     else paste0(source, " preset", if (!is.null(col.map)) " + user col.map" else ""),
            timezone = tz)
  # only worth reporting when the user overrode the auto-detection
  if (!is.null(datetime.format)) crit["datetime format"] <- datetime.format
  crit
}

# Report which canonical fields the mapping resolved. An unresolved OPTIONAL field is normal, so this
# is a note, not a warning: only the user can say whether a missing field matters to them.
.reportFieldMapping <- function(resolved, unresolved, verbose) {
  if (length(resolved) > 0)
    .mobyNote("Fields resolved: ", paste(resolved, collapse = ", "), verbose = verbose)
  if (length(unresolved) > 0)
    .mobyNote("Fields not found: ", paste(unresolved, collapse = ", "), verbose = verbose)
  invisible(NULL)
}


#' Canonical column names for the moby import functions
#'
#' @description The import functions [importDetections()], [importTags()] and [importDeployments()]
#' read a raw export and rename its columns onto a small, fixed set of *canonical* field names, so the
#' rest of moby always sees the same layout whatever the source system was. The `col.map` argument is
#' how you declare that mapping for a non-standard file: a named list of
#' `canonical_field = "your column name"`. For the built-in `source` presets (`"vue"`, `"glatos"`,
#' `"otn"`, `"etn"`, ...) the mapping is filled in for you, so `col.map` only needs the fields a preset
#' missed.
#'
#' Column matching is *tolerant*: names are compared ignoring case, spaces and punctuation, so
#' `"Station Name"`, `"station.name"` and `"STATION_NAME"` all match the canonical `station`. A field
#' may list several candidate names, and the first one present is used:
#' `list(lon = c("Longitude", "deploy_longitude"))`.
#'
#' Only a handful of fields are **required**; every other field is optional and is simply absent
#' downstream if you do not supply it (functions that need a missing field say so when you call them).
#'
#' @section Detection fields — `importDetections()`:
#' | Field | Description | Required |
#' |---|---|---|
#' | `datetime` | Detection timestamp, parsed to POSIXct (see `datetime.format`) | **yes** |
#' | `transmitter` | Full tag/transmitter code — the key used to join tag metadata | **yes** (or `ID`) |
#' | `ID` | Animal identifier | no — set from `transmitter` if absent |
#' | `transmitter_codespace` | Code space; joined as `codespace-transmitter` when both are present | no |
#' | `transmitter_name` | Human-readable tag name | no |
#' | `receiver` | Receiver serial / ID | no |
#' | `station` | Receiver station name | no |
#' | `lon`, `lat` | Coordinates (coerced to numeric) | no |
#' | `sensor_value`, `sensor_unit` | Sensor reading and its unit | no |
#'
#' If `ID` is absent it is initialised from `transmitter` (with a message); assign true animal IDs
#' afterwards with [assignAnimalIDs()].
#'
#' @section Tag fields — `importTags()`:
#' | Field | Description | Required |
#' |---|---|---|
#' | `transmitter` | Full tag code — the key joined to detections | **yes** |
#' | `ID` | Animal identifier | no — derived from `serial`/`transmitter` if absent |
#' | `serial` | Tag serial number | no |
#' | `tagging_date` | Release / tagging date (POSIXct) | no |
#' | `nominal_delay` | Transmitter mean delay, seconds (scales the min_lag filter) | no |
#' | `min_delay`, `max_delay` | Delay range; `nominal_delay` is set to their midpoint when it is absent | no |
#' | `species`, `sex`, `length`, `weight`, `tagging_location` | Biometrics / metadata | no |
#'
#' @section Deployment fields — `importDeployments()`:
#' | Field | Description | Required |
#' |---|---|---|
#' | `receiver` | Receiver serial / ID | **yes** |
#' | `station` | Station name | **yes** |
#' | `deploy` | Deployment date-time (POSIXct) | **yes** |
#' | `recover` | Recovery date-time (POSIXct); open (ongoing) deployments are handled downstream | no |
#' | `lon`, `lat` | Receiver coordinates (numeric) | no |
#' | `depth` | Station / bottom depth | no |
#'
#' @section What happens when an optional field is omitted:
#' \itemize{
#'   \item No `ID` — initialised from `transmitter`; attach real animal IDs with [assignAnimalIDs()].
#'   \item No `station` / `lon` / `lat`, but a `receiver` is present — [matchDeployments()] can
#'     back-fill them from the deployment log.
#'   \item No `receiver` — deployment matching ([matchDeployments()], [checkDeployments()]) is
#'     unavailable.
#'   \item No coordinates at all — the import still succeeds, but spatial functions (COAs, maps, step
#'     distances) error when you later call them.
#'   \item No `nominal_delay` (and no `min_delay`/`max_delay`) — the min_lag false-detection filter in
#'     [filterDetections()] stays off unless you set `min.lag.threshold` (see that function's
#'     \sQuote{Isolation-based filtering} section).
#' }
#'
#' @section Which function do I use?:
#' The import functions and [as_moby()] do different jobs and are normally used **in sequence**, not
#' as alternatives:
#'
#' | | [importDetections()] etc. | [as_moby()] |
#' |---|---|---|
#' | **Job** | *Reshape* a raw export into moby's canonical layout | *Declare* a tidy table as a `mobyData` |
#' | **Does** | renames columns, parses date-times, derives `ID` | records which column plays which role; attaches metadata |
#' | **Column names** | replaced by canonical names | kept as they are (`id.col = "animal_id"` is fine) |
#' | **Carries** | nothing — returns a plain data frame | `id.groups`, `tagging.dates`, `nominal.delay`, CRS, land layer |
#' | **Use when** | data comes from a receiver system, or its columns/dates need fixing | your table is already tidy, or you are attaching study metadata |
#'
#' Roughly: importing is to [as_moby()] what `read.csv()` is to `data.frame()` — a reader, then a
#' constructor. The usual pipeline runs both:
#'
#' ```r
#' det <- importDetections(file, source = "vue")      # 1. reshape  -> data frame
#' det <- assignAnimalIDs(det, tags)                  # 2. join tags -> mobyData (IDs, delays, dates)
#' det <- as_moby(det, id.groups = groups,            # 3. annotate  -> add the rest
#'                epsg.code = 3395, land.shape = coast)
#' ```
#'
#' If you **already have a tidy data frame in R** (your own column names, date-times already parsed),
#' skip the importers and call [as_moby()] directly — it keeps your names and can carry the full
#' metadata. Reach for `source = "generic"` only when you also want the reshaping: canonical renaming,
#' date parsing and `ID` derivation.
#'
#' @section Worked `col.map` examples:
#' ```r
#' # 1. A non-standard VUE export where the transmitter is split across two columns
#' importDetections(file, source = "vue",
#'                  col.map = list(transmitter_codespace = "Code.Space",
#'                                 transmitter = "ID", receiver = "Receiver.Name"))
#'
#' # 2. A GLATOS data frame already loaded in R — the preset handles it, no col.map needed
#' importDetections(glatos_df, source = "glatos")
#'
#' # 3. A plain generic CSV with your own column names
#' importDetections(file, source = "generic",
#'                  col.map = list(ID = "animal", datetime = "time_utc", station = "site",
#'                                 lon = "x", lat = "y", transmitter = "tag"))
#'
#' # 4. A tag database that stores a delay RANGE rather than a nominal delay
#' importTags(file, source = "generic",
#'            col.map = list(transmitter = "Codespace", tagging_date = "Tagdeployed",
#'                           min_delay = "Minoff_sec", max_delay = "Maxoff_sec"))
#' ```
#'
#' @seealso [importDetections()], [importTags()], [importDeployments()], [assignAnimalIDs()],
#' [matchDeployments()], [as_moby()]
#' @name moby_import_schema
NULL


#' Import and harmonise acoustic detection data
#'
#' @description Reads acoustic-telemetry detections from a range of common sources and
#' harmonises them into a single consistent schema - moby's canonical column names, with
#' date-times parsed - ready to be turned into a \code{\link{mobyData}} by
#' \code{\link{assignAnimalIDs}} / \code{\link{as_moby}}. Supported sources are Innovasea/VEMCO VUE
#' exports, Innovasea VDAT/Fathom `DET.csv` files, and detection extracts from the GLATOS,
#' OTN and ETN (`etn` package) systems. A `generic` mode plus a user-supplied `col.map`
#' handle non-standard layouts.
#'
#' @param x A path to a `.csv` (or `.xlsx`) detection file, or a data frame already loaded in
#' R (e.g. the output of `etn::get_acoustic_detections()` or `glatos::read_glatos_detections()`).
#' @param source One of `"vue"`, `"vdat"`, `"glatos"`, `"otn"`, `"etn"` or `"generic"`.
#' For `"generic"`, supply `col.map`.
#' @param tz Time zone used to parse date-times. Defaults to `"UTC"` (the convention for
#' GLATOS/OTN/ETN); set explicitly for VUE/VDAT exports recorded in another zone.
#' @param col.map Optional named list mapping canonical fields to the column name(s) in `x`, merged
#' over (and overriding) the chosen `source` preset. See [moby_import_schema] for the full list of
#' detection fields, which are required, and worked examples.
#' @param datetime.format Optional explicit `strptime` format for the datetime column; if
#' `NULL`, common layouts are auto-detected.
#' @param keep.extra Logical; retain source columns that were not mapped to a canonical field.
#' Defaults to `FALSE`.
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A data frame with harmonised columns (`ID`, `datetime`, `transmitter`, `receiver`,
#' `station`, `lon`, `lat`, ...), sorted by animal and time. When the source has no animal
#' identifier, `ID` is initialised from `transmitter` (assign true animal IDs afterwards with
#' \code{\link{assignAnimalIDs}}).
#'
#' This is a *harmonised table*, not yet a \code{\link{mobyData}}: importing reshapes columns, it
#' does not attach study metadata. Turn it into a `mobyData` with \code{\link{assignAnimalIDs}}
#' (animal IDs, tagging dates, nominal delays) and/or \code{\link{as_moby}} (`id.groups`, CRS, land
#' layer) - see the \sQuote{Which function do I use?} section of \code{\link{moby_import_schema}}.
#' Because the output already uses moby's canonical column names, those functions need no explicit
#' column arguments.
#'
#' @seealso \code{\link{moby_import_schema}} for the canonical field list and how importing relates
#' to \code{\link{as_moby}}; \code{\link{assignAnimalIDs}}, \code{\link{importDeployments}},
#' \code{\link{checkDeployments}}, \code{\link{as_moby}}
#' @examples
#' # harmonise a generic-layout detection CSV via an explicit column map
#' csv <- system.file("extdata", "rays_detections.csv", package = "moby")
#' det <- importDetections(csv, source = "generic",
#'                         col.map = list(ID = "animal_id", datetime = "timestamp",
#'                                        station = "station_name", lon = "deploy_longitude",
#'                                        lat = "deploy_latitude", receiver = "receiver_id",
#'                                        transmitter = "transmitter"))
#' head(det)
#'
#' @export

importDetections <- function(x,
                             source = c("vue", "vdat", "glatos", "otn", "etn", "generic"),
                             tz = "UTC",
                             col.map = NULL,
                             datetime.format = NULL,
                             keep.extra = FALSE,
                             verbose = getOption("moby.verbose", TRUE)) {

  source <- match.arg(source)

  # mapping assembled first: it needs no data, so an unusable 'col.map' or 'source' still errors
  # before anything is printed
  mapping <- if (source == "generic") list() else .detectionPresets()[[source]]
  if (!is.null(col.map)) {
    if (!is.list(col.map) || is.null(names(col.map))) stop("'col.map' must be a named list.", call. = FALSE)
    for (nm in names(col.map)) mapping[[nm]] <- col.map[[nm]]
  }
  if (length(mapping) == 0) stop("No column mapping available. For source='generic', supply 'col.map'.", call. = FALSE)

  # ---- header ---------------------------------------------------------------------------------
  # Printed BEFORE the file is read and parsed: on a large export those steps dominate the runtime,
  # and the user should see what is running rather than a silent console. The trade-off is that a
  # preset which turns out not to match the file now errors beneath the banner rather than in front
  # of it - the cheap argument checks above still error cleanly.
  .mobyHeader("importDetections()", "Reading and harmonising acoustic detections",
              input = .importSource(x), criteria = .importCriteria(source, col.map, tz, datetime.format),
              verbose = verbose)

  data <- .readSource(x)

  out <- .harmonise(data, mapping, datetime_fields = "datetime", tz = tz,
                    datetime.format = datetime.format, keep.extra = keep.extra)

  if (!"datetime" %in% names(out)) stop("Could not locate a datetime column. Check 'source' or provide 'col.map'.", call. = FALSE)

  # which canonical fields the mapping actually resolved, recorded by .harmonise() itself (see there:
  # names(out) cannot answer this under keep.extra)
  resolved <- attr(out, "mapped_fields")
  attr(out, "mapped_fields") <- NULL

  # combine GLATOS codespace + id into a full transmitter string when both present
  if (all(c("transmitter_codespace", "transmitter") %in% names(out))) {
    out$transmitter <- paste(out$transmitter_codespace, out$transmitter, sep = "-")
    out$transmitter_codespace <- NULL
    # merged away, so it is no longer a field the user can find in the result
    resolved <- setdiff(resolved, "transmitter_codespace")
  }

  # initialise animal ID from transmitter when the source carries none
  if (!"ID" %in% names(out)) {
    if (!"transmitter" %in% names(out)) stop("Neither an animal ID nor a transmitter column could be located.", call. = FALSE)
    out$ID <- as.character(out$transmitter)
    .mobyBlank(verbose)
    .mobyNote("No animal-ID column found; 'ID' initialised from 'transmitter'. ",
              "Join tag metadata to assign true animal IDs.", verbose = verbose)
  }
  out$ID <- factor(out$ID)

  # order canonical columns first
  pref <- c("ID", "datetime", "transmitter", "transmitter_name", "receiver", "station",
            "lon", "lat", "sensor_value", "sensor_unit")
  ord <- c(intersect(pref, names(out)), setdiff(names(out), pref))
  out <- out[, ord, drop = FALSE]
  out <- out[order(out$ID, out$datetime), , drop = FALSE]
  rownames(out) <- NULL

  # ---- outcome --------------------------------------------------------------------------------
  # The field-resolution lists are the high-value part: they are the only way for a user to tell
  # whether the chosen preset / col.map actually found their columns.
  .mobyBlank(verbose)
  # NA is not a transmitter: counting it would overstate the array
  n_tx <- if ("transmitter" %in% names(out)) length(unique(stats::na.omit(as.character(out$transmitter)))) else NULL
  .mobyOk(.fmtCount(nrow(out), "detection"), " imported",
          if (!is.null(n_tx)) paste0(" from ", .fmtCount(n_tx, "transmitter")) else "",
          verbose = verbose)
  .reportFieldMapping(resolved, setdiff(names(mapping), resolved), verbose)

  out
}


#' Import and harmonise receiver deployment metadata
#'
#' @description Reads a receiver-deployment / station log from a range of common sources and
#' harmonises it into a consistent schema for use with \code{\link{checkDeployments}} and the rest
#' of the `moby` workflow.
#'
#' @param x A path to a `.csv`/`.xlsx` deployment log, or a data frame (e.g. the output of
#' `etn::get_acoustic_deployments()`).
#' @param source One of `"vue"`, `"glatos"`, `"otn"`, `"etn"` or `"generic"`.
#' @param tz Time zone used to parse deploy/recover date-times. Defaults to `"UTC"`.
#' @param col.map Optional named list mapping canonical deployment fields to source column name(s),
#' merged over the `source` preset. See [moby_import_schema] for the full list of deployment fields
#' and which are required.
#' @param datetime.format Optional explicit `strptime` format for the deploy/recover columns.
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A data frame with columns `receiver`, `station`, `lon`, `lat`, `deploy` (POSIXct),
#' `recover` (POSIXct) and, where available, `depth`; sorted by receiver and deployment date.
#'
#' @seealso \code{\link{moby_import_schema}} for the canonical field list;
#' \code{\link{importDetections}}, \code{\link{checkDeployments}}
#' @examples
#' # read a raw ETN deployment export and harmonise it
#' csv <- system.file("extdata", "rays_deployments.csv", package = "moby")
#' deployments <- importDeployments(csv, source = "etn")
#' head(deployments)
#'
#' @export

importDeployments <- function(x,
                              source = c("vue", "glatos", "otn", "etn", "generic"),
                              tz = "UTC",
                              col.map = NULL,
                              datetime.format = NULL,
                              verbose = getOption("moby.verbose", TRUE)) {

  source <- match.arg(source)

  # mapping assembled first: it needs no data, so an unusable 'col.map' or 'source' still errors
  # before anything is printed
  mapping <- if (source == "generic") list() else .deploymentPresets()[[source]]
  if (!is.null(col.map)) {
    if (!is.list(col.map) || is.null(names(col.map))) stop("'col.map' must be a named list.", call. = FALSE)
    for (nm in names(col.map)) mapping[[nm]] <- col.map[[nm]]
  }
  if (length(mapping) == 0) stop("No column mapping available. For source='generic', supply 'col.map'.", call. = FALSE)

  # ---- header ---------------------------------------------------------------------------------
  # Printed BEFORE the file is read and parsed: on a large export those steps dominate the runtime,
  # and the user should see what is running rather than a silent console. The trade-off is that a
  # preset which turns out not to match the file now errors beneath the banner rather than in front
  # of it - the cheap argument checks above still error cleanly.
  .mobyHeader("importDeployments()", "Reading and harmonising the receiver deployment log",
              input = .importSource(x), criteria = .importCriteria(source, col.map, tz, datetime.format),
              verbose = verbose)

  data <- .readSource(x)

  out <- .harmonise(data, mapping, datetime_fields = c("deploy", "recover"), tz = tz,
                    datetime.format = datetime.format, keep.extra = FALSE)

  for (req in c("receiver", "station", "deploy")) {
    if (!req %in% names(out)) stop(paste0("Could not locate a '", req, "' column. Check 'source' or provide 'col.map'."), call. = FALSE)
  }

  # resolution record, kept by .harmonise() itself (names(out) cannot answer it under keep.extra)
  resolved <- attr(out, "mapped_fields")
  attr(out, "mapped_fields") <- NULL
  out$receiver <- as.character(out$receiver)
  out$station <- as.character(out$station)

  pref <- c("receiver", "station", "lon", "lat", "deploy", "recover", "depth")
  ord <- c(intersect(pref, names(out)), setdiff(names(out), pref))
  out <- out[, ord, drop = FALSE]
  out <- out[order(out$receiver, out$deploy), , drop = FALSE]
  rownames(out) <- NULL

  # ---- outcome --------------------------------------------------------------------------------
  .mobyBlank(verbose)
  .mobyOk(.fmtCount(nrow(out), "deployment record"), " imported across ",
          .fmtCount(length(unique(stats::na.omit(out$receiver))), "receiver"), verbose = verbose)
  .reportFieldMapping(resolved, setdiff(names(mapping), resolved), verbose)

  out
}


.tagPresets <- function() {
  list(
    vue = list(
      transmitter = c("Transmitter", "Tag ID", "Full ID", "Codespace"),
      ID = c("ID", "Animal", "Animal ID"),
      serial = c("Serialno", "Serial", "Serial Number"),
      tagging_date = c("Tagdeployed", "Tagging Date", "Tag Date", "tagging_date"),
      tagging_location = c("Location", "Tagging Location"),
      species = c("Common", "Species", "Scientific Name"),
      sex = c("Sex"), length = c("Tl_cm", "Length", "TL"),
      nominal_delay = c("Nominal Delay", "Nominal Delay (s)", "Delay", "Delay (s)", "nominal_delay"),
      min_delay = c("Min Delay", "Min Delay (s)", "min_delay"),
      max_delay = c("Max Delay", "Max Delay (s)", "max_delay")
    ),
    glatos = list(
      ID = c("animal_id"), transmitter = c("transmitter_id"),
      transmitter_codespace = c("transmitter_codespace"),
      serial = c("tag_serial_number"),
      tagging_date = c("utc_release_date_time", "release_date_time"),
      tagging_location = c("release_location"),
      species = c("common_name_e", "scientific_name"), sex = c("sex"),
      length = c("length"), weight = c("weight"),
      nominal_delay = c("tag_nominal_delay", "nominal_delay"),
      min_delay = c("tag_min_delay", "min_delay"), max_delay = c("tag_max_delay", "max_delay")
    ),
    otn = list(
      ID = c("animal_id", "catalognumber"), transmitter = c("tagname", "fieldnumber"),
      serial = c("tag_serial_number"),
      tagging_date = c("utc_release_date_time", "datereleasedtagger", "release_date"),
      species = c("scientificname", "commonname"), sex = c("sex"), length = c("length"),
      nominal_delay = c("nominal_delay", "tag_nominal_delay"),
      min_delay = c("min_delay"), max_delay = c("max_delay")
    ),
    etn = list(
      ID = c("animal_id"), transmitter = c("acoustic_tag_id"),
      serial = c("tag_serial_number"),
      tagging_date = c("release_date_time"), tagging_location = c("release_location"),
      species = c("scientific_name"), sex = c("sex"), length = c("length1", "length"),
      nominal_delay = c("nominal_delay"), min_delay = c("min_delay"), max_delay = c("max_delay")
    )
  )
}


#' Import and harmonise tag / animal metadata
#'
#' @description Reads a tag/animal metadata table and harmonises it into a consistent schema
#' (`transmitter`, `ID`, `tagging_date`, plus biometrics such as `species`, `sex`, `length`, and the
#' transmitter `nominal_delay` when the source provides one). Use together with
#' \code{\link{assignAnimalIDs}} to attach animal IDs, tagging dates and nominal delays to a
#' detection dataset.
#'
#' @details When the source specifies a delay RANGE (`min_delay`/`max_delay`) rather than a nominal
#' delay - as many tag-specification exports do - `nominal_delay` is derived as their midpoint. The
#' nominal delay is what \code{\link{filterDetections}} uses to scale its short-interval (min_lag)
#' false-detection filter, so carrying it through here means that filter can be enabled automatically.
#'
#' @param x A path to a `.csv`/`.xlsx` tag-metadata file, or a data frame (e.g. from
#' `etn::get_tags()` / `etn::get_animals()`).
#' @param source One of `"vue"`, `"glatos"`, `"otn"`, `"etn"` or `"generic"`.
#' @param tz Time zone used to parse the tagging date. Defaults to `"UTC"`.
#' @param col.map Optional named list mapping canonical tag fields to the column name(s) in `x`,
#' merged over (and extending) the `source` preset. See [moby_import_schema] for the full list of tag
#' fields and examples.
#' @param datetime.format Optional explicit `strptime` format for the tagging-date column.
#' @param keep.extra Logical; retain unmapped source columns. Defaults to `TRUE` so that
#' additional biometric fields are preserved.
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A data frame with at least `transmitter` and (when available) `ID`, `tagging_date`
#' (POSIXct), `nominal_delay` (seconds) and biometric columns.
#'
#' @seealso \code{\link{moby_import_schema}} for the canonical field list;
#' \code{\link{assignAnimalIDs}}, \code{\link{importDetections}}
#' @examples
#' # harmonise a tag-metadata table (here the bundled 'rays_tags' data frame)
#' tags <- importTags(rays_tags, source = "generic",
#'                    col.map = list(ID = "ID", transmitter = "transmitter",
#'                                   tagging_date = "tagging_date", species = "species"))
#' head(tags)
#'
#' @export

importTags <- function(x,
                       source = c("vue", "glatos", "otn", "etn", "generic"),
                       tz = "UTC",
                       col.map = NULL,
                       datetime.format = NULL,
                       keep.extra = TRUE,
                       verbose = getOption("moby.verbose", TRUE)) {

  source <- match.arg(source)

  # mapping assembled first: it needs no data, so an unusable 'col.map' or 'source' still errors
  # before anything is printed
  mapping <- if (source == "generic") list() else .tagPresets()[[source]]
  if (!is.null(col.map)) {
    if (!is.list(col.map) || is.null(names(col.map))) stop("'col.map' must be a named list.", call. = FALSE)
    for (nm in names(col.map)) mapping[[nm]] <- col.map[[nm]]
  }
  if (length(mapping) == 0) stop("No column mapping available. For source='generic', supply 'col.map'.", call. = FALSE)

  # ---- header ---------------------------------------------------------------------------------
  # Printed BEFORE the file is read and parsed: on a large export those steps dominate the runtime,
  # and the user should see what is running rather than a silent console. The trade-off is that a
  # preset which turns out not to match the file now errors beneath the banner rather than in front
  # of it - the cheap argument checks above still error cleanly.
  .mobyHeader("importTags()", "Reading and harmonising tag and animal metadata",
              input = .importSource(x), criteria = .importCriteria(source, col.map, tz, datetime.format),
              verbose = verbose)

  data <- .readSource(x)

  out <- .harmonise(data, mapping, datetime_fields = "tagging_date", tz = tz,
                    datetime.format = datetime.format, keep.extra = keep.extra)

  # resolution record, kept by .harmonise() itself (names(out) cannot answer it under keep.extra,
  # which is this function's default). Adjusted below as derived columns are added/merged.
  resolved <- attr(out, "mapped_fields")
  attr(out, "mapped_fields") <- NULL

  if (all(c("transmitter_codespace", "transmitter") %in% names(out))) {
    out$transmitter <- paste(out$transmitter_codespace, out$transmitter, sep = "-")
    out$transmitter_codespace <- NULL
    resolved <- setdiff(resolved, "transmitter_codespace")
  }
  if (!"transmitter" %in% names(out)) {
    stop("Could not locate a 'transmitter' column (the key used to join tags to detections). Provide 'col.map'.", call. = FALSE)
  }

  out$transmitter <- as.character(out$transmitter)
  if ("length" %in% names(out)) out$length <- suppressWarnings(as.numeric(out$length))
  for (dc in intersect(c("nominal_delay", "min_delay", "max_delay"), names(out)))
    out[[dc]] <- suppressWarnings(as.numeric(out[[dc]]))

  # many tag exports specify a delay RANGE rather than a nominal delay; the nominal (mean) delay is
  # the midpoint, which is what the short-interval false-detection filter is scaled to
  if (!"nominal_delay" %in% names(out) && all(c("min_delay", "max_delay") %in% names(out))) {
    out$nominal_delay <- (out$min_delay + out$max_delay) / 2
    # the result DOES carry a nominal_delay now, so it must not be listed as "not found" below - this
    # is the one field a user has to confirm is present (filterDetections reads it to scale min_lag)
    resolved <- union(resolved, "nominal_delay")
    .mobyBlank(verbose)
    .mobyNote("'nominal_delay' derived as the midpoint of 'min_delay' and 'max_delay'.", verbose = verbose)
  }

  pref <- c("ID", "transmitter", "serial", "tagging_date", "tagging_location", "species", "sex",
            "length", "nominal_delay", "min_delay", "max_delay")
  ord <- c(intersect(pref, names(out)), setdiff(names(out), pref))
  out <- out[, ord, drop = FALSE]

  # ---- outcome --------------------------------------------------------------------------------
  # The transmitter count is reported only when it differs from the row count, i.e. when several rows
  # share a code (a tag listed under more than one code space, or a re-tagged animal).
  .mobyBlank(verbose)
  n_tx <- length(unique(stats::na.omit(out$transmitter)))
  .mobyOk(.fmtCount(nrow(out), "tag"), " imported",
          if (n_tx != nrow(out)) paste0(" ", .mobyGlyph("mid"), " ",
                                        .fmtCount(n_tx, "distinct transmitter")) else "",
          verbose = verbose)
  .reportFieldMapping(resolved, setdiff(names(mapping), resolved), verbose)

  out
}


# Tolerant transmitter matching: exact string first, then the trailing numeric code.
#
# The numeric fallback exists for tag tables that store the bare code ("111") rather than the full
# string ("A69-1602-111"). It deliberately does NOT fire when both sides carry a code space and those
# code spaces DISAGREE: "A69-9999-111" and "A69-1602-111" are different transmitters, and matching
# them on the shared suffix would silently attribute another project's animal to one of yours.
# A tag that legitimately transmits on several code spaces (e.g. a Vemco V16P sensor tag) should list
# each of its codes as a row of the tag table, all mapped to the same animal ID; exact matching then
# resolves them, and assignAnimalIDs() reports any code that needs adding. Numeric codes shared by
# more than one tag are likewise left unmatched rather than resolved arbitrarily.
.matchTransmitter <- function(det_tx, tag_tx) {
  idx <- match(det_tx, tag_tx)
  clash <- character(0)
  na <- is.na(idx)
  if (any(na)) {
    num   <- function(x) sub(".*[-_ ]([0-9]+)$", "\\1", x)
    space <- function(x) ifelse(grepl("[-_ ][0-9]+$", x), sub("[-_ ][0-9]+$", "", x), NA_character_)
    det_num <- num(det_tx);   tag_num <- num(tag_tx)
    det_cs  <- space(det_tx); tag_cs  <- space(tag_tx)
    # only numeric codes unique within the tag table are eligible for the fallback
    ambiguous <- tag_num %in% tag_num[duplicated(tag_num)]
    cand <- match(det_num, replace(tag_num, ambiguous, NA_character_))
    # reject a candidate whose code space is known on both sides and differs
    bad <- !is.na(cand) & !is.na(det_cs) & !is.na(tag_cs[cand]) & det_cs != tag_cs[cand]
    clash <- unique(det_tx[na & bad])
    cand[bad] <- NA_integer_
    idx[na] <- cand[na]
  }
  attr(idx, "codespace_clash") <- clash
  idx
}


#' Assign animal IDs (and tagging dates) to detections from tag metadata
#'
#' @description Resolves the transmitter codes in a detection dataset to the animals that carried
#' them, using a tag table (see \code{\link{importTags}}). This is more than a column join: several
#' transmitter codes can map to the same animal, and the function also derives the per-animal
#' metadata - tagging dates and transmitter nominal delays - that later steps such as
#' \code{\link{filterDetections}} and \code{\link{summaryTable}} read automatically. Optional
#' biometric columns can be joined in as well.
#'
#' @details
#' \strong{Why several codes can mean one animal.} Most often this is because a single physical
#' transmitter emits on more than one code space: sensor tags (for example depth-sensing Vemco V16P
#' transmitters) report their identity and their sensor data under different code spaces, so one tag
#' appears in the detection file under two or more transmitter codes. Animals genuinely carrying
#' several transmitters at once are uncommon; the other routine case is an animal that was recaptured
#' and re-tagged, so its history spans consecutive tags. List each code as its own row of `tags`,
#' mapped to the same animal `ID`, and all of them resolve to that animal.
#'
#' \strong{Matching.} The full transmitter string is matched first. A trailing-numeric fallback then
#' covers tag tables that store the bare code (so `"A69-1602-111"` matches a tag stored as `"111"`).
#' The fallback is deliberately conservative: a detection whose code space differs from the candidate
#' tag's is left unmatched rather than assigned, because `"A69-9999-111"` and `"A69-1602-111"` are
#' different transmitters and guessing would attribute another project's animal to yours. Such codes
#' are reported, so you can add them to `tags` if they are further code spaces of your own tags.
#'
#' \strong{Where it sits in the pipeline.} This step is independent of
#' \code{\link{matchDeployments}}: one resolves transmitters to animals, the other resolves receivers
#' to deployment windows, and they key on different columns. Either order gives the same result.
#'
#' @param detections A detection dataset (`mobyData` or data frame) with a transmitter column.
#' @param tags A harmonised tag table from \code{\link{importTags}} (or a data frame with at
#' least a `transmitter` column).
#' @param transmitter.col Name of the transmitter column in `detections`. Defaults to
#' `"transmitter"`.
#' @param id.col Name of the animal-ID column to (re)create in `detections`. Resolved from the
#' `mobyData` metadata or `"ID"` when `NULL`.
#' @param keep.cols Optional character vector of additional `tags` columns (e.g. `"sex"`,
#' `"length"`, `"species"`) to join into the detections.
#' @param set.tagging.dates Logical; if `TRUE` (default) and `tags` has a `tagging_date` column,
#' attach per-individual tagging dates to the returned object's metadata.
#' @param set.nominal.delay Logical; if `TRUE` (default) and `tags` has a `nominal_delay` column
#' (see \code{\link{importTags}}), attach per-individual transmitter nominal delays (seconds) to the
#' returned object's metadata. \code{\link{filterDetections}} reads these automatically to scale its
#' short-interval (min_lag) false-detection filter, so arrays mixing tag families (e.g. 60 s and
#' 120 s tags) are handled per animal.
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A \code{\link{mobyData}} object with the `ID` column assigned (and, optionally,
#' tagging dates and biometric columns attached). Detections whose transmitter is absent from
#' `tags` keep `NA` IDs (with a warning).
#'
#' @seealso \code{\link{importTags}}, \code{\link{importDetections}}, \code{\link{as_moby}}
#' @examples
#' # join tag metadata to detections to assign animal IDs (and tagging dates)
#' tags <- importTags(rays_tags, source = "generic",
#'                    col.map = list(ID = "ID", transmitter = "transmitter",
#'                                   tagging_date = "tagging_date"))
#' # detections carrying the tagged rays' transmitters
#' det <- rays_detections[rays_detections$transmitter %in% rays_tags$transmitter, ]
#' det <- assignAnimalIDs(det, tags)
#' levels(det$ID)
#'
#' @export

assignAnimalIDs <- function(detections,
                            tags,
                            id.col = NULL,
                            transmitter.col = "transmitter",
                            keep.cols = NULL,
                            set.tagging.dates = TRUE,
                            set.nominal.delay = TRUE,
                            verbose = getOption("moby.verbose", TRUE)) {

  id.col <- .resolveArgs(detections, list(id.col = id.col))$id.col
  prev_meta <- attr(detections, "moby")
  det <- as.data.frame(detections)
  tg <- as.data.frame(tags)

  if (!transmitter.col %in% colnames(det)) {
    stop(paste0("Transmitter column ('", transmitter.col, "') not found in 'detections'."), call. = FALSE)
  }
  if (!"transmitter" %in% colnames(tg)) stop("'tags' must contain a 'transmitter' column (see importTags()).", call. = FALSE)

  # ---- header ---------------------------------------------------------------------------------
  # No criteria block: the join key is a column mapping, not a methodological choice, and the
  # framework reserves "-> Method" for settings that change what the result MEANS.
  .mobyHeader("assignAnimalIDs()", "Joining tag metadata to replace placeholder IDs with real animal IDs",
              input = paste0(.fmtCount(nrow(det), "detection"), " ", .mobyGlyph("mid"), " ",
                             .fmtCount(length(unique(stats::na.omit(as.character(det[[transmitter.col]])))),
                                       "transmitter")),
              verbose = verbose)

  # ensure the tag table has an animal-ID column
  if (!"ID" %in% colnames(tg)) {
    use_serial <- "serial" %in% colnames(tg) && !all(is.na(tg$serial))
    key <- if (use_serial) tg$serial else tg$transmitter
    tg$ID <- as.character(key)
    .mobyBlank(verbose)
    # report the column actually used: an all-NA 'serial' is present but unusable, and falls back
    .mobyNote("'tags' has no 'ID' column; animal IDs derived from ",
              if (use_serial) "'serial'" else "'transmitter'", ".", verbose = verbose)
  }
  tg$ID <- as.character(tg$ID)

  idx <- .matchTransmitter(as.character(det[[transmitter.col]]), as.character(tg$transmitter))
  clash <- attr(idx, "codespace_clash")
  n_unmatched <- sum(is.na(idx))
  if (n_unmatched > 0) {
    warning(paste0("- ", n_unmatched, " detection(s) had a transmitter not found in 'tags'; their ID is NA."), call. = FALSE)
  }
  # codes that share a tag's numeric part but sit in a different code space are left unmatched on
  # purpose: they are either another project's tags, or a further code space of one of your own tags
  # (common with sensor transmitters). Only the user can tell which, so name them.
  if (length(clash) > 0) {
    .mobyBlank(verbose)
    .mobyNote(.fmtCount(length(clash), "transmitter code"),
              if (length(clash) == 1) " matches" else " match", " a tag's number but a different code ",
              "space, so they were NOT assigned: ", paste(utils::head(clash, 5), collapse = ", "),
              if (length(clash) > 5) ", ..." else "",
              ". If these are further code spaces of your own transmitters (e.g. a sensor tag ",
              "transmitting on more than one code space), add each code as a row of 'tags' mapped to ",
              "the same animal ID.", verbose = verbose)
  }

  det[[id.col]] <- factor(tg$ID[idx])

  # optionally join additional biometric columns
  if (!is.null(keep.cols)) {
    miss <- setdiff(keep.cols, colnames(tg))
    if (length(miss) > 0) warning(paste0("- keep.cols not found in 'tags' and skipped: ", paste(miss, collapse = ", ")), call. = FALSE)
    for (cc in intersect(keep.cols, colnames(tg))) det[[cc]] <- tg[[cc]][idx]
  }

  # build per-individual tagging dates for the object metadata
  tagging.dates <- NULL
  if (set.tagging.dates && "tagging_date" %in% colnames(tg) && inherits(tg$tagging_date, "POSIXct")) {
    valid <- !is.na(tg$ID) & !is.na(tg$tagging_date)
    if (any(valid)) {
      epoch <- tapply(as.numeric(tg$tagging_date[valid]), tg$ID[valid], min)
      tagging.dates <- as.POSIXct(epoch, origin = "1970-01-01", tz = .dataTZ(tg$tagging_date))
      names(tagging.dates) <- names(epoch)
      # keep only IDs present in the detections
      tagging.dates <- tagging.dates[names(tagging.dates) %in% as.character(unique(det[[id.col]]))]
      if (length(tagging.dates) == 0) tagging.dates <- NULL
    }
  }

  # build per-individual nominal delays for the object metadata (scales filterDetections' min_lag
  # false-detection filter); the median guards against duplicate tag rows for one animal
  nominal.delay <- NULL
  if (set.nominal.delay && "nominal_delay" %in% colnames(tg)) {
    nd <- suppressWarnings(as.numeric(tg$nominal_delay))
    valid <- !is.na(tg$ID) & !is.na(nd) & nd > 0
    if (any(valid)) {
      agg <- tapply(nd[valid], tg$ID[valid], stats::median)
      nominal.delay <- stats::setNames(as.numeric(agg), names(agg))
      # keep only IDs present in the detections
      nominal.delay <- nominal.delay[names(nominal.delay) %in% as.character(unique(det[[id.col]]))]
      if (length(nominal.delay) == 0) nominal.delay <- NULL
    }
  }

  det <- det[order(det[[id.col]], det[[if ("datetime" %in% colnames(det)) "datetime" else id.col]]), , drop = FALSE]
  rownames(det) <- NULL

  # ---- outcome --------------------------------------------------------------------------------
  # Two of this function's effects are otherwise invisible: it writes tagging.dates and nominal.delay
  # into the returned object's metadata, which is what filterDetections() later reads. Report them, so
  # "did my tag table supply the delay the min_lag filter needs?" is answered here rather than by
  # inspecting attributes.
  .mobyBlank(verbose)
  # count detections that actually RECEIVED an ID, not merely those that found a tag row: a tag whose
  # own ID is NA (a blank 'serial', or a blank ID column) matches but assigns nothing, and reporting
  # the match rate would then claim every detection was assigned while half came back NA
  n_assigned <- sum(!is.na(det[[id.col]]))
  n_animals <- nlevels(droplevels(det[[id.col]]))
  .mobyOk(.fmtCount(n_assigned, "detection"), " matched to ", .fmtCount(n_animals, "animal"),
          verbose = verbose)
  meta_set <- character(0)
  if (!is.null(tagging.dates))
    meta_set <- c(meta_set, paste0("tagging.dates (", .fmtCount(length(tagging.dates), "individual"), ")"))
  if (!is.null(nominal.delay)) {
    # a single shared delay prints as one value; a mixed-tag-family array prints as a range. The
    # coverage is part of the answer: a delay attached for 2 of 8 animals scales the min_lag filter
    # for 2 of 8, and would otherwise be indistinguishable from full coverage.
    nd_rng <- range(nominal.delay)
    nd_txt <- if (isTRUE(all.equal(nd_rng[1], nd_rng[2]))) format(nd_rng[1])
              else paste0(format(nd_rng[1]), "-", format(nd_rng[2]))
    cover <- if (length(nominal.delay) < n_animals)
               paste0(", ", .fmtN(length(nominal.delay)), " of ", .fmtCount(n_animals, "individual")) else ""
    meta_set <- c(meta_set, paste0("nominal.delay (", nd_txt, " s", cover, ")"))
  }
  if (length(meta_set) > 0)
    .mobyNote("Metadata attached: ", paste(meta_set, collapse = ", "), verbose = verbose)
  # a requested metadata field that the tag table appeared to carry but yielded nothing is otherwise
  # discoverable only by inspecting attributes
  if (set.nominal.delay && "nominal_delay" %in% colnames(tg) && is.null(nominal.delay))
    .mobyNote("No usable 'nominal_delay' in 'tags'; none attached (filterDetections' min_lag filter ",
              "will need it supplied).", verbose = verbose)
  # NB: detections left unassigned are NOT reported here - the warning() above already states it, and
  # saying the same thing twice in two different styles is precisely the inconsistency to avoid.
  # Tag rows that matched but carry no ID have no such warning, so they get one of their own.
  n_blank <- sum(!is.na(idx) & is.na(tg$ID[idx]))
  if (n_blank > 0)
    warning(paste0("- ", n_blank, " detection(s) matched a tag row with no animal ID; their ID is NA."),
            call. = FALSE)

  # rebuild the mobyData, preserving the original metadata (column map, CRS, etc.) and updating the
  # ID column and (when available) the tagging dates and transmitter nominal delays
  base_meta <- if (!is.null(prev_meta)) prev_meta else list()
  base_meta$id.col <- id.col
  if (!is.null(tagging.dates)) base_meta$tagging.dates <- tagging.dates
  if (!is.null(nominal.delay)) base_meta$nominal.delay <- nominal.delay
  do.call(as_moby, c(list(det), base_meta, list(verbose = FALSE)))
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
