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
  for (field in names(mapping)) {
    src <- .matchColumn(mapping[[field]], data_cols)
    if (!is.na(src)) {
      out[[field]] <- data[[src]]
      used <- c(used, src)
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


#' Import and harmonise acoustic detection data
#'
#' @description Reads acoustic-telemetry detections from a range of common sources and
#' harmonises them into a single consistent schema, returning a \code{\link{mobyData}} object
#' ready for the rest of the `moby` workflow. Supported sources are Innovasea/VEMCO VUE
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
#' @param col.map Optional named list mapping canonical fields (`datetime`, `transmitter`,
#' `receiver`, `station`, `lon`, `lat`, `ID`, `sensor_value`, `sensor_unit`, ...) to the column
#' name(s) in `x`. Merged over (and overriding) the chosen `source` preset.
#' @param datetime.format Optional explicit `strptime` format for the datetime column; if
#' `NULL`, common layouts are auto-detected.
#' @param keep.extra Logical; retain source columns that were not mapped to a canonical field.
#' Defaults to `FALSE`.
#'
#' @return A \code{\link{mobyData}} object with harmonised columns
#' (`ID`, `datetime`, `transmitter`, `receiver`, `station`, `lon`, `lat`, ...). When the source
#' has no animal identifier, `ID` is initialised from `transmitter` (assign true animal IDs
#' later by joining tag metadata).
#'
#' @seealso \code{\link{importDeployments}}, \code{\link{checkDeployments}}, \code{\link{as_moby}}
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
                             keep.extra = FALSE) {

  source <- match.arg(source)
  data <- .readSource(x)

  mapping <- if (source == "generic") list() else .detectionPresets()[[source]]
  if (!is.null(col.map)) {
    if (!is.list(col.map) || is.null(names(col.map))) stop("'col.map' must be a named list.", call. = FALSE)
    for (nm in names(col.map)) mapping[[nm]] <- col.map[[nm]]
  }
  if (length(mapping) == 0) stop("No column mapping available. For source='generic', supply 'col.map'.", call. = FALSE)

  out <- .harmonise(data, mapping, datetime_fields = "datetime", tz = tz,
                    datetime.format = datetime.format, keep.extra = keep.extra)

  if (!"datetime" %in% names(out)) stop("Could not locate a datetime column. Check 'source' or provide 'col.map'.", call. = FALSE)

  # combine GLATOS codespace + id into a full transmitter string when both present
  if (all(c("transmitter_codespace", "transmitter") %in% names(out))) {
    out$transmitter <- paste(out$transmitter_codespace, out$transmitter, sep = "-")
    out$transmitter_codespace <- NULL
  }

  # initialise animal ID from transmitter when the source carries none
  if (!"ID" %in% names(out)) {
    if (!"transmitter" %in% names(out)) stop("Neither an animal ID nor a transmitter column could be located.", call. = FALSE)
    out$ID <- as.character(out$transmitter)
    message("Note: no animal-ID column found; 'ID' initialised from 'transmitter'. ",
            "Join tag metadata to assign true animal IDs.")
  }
  out$ID <- factor(out$ID)

  # order canonical columns first
  pref <- c("ID", "datetime", "transmitter", "transmitter_name", "receiver", "station",
            "lon", "lat", "sensor_value", "sensor_unit")
  ord <- c(intersect(pref, names(out)), setdiff(names(out), pref))
  out <- out[, ord, drop = FALSE]
  out <- out[order(out$ID, out$datetime), , drop = FALSE]
  rownames(out) <- NULL

  as_moby(out, id.col = "ID", datetime.col = "datetime", station.col = "station",
          lon.col = "lon", lat.col = "lat")
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
#' @param col.map Optional named list mapping canonical fields (`receiver`, `station`, `lon`,
#' `lat`, `deploy`, `recover`, `depth`) to source column name(s); merged over the `source` preset.
#' @param datetime.format Optional explicit `strptime` format for the deploy/recover columns.
#'
#' @return A data frame with columns `receiver`, `station`, `lon`, `lat`, `deploy` (POSIXct),
#' `recover` (POSIXct) and, where available, `depth`; sorted by receiver and deployment date.
#'
#' @seealso \code{\link{importDetections}}, \code{\link{checkDeployments}}
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
                              datetime.format = NULL) {

  source <- match.arg(source)
  data <- .readSource(x)

  mapping <- if (source == "generic") list() else .deploymentPresets()[[source]]
  if (!is.null(col.map)) {
    if (!is.list(col.map) || is.null(names(col.map))) stop("'col.map' must be a named list.", call. = FALSE)
    for (nm in names(col.map)) mapping[[nm]] <- col.map[[nm]]
  }
  if (length(mapping) == 0) stop("No column mapping available. For source='generic', supply 'col.map'.", call. = FALSE)

  out <- .harmonise(data, mapping, datetime_fields = c("deploy", "recover"), tz = tz,
                    datetime.format = datetime.format, keep.extra = FALSE)

  for (req in c("receiver", "station", "deploy")) {
    if (!req %in% names(out)) stop(paste0("Could not locate a '", req, "' column. Check 'source' or provide 'col.map'."), call. = FALSE)
  }
  out$receiver <- as.character(out$receiver)
  out$station <- as.character(out$station)

  pref <- c("receiver", "station", "lon", "lat", "deploy", "recover", "depth")
  ord <- c(intersect(pref, names(out)), setdiff(names(out), pref))
  out <- out[, ord, drop = FALSE]
  out <- out[order(out$receiver, out$deploy), , drop = FALSE]
  rownames(out) <- NULL
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
#' @param col.map Optional named list overriding/extending the `source` preset (e.g.
#' `list(transmitter = "Tag", tagging_date = "Deployed")`).
#' @param datetime.format Optional explicit `strptime` format for the tagging-date column.
#' @param keep.extra Logical; retain unmapped source columns. Defaults to `TRUE` so that
#' additional biometric fields are preserved.
#'
#' @return A data frame with at least `transmitter` and (when available) `ID`, `tagging_date`
#' (POSIXct), `nominal_delay` (seconds) and biometric columns.
#'
#' @seealso \code{\link{assignAnimalIDs}}, \code{\link{importDetections}}
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
                       keep.extra = TRUE) {

  source <- match.arg(source)
  data <- .readSource(x)

  mapping <- if (source == "generic") list() else .tagPresets()[[source]]
  if (!is.null(col.map)) {
    if (!is.list(col.map) || is.null(names(col.map))) stop("'col.map' must be a named list.", call. = FALSE)
    for (nm in names(col.map)) mapping[[nm]] <- col.map[[nm]]
  }
  if (length(mapping) == 0) stop("No column mapping available. For source='generic', supply 'col.map'.", call. = FALSE)

  out <- .harmonise(data, mapping, datetime_fields = "tagging_date", tz = tz,
                    datetime.format = datetime.format, keep.extra = keep.extra)

  if (all(c("transmitter_codespace", "transmitter") %in% names(out))) {
    out$transmitter <- paste(out$transmitter_codespace, out$transmitter, sep = "-")
    out$transmitter_codespace <- NULL
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
    message("Note: 'nominal_delay' derived as the midpoint of 'min_delay' and 'max_delay'.")
  }

  pref <- c("ID", "transmitter", "serial", "tagging_date", "tagging_location", "species", "sex",
            "length", "nominal_delay", "min_delay", "max_delay")
  ord <- c(intersect(pref, names(out)), setdiff(names(out), pref))
  out[, ord, drop = FALSE]
}


# tolerant transmitter matching: exact first, then by trailing numeric code
.matchTransmitter <- function(det_tx, tag_tx) {
  idx <- match(det_tx, tag_tx)
  na <- is.na(idx)
  if (any(na)) {
    det_num <- sub(".*[-_ ]([0-9]+)$", "\\1", det_tx)
    tag_num <- sub(".*[-_ ]([0-9]+)$", "\\1", tag_tx)
    idx2 <- match(det_num, tag_num)
    idx[na] <- idx2[na]
  }
  idx
}


#' Assign animal IDs (and tagging dates) to detections from tag metadata
#'
#' @description Joins tag metadata (see \code{\link{importTags}}) to a detection dataset on the
#' transmitter code, assigning each detection an animal `ID`. Matching is tolerant: it first
#' tries the full transmitter string and then falls back to the trailing numeric code (so
#' `"A69-1602-111"` matches a tag stored as `"111"`). When the tag table carries tagging dates,
#' these are attached to the returned \code{\link{mobyData}} object's metadata (and can flow
#' automatically into functions such as \code{\link{filterDetections}} and
#' \code{\link{summaryTable}}). Optional biometric columns can be joined in as well.
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
                            set.nominal.delay = TRUE) {

  id.col <- .resolveArgs(detections, list(id.col = id.col))$id.col
  prev_meta <- attr(detections, "moby")
  det <- as.data.frame(detections)
  tg <- as.data.frame(tags)

  if (!transmitter.col %in% colnames(det)) {
    stop(paste0("Transmitter column ('", transmitter.col, "') not found in 'detections'."), call. = FALSE)
  }
  if (!"transmitter" %in% colnames(tg)) stop("'tags' must contain a 'transmitter' column (see importTags()).", call. = FALSE)

  # ensure the tag table has an animal-ID column
  if (!"ID" %in% colnames(tg)) {
    key <- if ("serial" %in% colnames(tg) && !all(is.na(tg$serial))) tg$serial else tg$transmitter
    tg$ID <- as.character(key)
    message("Note: 'tags' has no 'ID' column; animal IDs derived from ",
            if ("serial" %in% colnames(tg)) "'serial'" else "'transmitter'", ".")
  }
  tg$ID <- as.character(tg$ID)

  idx <- .matchTransmitter(as.character(det[[transmitter.col]]), as.character(tg$transmitter))
  n_unmatched <- sum(is.na(idx))
  if (n_unmatched > 0) {
    warning(paste0("- ", n_unmatched, " detection(s) had a transmitter not found in 'tags'; their ID is NA."), call. = FALSE)
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

  # rebuild the mobyData, preserving the original metadata (column map, CRS, etc.) and updating the
  # ID column and (when available) the tagging dates and transmitter nominal delays
  base_meta <- if (!is.null(prev_meta)) prev_meta else list()
  base_meta$id.col <- id.col
  if (!is.null(tagging.dates)) base_meta$tagging.dates <- tagging.dates
  if (!is.null(nominal.delay)) base_meta$nominal.delay <- nominal.delay
  do.call(as_moby, c(list(det), base_meta))
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
