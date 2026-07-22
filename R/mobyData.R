#######################################################################################################
## mobyData object model ##############################################################################
#######################################################################################################

# Canonical fallback column names, used when a column argument is left NULL and no
# mobyData metadata is available (i.e. when a plain data.frame is supplied). This keeps
# moby's "works out-of-the-box with conventionally-named data frames" philosophy while
# allowing a mobyData object to override these per dataset.
.mobyDefaults <- c(
  id.col       = "ID",
  datetime.col = "datetime",
  timebin.col  = "timebin",
  station.col  = "station",
  lon.col      = "lon",
  lat.col      = "lat"
)


#' Create a moby telemetry dataset
#'
#' @description
#' Wraps a detection data frame into a `mobyData` object: a lightweight `data.frame`
#' subclass that carries dataset-level metadata (column mapping, coordinate reference
#' system, tagging dates, ID groups and an optional land layer) as attributes. Downstream
#' `moby` functions read this metadata automatically, so column names, EPSG codes and
#' tagging dates do not have to be re-specified on every call.
#'
#' Because a `mobyData` object **is** a `data.frame`, it can be used anywhere a data frame
#' is expected (including with base subsetting, `dplyr`, etc.), and plain data frames remain
#' fully supported by all `moby` functions. This replaces the former global
#' `setDefaults()`/`getDefaults()` mechanism with reproducible, object-carried metadata.
#'
#' Calling `as_moby()` on an existing `mobyData` object updates only the supplied fields and
#' inherits the rest, making it easy to add metadata (e.g. `tagging.dates`) after construction.
#'
#' @param data A data frame of detections (one row per detection, or per binned record).
#' @param id.col Name of the column containing animal IDs. Defaults to `"ID"`.
#' @param datetime.col Name of the column containing date-times in POSIXct format. Defaults to `"datetime"`.
#' @param timebin.col Name of the column containing time bins (in POSIXct format). Defaults to `"timebin"`.
#' @param station.col Name of the column containing station/receiver IDs. Defaults to `"station"`.
#' @param lon.col Name of the column containing longitude (or projected x) values. Defaults to `"lon"`.
#' @param lat.col Name of the column containing latitude (or projected y) values. Defaults to `"lat"`.
#' @param epsg.code Optional integer EPSG code of a **projected** (metre-based) coordinate
#' reference system, used when projecting coordinates or computing distances/areas.
#' @param tagging.dates Optional POSIXct vector of tagging/release dates. Either a single value
#' (applied to all individuals) or a named vector whose names match the animal IDs.
#' @param nominal.delay Optional transmitter nominal (mean) delay, in seconds. Either a single
#' value (applied to all individuals) or a named numeric vector whose names match the animal IDs
#' (for arrays mixing tag families, e.g. 60 s and 120 s tags). Stored in the metadata and read
#' automatically by \code{\link{filterDetections}} to scale its short-interval (min_lag)
#' false-detection filter. Usually populated for you by \code{\link{assignAnimalIDs}} when the tag
#' table carries a delay column.
#' @param id.groups Optional named list grouping IDs (e.g. by species, sex or life stage),
#' used by many functions to compute metrics or draw plots independently per group.
#' @param land.shape Optional `sf` (or `SpatialPolygons*`) object representing landmasses,
#' used by spatial functions (e.g. in-water distances, UD land clipping).
#'
#' @return A `mobyData` object (a `data.frame` with a `"moby"` metadata attribute).
#'
#' @aliases mobyData
#' @seealso \code{\link{mobyMeta}}, \code{\link{is_moby}}
#'
#' @examples
#' df <- data.frame(
#'   ID = c("A", "A", "B"),
#'   datetime = as.POSIXct(c("2023-01-01 00:00", "2023-01-01 01:00", "2023-01-01 00:00"),
#'                         tz = "UTC"),
#'   lon = c(-8.1, -8.2, -8.0),
#'   lat = c(37.0, 37.1, 37.0),
#'   station = c("R1", "R2", "R1")
#' )
#' md <- as_moby(df, tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
#' md
#'
#' # add or update metadata later
#' md <- as_moby(md, id.groups = list(grp1 = c("A", "B")))
#'
#' @export

as_moby <- function(data,
                    id.col       = .mobyDefaults[["id.col"]],
                    datetime.col = .mobyDefaults[["datetime.col"]],
                    timebin.col  = .mobyDefaults[["timebin.col"]],
                    station.col  = .mobyDefaults[["station.col"]],
                    lon.col      = .mobyDefaults[["lon.col"]],
                    lat.col      = .mobyDefaults[["lat.col"]],
                    id.groups    = NULL,
                    land.shape   = NULL,
                    epsg.code    = NULL,
                    tagging.dates = NULL,
                    nominal.delay = NULL) {

  if (!inherits(data, "data.frame")) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  # inherit existing metadata when re-wrapping a mobyData object: only override fields
  # that were explicitly supplied by the caller
  prev <- attr(data, "moby")
  supplied <- names(as.list(match.call())[-1])

  # coerce to a plain data.frame, then re-attach the class/metadata below
  data <- as.data.frame(data)

  pick <- function(name, value) {
    if (name %in% supplied) value
    else if (!is.null(prev) && !is.null(prev[[name]])) prev[[name]]
    else value
  }

  meta <- list(
    id.col       = pick("id.col", id.col),
    datetime.col = pick("datetime.col", datetime.col),
    timebin.col  = pick("timebin.col", timebin.col),
    station.col  = pick("station.col", station.col),
    lon.col      = pick("lon.col", lon.col),
    lat.col      = pick("lat.col", lat.col),
    epsg.code    = pick("epsg.code", epsg.code),
    tagging.dates = pick("tagging.dates", tagging.dates),
    nominal.delay = pick("nominal.delay", nominal.delay),
    id.groups    = pick("id.groups", id.groups),
    land.shape   = pick("land.shape", land.shape)
  )

  # ---- validation -------------------------------------------------------------
  errors <- c()

  # the ID column must exist; other columns are validated downstream when actually used
  if (!is.character(meta$id.col) || length(meta$id.col) != 1) {
    errors <- c(errors, "'id.col' must be a single character string (column name).")
  } else if (!meta$id.col %in% colnames(data)) {
    errors <- c(errors, paste0("ID column ('", meta$id.col, "') not found in the data."))
  }
  for (nm in c("datetime.col", "timebin.col", "station.col", "lon.col", "lat.col")) {
    if (!is.character(meta[[nm]]) || length(meta[[nm]]) != 1) {
      errors <- c(errors, paste0("'", nm, "' must be a single character string (column name)."))
    }
  }
  if (!is.null(meta$epsg.code) && !(is.numeric(meta$epsg.code) && length(meta$epsg.code) == 1)) {
    errors <- c(errors, "'epsg.code' must be a single numeric EPSG code.")
  }
  if (!is.null(meta$tagging.dates) && !inherits(meta$tagging.dates, "POSIXct")) {
    errors <- c(errors, "'tagging.dates' must be provided in POSIXct format (see as.POSIXct).")
  }
  if (!is.null(meta$nominal.delay) &&
      (!is.numeric(meta$nominal.delay) || any(!is.na(meta$nominal.delay) & meta$nominal.delay <= 0))) {
    errors <- c(errors, paste("'nominal.delay' must be a positive numeric value in seconds (a single value,",
                              "or a named numeric vector keyed by animal ID)."))
  }
  if (!is.null(meta$id.groups) && (!is.list(meta$id.groups) || is.null(names(meta$id.groups)))) {
    errors <- c(errors, "'id.groups' must be a named list of ID vectors.")
  }
  if (!is.null(meta$land.shape) && !inherits(meta$land.shape, c("sf", "SpatialPolygonsDataFrame", "SpatialPolygons"))) {
    errors <- c(errors, "'land.shape' must be an 'sf' or 'SpatialPolygons' object.")
  }
  if (length(errors) > 0) {
    stop(paste0("\n", paste0("- ", errors, collapse = "\n")), call. = FALSE)
  }

  # informational note about which mapped columns are present (helps catch typos early)
  mapped <- c(meta$datetime.col, meta$timebin.col, meta$station.col, meta$lon.col, meta$lat.col)
  missing_cols <- setdiff(mapped, colnames(data))
  if (length(missing_cols) > 0) {
    message("Note: the following mapped column(s) are not present in the data and will only ",
            "matter for functions that use them: ", paste(missing_cols, collapse = ", "), ".")
  }

  attr(data, "moby") <- meta
  class(data) <- unique(c("mobyData", "data.frame"))
  data
}


#' Check whether an object is a mobyData dataset
#'
#' @param x An object.
#' @return A logical value.
#' @seealso \code{\link{as_moby}}
#' @examples
#' data(rays)
#' is_moby(rays)                # TRUE
#' is_moby(as.data.frame(rays)) # FALSE (plain data frame)
#' @export
is_moby <- function(x) inherits(x, "mobyData")


#' Retrieve mobyData metadata
#'
#' @description Returns the metadata list stored on a `mobyData` object (column mapping,
#' EPSG code, tagging dates, ID groups and land layer), or `NULL` for a plain data frame.
#' @param x A `mobyData` object (or any object).
#' @return A named list of metadata, or `NULL`.
#' @seealso \code{\link{as_moby}}
#' @examples
#' data(rays)
#' meta <- mobyMeta(rays)
#' names(meta)
#' meta$epsg.code
#' meta$id.groups
#' @export
mobyMeta <- function(x) attr(x, "moby")


#' @export
print.mobyData <- function(x, ...) {
  meta <- attr(x, "moby")
  cat("<mobyData>", nrow(x), "records x", ncol(x), "columns\n")
  if (!is.null(meta)) {
    id <- meta$id.col
    n_ids <- if (!is.null(id) && id %in% colnames(x)) length(unique(x[[id]])) else NA
    cat("  individuals: ", if (is.na(n_ids)) "?" else n_ids,
        "  (id.col = '", meta$id.col, "')\n", sep = "")
    dt <- meta$datetime.col
    if (!is.null(dt) && dt %in% colnames(x) && inherits(x[[dt]], "POSIXct")) {
      rng <- range(x[[dt]], na.rm = TRUE)
      cat("  period: ", format(rng[1]), " to ", format(rng[2]),
          " (tz = ", .dataTZ(x[[dt]]), ")\n", sep = "")
    }
    mapped <- unlist(meta[c("datetime.col", "timebin.col", "station.col", "lon.col", "lat.col")])
    cat("  columns: ", paste(sprintf("%s='%s'", names(mapped), mapped), collapse = ", "), "\n", sep = "")
    extras <- c()
    if (!is.null(meta$epsg.code))    extras <- c(extras, paste0("epsg=", meta$epsg.code))
    if (!is.null(meta$tagging.dates)) extras <- c(extras, paste0("tagging.dates (", length(meta$tagging.dates), ")"))
    if (!is.null(meta$id.groups))    extras <- c(extras, paste0("id.groups (", length(meta$id.groups), ")"))
    if (!is.null(meta$land.shape))   extras <- c(extras, "land.shape")
    if (length(extras) > 0) cat("  metadata: ", paste(extras, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}


#' @export
`[.mobyData` <- function(x, ...) {
  meta <- attr(x, "moby")
  # subset as a plain data.frame, then re-attach metadata/class when the result is
  # still a (row/column) data frame
  out <- NextMethod()
  if (is.data.frame(out)) {
    attr(out, "moby") <- meta
    class(out) <- unique(c("mobyData", "data.frame"))
  }
  out
}


#######################################################################################################
## Internal: resolve column / metadata arguments #####################################################
#######################################################################################################

#' Resolve NULL column/parameter arguments from mobyData metadata or canonical defaults
#'
#' @description Internal helper. Given a data object and a set of argument values, returns the
#' resolved values: explicitly supplied values are kept; `NULL` column arguments are filled from
#' the object's `"moby"` metadata (if present) and otherwise from the canonical defaults
#' (\code{.mobyDefaults}); `NULL` `epsg.code`/`tagging.dates`/`id.groups`/`land.shape` are filled
#' from metadata only (no canonical fallback).
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.resolveArgs <- function(data, args) {
  meta <- attr(data, "moby")
  for (nm in names(.mobyDefaults)) {
    if (nm %in% names(args) && is.null(args[[nm]])) {
      args[[nm]] <- if (!is.null(meta) && !is.null(meta[[nm]])) meta[[nm]] else unname(.mobyDefaults[[nm]])
    }
  }
  for (nm in c("epsg.code", "tagging.dates", "id.groups", "land.shape")) {
    if (nm %in% names(args) && is.null(args[[nm]]) && !is.null(meta) && !is.null(meta[[nm]])) {
      args[[nm]] <- meta[[nm]]
    }
  }
  args
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
