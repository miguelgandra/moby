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
#' fully supported by all `moby` functions. Metadata travels with the object rather than
#' living in global state, so analyses stay reproducible and self-describing.
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
#' @seealso \code{\link{mobyMeta}}, \code{\link{is_moby}}. To read a raw export from a receiver
#' system (or a table whose columns need renaming and whose date-times need parsing) before declaring
#' it here, see \code{\link{importDetections}} and \code{\link{assignAnimalIDs}}; the
#' \sQuote{Which function do I use?} section of \code{\link{moby_import_schema}} contrasts importing
#' with `as_moby()`.
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

  # informational note about mapped columns that are absent (helps catch typos early) - but only for
  # roles the caller actually asked for, i.e. supplied explicitly or carried over with a non-default
  # name. A canonical default nobody requested is just moby's guess, and its absence is normal (a
  # dataset has no 'timebin' column until getTimeBins() has been run), so flagging it is pure noise.
  role_cols <- c("datetime.col", "timebin.col", "station.col", "lon.col", "lat.col")
  asked <- vapply(role_cols, function(nm)
    nm %in% supplied || !identical(unname(meta[[nm]]), unname(.mobyDefaults[[nm]])), logical(1))
  mapped <- unlist(meta[role_cols[asked]], use.names = FALSE)
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


# glyphs used by print.mobyData, with ASCII fallbacks for non-UTF-8 consoles (Windows in a legacy
# locale would otherwise render mojibake)
.mobyGlyphs <- function() {
  # \u escapes keep this file pure ASCII (so it behaves identically whatever locale R starts in)
  if (isTRUE(l10n_info()$`UTF-8`)) list(rule = "\u2500", sep = "\u00b7", arrow = "\u2192")
  else                             list(rule = "-",      sep = "|",      arrow = "->")
}

# thousands-separated integer, locale-independent
.fmtCount <- function(n) formatC(n, format = "d", big.mark = ",")

# pack items onto lines at item boundaries (never leaving a separator stranded at a line start)
.packItems <- function(items, sep, avail) {
  lines <- character(0); cur <- ""
  for (it in items) {
    cand <- if (nzchar(cur)) paste0(cur, " ", sep, " ", it) else it
    if (nchar(cand) > avail && nzchar(cur)) { lines <- c(lines, cur); cur <- it } else cur <- cand
  }
  if (nzchar(cur)) lines <- c(lines, cur)
  lines
}

# emit "label      value", with continuation lines aligned under the value column
.mobyField <- function(label, lines, indent = 11L) {
  if (length(lines) == 0) return(invisible(NULL))
  cat(formatC(label, width = -indent), lines[1], "\n", sep = "")
  for (l in lines[-1]) cat(strrep(" ", indent), l, "\n", sep = "")
  invisible(NULL)
}

# Keep the ID factor and the per-animal metadata in step with the rows that are actually present.
# Used after a step that DISCARDS detections: without it, an animal that loses all of its detections
# survives as an empty factor level (and a stale tagging.date / nominal.delay entry), which downstream
# per-individual loops then treat as a real animal with no data. Applying it also makes such steps
# commute - e.g. matching deployments before or after assigning animal IDs gives the same object.
.syncIDs <- function(data, id.col, meta) {
  if (!is.null(id.col) && id.col %in% colnames(data) && is.factor(data[[id.col]]))
    data[[id.col]] <- droplevels(data[[id.col]])
  present <- if (!is.null(id.col) && id.col %in% colnames(data))
    as.character(unique(data[[id.col]])) else character(0)
  for (nm in c("tagging.dates", "nominal.delay")) {
    v <- meta[[nm]]
    if (!is.null(v) && !is.null(names(v))) {
      keep <- names(v) %in% present
      meta[[nm]] <- if (any(keep)) v[keep] else NULL
    }
  }
  list(data = data, meta = meta)
}

# What does one row of this object represent? Decided from the columns, never from a stored label:
# moby already answers this question by column presence elsewhere (summaryTable() sums a 'detections'
# count column when present, else assumes one detection per row), and columns cannot go stale the way
# a declared type would when an object is transformed. Unrecognised shapes fall back to "records".
.mobyRowNoun <- function(x) {
  meta <- attr(x, "moby"); cn <- colnames(x)
  # raw/filtered detections: one row per decode, so a datetime column is present
  if (!is.null(meta$datetime.col) && meta$datetime.col %in% cn) return("detections")
  # aggregated positions (COAs, tracks): rows are time bins carrying a detection COUNT
  if ("detections" %in% cn && is.numeric(x[["detections"]])) return("positions")
  "records"
}

# human-readable span between two times
.fmtDuration <- function(from, to) {
  d <- as.numeric(difftime(to, from, units = "days"))
  if (!is.finite(d)) return(NA_character_)
  if (d >= 365)    sprintf("%.1f years", d / 365.25)
  else if (d >= 1) sprintf("%.0f day%s", d, if (round(d) == 1) "" else "s")
  else             sprintf("%.1f hours", d * 24)
}


#' @param x A `mobyData` object.
#' @param preview Number of leading rows shown under the summary. Defaults to 3; `0` prints the
#' summary only.
#' @param width Console width used to lay the summary out. Defaults to `getOption("width")`.
#' @param ... Ignored.
#' @rdname as_moby
#' @export
print.mobyData <- function(x, preview = 3L, width = getOption("width"), ...) {

  meta <- attr(x, "moby")
  g <- .mobyGlyphs()
  w <- max(40L, min(as.integer(width), 100L))

  # ---- header rule -----------------------------------------------------------------------------
  hdr <- paste0(strrep(g$rule, 2), " <mobyData> ")
  cat(hdr, strrep(g$rule, max(0L, w - nchar(hdr))), "\n\n", sep = "")

  # ---- overview --------------------------------------------------------------------------------
  noun <- .mobyRowNoun(x)
  id <- meta$id.col
  n_ids <- if (!is.null(id) && id %in% colnames(x)) length(unique(stats::na.omit(x[[id]]))) else NA_integer_
  # a position table knows how many detections it was aggregated from - report both, so the row count
  # is never mistaken for the detection count
  n_det <- if (noun == "positions") sum(x[["detections"]], na.rm = TRUE) else NA_real_
  overview <- c(paste(.fmtCount(nrow(x)), noun),
                if (is.finite(n_det) && n_det > 0) paste(.fmtCount(round(n_det)), "detections"),
                if (!is.na(n_ids)) paste(.fmtCount(n_ids), "individuals"),
                paste(.fmtCount(ncol(x)), "variables"))
  .mobyField("overview", .packItems(overview, g$sep, w - 11L))

  # ---- period ----------------------------------------------------------------------------------
  dt <- meta$datetime.col
  if (!is.null(dt) && dt %in% colnames(x) && inherits(x[[dt]], "POSIXct") &&
      any(!is.na(x[[dt]]))) {
    rng <- range(x[[dt]], na.rm = TRUE)
    dur <- .fmtDuration(rng[1], rng[2])
    span <- paste0(format(rng[1], "%Y-%m-%d"), " ", g$arrow, " ", format(rng[2], "%Y-%m-%d"))
    qual <- paste(c(dur, paste0("tz = ", .dataTZ(x[[dt]]))), collapse = ", ")
    .mobyField("period", paste0(span, " (", qual, ")"))
  }

  # ---- space -----------------------------------------------------------------------------------
  st <- meta$station.col
  space <- character(0)
  if (!is.null(st) && st %in% colnames(x))
    space <- c(space, paste(.fmtCount(length(unique(stats::na.omit(x[[st]])))), "stations"))
  for (ax in c("lon", "lat")) {
    cc <- meta[[paste0(ax, ".col")]]
    if (!is.null(cc) && cc %in% colnames(x) && is.numeric(x[[cc]]) && any(!is.na(x[[cc]]))) {
      r <- range(x[[cc]], na.rm = TRUE)
      space <- c(space, sprintf("%s [%.2f, %.2f]", ax, r[1], r[2]))
    }
  }
  .mobyField("space", .packItems(space, g$sep, w - 11L))

  # ---- metadata (only what is actually attached) -----------------------------------------------
  extras <- c(
    if (!is.null(meta$tagging.dates)) paste0("tagging.dates (", length(meta$tagging.dates), ")"),
    if (!is.null(meta$nominal.delay)) paste0("nominal.delay (", length(meta$nominal.delay), ")"),
    if (!is.null(meta$id.groups))     paste0("id.groups (", length(meta$id.groups), ")"),
    if (!is.null(meta$epsg.code))     paste0("EPSG: ", meta$epsg.code),
    if (!is.null(meta$land.shape))    "land.shape")
  .mobyField("metadata", .packItems(extras, g$sep, w - 11L))

  # ---- preview ---------------------------------------------------------------------------------
  preview <- as.integer(preview)
  if (!is.na(preview) && preview > 0L && nrow(x) > 0L && ncol(x) > 0L) {
    n <- min(preview, nrow(x))
    cat("\nPreview (first ", n, if (n == 1L) " row)" else " rows)", "\n", sep = "")
    # print as a plain data.frame: shows every variable and wraps to the console width, so nothing
    # is silently hidden
    print(utils::head(as.data.frame(x), n))
  }

  cat("\n", strrep(g$rule, w), "\n", sep = "")
  invisible(x)
}


#' Return the first or last rows of a `mobyData`
#'
#' @description `head()` and `tail()` behave as they do for any data frame: they return the first
#' (or last) `n` rows of the detection table for inspection. They return a plain `data.frame`, so
#' the rows print as rows - the study metadata summary is what \code{print()} on the `mobyData`
#' itself is for.
#'
#' To subset a dataset and *keep* it a `mobyData` (metadata and all), use `[` instead:
#' `x[1:100, ]`.
#'
#' @param x A \code{\link{mobyData}} object.
#' @param n Number of rows. Defaults to 6, as for \code{utils::head}.
#' @param ... Further arguments passed to \code{utils::head} / \code{utils::tail}.
#'
#' @return A `data.frame` of the first (or last) `n` rows.
#'
#' @seealso \code{\link{as_moby}}, \code{\link{mobyMeta}}
#' @examples
#' data(rays)
#' head(rays, 3)      # rows, like any data frame
#' dim(rays[1:100, ]) # '[' keeps the mobyData class and its metadata
#'
#' @export
head.mobyData <- function(x, n = 6L, ...) utils::head(as.data.frame(x), n = n, ...)

#' @rdname head.mobyData
#' @export
tail.mobyData <- function(x, n = 6L, ...) utils::tail(as.data.frame(x), n = n, ...)


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
