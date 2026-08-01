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
#' Unlike the analytical functions, `as_moby()` prints no banner and no completion summary: the
#' object it returns already describes itself when printed (see the `print()` method below). Its only
#' console output is a note listing mapped columns that are absent from the data, which `verbose`
#' (or `options(moby.verbose = FALSE)`) silences.
#'
#' @param data A data frame of detections (one row per detection, or per binned record).
#' @param tags Optional tag table (see \code{\link{importTags}}) from which to derive the
#' per-animal `tagging.dates` and `nominal.delay` metadata, keyed by the IDs present in `data`.
#' This is where tag-derived metadata enters a dataset: \code{\link{matchTags}} joins the tag
#' table to the detections but attaches nothing, so nothing is added to your object without an
#' explicit `as_moby()` call. Values passed directly to `tagging.dates` or `nominal.delay` take
#' precedence over anything derived here.
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
#' false-detection filter. Usually populated by passing `tags` to `as_moby()`, when the tag
#' table carries a delay column.
#' @param id.groups Optional named list grouping IDs (e.g. by species, sex or life stage),
#' used by many functions to compute metrics or draw plots independently per group.
#' @param land.shape Optional `sf` (or `SpatialPolygons*`) object representing landmasses,
#' used by spatial functions (e.g. in-water distances, UD land clipping).
#' @param verbose Logical; print informational notes about the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A `mobyData` object (a `data.frame` with a `"moby"` metadata attribute).
#'
#' @aliases mobyData
#' @seealso \code{\link{mobyMeta}}, \code{\link{is_moby}}. To read a raw export from a receiver
#' system (or a table whose columns need renaming and whose date-times need parsing) before declaring
#' it here, see \code{\link{importDetections}} and \code{\link{matchTags}}; the
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
                    tags         = NULL,
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
                    nominal.delay = NULL,
                    verbose = getOption("moby.verbose", TRUE)) {

  if (!inherits(data, "data.frame")) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  # inherit existing metadata when re-wrapping a mobyData object: only override fields
  # that were explicitly supplied by the caller
  prev <- attr(data, "moby")
  # `supplied` records which arguments the caller named, and is only ever queried for the eleven
  # metadata field names (via pick(), below) and the five role columns (via `asked`, further down).
  # Non-metadata arguments such as `verbose` therefore land in this vector harmlessly: nothing
  # consults them, so passing verbose = FALSE cannot alter which metadata is inherited or reported.
  supplied <- names(as.list(match.call())[-1])

  land_expr <- substitute(land.shape)
  land_name <- if (is.symbol(land_expr)) as.character(land_expr) else NULL

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

  # ---- metadata derived from a tag table --------------------------------------
  # as_moby() is the package's only constructor, so per-animal metadata enters the object HERE
  # rather than being attached invisibly by the join step that happens to have seen the tag table.
  # Precedence: an explicitly supplied tagging.dates/nominal.delay always wins, then whatever the
  # tag table yields, then metadata inherited from a mobyData being re-declared.
  if (!is.null(tags)) {
    tg <- as.data.frame(tags)
    if (!any(c("ID", "serial", "transmitter") %in% colnames(tg))) {
      stop("'tags' must contain an 'ID', 'serial' or 'transmitter' column (see importTags()).",
           call. = FALSE)
    }
    tg$ID <- .tagIDs(tg)
    # Key on the DECLARED roster, not the observed values: a factor ID column states which animals
    # the study covers, and an animal that was tagged but never detected is a legitimate level (it
    # keeps per-individual summaries honest). filterDetections() validates tagging.dates against
    # those levels, so deriving from unique() instead would attach 24 dates for a 26-animal roster
    # and then reject the dataset. Fall back to the observed values for a plain character column,
    # which declares nothing. A missing id.col is left to the validation below to report.
    if (is.character(meta$id.col) && length(meta$id.col) == 1 && meta$id.col %in% colnames(data)) {
      id_values <- data[[meta$id.col]]
      ids <- if (is.factor(id_values)) levels(id_values) else as.character(unique(id_values))
      if (!("tagging.dates" %in% supplied)) {
        derived <- .deriveTaggingDates(tg, ids)
        if (!is.null(derived)) meta$tagging.dates <- derived
      }
      if (!("nominal.delay" %in% supplied)) {
        derived <- .deriveNominalDelay(tg, ids)
        if (!is.null(derived)) meta$nominal.delay <- derived
      }
    }
  }

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
    .mobyNote("The following mapped column(s) are not present in the data and will only ",
              "matter for functions that use them: ", paste(missing_cols, collapse = ", "), ".",
              verbose = verbose)
  }

  attr(data, "moby") <- meta
  # The name the caller gave the land layer, kept for print() only. Deliberately NOT part of `meta`:
  # the internal callers rebuild a dataset with do.call(as_moby, meta), so anything in that list has
  # to be a formal argument, and a display label has no business in the user-facing signature.
  if ("land.shape" %in% supplied) attr(data, "moby.land.name") <- land_name
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


#' Demote a dataset to a plain data frame.
#'
#' `as.data.frame()` drops `oldClass()` but keeps every other attribute, so coercing a `mobyData`
#' left the metadata attached to an object that `is_moby()` reported as plain - and every metadata
#' read in the package goes to `attr(x, "moby")` rather than checking the class, so that orphaned
#' metadata stayed *live*. One function's "plain data frame" was another's fully-configured dataset.
#'
#' Class and metadata are two halves of one fact, so they are removed together. Use this wherever a
#' function coerces its input for plain-data-frame indexing, and note that anything needing the
#' metadata must read it BEFORE coercing (which is why `.validateArguments()` resolves column roles
#' first and demotes second).
#' @keywords internal
#' @noRd
.stripMoby <- function(x) {
  # drop the subclass BEFORE coercing, so as.data.frame() dispatches to the data.frame (or tibble /
  # data.table) method rather than recursing into as.data.frame.mobyData() below
  if (inherits(x, "mobyData")) class(x) <- setdiff(class(x), "mobyData")
  x <- as.data.frame(x)
  attr(x, "moby") <- NULL
  attr(x, "moby.land.name") <- NULL
  x
}

#' Demote a mobyData to a plain data frame
#'
#' @description Drops the `mobyData` class **and** its metadata, returning an ordinary data frame.
#' Without this method `as.data.frame()` would remove the class but keep the metadata attribute,
#' leaving an object that `is_moby()` calls plain while moby's column resolution still reads its
#' stored roles.
#'
#' @param x A `mobyData` object.
#' @param ... Ignored.
#' @return A plain `data.frame`, carrying no moby metadata.
#' @seealso \code{\link{as_moby}}, \code{\link{is_moby}}, \code{\link{mobyMeta}}
#' @examples
#' data(rays)
#' plain <- as.data.frame(rays)
#' is_moby(plain)         # FALSE
#' mobyMeta(plain)        # NULL - the metadata is gone, not merely hidden
#' @export
as.data.frame.mobyData <- function(x, ...) .stripMoby(x)

#' Re-attach mobyData metadata to a result, when the input carried it.
#'
#' The inverse of `.stripMoby()`, and the package's single convention for class handling: an
#' operation PRESERVES its input's class (mobyData in, mobyData out; plain in, plain out), while
#' only `as_moby()` CONSTRUCTS. Passing `meta = NULL` - which is what `attr(data, "moby")` yields
#' for a plain input - is therefore not an edge case but the ordinary plain-data-frame path.
#'
#' Deliberately re-attaches rather than calling `as_moby()`: the constructor fills unsupplied roles
#' from `.mobyDefaults`, so rebuilding through it silently rewrites a caller's column map to the
#' canonical names (see the bug this replaced in `matchDeployments()`).
#' @keywords internal
#' @noRd
.restoreClass <- function(x, meta) {
  if (is.null(meta) || !is.data.frame(x)) return(x)
  attr(x, "moby") <- meta
  class(x) <- unique(c("mobyData", "data.frame"))
  x
}


# (thousands-separated integers come from .fmtN() in helpers-messaging.R; this file used to carry a
# byte-identical private copy, which collided with the count+noun .fmtCount() helper)

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
# Title-case the first letter, leaving the rest alone (so "min_lag" style names keep their shape)
.capitalise <- function(x) sub("^(.)", "\\U\\1", x, perl = TRUE)

# the CRS may be a bare code or an sf crs object (whose $input already reads "EPSG:32629")
.epsgLabel <- function(e) {
  if (inherits(e, "crs")) sub("^EPSG:", "", as.character(e$input)[1]) else as.character(e)[1]
}

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
#' @param width Console width used to lay the summary out. Defaults to `getOption("width")`.
#' @param ... Ignored.
#' @rdname as_moby
#' @export
print.mobyData <- function(x, width = getOption("width"), ...) {

  meta <- attr(x, "moby")
  g <- .mobyGlyphs()
  w <- .printWidth(width)
  .printOpen("mobyData", w)

  # ---- Input: the scale of what is held ----------------------------------------------------------
  noun <- .mobyRowNoun(x)
  n_det <- if (noun == "positions") sum(x[["detections"]], na.rm = TRUE) else NA_real_
  scale <- c(.fmtCount(nrow(x), sub("s$", "", noun)),
             if (is.finite(n_det) && n_det > 0) .fmtCount(round(n_det), "detection"))
  id <- meta$id.col
  if (!is.null(id) && id %in% colnames(x)) {
    # TAGGED is the declared roster (the factor's levels), DETECTED the animals with data. They differ
    # whenever an animal was tagged and never heard from, which is worth seeing rather than inferring.
    n_obs <- .nObserved(x[[id]])
    n_tagged <- if (is.factor(x[[id]])) nlevels(x[[id]]) else n_obs
    scale <- c(scale,
               if (n_tagged > n_obs) paste0(.fmtCount(n_tagged, "individual"), " tagged"),
               paste0(.fmtCount(n_obs, "individual"), if (n_tagged > n_obs) " detected" else ""))
  }
  st <- meta$station.col
  if (!is.null(st) && st %in% colnames(x))
    scale <- c(scale, .fmtCount(length(unique(stats::na.omit(x[[st]]))), "station"))
  scale <- c(scale, .fmtCount(ncol(x), "variable"))
  .printSection("Summary", stats::setNames(rep("", length(scale)), scale), blank.before = FALSE)

  # ---- Coverage: the extent the data span --------------------------------------------------------
  dt <- meta$datetime.col
  has_time <- !is.null(dt) && dt %in% colnames(x) && inherits(x[[dt]], "POSIXct") && any(!is.na(x[[dt]]))
  axes <- Filter(function(a) {
    cc <- meta[[paste0(a, ".col")]]
    !is.null(cc) && cc %in% colnames(x) && is.numeric(x[[cc]]) && any(!is.na(x[[cc]]))
  }, c("lon", "lat"))
  cov <- character(0)
  if (has_time) {
    rng <- range(x[[dt]], na.rm = TRUE)
    cov["Time"] <- paste0(format(rng[1], "%Y-%m-%d"), " ", g$arrow, " ", format(rng[2], "%Y-%m-%d"),
                          " (", .fmtDuration(rng[1], rng[2]), ")")
  }
  for (a in axes) {
    r <- range(x[[meta[[paste0(a, ".col")]]]], na.rm = TRUE)
    cov[c(lon = "Longitude", lat = "Latitude")[[a]]] <- sprintf("%.2f %s %.2f", r[1], g$arrow, r[2])
  }
  if (has_time) cov["Time zone"] <- .dataTZ(x[[dt]])
  .printSection("Coverage", cov)

  # ---- Metadata: every supported field, present or not -------------------------------------------
  # Listing the absent ones too is the point: it says at a glance what this dataset can and cannot do
  # downstream (no tagging dates means no residency; no CRS means no distances).
  n_roster <- if (!is.null(id) && id %in% colnames(x) && is.factor(x[[id]])) nlevels(x[[id]]) else NA_integer_
  covered <- function(v) {
    if (is.na(n_roster) || is.null(names(v))) paste0("(", length(v), ")")
    else paste0("(", sum(levels(x[[id]]) %in% names(v)), "/", n_roster, ")")
  }
  entries <- list(
    "Tagging dates"  = if (!is.null(meta$tagging.dates)) covered(meta$tagging.dates),
    "Nominal delays" = if (!is.null(meta$nominal.delay)) covered(meta$nominal.delay),
    "ID groups"      = if (!is.null(meta$id.groups)) paste0("(", length(meta$id.groups), ")"),
    "Land polygon"   = if (!is.null(meta$land.shape))
                         if (!is.null(attr(x, "moby.land.name")))
                           paste0("(", attr(x, "moby.land.name"), ")") else "",
    "Coordinate CRS" = if (!is.null(meta$epsg.code)) paste0("EPSG:", .epsgLabel(meta$epsg.code)))
  present <- !vapply(entries, is.null, logical(1))
  entries[!present] <- ""
  .printSection("Metadata", unlist(entries), marker = ifelse(present, "tick", "cross"))

  .printClose(w)
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
head.mobyData <- function(x, n = 6L, ...) utils::head(.stripMoby(x), n = n, ...)

#' @rdname head.mobyData
#' @export
tail.mobyData <- function(x, n = 6L, ...) utils::tail(.stripMoby(x), n = n, ...)


#' @export
`[.mobyData` <- function(x, ...) {
  meta <- attr(x, "moby")
  # subset as a plain data.frame, then re-attach metadata/class when the result is
  # still a (row/column) data frame
  .restoreClass(NextMethod(), meta)
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
