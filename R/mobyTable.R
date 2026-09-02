#######################################################################################################
## mobyTable: one typed summary table, and one place that formats it ##################################
##
## moby's summary tables used to be built already formatted: every column a character string, the
## mean row and the group headings baked in as data rows, missing values written as "-". That is fine
## to look at once and useless afterwards - mean(tbl$`N Detect`) fails, and writing the table to CSV
## exports the blank group-heading rows as records.
##
## So the tables now return their NUMBERS, and everything about presentation happens on the way out:
##
##   summaryTable()/movementTable()/transitionsTable()  ->  typed mobyTable (numeric, POSIXct)
##   format(x)                                          ->  character table, ready for write.csv()
##   print(x)                                           ->  the console rendering
##
## This is the same split the rest of the package already uses (mobyData, mobyFilter, mobyQC all
## print through helpers-print.R); the tables were the only things still formatting themselves.
#######################################################################################################


#' Construct a mobyTable.
#'
#' @param df The typed table: numeric columns numeric, dates POSIXct/Date.
#' @param kind Which table this is ("summary", "movement", "transitions"); selects the header
#'   dictionary and the display precision.
#' @param error.stat "sd" or "se" - which statistic the mean row reports.
#' @param label.col The column that names each row (the ID column, or the transition). The mean row
#'   writes its label there.
#' @param units Optional named list of runtime-decided units, e.g. `list(rom = "km/h")`, substituted
#'   into the headers at format time. Kept out of the column NAMES so an override keyed on a column
#'   does not stop working when the data changes the unit.
#' @param extra Any further attributes to carry (residency.index, id.groups, ...).
#' @keywords internal
#' @noRd
.newMobyTable <- function(df, kind, error.stat = "sd", label.col = "ID", units = list(),
                          extra = list()) {
  rownames(df) <- NULL
  attr(df, "moby.table") <- list(kind = kind, error.stat = error.stat, label.col = label.col,
                                 units = units)
  for (nm in names(extra)) attr(df, nm) <- extra[[nm]]
  class(df) <- unique(c("mobyTable", "data.frame"))
  df
}

#' The metadata a mobyTable carries, with defaults for a table that lost it (e.g. through `[`).
#' @keywords internal
#' @noRd
.tableMeta <- function(x) {
  m <- attr(x, "moby.table")
  if (is.null(m)) m <- list()
  if (is.null(m$kind)) m$kind <- "summary"
  if (is.null(m$error.stat)) m$error.stat <- "sd"
  if (is.null(m$label.col)) m$label.col <- names(x)[1]
  if (is.null(m$units)) m$units <- list()
  m
}


#' Display precision per column, by table kind.
#'
#' Fixed, not derived from the data. The tables used to take the maximum number of decimals any value
#' in a column happened to be written with, so a single 0.123456 dragged the whole column to six
#' decimal places and the same metric printed differently for different datasets.
#' @keywords internal
#' @noRd
.tablePrecision <- function(kind) {
  switch(kind,
    summary = c(n_detections = 0, n_receivers = 0, monitoring_duration_d = 0, detection_span_d = 0,
                n_days_detected = 0, IR1 = 2, IR2 = 2, IWR = 2, `IR2/IR1` = 2),
    movement = c(distance_km = 1, rom = 1, rom_max = 1, linearity_index = 2),
    transitions = c(n_movements = 0, n_individuals = 0, pct_individuals = 0,
                    mean_duration = 1, error_duration = 1),
    numeric(0))
}

#' Columns a cross-individual mean must not be taken over.
#'
#' A percentage that is already a share of the group, and anything that is an identifier rather than
#' a measurement. Averaging them produces a number, just not one the row's label claims.
#' @keywords internal
#' @noRd
.tableNoMeanCols <- function(kind) {
  switch(kind,
    transitions = c("pct_individuals"),
    character(0))
}


#' Human-readable headers, per kind and style.
#'
#' Internal snake_case names are what the API guarantees - `decimals` and `group.by` are keyed on
#' them, and they stay put when the display style changes. These dictionaries are the display layer.
#' @keywords internal
#' @noRd
.tableHeaders <- function(cols, meta, style = "report") {
  u <- function(nm, dflt) if (!is.null(meta$units[[nm]])) meta$units[[nm]] else dflt
  err <- toupper(meta$error.stat)
  report <- switch(meta$kind,
    summary = c(
      tagging_date = "Tagging date", first_detection = "First detection",
      last_detection = "Last detection", n_detections = "N Detect", n_receivers = "N Receiv",
      monitoring_duration_d = "Monitoring duration (d)", detection_span_d = "Detection span (d)",
      n_days_detected = "N days detected", group = "Group"),
    movement = c(
      distance_km = "Distance (km)", rom = paste0("ROM (", u("rom", "m/h"), ")"),
      rom_max = paste0("Max ROM (", u("rom", "m/h"), ")"),
      linearity_index = "Linearity index", group = "Group"),
    transitions = c(
      transition = "Transition", n_movements = "Movements", n_individuals = "Individuals (n)",
      pct_individuals = "Individuals (%)",
      mean_duration = paste0("Mean duration (", u("duration", "h"), ")"),
      error_duration = paste0(err, " duration (", u("duration", "h"), ")"), group = "Group"),
    character(0))
  concise <- switch(meta$kind,
    summary = c(
      tagging_date = "Tagged", first_detection = "First det.", last_detection = "Last det.",
      n_detections = "N det.", n_receivers = "N rec.",
      monitoring_duration_d = "Monit. (d)", detection_span_d = "Span (d)",
      n_days_detected = "Days det.", group = "Group"),
    movement = c(
      distance_km = "Dist. (km)", rom = paste0("ROM (", u("rom", "m/h"), ")"),
      rom_max = paste0("Max ROM (", u("rom", "m/h"), ")"),
      linearity_index = "LI", group = "Group"),
    transitions = c(
      transition = "Type", n_movements = "Moves", n_individuals = "Ind. (n)",
      pct_individuals = "Ind. (%)",
      mean_duration = paste0("Mean dur. (", u("duration", "h"), ")"),
      error_duration = paste0(err, " dur. (", u("duration", "h"), ")"), group = "Group"),
    character(0))

  dict <- if (identical(style, "concise")) concise else report
  out <- unname(dict[cols])
  # Anything the dictionary does not know - the ID column under whatever name the caller gave it, a
  # residency index, a sensor statistic, a user metadata column - keeps its own name, lightly
  # prettified. Renaming those would break the link back to where they came from.
  miss <- is.na(out)
  if (any(miss)) out[miss] <- vapply(cols[miss], function(s) {
    if (grepl("^(IR|IWR)", s) || s == meta$label.col) return(s)
    s <- gsub("_", " ", s)
    paste0(toupper(substring(s, 1, 1)), substring(s, 2))
  }, character(1))
  out
}

#' ASCII spellings for the typographic symbols a header may carry.
#'
#' Applied on the way out of format(), so the dictionaries stay readable and one function owns every
#' substitution.
#' @keywords internal
#' @noRd
.foldTableSymbols <- function(s) {
  s <- gsub("\u00b0C", "deg C", s, fixed = TRUE)
  s <- gsub("\u00b0", " deg", s, fixed = TRUE)
  s <- gsub("\u00b1", "+/-", s, fixed = TRUE)
  s
}

#' The grouping column as a factor, missing values kept as their own trailing group.
#'
#' A row whose group is unknown is a fact about the data, not a row to drop.
#' @keywords internal
#' @noRd
.tableGroupFactor <- function(v, missing.label = "(missing)") {
  na <- is.na(v)
  lv <- if (is.factor(v)) levels(droplevels(v[!na])) else sort(unique(as.character(v[!na])))
  ch <- as.character(v); ch[na] <- missing.label
  factor(ch, levels = c(lv, if (any(na)) missing.label))
}

#' One `mean +/- error` row over a set of rows.
#'
#' Shared by the ungrouped footer and the per-group ones, so the two cannot drift.
#' @keywords internal
#' @noRd
.tableFooterRow <- function(sub, agg_cols, prec_of, errfun, pm, err_stat, cols, label.col) {
  foot <- stats::setNames(rep(NA_character_, length(cols)), cols)
  for (nm in agg_cols) {
    m <- mean(sub[[nm]], na.rm = TRUE)
    if (is.finite(m)) {
      e <- errfun(sub[[nm]])
      foot[nm] <- if (is.finite(e))
        sprintf(paste0("%.", prec_of(nm), "f ", pm, " %.", prec_of(nm), "f"), m, e)
      else sprintf(paste0("%.", prec_of(nm), "f"), m)
    }
  }
  foot[[label.col]] <- paste0("mean ", pm, " ", err_stat)
  as.data.frame(as.list(foot), stringsAsFactors = FALSE, check.names = FALSE)
}

#' Validate a `decimals` override.
#' @keywords internal
#' @noRd
.assertTableDecimals <- function(decimals, df, meta) {
  if (is.null(decimals)) return(stats::setNames(integer(0), character(0)))
  if (!is.numeric(decimals) || !length(decimals) || is.null(names(decimals)) ||
      any(!nzchar(names(decimals))))
    stop("'decimals' must be a named numeric vector, e.g. c(distance_km = 2). Use the column names ",
         "format(x, style = \"internal\") shows.", call. = FALSE)
  if (anyNA(decimals) || any(!is.finite(decimals)) || any(decimals < 0) ||
      any(decimals != trunc(decimals)))
    stop("'decimals' must be whole numbers of decimal places, zero or more.", call. = FALSE)
  unknown <- setdiff(names(decimals), names(df))
  if (length(unknown)) {
    # a display header the caller read off the formatted table, mapped back to its column
    hits <- unlist(lapply(c("report", "concise"), function(st) {
      h <- .foldTableSymbols(.tableHeaders(names(df), meta, st))
      stats::setNames(names(df), h)[intersect(unknown, h)]
    }))
    hits <- hits[!duplicated(names(hits))]
    hint <- if (length(hits))
      paste0(" '", names(hits)[1], "' is a display header; use '", unname(hits)[1], "'.")
    else " Column names are the ones format(x, style = \"internal\") shows."
    stop("'decimals' names column(s) not in this table: ", paste(unknown, collapse = ", "), ".",
         hint, call. = FALSE)
  }
  bad <- names(decimals)[!vapply(df[names(decimals)], is.numeric, logical(1))]
  if (length(bad))
    stop("'decimals' sets a decimal place on the non-numeric column(s): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  stats::setNames(as.integer(decimals), names(decimals))
}


#' Format a moby summary table for display or export
#'
#' @description Renders a `mobyTable` as a character data frame: fixed per-metric precision, missing
#' values as `"-"`, and (where a group has more than one row) a display-only `mean +/- error` row.
#' This is the export route - `write.csv(format(x), "table.csv", row.names = FALSE)` - and the same
#' rendering `print()` shows.
#'
#' The returned table stays rectangular. Group headings and the blank lines between groups belong to
#' the console rendering, not to the exported object, so a CSV never carries an empty record.
#'
#' @param x A `mobyTable` (from \code{\link{summaryTable}}, \code{\link{movementTable}} or
#'   \code{\link{transitionsTable}}).
#' @param style Column-name style. `"internal"` (default) keeps the snake_case names the API
#'   guarantees; `"report"` uses publication-ready headers (`n_detections` -> "N Detect");
#'   `"concise"` uses the same headers abbreviated for narrow tables. Only the names differ.
#' @param symbols Whether the rendered table may use typographic symbols: `"ascii"` (default) writes
#'   `+/-`, `"unicode"` the plus-minus sign. ASCII is the default because this table is usually
#'   written to a file, and a spreadsheet opening a UTF-8 CSV with no byte-order mark guesses the
#'   encoding. `print()` picks the right one for the terminal on its own.
#' @param decimals Optional per-column override of the display precision, as a named numeric vector
#'   of decimal places - `c(distance_km = 2)`. Merged OVER the built-in precision, so naming one
#'   column leaves the rest untouched, and both the values and the `mean +/- error` row follow it.
#'   Named by the INTERNAL column names, which do not change with `style`.
#' @param group.by Column to group the rendering by. Defaults to the table's own `group` column when
#'   it has one (the `id.groups` split), or to no grouping. Each group gets its own
#'   `mean +/- error` row. Pass `FALSE` to render ungrouped.
#' @param include.summary.row Logical; append the display-only `mean +/- error` row(s). Default
#'   `TRUE`. A group of ONE row never gets one: the "mean" would just restate that row.
#' @param datetime.format `strftime` format for date/time columns. Default `"%d/%m/%Y"`.
#' @param ... Unused.
#' @return A character `data.frame`.
#' @seealso \code{\link{summaryTable}}, \code{\link{movementTable}}, \code{\link{transitionsTable}}
#' @examples
#' data(rays)
#' tbl <- summaryTable(rays, last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"))
#' # the object itself is typed - you can compute on it
#' mean(tbl$n_detections)
#' # ...and format() renders the version you export
#' head(format(tbl, style = "report"))
#' @exportS3Method format mobyTable
format.mobyTable <- function(x, style = c("internal", "report", "concise"),
                             symbols = c("ascii", "unicode"), decimals = NULL,
                             group.by = NULL, include.summary.row = TRUE,
                             datetime.format = "%d/%m/%Y", ...) {
  style <- match.arg(style)
  symbols <- match.arg(symbols)
  meta <- .tableMeta(x)
  df <- as.data.frame(x)
  attr(df, "moby.table") <- NULL
  if (nrow(df) == 0 || ncol(df) == 0) return(data.frame())

  # default to the table's own grouping; FALSE renders ungrouped
  if (is.null(group.by)) group.by <- if ("group" %in% names(df)) "group" else FALSE
  grp <- NULL
  if (!isFALSE(group.by)) {
    if (!is.character(group.by) || length(group.by) != 1 || !group.by %in% names(df))
      stop("'group.by' must name a column of the table. Available: ",
           paste(names(df), collapse = ", "), call. = FALSE)
    grp <- .tableGroupFactor(df[[group.by]])
    ord <- order(as.integer(grp), seq_along(grp))    # stable: incoming order kept within a group
    df <- df[ord, , drop = FALSE]; grp <- grp[ord]
  }

  dec <- .assertTableDecimals(decimals, df, meta)
  err_stat <- meta$error.stat
  pm <- if (identical(symbols, "unicode")) "\u00b1" else "+/-"
  errfun <- if (err_stat == "se") function(v) .stdError(v) else function(v) stats::sd(v, na.rm = TRUE)

  prec_map <- .tablePrecision(meta$kind)
  prec_map[names(dec)] <- dec                        # user override merges over the defaults
  prec_of <- function(nm) { p <- unname(prec_map[nm]); if (is.na(p)) 2L else as.integer(p) }

  num_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  agg_cols <- setdiff(num_cols, .tableNoMeanCols(meta$kind))

  disp <- as.data.frame(lapply(names(df), function(nm) {
    col <- df[[nm]]
    if (inherits(col, c("POSIXt", "Date"))) format(col, datetime.format)
    else if (nm %in% num_cols) ifelse(is.na(col), NA_character_,
                                      sprintf(paste0("%.", prec_of(nm), "f"), col))
    else as.character(col)
  }), stringsAsFactors = FALSE)
  names(disp) <- names(df)

  # Display-only mean rows. A group of one never gets one: its "mean" is the row itself, so the line
  # would add a restatement rather than a summary.
  row_grp <- if (is.null(grp)) NULL else as.character(grp)
  if (isTRUE(include.summary.row)) {
    if (is.null(grp)) {
      if (nrow(df) > 1)
        disp <- rbind(disp, .tableFooterRow(df, agg_cols, prec_of, errfun, pm, err_stat,
                                            names(df), meta$label.col))
    } else {
      out <- disp[0, , drop = FALSE]; tag <- character(0)
      for (lv in levels(grp)) {
        idx <- which(grp == lv)
        if (!length(idx)) next
        block <- disp[idx, , drop = FALSE]
        n_new <- length(idx)
        if (length(idx) > 1) {
          f <- .tableFooterRow(df[idx, , drop = FALSE], agg_cols, prec_of, errfun, pm, err_stat,
                               names(df), meta$label.col)
          # the footer carries its own group value, so an exported table stays self-describing
          if (!identical(group.by, meta$label.col)) f[[group.by]] <- lv
          block <- rbind(block, f); n_new <- n_new + 1L
        }
        out <- rbind(out, block); tag <- c(tag, rep(lv, n_new))
      }
      disp <- out; row_grp <- tag
    }
  }
  rownames(disp) <- NULL
  disp[is.na(disp)] <- "-"
  if (style != "internal") names(disp) <- .tableHeaders(names(disp), meta, style)
  if (identical(symbols, "ascii")) names(disp) <- .foldTableSymbols(names(disp))
  # Which group each OUTPUT row belongs to, for print() to break on. An attribute rather than blank
  # separator rows: format() is the export route, and an empty record reads as a malformed row.
  if (!is.null(row_grp)) attr(disp, "table.groups") <- row_grp
  disp
}


#' Print a moby summary table
#'
#' @description Renders the formatted table (see \code{\link[=format.mobyTable]{format}}) beneath a
#' one-line banner, with grouped tables broken by a blank line and a group heading. The object itself
#' stays typed - this affects only what is shown.
#' @param x A `mobyTable`.
#' @param ... Passed to \code{\link[=format.mobyTable]{format}}, so the display can be tuned in
#'   place, e.g. `print(x, style = "report", decimals = c(distance_km = 2))`.
#' @return `x`, invisibly.
#' @exportS3Method print mobyTable
print.mobyTable <- function(x, ...) {
  meta <- .tableMeta(x)
  df <- as.data.frame(x)
  if (nrow(df) == 0) { cat("<mobyTable> 0 rows\n"); return(invisible(x)) }
  # The console IS the right place to ask what the terminal can render - and to ask it once, for the
  # banner and the table together.
  uni <- cli::is_utf8_output()
  pm <- if (uni) "\u00b1" else "+/-"
  dots <- list(...)
  fmt <- format(x, symbols = if (uni) "unicode" else "ascii", ...)
  rg <- attr(fmt, "table.groups")

  banner <- if (!is.null(rg) && length(unique(rg)) > 1)
    sprintf(" (%d groups; one mean %s %s row per group)", length(unique(rg)), pm, meta$error.stat)
  else if (nrow(df) > 1) sprintf(" (final row: mean %s %s)", pm, meta$error.stat) else ""
  cat(sprintf("<mobyTable: %s> %d row%s%s\n", meta$kind, nrow(df),
              if (nrow(df) != 1) "s" else "", banner))

  if (!is.null(rg) && length(rg) == nrow(fmt) && length(unique(rg)) > 1L) {
    # one block per group, each under its own heading: the group name is a label, not a data row
    g <- .mobyGlyphs()
    # the grouping column is dropped from the blocks: the heading already carries it, and repeating
    # it on every row costs the width that makes a wide table wrap. format() keeps it, because an
    # exported table has no headings to carry it.
    gcol <- if (is.null(dots$group.by) || isTRUE(is.character(dots$group.by))) {
      nm <- if (is.null(dots$group.by)) "group" else dots$group.by
      hdr <- .tableHeaders(nm, meta, "internal")
      intersect(unique(c(nm, hdr)), names(fmt))
    } else character(0)
    body <- if (length(gcol)) fmt[, setdiff(names(fmt), gcol), drop = FALSE] else fmt
    for (lv in unique(rg)) {
      cat("\n", g$rule, " ", lv, "\n", sep = "")
      print(body[rg == lv, , drop = FALSE], row.names = FALSE)
    }
  } else {
    print(fmt, row.names = FALSE)
  }
  invisible(x)
}
