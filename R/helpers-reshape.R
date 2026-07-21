##################################################################################################
## Reshape / join helpers  #######################################################################
##################################################################################################
## Base-R replacements for the reshape2 and plyr functions moby relied on, so that neither
## package (both author-retired: reshape2 -> tidyr, plyr -> dplyr) is required at install time.
## Each helper reproduces the exact behaviour moby depended on; the equivalence with the original
## functions is pinned in tests/testthat/test-helpers-reshape.R (in particular the .castWide
## fill / drop = FALSE / factor-level-preservation grid, on which the association, UD and network
## matrices depend numerically).


#' Left join preserving the row order and column layout of `x`
#'
#' @description Base-R replacement for `.joinKeep(x, y, by, type = "left")`. Unlike
#'   `merge()`, which does not guarantee the order of its output, this keeps every row of
#'   `x` in its original order and appends `y`'s non-key columns after `x`'s columns.
#' @param x,y Data frames to join.
#' @param by Character vector of key column name(s) common to `x` and `y`.
#' @param type Only `"left"` is supported (kept for call-site compatibility).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.joinKeep <- function(x, y, by, type = "left") {
  if (!identical(type, "left")) stop(".joinKeep only supports type = 'left'.", call. = FALSE)
  ycols <- setdiff(names(y), by)
  keyx <- if (length(by) == 1L) x[[by]] else do.call(paste, c(x[by], list(sep = "\r")))
  keyy <- if (length(by) == 1L) y[[by]] else do.call(paste, c(y[by], list(sep = "\r")))
  idx  <- match(as.character(keyx), as.character(keyy))     # first match; moby's y tables have unique keys
  # Keep x exactly (its columns, row order, and any list-columns — shadeSeasons feeds a list-column from
  # aggregate(simplify = FALSE)), then append y's non-key columns by POSITION. Positional append both
  # allows duplicate names (e.g. both frames carry an aggregate column "x"; the call sites rename the
  # result positionally) and preserves list-columns, which as.data.frame(as.list(x)) would expand. This
  # matches plyr::join, which returned a plain data.frame.
  out <- x
  class(out) <- "data.frame"
  attr(out, "moby") <- NULL
  for (col in ycols) out[[length(out) + 1L]] <- y[[col]][idx]
  names(out) <- c(names(x), ycols)
  rownames(out) <- NULL
  out
}


#' Row-bind data frames, filling absent columns with NA
#'
#' @description Base-R replacement for `.rbindFill()`. Accepts either several
#'   data-frame arguments or a single list of them. `NULL` inputs are dropped (returning
#'   `NULL` only if every input is `NULL`); empty (0-row) frames are retained so their
#'   columns still contribute to the union, matching `.rbindFill()`.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.rbindFill <- function(...) {
  dfs <- list(...)
  if (length(dfs) == 1L && is.list(dfs[[1]]) && !is.data.frame(dfs[[1]])) dfs <- dfs[[1]]
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (!length(dfs)) return(NULL)
  cols <- unique(unlist(lapply(dfs, names)))
  # Rebuild each frame with the full column union: existing columns kept, absent ones filled with an
  # NA vector of the correct length (`rep(NA, nrow(d))` — works for 0-row frames, where `d[miss] <- NA`
  # would error). Rebuilding via as.data.frame(list) also yields a clean, plain data.frame, stripping
  # any subclass such as 'mobyData' (whose `[` method would otherwise re-attach metadata), exactly as
  # plyr::rbind.fill returned a plain data.frame.
  filled <- lapply(dfs, function(d) {
    n <- nrow(d)
    cells <- lapply(cols, function(cc) if (cc %in% names(d)) d[[cc]] else rep(NA, n))
    as.data.frame(stats::setNames(cells, cols), check.names = FALSE, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, filled)
  rownames(out) <- NULL
  out
}


#' Replace values found in `from` with the matching `to`
#'
#' @description Base-R replacement for `.mapValues(x, from, to, warn_missing = FALSE)`.
#'   Values of `x` present in `from` are replaced by the corresponding element of `to`;
#'   unmatched values are left unchanged. Returns a character vector (as the moby call sites
#'   expect). `...` absorbs the `warn_missing` argument carried over from the plyr calls.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.mapValues <- function(x, from, to, ...) {
  i <- match(x, from)
  x <- as.character(x)
  hit <- !is.na(i)
  x[hit] <- as.character(to)[i[hit]]
  x
}


#' Melt a named list into a long data frame
#'
#' @description Base-R replacement for `reshape2::melt()` applied to a named list: returns
#'   a data frame with columns `value` (the concatenated elements) and `L1` (the list name
#'   each element came from).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.meltList <- function(x) {
  data.frame(value = unlist(x, use.names = FALSE),
             L1 = rep(names(x), lengths(x)),
             stringsAsFactors = FALSE)
}


#' Spread long data into a wide (row x col) table
#'
#' @description Base-R replacement for
#'   `reshape2::dcast(data, row.col ~ col.col, value.var, fun.aggregate, fill, drop = FALSE)`.
#'   The value column is aggregated within each (row, col) cell by `fun.aggregate`, and only
#'   structurally-empty cells (no underlying rows) are set to `fill` — a cell whose aggregator
#'   legitimately returns `NA` keeps that `NA`, exactly as reshape2 does. Factor row/column
#'   variables keep ALL their levels (the padded square grid the network/UD code relies on);
#'   non-factor variables use their sorted unique observed values and preserve class and
#'   attributes (e.g. a POSIXct time-bin column keeps its class and `tzone`).
#' @details `fun.aggregate` must return a single value per cell (every moby aggregator does:
#'   `sum`, `max`, and the custom most-frequent-value / tie helpers). A vectorised aggregator
#'   would make `tapply()` return a list-matrix and fail loudly rather than silently.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.castWide <- function(data, row.col, col.col, value.col, fun.aggregate, fill = NA) {
  rv <- data[[row.col]]
  cv <- data[[col.col]]

  # Column keys: all levels for a factor (padded grid, drop = FALSE), else sorted unique observed.
  ckeys <- if (is.factor(cv)) levels(cv) else sort(unique(cv))
  ccode <- if (is.factor(cv)) as.integer(cv) else match(cv, ckeys)
  cf <- factor(ccode, levels = seq_along(ckeys))

  # Row keys, same rule. Integer-code via match() rather than factor(rv, levels = ...): passing a
  # POSIXct/numeric key column as factor levels silently mis-codes (factor coerces x but not levels),
  # so the class-preserving first column is kept separately in `firstcol`.
  if (is.factor(rv)) {
    rkeys <- levels(rv); rcode <- as.integer(rv); firstcol <- factor(rkeys, levels = rkeys)
  } else {
    rkeys <- sort(unique(rv)); rcode <- match(rv, rkeys); firstcol <- rkeys
  }
  rf <- factor(rcode, levels = seq_along(rkeys))

  # Enforce a single value per cell so `m` is always an atomic matrix. A zero-length return (e.g. the
  # most-frequent aggregator names(which.max(table(x))) on a cell whose values are all NA) becomes NA
  # rather than a silent list-column NULL; a genuine multi-value return fails loudly.
  agg1 <- function(z) {
    r <- fun.aggregate(z)
    if (length(r) == 1L) r
    else if (length(r) == 0L) NA
    else stop(".castWide: fun.aggregate must return a single value per cell.", call. = FALSE)
  }
  m   <- tapply(data[[value.col]], list(rf, cf), FUN = agg1)
  cnt <- tapply(data[[value.col]], list(rf, cf), FUN = length)
  m[is.na(cnt)] <- fill                                   # fill only structurally-empty cells
  colnames(m) <- as.character(ckeys)

  out <- data.frame(firstcol, m, check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)
  names(out)[1] <- row.col
  out
}
