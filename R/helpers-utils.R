#######################################################################################################
## Internal helpers: utilities #######################################################################
## Console output, datetime/timezone, data-frame manipulation, rescaling and numeric formatting.
#######################################################################################################

##################################################################################################
## Import namespaces   ###########################################################################

#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
NULL


##################################################################################################
## Print to console   ############################################################################

#' Print to console
#'
#' @description Prints a string to the console with a specific formatting.
#' @param string A character string to be printed to the console.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.printConsole <- function(string, verbose = getOption("moby.verbose", TRUE)){
  if (!isTRUE(verbose)) return(invisible(NULL))
  wrapped_text <- paste(strwrap(string, width=getOption("width")*1.2), collapse="\n")
  # Routed through message() (stderr) so section headers are suppressible and never contaminate
  # captured stdout. ANSI bold is emitted only on interactive colour-capable terminals, so log
  # files and non-ANSI front-ends stay clean.
  if (interactive() && isTRUE(getOption("crayon.enabled", interactive()))) {
    message("\033[0;1m", wrapped_text, "\033[0m")
  } else {
    message(wrapped_text)
  }
  invisible(NULL)
}


##################################################################################################
## Plot summary cards   ##########################################################################

#' Print a plot summary "card" to the console
#'
#' @description Shared formatter for the summary blocks the plotting functions print after rendering
#' (e.g. "Abacus plot", "Movement network"). Centralises the rule width and key-column padding so
#' every card looks identical, replacing the per-function \code{hr <- strrep("-", NN)} /
#' \code{kv <- function(k, v) cat(sprintf("  %-Ns ...))} definitions that had drifted (rule widths of
#' 50/52/54 and key paddings of 11/12/14/16). These are object-description displays, so \code{cat()}
#' (stdout) is intentional here - unlike status output, which goes through \code{.mobyInform()}.
#'
#' \code{.summaryOpen()} prints a blank line, the title, and the top rule; \code{.summaryClose()}
#' prints the bottom rule; \code{.kv()} prints one aligned \code{key: value} row; \code{.kvCont()}
#' prints a continuation row (blank key column, aligned under the values).
#'
#' @param title Card title (e.g. "Abacus plot").
#' @param key,value A key/value pair (the colon is added automatically).
#' @note These functions are intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.summaryRule  <- function() strrep("-", 54L)               # single, consistent card width
.summaryOpen  <- function(title) cat("\n", title, "\n", .summaryRule(), "\n", sep = "")
.summaryClose <- function() cat(.summaryRule(), "\n", sep = "")
.kv           <- function(key, value) cat(sprintf("  %-16s %s\n", paste0(key, ":"), value))
.kvCont       <- function(value) cat(sprintf("  %-16s %s\n", "", value))


##################################################################################################
## Timezone helper   #############################################################################

#' Retrieve the timezone of a date-time vector
#'
#' @description Returns the timezone (`tzone` attribute) of a POSIXct vector, falling back
#' to "UTC" when the attribute is unset or empty. Used throughout the package to preserve
#' the timezone of user-supplied timestamps when round-tripping through numeric/epoch
#' representations (e.g. after \code{tapply} or \code{aggregate}), instead of silently
#' relabelling them as UTC.
#' @param x A POSIXct vector (or any object carrying a `tzone` attribute).
#' @return A single character string with the timezone name.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.dataTZ <- function(x){
  tz <- attr(x, "tzone")
  if (is.null(tz) || length(tz) == 0 || tz[1] == "") "UTC" else tz[1]
}


##################################################################################################
## Safe column drop   ############################################################################

#' Drop columns from a data frame by name (safely)
#'
#' @description Removes the named columns from a data frame. Unlike the
#' \code{df[, -which(colnames(df) \%in\% cols)]} idiom, this helper returns the data frame
#' unchanged when none of the requested columns are present (rather than silently dropping
#' every column, which is what negative indexing with \code{integer(0)} does).
#' @param df A data frame.
#' @param cols A character vector of column names to remove.
#' @return The data frame without the specified columns (always with `drop = FALSE`).
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.dropCols <- function(df, cols){
  keep <- setdiff(colnames(df), cols)
  df[, keep, drop = FALSE]
}


##################################################################################################
## Lag / lead   ##################################################################################

#' Shift a vector by n positions (base-R replacements for dplyr::lag / dplyr::lead)
#'
#' @description Shift `x` forward (`.lag`) or backward (`.lead`) by `n`, padding the
#'   vacated positions with `NA`. Uses `[`-indexing (not `c(NA, head(x, -n))`) so the
#'   class and attributes of `x` are preserved exactly — in particular the `tzone` of a
#'   POSIXct vector, which matters throughout moby's time-interval calculations.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.lag <- function(x, n = 1L) {
  len <- length(x)
  x[c(rep(NA_integer_, min(n, len)), seq_len(max(0L, len - n)))]
}

.lead <- function(x, n = 1L) {
  len <- length(x)
  x[c(seq_len(max(0L, len - n)) + n, rep(NA_integer_, min(n, len)))]
}


##################################################################################################
## Standard error / axis break  ##################################################################
## Base-R replacements for plotrix::std.error() and plotrix::axis.break()

#' Standard error of the mean
#'
#' @description Base-R replacement for `plotrix::std.error()`: the sample standard
#'   deviation divided by the square root of the number of non-missing values. Matches
#'   plotrix on a plain numeric vector (the only form used within moby).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.stdError <- function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

#' Draw a slash-style axis break on the x-axis
#'
#' @description Reproduces the `axis = 1`, `style = "slash"` branch of
#'   `plotrix::axis.break()`: a background-coloured masking rectangle over the axis at
#'   `breakpos`, overdrawn with two parallel slashes. Adapted from
#'   `plotrix::axis.break()` (plotrix, GPL (>= 2); authors Jim Lemon and Duncan Murdoch).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.axisBreak <- function(breakpos, brw = 0.02, bgcol = "white", breakcol = "black") {
  usr <- graphics::par("usr")
  xw <- (usr[2] - usr[1]) * brw
  yw <- (usr[4] - usr[3]) * brw
  old.xpd <- graphics::par("xpd"); on.exit(graphics::par(xpd = old.xpd))
  graphics::par(xpd = TRUE)
  graphics::rect(breakpos - xw / 2, usr[3] - yw / 2,
                 breakpos + xw / 2, usr[3] + yw / 2, col = bgcol, border = bgcol)
  graphics::segments(c(breakpos - xw, breakpos), c(usr[3] - yw / 2, usr[3] - yw / 2),
                     c(breakpos, breakpos + xw), c(usr[3] + yw / 2, usr[3] + yw / 2),
                     col = breakcol, lty = 1)
}


##################################################################################################
## NA fill: last-obs-carried-forward / linear interpolation  #####################################
## Base-R replacements for zoo::na.locf() and zoo::na.approx()

#' Last observation carried forward
#'
#' @description Base-R replacement for `zoo::na.locf(x, na.rm = FALSE)`: each `NA` takes
#'   the value of the most recent non-`NA` element; leading `NA`s (with no prior value to
#'   carry) are kept as `NA`, so length is always preserved. Works for numeric, character
#'   and factor input (`x[idx]` preserves type and levels).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.naLocf <- function(x) {
  idx <- seq_along(x)
  idx[is.na(x)] <- 0L
  idx <- cummax(idx)
  idx[idx == 0L] <- NA_integer_          # unfillable leading NAs stay NA
  x[idx]
}

#' Linear interpolation of interior NAs
#'
#' @description Base-R replacement for `zoo::na.approx(v, na.rm = FALSE, rule = 2)`:
#'   interior `NA`s are filled by linear interpolation over the integer index; with
#'   `rule = 2`, leading/trailing `NA`s take the nearest end value. With fewer than two
#'   non-missing anchors the input is returned unchanged.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.naApprox <- function(v) {
  na <- is.na(v)
  if (!any(na) || sum(!na) < 2L) return(v)
  i <- seq_along(v)
  v[na] <- stats::approx(x = i[!na], y = v[!na], xout = i[na], rule = 2)$y
  v
}


##################################################################################################
## Rescale function  #############################################################################
## Sourced from 'scales' package (https://CRAN.R-project.org/package=scales)

#' Rescale
#'
#' @description Rescales a numeric vector to a specified range.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.rescale <- function (x, to=c(0, 1), from=range(x, na.rm=TRUE, finite=TRUE), ...) {
  if (.zero_range(from) || .zero_range(to)) return(ifelse(is.na(x), NA, mean(to)))
  (x - from[1])/diff(from) * diff(to) + to[1]
}

.zero_range <- function(x, tol=1000*.Machine$double.eps) {
  if (length(x) == 1) return(TRUE)
  if (length(x) != 2) stop("'x' must be length 1 or 2", call.=FALSE)
  if (any(is.na(x))) return(NA)
  if (x[1] == x[2]) return(TRUE)
  if (all(is.infinite(x))) return(FALSE)
  m <- min(abs(x))
  if (m == 0) return(FALSE)
  abs((x[1] - x[2]) / m) < tol
}



##################################################################################################
## Decimal Places   ##############################################################################
## Sourced from https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r

#' Decimal Places
#'
#' @description Determines the number of decimal places in a numeric value.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.decimalPlaces <- function(x) {
  if(is.na(x)){return(NA)}
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', format(x, scientific=FALSE)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

# vectorised form (handles vectors of values element-wise)
.decimalPlaces <- Vectorize(.decimalPlaces)


##################################################################################################
## Time-bin width (seconds) ######################################################################

#' Resolve a time-bin specification to a width in seconds
#'
#' @description Converts a bin specification - a unit string (`"hour"`, `"day"`, `"week"`,
#' `"month"`), a number of minutes, or `"auto"` - into a bin width in seconds. `"auto"` targets
#' roughly 150 bins across the temporal span, snapped to the nearest natural unit.
#' @param bin A unit string, a numeric number of minutes, or `"auto"`.
#' @param times A POSIXct vector (used to size the `"auto"` bin from the data span).
#' @return A single numeric bin width, in seconds.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.binSeconds <- function(bin, times) {
  units_sec <- c(hour = 3600, day = 86400, week = 604800, month = 2629800)
  if (is.numeric(bin)) return(bin * 60)                         # minutes -> seconds
  if (identical(bin, "auto")) {
    span_sec <- as.numeric(diff(range(times, na.rm = TRUE)), units = "secs")
    raw <- span_sec / 150
    return(unname(units_sec[which.min(abs(log(units_sec / raw)))]))
  }
  bin <- match.arg(bin, names(units_sec))
  unname(units_sec[bin])
}


##################################################################################################
## ID-group comparison series ####################################################################

#' Build the list of ID series for group comparisons
#'
#' @description Given a named list of ID groups and a comparison mode, returns the ordered set of
#' "series" to analyse - one entry per within-group ("all"/"within") and/or per between-group pair
#' ("all"/"between"). Each element is the vector of IDs in that series; pair series are named
#' `"A <-> B"` (the `<->` marker lets callers detect a between-group comparison). Shared by the
#' co-occurrence-aware plotting functions so the within/between/all combinatorics live in one place.
#' @param id.groups A named list of ID groups.
#' @param group.comparisons One of "within", "between" or "all".
#' @return A named list of ID vectors (names are group names and/or `"A <-> B"` pair labels).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.buildGroupComparisonList <- function(id.groups, group.comparisons) {
  ids_table <- .meltList(id.groups)
  colnames(ids_table) <- c("ID", "group")
  out <- list()
  if (group.comparisons %in% c("within", "all")) {
    for (g in names(id.groups)) out[[g]] <- id.groups[[g]]
  }
  if (group.comparisons %in% c("between", "all")) {
    pairs <- utils::combn(names(id.groups), 2, simplify = FALSE)
    for (p in pairs) out[[paste(p, collapse = " <-> ")]] <- ids_table$ID[ids_table$group %in% p]
  }
  out
}


##################################################################################################
## Co-occurrence within a single time-bin ########################################################

#' Detect co-occurrences within one time-bin row
#'
#' @description Given one row of a wide `time-bin x individual` table (values = the station each
#' individual was at, `NA` if absent) and the group label of each individual, returns the co-occurring
#' clusters - stations where 2+ individuals were detected together. `comparison` selects which clusters
#' count: `"within"` keeps clusters containing 2+ individuals of the *same* group; `"between"` keeps
#' clusters spanning *2+* groups (checked **per cluster** - the authoritative definition, shared by
#' the co-occurrence plots). `return` selects the output: the station name of each valid cluster
#' (`"stations"`, for per-location tallies) or its size (`"sizes"`, for group-size distributions).
#' @param row A character/atomic vector: one wide-table row (station per individual, `NA` if absent).
#' @param comparison One of `"within"` or `"between"`.
#' @param groups Group label of each individual (same length/order as `row`).
#' @param return `"stations"` (default) or `"sizes"`.
#' @return A vector of station names or cluster sizes for the valid clusters, or `NULL` if none.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.countCooccurrences <- function(row, comparison, groups, return = c("stations", "sizes")) {
  return <- match.arg(return)
  row <- as.character(row); names(row) <- groups
  counts <- table(row)
  if (!any(counts > 1)) return(NULL)
  clusters <- names(counts)[counts > 1]
  members  <- sapply(clusters, function(x) names(row)[which(row == x)], simplify = FALSE)
  valid <- if (comparison == "between") vapply(members, function(x) length(unique(x)) > 1, logical(1))
           else                          vapply(members, function(x) any(table(x) > 1), logical(1))
  if (return == "sizes") as.numeric(counts[clusters][valid]) else clusters[valid]
}
