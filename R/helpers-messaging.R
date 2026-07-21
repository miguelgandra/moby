#######################################################################################################
## Internal helpers: console messaging ################################################################
## A single chokepoint for user-facing status, warnings, errors, progress bars and runtime reporting.
## Routing everything through these wrappers keeps tone/behaviour consistent (status -> stderr and
## suppressible; every condition hides the call frame) and means a future migration to a richer
## back-end (e.g. cli) is a change to this one file rather than a package-wide sweep.
##
## Conventions:
##   * .mobyInform() - process status / progress notes. Goes to stderr via message() (so it honours
##     suppressMessages()), and is gated on `verbose`. NEVER use cat() for status: cat() writes to
##     stdout, cannot be silenced, and contaminates captured output in reproducible pipelines. Reserve
##     cat() for S3 print.* methods and the plot summary cards, which display an object.
##   * .mobyWarn() / .mobyAbort() - conditions, always with call. = FALSE (users never see the frame).
##   * verbose defaults read getOption("moby.verbose", TRUE), so setOption(moby.verbose = FALSE) is a
##     single global silence switch for the whole package.
#######################################################################################################


#' Emit a status/progress message (verbose-gated, to stderr)
#'
#' @description Package-wide wrapper for user-facing status output. Concatenates its arguments like
#' \code{message()} and prints them to stderr (so \code{suppressMessages()} works) only when
#' \code{verbose} is TRUE.
#' @param ... Objects pasted into the message (as in \code{message()}).
#' @param verbose Logical; print the message. Defaults to \code{getOption("moby.verbose", TRUE)}.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.mobyInform <- function(..., verbose = getOption("moby.verbose", TRUE)) {
  if (isTRUE(verbose)) message(...)
  invisible(NULL)
}


#' Emit a warning (always hiding the call frame)
#'
#' @description Package-wide wrapper for \code{warning(..., call. = FALSE)}, so every moby warning
#' follows the house rule of not exposing the internal call. Warnings are not gated on \code{verbose}
#' (a warning always signals something the user should see).
#' @param ... Objects pasted into the warning message.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.mobyWarn <- function(...) {
  warning(..., call. = FALSE)
  invisible(NULL)
}


#' Raise an error (always hiding the call frame)
#'
#' @description Package-wide wrapper for \code{stop(..., call. = FALSE)}, so every moby error follows
#' the house rule of not exposing the internal call frame.
#' @param ... Objects pasted into the error message.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.mobyAbort <- function(...) {
  stop(..., call. = FALSE)
}


#' Create a progress bar (verbose- and interactive-gated)
#'
#' @description Returns a \code{utils::txtProgressBar} only when \code{verbose} is TRUE and the session
#' is interactive, otherwise \code{NULL}. Gating on \code{interactive()} keeps progress bars out of log
#' files, knitted documents and non-interactive CRAN/batch runs. Pair with \code{.progressSet()} /
#' \code{.progressEnd()}, which no-op on a \code{NULL} bar so call sites need no guards.
#' @param max Maximum value of the bar.
#' @param verbose Logical; whether progress is wanted. Defaults to \code{getOption("moby.verbose", TRUE)}.
#' @param min Minimum value of the bar (default 0).
#' @return A \code{txtProgressBar}, or \code{NULL}.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.progressBar <- function(max, verbose = getOption("moby.verbose", TRUE), min = 0) {
  if (!isTRUE(verbose) || !interactive() || !is.finite(max) || max <= min) return(NULL)
  utils::txtProgressBar(min = min, max = max, initial = min, style = 3)
}

#' Advance a progress bar (no-op on NULL)
#' @keywords internal
#' @noRd
.progressSet <- function(pb, value) {
  if (!is.null(pb)) utils::setTxtProgressBar(pb, value)
  invisible(NULL)
}

#' Close a progress bar (no-op on NULL)
#' @keywords internal
#' @noRd
.progressEnd <- function(pb) {
  if (!is.null(pb)) { close(pb); .mobyInform("") }  # newline after the bar
  invisible(NULL)
}


#' Report elapsed run time (verbose-gated, to stderr)
#'
#' @description Prints a single "Total execution time: X units" line via \code{.mobyInform()},
#' replacing the block of hand-duplicated timing \code{cat()} lines scattered across the compute
#' functions.
#' @param start.time A \code{POSIXct} start stamp (typically \code{Sys.time()} captured on entry).
#' @param verbose Logical; print the line. Defaults to \code{getOption("moby.verbose", TRUE)}.
#' @param prefix Optional leading label (default "Total execution time:").
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.reportRuntime <- function(start.time, verbose = getOption("moby.verbose", TRUE),
                           prefix = "Total execution time:") {
  if (!isTRUE(verbose)) return(invisible(NULL))
  elapsed <- difftime(Sys.time(), start.time)
  .mobyInform(sprintf("%s %.2f %s", prefix, as.numeric(elapsed), attr(elapsed, "units")), verbose = TRUE)
  invisible(NULL)
}
