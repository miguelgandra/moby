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
.progressBar <- function(max, verbose = getOption("moby.verbose", TRUE), min = 0,
                         name = "Processing", .envir = parent.frame()) {
  if (!isTRUE(verbose) || !interactive() || !is.finite(max) || max <= min) return(NULL)
  # Separated from whatever preceded it: a bar drawn flush against the header's criteria bullets reads
  # as part of that block rather than as the next step of the run. Emitted here rather than at the
  # call sites so every bar in the package is spaced the same way, and only when a bar is actually
  # created (a quiet or non-interactive run returns above, and leaves no stray blank line).
  .mobyBlank(verbose)
  # cli bar: tied to the caller's frame (auto-cleaned on early exit) and auto-suppressed for loops
  # that finish faster than getOption("cli.progress_show_after") (~2 s), so quick runs stay clean.
  cli::cli_progress_bar(name, total = as.integer(max - min), .envir = .envir)
}

#' Advance a progress bar (no-op on NULL)
#' @keywords internal
#' @noRd
.progressSet <- function(pb, value) {
  if (!is.null(pb)) cli::cli_progress_update(id = pb, set = value)
  invisible(NULL)
}

#' Close a progress bar (no-op on NULL)
#' @keywords internal
#' @noRd
.progressEnd <- function(pb) {
  if (!is.null(pb)) cli::cli_progress_done(id = pb)   # cli clears the bar itself; no manual newline
  invisible(NULL)
}


#######################################################################################################
## Structured console framework (cli backend) #########################################################
## The vocabulary shared by every user-facing moby function, so output looks the same wherever the
## user lands. Built on cli for consistent symbols, alignment and automatic ASCII fallback on consoles
## that cannot render Unicode (the previous hand-written bullet byte-escaped to "<U+2022>" there).
##
## Verbosity is a single logical, honoured package-wide via getOption("moby.verbose", TRUE). There is
## deliberately no second "detailed" level: moby has no per-item diagnostic stream that would need one,
## and the material that looks like detail (per-stage removals, per-station QC rows) already has a
## better home in the returned object's print method.
##
## Tiers (see ?moby_pipeline for which function is which):
##   Silent    no output; accessors, predicates, inline classifiers
##   Summary   .mobyHeader() + a short outcome
##   Progress  .mobyHeader() + a progress bar + an outcome (+ runtime when the run is long)
##
## Nesting: when one user-facing function calls another, the PARENT passes verbose = FALSE to the
## child. The header depth guard below is a safety net for the case where someone forgets, not the
## primary mechanism.
#######################################################################################################

# internal state: the frame that opened the current header (nesting guard). Storing the FRAME rather
# than a flag makes the guard self-expiring: once that function returns, its environment is no longer
# among sys.frames(), so no on.exit() bookkeeping is needed (and none can be forgotten).
.moby_state <- new.env(parent = emptyenv())
.moby_state$opener <- NULL

#' Is a moby header already open further up the call stack?
#' @keywords internal
#' @noRd
.mobyHeaderOpen <- function() {
  op <- .moby_state$opener
  if (is.null(op)) return(FALSE)
  any(vapply(sys.frames(), identical, logical(1), op))
}

#' The standard function header
#'
#' @description Opens a user-facing function's output with the package's shared banner: a rule naming
#' the function, an `i` line saying what is being done to the data, an optional `Input:` scale line, and
#' an optional labelled block of the methodological choices that define the run.
#'
#' Only choices that change what the analysis MEANS belong in `criteria` (a bandwidth, a speed limit, a
#' threshold) - never cosmetic arguments, and never values already recoverable from the returned
#' object, which the print method reports.
#'
#' @param title Function name, e.g. `"filterDetections()"`.
#' @param intro One line: what the function is doing with the user's data.
#' @param input Optional one-line scale of the input (e.g. `"1,643 detections | 8 individuals"`).
#' @param facts Optional named character vector of facts about the data being handled (counts, period,
#' what was dropped). Rendered as aligned bullets directly under `intro`, with no heading - use this
#' for "what am I looking at", and `criteria` for "how is it being analysed".
#' @param criteria Optional named character vector of methodological choices; names are padded so the
#' values align.
#' @param criteria.label Heading for the `criteria` block. Defaults to `"Method"`.
#' @param verbose Logical; print. Defaults to \code{getOption("moby.verbose", TRUE)}.
#' @param .envir Frame used to register the guard reset; defaults to the caller.
#' @keywords internal
#' @noRd
.mobyHeader <- function(title, intro, input = NULL, facts = NULL, criteria = NULL,
                        criteria.label = "Method",
                        verbose = getOption("moby.verbose", TRUE), .envir = parent.frame()) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  # A header is already open further up the stack: stay quiet rather than nest banners. This is a
  # SAFETY NET only - the primary mechanism is that a parent passes verbose = FALSE to child calls.
  if (.mobyHeaderOpen()) return(invisible(NULL))
  .moby_state$opener <- .envir

  cli::cli_rule(left = "{.strong {title}}", right = "moby")
  cli::cli_text("")
  cli::cli_alert_info("{intro}")
  if (!is.null(input)) cli::cli_text("{cli::symbol$bullet} Input: {input}")
  # cli_verbatim: emitted exactly as given, so the padding that aligns the values survives
  # (cli_text/cli_bullets run the string through glue and collapse repeated spaces).
  .kvBlock <- function(x, indent) {
    w <- max(nchar(names(x)))
    for (nm in names(x))
      cli::cli_verbatim(paste0(indent, cli::symbol$bullet, " ", formatC(paste0(nm, ":"), width = -(w + 1)),
                               " ", x[[nm]]))
  }
  if (length(facts) > 0) .kvBlock(facts, "")
  if (length(criteria) > 0) .mobyBlock(criteria.label, criteria, verbose = TRUE)
  invisible(NULL)
}

#' A labelled block of aligned bullets ("-> Label" followed by indented entries).
#'
#' The shape the header uses for its criteria, available on its own so a completion summary can group
#' related lines the same way instead of emitting a loose run of bullets. An entry whose value is
#' empty prints as a bare name, which is how a filter with no tunable parameter is listed.
#' @keywords internal
#' @noRd
.mobyBlock <- function(label, items, verbose = getOption("moby.verbose", TRUE)) {
  if (!isTRUE(verbose) || length(items) == 0) return(invisible(NULL))
  nms <- names(items)
  vals <- as.character(unlist(items))
  has_val <- !is.na(vals) & nzchar(vals)
  cli::cli_text("")
  cli::cli_text("{cli::symbol$arrow_right} {label}")
  # width from the named entries only, so a bare name does not pad the others
  w <- if (any(has_val)) max(nchar(nms[has_val])) else 0L
  for (i in seq_along(items)) {
    # cli_verbatim: emitted exactly as given, so the padding that aligns the values survives
    cli::cli_verbatim(if (has_val[i])
      paste0("  ", cli::symbol$bullet, " ", formatC(nms[i], width = -w), "  ", vals[i])
      else paste0("  ", cli::symbol$bullet, " ", nms[i]))
  }
  invisible(NULL)
}

#' A blank separator before an outcome block
#' @keywords internal
#' @noRd
.mobyBlank <- function(verbose = getOption("moby.verbose", TRUE)) {
  if (isTRUE(verbose)) cli::cli_text("")
  invisible(NULL)
}

#' A success / outcome line ("this is what came out")
#' @keywords internal
#' @noRd
.mobyOk <- function(..., verbose = getOption("moby.verbose", TRUE)) {
  if (isTRUE(verbose)) cli::cli_alert_success("{paste0(...)}")
  invisible(NULL)
}

#' An informational line ("context you may want")
#' @keywords internal
#' @noRd
.mobyNote <- function(..., verbose = getOption("moby.verbose", TRUE)) {
  if (isTRUE(verbose)) cli::cli_alert_info("{paste0(...)}")
  invisible(NULL)
}

#' An attention line ("worth a second look") - NOT a warning; the run succeeded
#'
#' @description Reserved for outcomes the user should actually look at (records held back for review,
#' individuals lost entirely). A non-zero count on its own is not attention-worthy: an expected removal
#' belongs on the \code{.mobyOk()} line.
#' @keywords internal
#' @noRd
.mobyAttention <- function(..., verbose = getOption("moby.verbose", TRUE)) {
  if (isTRUE(verbose)) cli::cli_alert_warning("{paste0(...)}")
  invisible(NULL)
}

#' A pointer line ("where the output went", "what was attached")
#' @keywords internal
#' @noRd
.mobyArrow <- function(..., verbose = getOption("moby.verbose", TRUE)) {
  if (isTRUE(verbose)) cli::cli_text("{cli::symbol$arrow_right} {paste0(...)}")
  invisible(NULL)
}

#' Format an elapsed time compactly ("0.4s", "3m 12s", "1h 04m")
#' @keywords internal
#' @noRd
.fmtElapsed <- function(start.time) {
  secs <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
  if (!is.finite(secs)) return("?")
  if (secs < 60) return(sprintf("%.1fs", secs))
  m <- floor(secs / 60); s <- round(secs - 60 * m)
  if (m < 60) return(sprintf("%dm %02ds", m, s))
  h <- floor(m / 60); m <- m - 60 * h
  sprintf("%dh %02dm", h, m)
}

#' Runtime footer, for genuinely long-running functions only
#'
#' @description Prints a stopwatch line. Deliberately NOT used by fast functions - a runtime line on a
#' millisecond call is noise. `min.secs` suppresses it when the run turned out to be trivial.
#' @keywords internal
#' @noRd
.mobyRuntime <- function(start.time, verbose = getOption("moby.verbose", TRUE), min.secs = 0) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  secs <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
  if (is.finite(secs) && secs < min.secs) return(invisible(NULL))
  sw <- .mobyGlyph("stopwatch")
  rt <- .fmtElapsed(start.time)
  cli::cli_text("{sw} runtime: {rt}")
  invisible(NULL)
}

#' Locale-aware glyphs for message CONTENT
#'
#' @description cli maps its own symbols to ASCII automatically, but a literal Unicode character placed
#' in message text is byte-escaped on a non-UTF-8 console (a raw "\u00b7" prints as "<U+00B7>"). Route
#' any such character through here so content degrades as gracefully as cli's symbols do.
#' @param name One of "mid" (middle dot separator), "times" (multiplication sign), "stopwatch".
#' @keywords internal
#' @noRd
.mobyGlyph <- function(name) {
  utf8 <- cli::is_utf8_output()
  switch(name,
    mid       = if (utf8) "\u00b7" else "|",
    times     = if (utf8) "\u00d7" else "x",
    stopwatch = if (utf8) "\u23f1" else cli::symbol$bullet,
    dash      = if (utf8) "\u2014" else "-",
    stop("unknown glyph: ", name))
}

#' Format an integer with thousands separators, for message text
#' @keywords internal
#' @noRd
.fmtN <- function(n) {
  # Counts are integers, but this also formats user-supplied thresholds (a max.gap of 0.5 hours, say),
  # and format = "d" TRUNCATES those - reporting "0 hours" for a 30-minute rule. Keep the integer
  # path (thousands separators, no decimals) and fall back to a faithful rendering otherwise.
  # a non-finite value has no integer representation at all (format = "d" yields NA with a coercion
  # warning), so render it as itself
  if (is.numeric(n) && length(n) == 1L && !is.finite(n)) return(as.character(n))
  if (is.numeric(n) && length(n) == 1L && n != round(n))
    return(format(n, big.mark = ",", trim = TRUE, scientific = FALSE))
  formatC(n, format = "d", big.mark = ",")
}

#' Number of individuals OBSERVED in the data (as opposed to DECLARED by the factor levels).
#'
#' moby's ID column is always a factor, so two different quantities can be read off it and they are
#' not interchangeable:
#'   OBSERVED  individuals that actually contributed data  -> this helper
#'   DECLARED  the roster carried by the factor levels     -> nlevels()
#' They diverge whenever `id.groups` names a tagged animal that was never detected (validateArguments
#' re-levels the factor to the roster), whenever a caller subsets without droplevels(), or whenever a
#' pre-factored column arrives with spare levels.
#'
#' Use OBSERVED to answer "how much data is there?" - console input lines, plot summary cards, and any
#' denominator describing the animals that contributed. Use DECLARED where the roster IS the subject:
#' a per-individual results table is deliberately roster-shaped, so that a tagged-but-never-detected
#' animal still gets a row and positionally-supplied tagging.dates stay aligned.
#' @keywords internal
#' @noRd
.nObserved <- function(x) length(unique(as.character(stats::na.omit(x))))

#' Note the tagged individuals that contributed no data, when the roster is larger than the data.
#'
#' The companion to `.nObserved()`. A header's input line reports the OBSERVED count, because that is
#' the plain answer to "how much data is there"; this states the difference separately, so a reader
#' does not have to reconcile two numbers inside one sentence. Silent when they agree, so an ordinary
#' run gains nothing. `detail` says what became of the animal, which differs by function (a
#' roster-shaped table keeps its row; a wide table keeps its column).
#' @keywords internal
#' @noRd
.noteUndetected <- function(x, verbose, detail = NULL) {
  n <- nlevels(x) - .nObserved(x)
  if (n > 0) {
    # separated from the header block: without it the note lands flush under the criteria bullets and
    # reads as a further "Method" entry rather than as a fact about the data
    .mobyBlank(verbose)
    .mobyNote(.fmtCount(n, "tagged individual"), " with no detections",
              if (!is.null(detail)) paste0(" (", detail, ")") else "", verbose = verbose)
  }
  invisible(NULL)
}

#' Detection count for a header's input line.
#'
#' A row is not always a detection: an aggregated dataset (COAs, binned records) carries a numeric
#' `detections` column holding how many decodes each row stands for, so nrow() would report records
#' instead. Shared so that a parent and its child never state different totals for the same input.
#' @keywords internal
#' @noRd
.nDetections <- function(data) {
  if ("detections" %in% colnames(data) && is.numeric(data[["detections"]]))
    sum(data[["detections"]], na.rm = TRUE)
  else nrow(data)
}

#' Format a count with its noun, agreeing in number ("1 individual", "8 individuals").
#'
#' Console text is assembled as plain strings before reaching cli, so cli's own `{?s}` pluralisation
#' (which needs a glue expression) is unavailable - hence this helper. Irregular plurals are passed
#' explicitly via `plural`.
#' @keywords internal
#' @noRd
.fmtCount <- function(n, singular, plural = paste0(singular, "s")) {
  paste0(.fmtN(n), " ", if (isTRUE(n == 1)) singular else plural)
}
