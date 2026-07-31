#######################################################################################################
## Internal helpers: the shared layout of moby's print methods ########################################
##
## Every moby object prints to the same shape, so a user reading one has already learned the others:
##
##   -- <class> ----------------------------------------------------
##   Input      what was analysed
##   Results    the overall outcome
##   <Summary>  a breakdown by category (named for the object: "Metadata", "Filtering", ...)
##   Details    where the complete diagnostics live
##   ---------------------------------------------------------------
##
## Not every object needs every section - a mobyData has nothing to "result" - but the ones it does
## use appear in that order, with the same glyphs, indentation and alignment. These helpers own that
## layout so a new class gets it by calling them rather than by copying a previous method.
##
## print methods legitimately use cat(): they DISPLAY an object rather than report progress, so they
## write to stdout and are not gated on `verbose` (see helpers-messaging.R for the status vocabulary).
#######################################################################################################


#' Glyphs for the print methods, with ASCII fallbacks off UTF-8.
#'
#' Separate from `.mobyGlyph()`, which serves the cli-backed status output: cli picks its own fallbacks
#' from the console's encoding, whereas these are written by cat() and must choose for themselves.
#' @keywords internal
#' @noRd
.mobyGlyphs <- function() {
  # \u escapes keep this file pure ASCII (so it behaves identically whatever locale R starts in)
  if (isTRUE(l10n_info()$`UTF-8`))
    list(rule = "\u2500", sep = "\u00b7", arrow = "\u2192", tick = "\u2713",
         cross = "\u2717", bullet = "\u2022", dash = "\u2014", info = "\u2139")
  else
    list(rule = "-", sep = "|", arrow = "->", tick = "v", cross = "x", bullet = "*",
         dash = "-", info = "i")
}

#' Console width to lay a printed object out in.
#' @keywords internal
#' @noRd
.printWidth <- function(width = getOption("width")) max(40L, min(as.integer(width), 100L))

#' The opening rule: `-- <class> --------------------`.
#' @keywords internal
#' @noRd
.printOpen <- function(class, width = .printWidth()) {
  g <- .mobyGlyphs()
  hdr <- paste0(strrep(g$rule, 2), " <", class, "> ")
  cat(hdr, strrep(g$rule, max(0L, width - nchar(hdr))), "\n", sep = "")
  invisible(NULL)
}

#' The closing rule.
#' @keywords internal
#' @noRd
.printClose <- function(width = .printWidth()) {
  cat(strrep(.mobyGlyphs()$rule, width), "\n", sep = "")
  invisible(NULL)
}

#' One section: a heading, then its entries.
#'
#' `items` is a named character vector: the name is the left-hand label, the value the right-hand
#' text. An entry whose value is empty prints as a bare label, and does not pad the others - so a
#' list of plain statements ("551 deployment records") and a list of label/value pairs ("Time  2010
#' -> 2022") both look right without the caller doing any spacing.
#'
#' `marker` is `"none"`, `"bullet"`, `"tick"`, or a character vector parallel to `items` (which is how
#' a present/absent list mixes ticks and crosses).
#' @keywords internal
#' @noRd
.printSection <- function(title, items, marker = "none", indent = 2L, blank.before = TRUE) {
  if (length(items) == 0) return(invisible(NULL))
  g <- .mobyGlyphs()
  if (blank.before) cat("\n")
  if (!is.null(title)) cat(title, "\n", sep = "")

  nms <- names(items)
  if (is.null(nms)) nms <- rep("", length(items))
  vals <- as.character(unlist(items))
  vals[is.na(vals)] <- ""
  has_val <- nzchar(vals)

  mk <- if (length(marker) == 1L) rep(marker, length(items)) else marker
  mk <- vapply(mk, function(m) switch(m, none = "", bullet = g$bullet, tick = g$tick,
                                      cross = g$cross, m), character(1))

  # align the values on the widest LABEL THAT HAS ONE, so a bare statement cannot stretch the column
  w <- if (any(has_val)) max(nchar(nms[has_val])) else 0L
  pad <- strrep(" ", indent)
  for (i in seq_along(items)) {
    lead <- if (nzchar(mk[i])) paste0(mk[i], " ") else ""
    cat(pad, lead,
        if (has_val[i]) paste0(formatC(nms[i], width = -w), "  ", vals[i]) else nms[i],
        "\n", sep = "")
  }
  invisible(NULL)
}

#' A closing pointer at where the full detail lives.
#' @keywords internal
#' @noRd
.printNote <- function(...) {
  cat("\n", .mobyGlyphs()$info, " ", paste0(...), "\n", sep = "")
  invisible(NULL)
}
