# Fixtures for the OBSERVED-vs-DECLARED individual-count contract.
#
# moby's ID column is always a factor, so two different quantities can be read off it:
#   OBSERVED  individuals actually present in the data   length(unique(as.character(x)))
#   DECLARED  the roster carried by the factor levels    nlevels(x) / levels(x)
#
# Both are legitimate: a per-individual results table is deliberately roster-shaped (so a tagged but
# never-detected animal still gets a row, and so positionally-supplied tagging.dates stay aligned),
# while a count of "how much data do I have" must be observed. The bugs live where one is used in
# place of the other.
#
# The two quantities diverge by three routes, and the fixtures below cover all three because they are
# NOT interchangeable - route 1 is a genuine user DECLARATION, routes 2 and 3 are accidents of how the
# data frame was built, and several call sites are correct under one and wrong under another:
#
#   route 1  id.groups (passed explicitly, or auto-resolved from mobyData metadata) makes
#            validateArguments() filter rows to the roster AND re-level the factor to it
#   route 2  subsetting without droplevels(), leaving stale levels behind
#   route 3  a pre-factored ID column that carries levels the data never had


# Capture moby console output. cli headers, notes and the plot summary cards are emitted as MESSAGE
# conditions on stderr, so expect_output() cannot see them and the suite's suppressMessages() wrappers
# swallow them; only S3 print methods (print.mobyData, print.mobyFilter) use cat() to stdout. Grabbing
# both streams is what stops a console assertion from silently testing an empty string.
console_text <- function(expr) {
  out <- character(0)
  err <- capture.output(
    out <- capture.output({ .res <- force(expr); invisible(NULL) }, type = "output"),
    type = "message")
  paste(c(out, err), collapse = "\n")
}

# ROUTE 3 - a pre-factored ID column carrying one level the data never had.
# The universal fixture: it needs no id.groups formal, so it reaches the functions route 1 cannot
# (calculateResidency, calculateROM, calculateLinearityIndex, createWideTable, filterDetections,
# plotChronogram, plotContours). `pos` places the phantom first/middle/last in the level order, which
# is what proves POSITIONAL alignment rather than mere tolerance of a trailing NA.
ghost_levels_df <- function(pos = c("last", "first", "middle"), ghost = "GHOST") {
  pos <- match.arg(pos)
  d <- make_detections(ids = c("A", "B", "C"), n_per_id = 12)
  d$timebin <- d$datetime          # make_detections emits hourly stamps; they are the bins
  real <- c("A", "B", "C")
  lv <- switch(pos,
               last   = c(real, ghost),
               first  = c(ghost, real),
               middle = c("A", ghost, "B", "C"))
  d$ID <- factor(d$ID, levels = lv)
  d
}

# Tagging dates covering every LEVEL of ghost_levels_df(), including the phantom. Needed because
# .checkAnimalParams() errors when a level has no date, and because a level-ordered unnamed vector is
# how a positional-alignment bug would show up.
ghost_tagging_dates <- function(d, id.col = "ID") {
  stats::setNames(rep(min(d$datetime) - 86400, nlevels(d[[id.col]])), levels(d[[id.col]]))
}

# ROUTE 2 - subset without droplevels, on a PLAIN data frame with no moby metadata at all.
# The `moby` attribute must be absent: .resolveArgs() injects id.groups from metadata whenever it can,
# which sends calculateTransitions down its roster branch and hides the NULL-id.groups bug entirely.
stale_levels_df <- function(keep = c("A", "B")) {
  d <- make_detections(ids = c("A", "B", "C", "D"), n_per_id = 12)
  d$timebin <- d$datetime
  d$ID <- factor(d$ID)
  out <- d[d$ID %in% keep, ]              # rows dropped, levels deliberately retained
  attr(out, "moby") <- NULL
  rownames(out) <- NULL
  out
}

# ROUTE 2 through `[.mobyData` - the real-world case, needing no contrivance at all.
# `[.mobyData` re-attaches metadata verbatim and never calls .syncIDs(), so an ordinary subset of the
# shipped dataset keeps the full 8-ID roster; .resolveArgs() then feeds it back as id.groups and
# validateArguments() re-levels the subset to 8.
stale_moby_subset <- function(keep = c("R01", "R02")) {
  data(rays, envir = environment())
  rays[as.character(rays$ID) %in% keep, ]
}

# A wide (time bin x individual) table carrying a phantom all-NA column, for the association tests.
# The two date regimes are BOTH required: .pairwiseOverlap()'s "no detections" guard tests
# is.na(end_dates[id]), which is only a valid proxy when createWideTable DERIVED the dates from the
# data. A user-supplied scalar date gives the phantom a real monitoring window and defeats the guard.
ghost_wide <- function(dates = c("derived", "scalar")) {
  dates <- match.arg(dates)
  d <- ghost_levels_df()
  if (dates == "derived")
    createWideTable(d, value.col = "station", verbose = FALSE)
  else
    createWideTable(d, value.col = "station", verbose = FALSE,
                    start.dates = min(d$timebin), end.dates = max(d$timebin))
}
