
mobyEnv <- NULL

# this will run automatically when package is loaded
.onLoad <- function(libname, pkgname){

  # internal environment used to track session state (e.g. one-off warning flags).
  # Dataset-level defaults are no longer stored here: they now travel with the data via
  # the mobyData object (see as_moby()).
  mobyEnv <<- new.env()
}


# Declare symbols that are bound dynamically at run time (e.g. via list2env(),
# foreach iterators and non-standard evaluation) so that R CMD check does not
# report them as undefined global variables.
utils::globalVariables(c(
  "pairwise_combinations", "complete_ids", "id.groups", "group.comparisons",
  "start_dates", "end_dates", "timebin.col", "metric", "data", "p", "reason"
))


# this will run automatically when package is attached
.onAttach <- function(libname, pkgname){

  # Startup banner in the shared "house style" used across the telemetry package family (a borderless
  # =/- frame with three zones: identity, capabilities, pointers). The name is bold and the frame
  # rules + pointer lines are dim, but ONLY on interactive colour-capable terminals; everywhere else
  # (log files, knitr, non-ANSI front-ends) it degrades to plain text, so nothing depends on styling.
  W    <- 54L
  rule <- function(ch) strrep(ch, W)
  ansi <- interactive() && isTRUE(getOption("crayon.enabled", interactive()))
  bold <- function(x) if (ansi) paste0("\033[1m", x, "\033[0m") else x
  dim2 <- function(x) if (ansi) paste0("\033[2m", x, "\033[0m") else x
  ptr  <- function(label, target) dim2(paste0("  ", formatC(label, width = -17L), " ", target))

  packageStartupMessage(paste(c(
    dim2(rule("=")),
    bold(sprintf(" moby v%s", utils::packageVersion("moby"))),
    dim2(rule("-")),
    "  - filter & quality-control acoustic detections",
    "  - estimate residency, home range & movement networks",
    dim2(rule("-")),
    ptr("dive in:",          "help(package = \"moby\")"),
    ptr("chart the depths:", "https://github.com/miguelgandra/moby"),
    ptr("cite your voyage:", "citation(\"moby\")"),
    dim2(rule("="))
  ), collapse = "\n"))
}
