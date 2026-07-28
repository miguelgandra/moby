# Console regression tests for the OBSERVED-vs-DECLARED contract (Tier B of the ID-counting audit).
#
# These assert PRINTED TEXT, so they must go through console_text() (helper-id-divergence.R): moby's
# headers and plot summary cards are message conditions on stderr, which expect_output() cannot see.
#
# The rule under test: a count describing "how much data is there" reports OBSERVED individuals.
# Counts where the roster IS the subject stay DECLARED and are covered elsewhere.

test_that("calculateVisits and calculateTransitions report observed individuals, not stale levels", {
  d <- stale_levels_df()                    # 4 levels retained by the subset, 2 animals present
  expect_equal(nlevels(d$ID), 4L)
  expect_equal(length(unique(as.character(d$ID))), 2L)

  # before the fix both announced "4 individuals" for a run covering 2
  expect_match(console_text(calculateVisits(d, spatial.col = "station")), "2 individuals")
  expect_match(console_text(calculateTransitions(d, spatial.col = "station")), "2 individuals")
  expect_false(grepl("4 individuals", console_text(calculateVisits(d, spatial.col = "station"))))

  # the two share .residenceRuns() and their headers are deliberately mirrored: the individual count
  # must be worded identically or they will describe the same dataset differently
  vis <- grep("Input:", strsplit(console_text(calculateVisits(d, spatial.col = "station")), "\n")[[1]],
              value = TRUE)
  tra <- grep("Input:", strsplit(console_text(calculateTransitions(d, spatial.col = "station")), "\n")[[1]],
              value = TRUE)
  expect_equal(sub(".*(\\d+ individuals?).*", "\\1", vis), sub(".*(\\d+ individuals?).*", "\\1", tra))
})


test_that("a phantom roster level does not inflate the reported individual count", {
  d <- ghost_levels_df()                    # 4 levels, 3 observed
  expect_match(console_text(calculateVisits(d, spatial.col = "station")), "3 individuals")
  expect_match(console_text(calculateTransitions(d, spatial.col = "station")), "3 individuals")
})


# The plot cards need a series long enough to bin (the 12-hour fixture yields a single date bin), so
# these build a 10-day version of the same route-3 divergence: 4 declared levels, 3 observed.
ghost_plot_df <- function(days = 10) {
  d <- make_detections(ids = c("A", "B", "C"), n_per_id = 24 * days)
  d$timebin <- d$datetime
  d$ID <- factor(d$ID, levels = c("A", "B", "C", "GHOST"))
  d
}

test_that("plotContours' card matches the sample size printed inside the figure", {
  # plotContours bins by date, so it needs a long series (mirrors test-plotContours.R's fixture);
  # the divergence is the same route 3, with a fourth declared level nothing was recorded against.
  set.seed(1)
  n <- 1500
  dt <- as.POSIXct("2022-11-01", tz = "UTC") + sort(stats::runif(n, 0, 300 * 86400))
  d <- data.frame(ID = factor(sample(c("A", "B", "C"), n, replace = TRUE),
                              levels = c("A", "B", "C", "GHOST")),
                  datetime = dt,
                  depth = 15 + stats::rnorm(n, 0, 2), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".png"); on.exit(unlink(f), add = TRUE)

  # plotContours already prints the observed count INSIDE the panel as "(n=...)", so a declared card
  # made the function contradict itself within a single call
  txt <- console_text(plotContours(d, variables = "depth", id.col = "ID",
                                   datetime.col = "datetime", diel.lines = 0, file = f))
  expect_match(txt, "Individuals:\\s+3")
  expect_false(grepl("Individuals:\\s+4", txt))
})


test_that("plotChronogram's card counts the individuals actually drawn", {
  d <- ghost_plot_df()
  f <- tempfile(fileext = ".png"); on.exit(unlink(f), add = TRUE)
  txt <- console_text(plotChronogram(d, id.col = "ID", timebin.col = "timebin",
                                     station.col = "station", diel.lines = 0, file = f))
  expect_match(txt, "Individuals:\\s+3")
  expect_false(grepl("Individuals:\\s+4", txt))
})


test_that("plotMaps' card arithmetic is internally consistent", {
  set.seed(42)
  d <- ghost_plot_df(3)
  d$lon <- -8 + stats::runif(nrow(d), 0, 0.05)
  d$lat <- 37 + stats::runif(nrow(d), 0, 0.05)
  f <- tempfile(fileext = ".png"); on.exit(unlink(f), add = TRUE)

  txt <- console_text(plotMaps(d, epsg.code = 32629, file = f))
  # "Individuals: n_ids of n_total" must reconcile with the lines beneath it. Before the fix the
  # denominator was re-factored (dropping the phantom) while the discard count was not, so the card
  # printed "3 of 3" directly above a line reporting one individual absent.
  m <- regmatches(txt, regexpr("Individuals:\\s+(\\d+) of (\\d+)", txt))
  expect_length(m, 1L)
  n <- as.integer(regmatches(m, gregexpr("\\d+", m))[[1]])
  n_ids <- n[1]; n_total <- n[2]
  expect_equal(n_total, 4L)                 # the roster the discard count is subtracted FROM
  expect_lt(n_ids, n_total)

  # and the absence is attributed to the right cause: a never-detected animal had zero detections,
  # so describing it as "< 5 detections" was a false statement
  expect_match(txt, "with no detections \\(not mapped\\)")
})


test_that("roster-shaped tables report observed individuals and disclose the rest separately", {
  d <- ghost_levels_df(); d$dist_m <- stats::runif(nrow(d), 0, 500)   # 4 declared, 3 observed
  NOTE <- "1 tagged individual with no detections"

  # The header answers "how much data is there" plainly, and the roster difference is stated on its
  # own line rather than folded into the count ("3 of 4 individuals"), which is harder to scan.
  rom <- console_text(calculateROM(d))
  expect_match(rom, "3 individuals")
  expect_false(grepl("4 individuals", rom))
  expect_match(rom, NOTE, fixed = TRUE)

  # INVARIANT: the returned table stays roster-shaped. The header count and the row count differ ON
  # PURPOSE - that is what the note explains - and the row is what keeps positionally-supplied
  # tagging.dates aligned.
  expect_equal(nrow(suppressMessages(calculateROM(d, verbose = FALSE))), nlevels(d$ID))

  wide <- console_text(createWideTable(d, value.col = "station"))
  expect_match(wide, "3 individuals")
  expect_match(wide, "column retained, all NA", fixed = TRUE)
  # the wide table keeps one column per declared individual
  expect_equal(ncol(suppressMessages(createWideTable(d, value.col = "station", verbose = FALSE))),
               nlevels(d$ID) + 1L)
})


test_that("a clean run gains nothing from the disclosure note", {
  d <- droplevels(ghost_levels_df()); d$dist_m <- stats::runif(nrow(d), 0, 500)
  expect_equal(nlevels(d$ID), 3L)
  txt <- console_text(calculateROM(d))
  expect_match(txt, "3 individuals")
  expect_false(grepl("no detections", txt))
  expect_false(grepl("tagged individual", txt))
})
