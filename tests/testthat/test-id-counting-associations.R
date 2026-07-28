# OBSERVED vs DECLARED individuals in calculateAssociations().
#
# The wide table carries a column per DECLARED individual (one per factor level), so a tagged but
# never-detected animal gets an all-NA column and a node in the network. That is correct and must
# stay: the roster is the node set. What must NOT happen is for that animal to be handed a numeric
# association with everybody else.
#
# THE SCIENTIFIC POINT. An association index of 0 is a SAMPLING zero: both animals were monitored
# over a shared window and were never once recorded in the same receiver/time bin - that is real
# evidence of non-association. An animal that was never detected at all yields a STRUCTURAL zero:
# there is no evidence either way, and the dyad is unmeasurable. Conflating the two deflates the
# mean association, inflates apparent network sparsity, and - because randomizeAssociations()
# permutes the same wide table - carries the fabricated zeros into BOTH the observed and the null
# distributions, so the null shifts to meet the fabricated observation and the deflation is invisible
# in the p-values.
#
# .pairwiseOverlap() already INTENDS the right thing:
#     # if any of the individuals doesn't have detections, jump to next pair
#     if (any(is.na(end_dates[c(id1, id2)]))) return(result)
# but it tests is.na(end_dates[id]) as a proxy for "has no detections". That proxy only holds when
# createWideTable() DERIVED the monitoring window from the data (no detections -> no last detection
# -> NA). A user-supplied SCALAR start.dates/end.dates gives every individual - phantom included - a
# real monitoring window, the proxy returns FALSE, the guard never fires, and the phantom is scored
# 0 against everyone.
#
# The two date regimes are therefore pinned SIDE BY SIDE throughout this file. Identical detections,
# identical roster, two different answers for the same dyads: that inconsistency is what makes this a
# bug rather than a defensible convention.

# ghost_wide() covers the phantom-last case, which is all the two-regime comparison needs. This
# variant (declared, as required) additionally routes the phantom through the FIRST and MIDDLE
# positions of the level order, so the invariant below proves POSITIONAL alignment of the real dyads
# rather than mere tolerance of a trailing column.
ghost_wide_at <- function(pos, dates = c("derived", "scalar")) {
  dates <- match.arg(dates)
  d <- ghost_levels_df(pos)
  if (dates == "derived") {
    createWideTable(d, value.col = "station", verbose = FALSE)
  } else {
    createWideTable(d, value.col = "station", verbose = FALSE,
                    start.dates = min(d$timebin), end.dates = max(d$timebin))
  }
}

assoc <- function(w) networkEdges(suppressMessages(calculateAssociations(w, verbose = FALSE)))

is_ghost <- function(e, ghost = "GHOST") e$id1 == ghost | e$id2 == ghost

EDGE_COLS <- c("id1", "id2", "association", "co_occurrences", "shared_monit_days",
               "start", "end", "mean_consec_overlap", "max_consec_overlap")

# Make two edge tables comparable across runs that differ only in whether a phantom was present.
# Two incidental differences have to be normalised away, neither of them this bug:
#   * networkEdges() carries the run's `ids` attribute, which legitimately differs (4 vs 3);
#   * when the FIRST result row is a skipped dyad (phantom at pos = "first"/"middle"), its all-NA
#     start/end costs the whole column its POSIXct class in .rbindFill(), so the real dyads' dates
#     come back as raw epoch seconds. Separate latent quirk, orthogonal to the guard under audit -
#     and one the fix will ALSO trigger in the scalar regime once phantom dyads start being skipped
#     there too, so comparing the instants rather than the class is what keeps this invariant honest.
norm_edges <- function(e) {
  e <- e[, EDGE_COLS]
  out <- data.frame(e[, c("id1", "id2")], stringsAsFactors = FALSE)
  for (v in c("association", "co_occurrences", "shared_monit_days",
              "mean_consec_overlap", "max_consec_overlap")) out[[v]] <- e[[v]]
  for (v in c("start", "end")) out[[v]] <- as.numeric(e[[v]])   # epoch seconds: class-agnostic
  rownames(out) <- NULL
  out[, EDGE_COLS]
}

# The three real dyads of ghost_levels_df(), computed once from a droplevels() run so that every
# assertion below compares against an independently derived truth rather than a magic number. The
# phantom cannot influence these: it contributes no detections and no shared window.
real_truth <- function() {
  d <- ghost_levels_df()
  d$ID <- droplevels(d$ID)                       # no phantom level at all -> 3 individuals, 3 dyads
  assoc(createWideTable(d, value.col = "station", verbose = FALSE))
}


test_that("derived monitoring window: a never-detected individual yields NA, not zero", {
  # The correct regime, pinned so that the buggy one below has something to be inconsistent with.
  e <- assoc(ghost_wide("derived"))
  truth <- real_truth()

  # INVARIANT: must be identical before and after the fix - the roster is the node set, so the
  # phantom still gets its three dyads, it just has no measurable association.
  expect_equal(nrow(e), choose(4, 2))
  expect_equal(sum(is_ghost(e)), 3L)

  # INVARIANT: a structural zero is reported as "unmeasurable", across every column that would
  # otherwise imply the pair had been observed together.
  g <- e[is_ghost(e), ]
  expect_true(all(is.na(g$association)))
  expect_true(all(is.na(g$co_occurrences)))
  expect_equal(g$shared_monit_days, rep(0, 3))
  expect_true(all(is.na(g$start)) && all(is.na(g$end)))

  # INVARIANT: the column set is the documented edge table, and the real dyads match the no-phantom
  # run exactly.
  expect_named(e, EDGE_COLS)
  expect_equal(norm_edges(e[!is_ghost(e), ]), norm_edges(truth))

  # INVARIANT: with the phantom's dyads excluded rather than zeroed, the headline statistic is the
  # mean over the three measurable dyads.
  expect_equal(mean(e$association, na.rm = TRUE), mean(truth$association))
})


test_that("scalar monitoring window fabricates association = 0 for a never-detected individual", {
  e <- assoc(ghost_wide("scalar"))
  truth <- real_truth()
  g <- e[is_ghost(e), ]

  # REGRESSION: NA, not 0, for every dyad involving the never-detected GHOST - identical to what the
  # derived-date run above produces from the very same detections. An SRI of 0 is a SAMPLING zero
  # ("both monitored, never co-occurred"); a never-detected animal is a STRUCTURAL zero ("no evidence
  # either way"). Before the fix this regime reported 0, purely because a scalar date defeated the
  # guard's is.na(end_dates) proxy.
  expect_true(all(is.na(g$association)))
  expect_true(all(is.na(g$co_occurrences)))
  # and it is no longer credited with the monitoring window the scalar dates handed it
  expect_equal(g$shared_monit_days, rep(0, 3))
  expect_true(all(is.na(g$start)))
  expect_true(all(is.na(g$end)))

  # REGRESSION: the headline statistic no longer depends on how the monitoring window was specified.
  # Before the fix the three fabricated zeros halved it (18.06 vs 36.11); now both date regimes agree
  # on the mean over the measurable dyads.
  expect_equal(mean(e$association, na.rm = TRUE), mean(truth$association))

  # INVARIANT: must be identical before and after the fix. Only the phantom's dyads are wrong; the
  # real animals' associations, and the shape of the table, must not move.
  expect_equal(nrow(e), choose(4, 2))
  expect_named(e, EDGE_COLS)
  expect_equal(norm_edges(e[!is_ghost(e), ]), norm_edges(truth))
})


test_that("real dyads are identical across both date regimes and every phantom position", {
  # INVARIANT: must be identical before and after the fix. This is the surgical-fix proof: whatever
  # the fix does to the phantom's dyads, the six combinations below must keep returning the same
  # three real dyads, in the same order, with the same values. Varying the phantom's position in the
  # level order (first / middle / last) is what makes this an alignment test and not just a "the
  # numbers happen to match" test - with pos = "first" or "middle" the phantom's dyads are
  # interleaved with the real ones rather than trailing them.
  truth <- real_truth()
  for (pos in c("last", "first", "middle")) {
    for (dates in c("derived", "scalar")) {
      e <- assoc(ghost_wide_at(pos, dates))
      expect_equal(nrow(e), choose(4, 2), info = paste(pos, dates))
      expect_equal(norm_edges(e[!is_ghost(e), ]), norm_edges(truth), info = paste(pos, dates))
    }
  }
})


test_that("the 'no detections' console note is silenced by a scalar monitoring window", {
  # calculateAssociations() reports missing individuals from the same proxy the dyad guard uses
  # (missing_individuals <- names(which(is.na(end_dates)))), so the user-facing warning about the
  # phantom disappears in lockstep with the fabricated zeros - the one regime that produces wrong
  # numbers is also the one that says nothing about it.
  # (console_text() because moby's notes are message conditions on stderr.)

  # INVARIANT: must be identical before and after the fix.
  expect_match(console_text(calculateAssociations(ghost_wide("derived"), verbose = TRUE)),
               "1 individual with no detections")

  # CHARACTERISE-BEFORE (bug): no note is emitted, because end_dates has no NA to count. The correct
  # behaviour is the same note as the derived run - the animal has no detections under either date
  # regime. This assertion FLIPS when the fix lands (expect_match, same pattern).
  expect_no_match(console_text(calculateAssociations(ghost_wide("scalar"), verbose = TRUE)),
                  "with no detections")
})


test_that("fabricated zeros are carried into randomizeAssociations", {
  rz <- function(w) suppressWarnings(suppressMessages(
    randomizeAssociations(w, suppressMessages(calculateAssociations(w, verbose = FALSE)),
                          iterations = 20, random.seed = 1, verbose = FALSE)))$pairwise_results

  pd <- rz(ghost_wide("derived"))
  ps <- rz(ghost_wide("scalar"))
  gd <- grepl("GHOST", pd$pair); gs <- grepl("GHOST", ps$pair)

  # INVARIANT: must be identical before and after the fix. An unmeasurable dyad has no null
  # distribution and no test: NaN mean, no p-value, no significance label.
  expect_true(all(is.nan(pd$mean_null[gd])))
  expect_true(all(is.na(pd$p_value[gd])))
  expect_equal(sum(!is.na(pd$p_value)), 3L)          # only the 3 measurable dyads are tested

  # REGRESSION: an unmeasurable dyad gets no null distribution and no test, under EITHER date
  # regime. Before the fix the scalar regime manufactured a formal result from no data - a null of
  # fabricated zeros (permuting an all-NA column can only yield 0), an observed 0 matching it
  # exactly, and p = 1.000 labelled "non-significant".
  expect_true(all(is.na(ps$p_value[gs])))
  expect_true(all(is.na(ps$significance[gs])))

  # REGRESSION: only the 3 measurable dyads enter the multiple-testing correction. Before the fix
  # all 6 did, so the fabricated tests diluted the adjusted p-values of the real ones.
  expect_equal(sum(!is.na(ps$p_adjusted)), 3L)
})
