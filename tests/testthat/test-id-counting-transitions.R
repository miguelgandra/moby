# CHARACTERISE-BEFORE for R/calculateTransitions.R:122
#
#   group_sizes <- stats::setNames(nlevels(data[, id.col]), "all")
#
# This is the NULL-id.groups branch. It stores the DECLARED count (the factor roster) as the
# network's "group.sizes" attribute, and transitionsTable() divides by that to produce the
# percentage in its "Individuals" column. When the roster carries levels the data never had, every
# percentage is deflated and nothing warns.
#
# The sibling branch at :127 (group_sizes <- lengths(id.groups)) is CORRECT to stay declared: an
# explicit id.groups is a real user declaration, and "3 of the 5 tagged (60%)" is the intended
# reading. The final test in this file pins that, so the fix cannot be applied to both branches.


# The same call on the same data, with the never-observed levels dropped. This is the answer the
# fixed code must produce; deriving it here rather than hard-coding percentages keeps the test a
# statement about the contract instead of a record of the symptom.
observed_truth <- function(d) {
  net <- calculateTransitions(droplevels(d), spatial.col = "station", verbose = FALSE)
  list(net = net, table = transitionsTable(net, verbose = FALSE))
}

# networkEdges() hands back the network object itself, so its attributes (group.sizes,
# processing.date) ride along and identical() would compare those too. Comparing the *network* of
# two runs means comparing the bare edge columns.
edge_core <- function(net) {
  e <- networkEdges(net)
  data.frame(e[, c("group", "from", "to", "n_movements", "n_individuals", "mean_duration_h")],
             row.names = NULL, check.names = FALSE)
}


test_that("group.sizes on a NULL-id.groups network counts declared levels, not observed animals", {
  d <- stale_levels_df()                       # 4 levels retained by the subset, 2 animals present
  # the fixture only bites if the NULL-id.groups branch is actually reached: no moby metadata means
  # .resolveArgs() has no id.groups to inject
  expect_null(attr(d, "moby"))
  expect_equal(nlevels(d$ID), 4L)
  expect_equal(length(unique(as.character(d$ID))), 2L)

  net <- calculateTransitions(d, spatial.col = "station", verbose = FALSE)
  truth <- observed_truth(d)

  # REGRESSION: counts the animals that contributed data, not the 4 retained factor levels.
  expect_equal(unname(attr(net, "group.sizes")), 2L)

  # the truth, computed from the same data with the stale levels dropped
  expect_equal(unname(attr(truth$net, "group.sizes")), 2L)

  # REGRESSION: the stale and dropped runs now agree, as they must - they describe the same 2
  # animals moving between the same stations. Stale levels no longer reach the denominator.
  expect_equal(attr(net, "group.sizes"), attr(truth$net, "group.sizes"))

  # the correction is silent by design: stale levels are a normal consequence of subsetting, and
  # once the count ignores them there is nothing for the user to act on
  expect_no_warning(calculateTransitions(d, spatial.col = "station", verbose = FALSE))

  # INVARIANT: must be identical before and after the fix. The attribute stays a single named
  # integer under the "all" pseudo-group, and no id.groups is invented.
  expect_named(attr(net, "group.sizes"), "all")
  expect_type(attr(net, "group.sizes"), "integer")
  expect_null(attr(net, "id.groups"))
})


test_that("transitionsTable percentages divide by the declared roster", {
  d <- stale_levels_df()
  net <- calculateTransitions(d, spatial.col = "station", verbose = FALSE)
  tt <- transitionsTable(net, verbose = FALSE)
  truth <- observed_truth(d)

  # REGRESSION: both animals made the R1 --> R2 move, so the denominator is 2, not the 4 declared
  # levels. Before the fix this read "2 (50%)".
  expect_equal(tt$Individuals[tt$Type == "R1 --> R2"], "2 (100%)")
  # one of the two made R2 --> R1. Before the fix this read "1 (25%)".
  expect_equal(tt$Individuals[tt$Type == "R2 --> R1"], "1 (50%)")

  # the whole column, against the same column computed from the droplevels() run
  expect_equal(truth$table$Individuals[truth$table$Type == "R1 --> R2"], "2 (100%)")
  expect_equal(truth$table$Individuals[truth$table$Type == "R2 --> R1"], "1 (50%)")
  # REGRESSION: every percentage now matches the droplevels() run. This is the headline assertion -
  # a published transitions table no longer depends on whether the input carried stale levels.
  expect_equal(tt$Individuals, truth$table$Individuals)

  # INVARIANT: must be identical before and after the fix. Only the denominator is wrong - the
  # numerator (how many distinct animals made each move) is already counted from the data.
  expect_equal(sub(" .*$", "", tt$Individuals), sub(" .*$", "", truth$table$Individuals))
  expect_equal(tt$Type, truth$table$Type)
  expect_equal(tt$Movements, truth$table$Movements)
  expect_equal(colnames(tt), c("Type", "Movements", "Individuals", "Mean duration (h)"))
  expect_equal(nrow(tt), 6L)                   # 3 stations, both directions, all realised
})


test_that("stale ID levels do not move the network itself", {
  d <- stale_levels_df()
  net <- calculateTransitions(d, spatial.col = "station", verbose = FALSE)
  truth <- observed_truth(d)

  # INVARIANT: must be identical before and after the fix. The edge and node tables are built from
  # the rows, so the phantom levels never reach them; this is what proves the fix is confined to the
  # group.sizes denominator.
  expect_equal(edge_core(net), edge_core(truth$net))
  expect_equal(networkNodes(net), networkNodes(truth$net))

  e <- networkEdges(net)
  expect_equal(nrow(e), 6L)
  expect_equal(paste(e$from, e$to), c("R1 R2", "R1 R3", "R2 R1", "R2 R3", "R3 R1", "R3 R2"))
  expect_equal(e$n_movements, c(2L, 3L, 1L, 2L, 4L, 1L))
  expect_equal(e$n_individuals, c(2L, 2L, 1L, 2L, 2L, 1L))
  expect_equal(attr(net, "ordered.sites"), c("R1", "R2", "R3"))
  # the per-transition records are the raw material of the edges; 13 moves across the 2 animals
  expect_equal(nrow(attr(net, "transition_records")[["all"]]), 13L)
  expect_setequal(unique(attr(net, "transition_records")[["all"]]$id), c("A", "B"))
})


test_that("a phantom factor level deflates the percentages the same way", {
  # route 3 rather than route 2: the ID column was pre-factored with a level the data never had, so
  # no subsetting is involved. Same bug, reached without touching the rows.
  d <- ghost_levels_df("first")                # phantom first in the level order: 4 levels, 3 animals
  expect_equal(nlevels(d$ID), 4L)
  expect_equal(length(unique(as.character(d$ID))), 3L)

  net <- calculateTransitions(d, spatial.col = "station", verbose = FALSE)
  tt <- transitionsTable(net, verbose = FALSE)
  truth <- observed_truth(d)

  # REGRESSION: route 3 (a pre-factored phantom level, no subsetting) is corrected too - the fix is
  # to the count itself, not to any subsetting path. Before the fix this read 4L.
  expect_equal(unname(attr(net, "group.sizes")), 3L)
  expect_equal(unname(attr(truth$net, "group.sizes")), 3L)

  # REGRESSION: all three animals made R1 --> R2. Before the fix this read "3 (75%)".
  expect_equal(tt$Individuals[tt$Type == "R1 --> R2"], "3 (100%)")
  expect_equal(truth$table$Individuals[truth$table$Type == "R1 --> R2"], "3 (100%)")

  # INVARIANT: must be identical before and after the fix.
  expect_equal(edge_core(net), edge_core(truth$net))
  expect_equal(sub(" .*$", "", tt$Individuals), sub(" .*$", "", truth$table$Individuals))
})


test_that("an explicit id.groups roster stays the denominator (the branch that must NOT change)", {
  # INVARIANT: must be identical before and after the fix. calculateTransitions() has two branches;
  # only the NULL-id.groups one is wrong. Here the user declares a 4-animal roster and only 3 were
  # detected, and "of the 4 tagged" is the intended reading - lengths(id.groups) must survive the fix.
  sub <- stale_moby_subset(c("R01", "R02", "R03"))
  net <- suppressWarnings(
    calculateTransitions(sub, spatial.col = "station", verbose = FALSE,
                         id.groups = list(clavata = c("R01", "R02", "R03", "R04"))))

  expect_equal(attr(net, "group.sizes"), c(clavata = 4L))
  expect_equal(names(attr(net, "id.groups")), "clavata")

  tt <- suppressWarnings(transitionsTable(net, verbose = FALSE))
  # 1 of the 4 declared animals = 25%, and 2 of 4 = 50%: percentages against the declared roster,
  # which is correct here precisely because the roster was declared.
  expect_true(all(grepl("^[0-9]+ \\([0-9]+%\\)$", tt$Individuals[nzchar(tt$Individuals)])))
  expect_equal(tt$Individuals[tt$Type == "ST01 --> ST05"], "1 (25%)")
  expect_equal(tt$Individuals[tt$Type == "ST01 --> ST04"], "2 (50%)")
})
