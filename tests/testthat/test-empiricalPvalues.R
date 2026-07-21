# Regression pins for the shared permutation-inference helper. The one-tailed tails were previously
# inverted in randomizeAssociations (a strong positive association returned p ~ 1, "non-significant");
# these tests lock the correct direction so it cannot silently regress again.

test_that(".empiricalPvalues uses the correct one-tailed tails", {
  set.seed(1)
  null <- matrix(rnorm(1000, 0, 1), nrow = 1)          # null centred at 0
  obs_high <- 10                                        # observed clearly ABOVE the null

  g <- .empiricalPvalues(obs_high, null, alternative = "greater")
  expect_lt(g$p_value, 0.02)                            # right tail: obs >> null -> significant
  expect_equal(g$association, "more")

  l <- .empiricalPvalues(obs_high, null, alternative = "less")
  expect_gt(l$p_value, 0.9)                             # left tail: obs >> null -> NOT less -> ns
  expect_equal(l$association, "non-significant")

  # a low observation flips the picture
  obs_low <- -10
  expect_gt(.empiricalPvalues(obs_low, null, alternative = "greater")$p_value, 0.9)
  expect_lt(.empiricalPvalues(obs_low, null, alternative = "less")$p_value, 0.02)
})

test_that(".empiricalPvalues caps the two-sided p-value at 1 and honours custom labels", {
  set.seed(2)
  null <- matrix(rnorm(500, 5, 1), nrow = 1)
  ts <- .empiricalPvalues(5, null, alternative = "two.sided")   # obs at the null centre
  expect_lte(ts$p_value, 1)                                     # previously could exceed 1

  lab <- .empiricalPvalues(10, matrix(rnorm(500), nrow = 1), alternative = "greater",
                           labels = c(more = "positive", less = "negative"))
  expect_equal(lab$association, "positive")                     # association-plot labels
})

test_that(".empiricalPvalues is NA-robust and supports method='none'", {
  null <- matrix(c(rnorm(100), rep(NA, 100)), nrow = 2, byrow = TRUE)  # row 2 all-NA
  obs <- c(3, NA)
  res <- .empiricalPvalues(obs, null, alternative = "greater", p.adjust.method = "none")
  expect_true(is.na(res$p_value[2]) && is.na(res$association[2]))      # undefined row -> NA
  expect_false(is.na(res$p_value[1]))
  expect_equal(res$p_value, res$p_adjusted)                            # method='none' -> unchanged
})
