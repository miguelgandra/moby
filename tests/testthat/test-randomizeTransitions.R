rt_network <- function(seed = 1) {
  set.seed(seed)
  mk <- function(id, sites) data.frame(
    ID = id, datetime = as.POSIXct("2023-01-01", tz = "UTC") + seq_along(sites) * 3600, site = sites)
  d <- rbind(mk("A", rep(c("S1", "S2"), 5)), mk("B", rep(c("S1", "S2"), 5)),
             mk("C", rep(c("S2", "S3"), 5)), mk("D", c("S1", "S2", "S3", "S1", "S2", "S3")))
  d$ID <- factor(d$ID)
  md <- as_moby(d, station.col = "site", tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  calculateTransitions(md, spatial.col = "site")
}

test_that("randomizeTransitions returns edge and network tests with adjusted p-values", {
  rt <- randomizeTransitions(rt_network(), iterations = 200, random.seed = 1)
  expect_true(all(c("edges", "network") %in% names(rt)))
  expect_true(all(c("n_movements", "mean_null", "p_value", "p_adjusted", "association") %in% colnames(rt$edges)))
  expect_true(all(rt$edges$p_adjusted >= rt$edges$p_value - 1e-9))   # BH never decreases p
  expect_true(all(c("n_transitions", "mean_null", "p_value") %in% colnames(rt$network)))
})

test_that("randomizeTransitions detects strongly-used transitions", {
  rt <- randomizeTransitions(rt_network(), iterations = 300, random.seed = 1)
  s1s2 <- rt$edges[rt$edges$from == "S1" & rt$edges$to == "S2", ]
  expect_lt(s1s2$p_adjusted, 0.05)          # a dominant transition is non-random
  expect_equal(s1s2$association, "more")
})

test_that("randomizeTransitions is reproducible and supports p.adjust.method='none'", {
  a <- randomizeTransitions(rt_network(), iterations = 150, random.seed = 7)
  b <- randomizeTransitions(rt_network(), iterations = 150, random.seed = 7)
  expect_equal(a$edges$p_value, b$edges$p_value)
  n <- randomizeTransitions(rt_network(), iterations = 150, random.seed = 7, p.adjust.method = "none")
  expect_equal(n$edges$p_value, n$edges$p_adjusted)
})

test_that("randomizeTransitions validates its input", {
  expect_error(randomizeTransitions(data.frame(a = 1)), "movement")
  expect_error(randomizeTransitions(rt_network(), p.adjust.method = "bogus"), "p.adjust.method")
})
