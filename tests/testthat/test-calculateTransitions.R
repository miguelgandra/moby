trans_dataset <- function() {
  d <- data.frame(
    ID = factor(c(rep("A", 6), rep("B", 4))),
    datetime = as.POSIXct("2023-01-01", tz = "UTC") + c(0, 1, 2, 3, 4, 5, 0, 1, 2, 3) * 3600,
    site = c("S1", "S1", "S2", "S2", "S3", "S1", "S2", "S3", "S3", "S2"),
    lon = c(-8, -8, -7, -7, -6, -8, -7, -6, -6, -7),
    lat = c(37, 37, 38, 38, 39, 37, 38, 39, 39, 38))
  as_moby(d, station.col = "site", lon.col = "lon", lat.col = "lat",
          tagging.dates = as.POSIXct(c(A = "2023-01-01", B = "2023-01-01"), tz = "UTC"))
}

test_that("calculateTransitions builds a directed movement network with correct edge weights", {
  md <- trans_dataset()
  net <- calculateTransitions(md, spatial.col = "site")
  expect_true(is_mobyNetwork(net))
  expect_equal(networkType(net), "movement")
  e <- networkEdges(net); n <- networkNodes(net)
  expect_true(all(e$from != e$to))                            # no self-loops
  expect_true(all(c("n_movements", "n_individuals", "mean_duration_h") %in% colnames(e)))

  s2s3 <- e[e$from == "S2" & e$to == "S3", ]
  expect_equal(s2s3$n_movements, 2L)                          # A and B both did S2->S3
  expect_equal(s2s3$n_individuals, 2L)
  s1s2 <- e[e$from == "S1" & e$to == "S2", ]
  expect_equal(s1s2$n_movements, 1L)                          # only A

  expect_true(all(c("site", "n_detections", "n_individuals", "n_residence",
                    "mean_residence_h", "lon", "lat") %in% colnames(n)))
  # n_residence = number of residence events (visits): A visits S1 twice; S2 = A(1)+B(2); S3 = A(1)+B(1)
  expect_equal(n$n_residence[n$site == "S1"], 2L)
  expect_equal(n$n_residence[n$site == "S2"], 3L)
  expect_equal(n$n_residence[n$site == "S3"], 2L)
  expect_true(all(n$mean_residence_h >= 0, na.rm = TRUE))
  recs <- attr(net, "transition_records")[["all"]]
  expect_equal(nrow(recs), 5L)                                # 3 (A) + 2 (B) transitions (1 h gaps, no split)
  expect_false(any(recs$crossed_gap))                         # nothing exceeds the 48 h default
})

test_that("cross-site gaps are flagged and excluded from mean transit time", {
  d <- data.frame(ID = factor("A"),
                  datetime = as.POSIXct("2023-01-01", tz = "UTC") + c(0, 1, 100, 101) * 3600,  # 99 h S1->S2
                  site = c("S1", "S1", "S2", "S2"), lon = c(-8, -8, -7, -7), lat = c(37, 37, 38, 38))
  md  <- as_moby(d, station.col = "site", lon.col = "lon", lat.col = "lat")
  net <- suppressMessages(calculateTransitions(md, spatial.col = "site"))
  rec <- attr(net, "transition_records")[["all"]]; e <- networkEdges(net)
  expect_true(rec$crossed_gap[rec$from == "S1" & rec$to == "S2"])
  expect_true(is.na(e$mean_duration_h[e$from == "S1" & e$to == "S2"]))     # transit time excluded
  expect_equal(e$n_movements[e$from == "S1" & e$to == "S2"], 1L)           # but the edge is kept

  # max.gap = Inf recovers the legacy behaviour: nothing flagged, transit time retained
  net2 <- suppressMessages(calculateTransitions(md, spatial.col = "site", max.gap = Inf))
  rec2 <- attr(net2, "transition_records")[["all"]]; e2 <- networkEdges(net2)
  expect_false(any(rec2$crossed_gap))
  expect_false(is.na(e2$mean_duration_h[e2$from == "S1" & e2$to == "S2"]))
})

test_that("gap-split revisits count as separate visits and emit no self-edge", {
  d <- data.frame(ID = factor("A"),
                  datetime = as.POSIXct("2023-01-01", tz = "UTC") + c(0, 1, 100, 101) * 3600,  # S1, 99 h, S1
                  site = c("S1", "S1", "S1", "S1"), lon = -8, lat = 37)
  md  <- as_moby(d, station.col = "site", lon.col = "lon", lat.col = "lat")
  net <- suppressMessages(calculateTransitions(md, spatial.col = "site"))
  n <- networkNodes(net); e <- networkEdges(net)
  expect_equal(n$n_residence[n$site == "S1"], 2L)             # two visits (gap-split)
  expect_equal(nrow(e[e$from == "S1" & e$to == "S1", ]), 0L)  # no S1->S1 self-edge
  net_inf <- suppressMessages(calculateTransitions(md, spatial.col = "site", max.gap = Inf))
  ninf <- networkNodes(net_inf)
  expect_equal(ninf$n_residence[ninf$site == "S1"], 1L)       # one continuous visit
})

test_that("transition counts agree with transitionsTable", {
  md <- trans_dataset()
  net <- calculateTransitions(md, spatial.col = "site")
  tt <- transitionsTable(net)
  expect_equal(sum(networkEdges(net)$n_movements), sum(as.numeric(tt$Movements)))
})

test_that("calculateTransitions supports id.groups", {
  md <- trans_dataset()
  net <- calculateTransitions(md, spatial.col = "site",
                              id.groups = list(grp1 = "A", grp2 = "B"))
  expect_setequal(unique(networkEdges(net)$group), c("grp1", "grp2"))
  expect_setequal(unique(networkNodes(net)$group), c("grp1", "grp2"))
})

test_that("plot.mobyNetwork draws a movement network", {
  skip_if_not_installed("igraph")
  md <- trans_dataset()
  net <- calculateTransitions(md, spatial.col = "site")
  tmp <- tempfile(fileext = ".pdf"); grDevices::pdf(tmp); on.exit({ grDevices::dev.off(); unlink(tmp) })
  expect_no_error(plot(net))
})

test_that("spatial.col defaults to the mobyData station column when omitted", {
  # a mobyData whose station column is 'station' -> calculateTransitions(md) needs no spatial.col
  data(rays)
  net_default  <- calculateTransitions(rays)                       # resolved from metadata
  net_explicit <- calculateTransitions(rays, spatial.col = mobyMeta(rays)$station.col)
  expect_true(is_mobyNetwork(net_default))
  expect_equal(networkNodes(net_default), networkNodes(net_explicit))
})
