# Residence events (visits). Segmentation is deterministic (no model fitting), so assertions are exact.

visit_data <- function(times_h, sites, id = "A", tz = "UTC") {
  d <- data.frame(ID = factor(id),
                  datetime = as.POSIXct("2023-01-01", tz = tz) + times_h * 3600,
                  site = sites, lon = -8, lat = 37)
  as_moby(d, station.col = "site", lon.col = "lon", lat.col = "lat")
}

test_that("visits table has the expected schema and valid values", {
  md <- visit_data(c(0, 1, 2), c("S1", "S1", "S2"))
  v  <- suppressMessages(calculateVisits(md, spatial.col = "site"))
  expect_s3_class(v, "data.frame")
  expect_equal(names(v), c("group", "id", "site", "arrival", "departure",
                           "n_detections", "residence_h"))
  expect_true(all(v$residence_h >= 0))
  expect_true(all(v$arrival <= v$departure))
  expect_equal(nrow(v), 2L)                      # S1 (2 detections) then S2 (1 detection)
  expect_equal(v$n_detections, c(2L, 1L))
  expect_equal(v$residence_h, c(1, 0))           # S1 spans 1 h; S2 is a single detection
})

test_that("a gap longer than max.gap splits one stay into two visits", {
  md <- visit_data(c(0, 1, 100, 101), rep("S1", 4))            # same site, 99 h gap in the middle
  expect_equal(nrow(suppressMessages(calculateVisits(md, spatial.col = "site"))), 2L)                 # default 48 h
  expect_equal(nrow(suppressMessages(calculateVisits(md, spatial.col = "site", max.gap = 200))), 1L)  # tolerant
  expect_equal(nrow(suppressMessages(calculateVisits(md, spatial.col = "site", max.gap = Inf))), 1L)  # disabled
})

test_that("max.gap.units are honoured", {
  md <- visit_data(c(0, 1, 30, 31), rep("S1", 4))              # 29 h gap
  expect_equal(nrow(suppressMessages(calculateVisits(md, spatial.col = "site",
                    max.gap = 1, max.gap.units = "days"))), 2L)          # 24 h < 29 h -> split
  expect_equal(nrow(suppressMessages(calculateVisits(md, spatial.col = "site",
                    max.gap = 2, max.gap.units = "days"))), 1L)          # 48 h > 29 h -> one visit
})

test_that("visits are computed within id.groups", {
  d <- data.frame(ID = factor(c("A", "A", "B", "B")),
                  datetime = as.POSIXct("2023-01-01", tz = "UTC") + c(0, 1, 0, 1) * 3600,
                  site = c("S1", "S2", "S3", "S3"), lon = -8, lat = 37)
  md <- as_moby(d, station.col = "site", lon.col = "lon", lat.col = "lat")
  v  <- suppressMessages(calculateVisits(md, spatial.col = "site",
          id.groups = list(g1 = "A", g2 = "B")))
  expect_setequal(unique(v$group), c("g1", "g2"))
  expect_equal(nrow(v[v$id == "A", ]), 2L)       # S1 then S2
  expect_equal(nrow(v[v$id == "B", ]), 1L)       # S3 (2 detections, one visit)
})

test_that("bad max.gap is rejected", {
  md <- visit_data(c(0, 1), c("S1", "S1"))
  expect_error(suppressMessages(calculateVisits(md, spatial.col = "site", max.gap = -1)), "max.gap")
  expect_error(suppressMessages(calculateVisits(md, spatial.col = "site", max.gap = c(1, 2))), "max.gap")
})
