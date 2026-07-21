# Tests for the refactored filterDetections(): opt-in filters, the per-receiver min_lag false-detection
# filter, the corroboration-shielded speed filter, duplicate removal, per-ID parameters, and the
# mobyFilter return object. Speed fixtures place a far-off station ~17.8 km east so an A->B step in a
# few minutes is impossible while within-station steps are slow.

tg <- as.POSIXct(c(F = "2023-01-01"), tz = "UTC")
mk <- function(id, hours, lon, lat, station)
  data.frame(ID = id, datetime = as.POSIXct("2023-06-01 00:00", tz = "UTC") + hours * 3600,
             lon = lon, lat = lat, station = station, stringsAsFactors = FALSE)
Alon <- -8.0; Blon <- -7.8

test_that("a single detection before the tagging date is removed (off-by-one regression)", {
  tag <- as.POSIXct("2020-01-10", tz = "UTC")
  df <- data.frame(ID = rep(c("A", "B"), each = 5),
                   datetime = rep(as.POSIXct("2020-01-10 12:00:00", tz = "UTC") +
                                    c(-86400, 0, 3600, 7200, 10800), 2), lon = -8, lat = 37)
  res <- suppressWarnings(suppressMessages(filterDetections(df, tagging.dates = tag)))
  expect_equal(nrow(res$data), 8L)
  expect_equal(nrow(res$data_discarded), 2L)
})

test_that("filterDetections returns a printable mobyFilter with qc_flag and a parameters attribute", {
  data(rays)
  res <- suppressWarnings(suppressMessages(filterDetections(rays)))
  expect_s3_class(res, "mobyFilter")
  expect_s3_class(res$data, "mobyData")
  expect_true("qc_flag" %in% names(res$data))
  expect_true(all(res$data$qc_flag %in% c("valid", "overspeed_review")))
  expect_false(is.null(attr(res, "parameters")))
  expect_output(print(res), "mobyFilter")
  # default is conservative: nothing but the temporal bounds runs (isolation + min_lag are off)
  expect_true("reason" %in% names(res$data_discarded))
})

test_that("the per-individual summary total counts pre-tagging and cut-off removals", {
  tag <- as.POSIXct("2020-01-10", tz = "UTC"); cutoff <- as.POSIXct("2020-01-10 23:00:00", tz = "UTC")
  df <- data.frame(ID = rep("A", 4),
                   datetime = as.POSIXct(c("2020-01-09 12:00:00", "2020-01-10 12:00:00",
                                           "2020-01-10 13:00:00", "2020-01-11 12:00:00"), tz = "UTC"),
                   lon = -8, lat = 37)
  res <- suppressWarnings(suppressMessages(filterDetections(df, tagging.dates = tag, cutoff.dates = cutoff)))
  expect_equal(nrow(res$data_discarded), 2L)
})

test_that("duplicate removal drops exact duplicates and is toggleable", {
  d <- mk("F", c(0, 1, 2), Alon, 37, "SA"); d <- rbind(d, d[3, ]); d$ID <- factor(d$ID)
  r <- suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg)))
  expect_equal(nrow(r$data_discarded), 1L)
  expect_equal(r$data_discarded$reason, "duplicate detection")
  r0 <- suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg, remove.duplicates = FALSE)))
  expect_equal(nrow(r0$data_discarded), 0L)
})

test_that("min_lag removes an uncorroborated lone decode and is off without nominal.delay", {
  d <- rbind(
    data.frame(ID = "F", datetime = as.POSIXct("2023-06-01 00:00", tz = "UTC") + (0:5) * 60,
               lon = Alon, lat = 37, station = "SA"),                       # corroborated burst at SA
    data.frame(ID = "F", datetime = as.POSIXct("2023-06-01 01:00", tz = "UTC"),
               lon = Blon, lat = 37, station = "SB"))                       # lone decode at SB
  d$ID <- factor(d$ID)
  r <- suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg, nominal.delay = 120)))
  expect_equal(nrow(r$data), 6L)
  expect_true(all(r$data$station == "SA"))
  expect_match(r$data_discarded$reason, "min_lag")
  # off without a nominal delay
  expect_equal(nrow(suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg)))$data_discarded), 0L)
  # resolved from mobyData metadata when present
  md <- as_moby(d, tagging.dates = tg); attr(md, "moby")$nominal.delay <- 120
  expect_equal(nrow(suppressWarnings(suppressMessages(filterDetections(md)))$data_discarded), 1L)
})

test_that("the isolation filter is off by default and removes both-sided isolated detections when set", {
  d <- data.frame(ID = "F", datetime = as.POSIXct("2023-06-01", tz = "UTC") + c(0, 1, 49, 97, 98) * 3600,
                  lon = Alon, lat = 37, station = "SA"); d$ID <- factor(d$ID)
  expect_equal(nrow(suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg)))$data_discarded), 0L)
  r <- suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg, isolation.window = 24)))
  expect_equal(nrow(r$data_discarded), 1L)                                  # the isolated middle detection
  # FALSE is an accepted synonym for the off-switch
  expect_equal(nrow(suppressWarnings(suppressMessages(
    filterDetections(d, tagging.dates = tg, isolation.window = FALSE)))$data_discarded), 0L)
})

test_that("speed filter removes an isolated spatial outlier (spike)", {
  spike <- rbind(mk("F", c(0, 1, 2, 3, 4), Alon, 37, "SA"),
                 mk("F", 4 + 10 / 60, Blon, 37, "SB"),        # impossible A->B jump
                 mk("F", 5, Alon, 37, "SA"))
  spike$ID <- factor(spike$ID)
  r <- suppressWarnings(suppressMessages(filterDetections(spike, tagging.dates = tg, max.speed = 5)))
  expect_equal(nrow(r$data_discarded), 1L)
  expect_equal(r$data_discarded$station, "SB")                # the outlier, not a neighbour
  expect_match(r$data_discarded$reason, "max speed")
})

test_that("corroboration shield: a multi-detection relocation is flagged, not cascaded away", {
  casc <- rbind(mk("F", c(0, 1, 2, 3, 4), Alon, 37, "SA"),
                mk("F", c(4 + 10 / 60, 4.5, 5.0), Blon, 37, "SB"))          # 3 corroborating B detections
  casc$ID <- factor(casc$ID)
  r <- suppressWarnings(suppressMessages(filterDetections(casc, tagging.dates = tg, max.speed = 5)))
  expect_equal(nrow(r$data), 8L)                                            # nothing deleted
  expect_equal(nrow(r$data_discarded), 0L)
  expect_equal(sum(r$data$qc_flag == "overspeed_review"), 3L)               # the B group flagged for review
  expect_true(all(r$data$qc_flag[r$data$station == "SB"] == "overspeed_review"))
  # min.corroboration is monotonic: raising it above the cluster size lifts the shield (deletes)...
  r4 <- suppressWarnings(suppressMessages(filterDetections(casc, tagging.dates = tg, max.speed = 5, min.corroboration = 4L)))
  expect_true(nrow(r4$data_discarded) > 0)
  # ...and min.corroboration = 1 shields everything (flag-only, deletes nothing)
  r1 <- suppressWarnings(suppressMessages(filterDetections(casc, tagging.dates = tg, max.speed = 5, min.corroboration = 1L)))
  expect_equal(nrow(r1$data_discarded), 0L)
})

test_that("shield wins over removal on a multi-leg transit (no split-cluster deletion)", {
  # A A A -> B B -> C C, all far apart: two overspeed legs. A per-window singleton (a co-located
  # detection whose twin is outside that +/-3 window) must NOT be deleted once any window shields it.
  d <- data.frame(ID = "F", datetime = as.POSIXct("2023-06-01 00:00", tz = "UTC") + (0:6) * 60, lon = -8.0,
                  lat = c(37.0, 37.0, 37.0, 37.36, 37.36, 37.54, 37.54),
                  station = c("A", "A", "A", "B", "B", "C", "C"), stringsAsFactors = FALSE)
  d$ID <- factor(d$ID)
  r <- suppressWarnings(suppressMessages(
    filterDetections(d, tagging.dates = tg, max.speed = 1, speed.unit = "km/h",
                     acoustic.range = 0, remove.duplicates = FALSE)))
  expect_equal(nrow(r$data), 7L)                 # nothing deleted
  expect_equal(nrow(r$data_discarded), 0L)
  expect_false("iteration" %in% names(r$data_discarded))   # internal column does not leak
})

test_that("max.speed and nominal.delay accept per-ID named vectors", {
  casc <- rbind(mk("F", c(0, 1, 2, 3, 4), Alon, 37, "SA"), mk("F", c(4 + 10 / 60, 4.5, 5.0), Blon, 37, "SB"))
  casc$ID <- factor(casc$ID)
  expect_no_error(suppressWarnings(suppressMessages(
    filterDetections(casc, tagging.dates = tg, max.speed = c(F = 5), speed.unit = "m/s"))))
})

test_that("min.detections and min.days drop under-sampled individuals", {
  d <- mk("F", c(0, 1, 2), Alon, 37, "SA"); d$ID <- factor(d$ID)
  r_det <- suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg, min.detections = 10)))
  expect_equal(nrow(r_det$data), 0L)
  expect_match(r_det$data_discarded$reason[1], "fewer than 10 detections")
  r_day <- suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg, min.days = 5)))
  expect_equal(nrow(r_day$data), 0L)
  expect_match(r_day$data_discarded$reason[1], "fewer than 5 days")        # note: the word 'days' is present
})

test_that("input validation: NA datetimes error, and 0 no longer silently disables a filter", {
  d <- mk("F", c(0, 1, 2), Alon, 37, "SA"); d$ID <- factor(d$ID); d$datetime[2] <- NA
  expect_error(suppressWarnings(suppressMessages(filterDetections(d, tagging.dates = tg))), "missing value")
  d2 <- data.frame(ID = "F", datetime = as.POSIXct("2023-06-01", tz = "UTC") + c(0, 48, 96) * 3600,
                   lon = Alon, lat = 37, station = "SA"); d2$ID <- factor(d2$ID)
  # isolation.window = 0 is a real (zero-hour) threshold, not an accidental off-switch
  expect_no_error(suppressWarnings(suppressMessages(filterDetections(d2, tagging.dates = tg, isolation.window = 0))))
})
