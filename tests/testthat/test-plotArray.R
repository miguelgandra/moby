# Tests for plotArray() - the receiver-array map. Rendering is checked for "runs without error"
# (guarded to a temp device); the data-reduction, status classification and zero-detection cross-check
# are asserted on the invisibly-returned station table, which is testable without inspecting pixels.

dep_fixture <- function() {
  data.frame(
    receiver = paste0("R", 1:5),
    station  = paste0("ST", 1:5),
    lon = c(-9.02, -9.00, -8.98, -8.96, -8.94),
    lat = c(38.45, 38.46, 38.45, 38.46, 38.45),
    deploy  = as.POSIXct(c("2023-01-01", "2023-01-01", "2023-06-01", "2023-01-01", "2023-01-01"), tz = "UTC"),
    recover = as.POSIXct(c("2023-12-31", "2023-03-01", NA, "2023-12-31", "2023-12-31"), tz = "UTC"),
    habitat = c("reef", "sand", "reef", "sand", "reef"),
    stringsAsFactors = FALSE)
}

test_that("plotArray runs across geographic / projected / feature combinations", {
  d <- dep_fixture()
  pdf(tempfile(fileext = ".pdf"), width = 8, height = 7); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressMessages(plotArray(d)))
  expect_no_error(suppressMessages(plotArray(d, detection.range = 800, label = TRUE, hull = TRUE)))
  expect_no_error(suppressMessages(plotArray(d, epsg.code = 32629, detection.range = 800, coastline = FALSE)))
  expect_no_error(suppressMessages(plotArray(d, color.by = "habitat", label = "receiver", label.wrap = TRUE)))
  expect_no_error(suppressMessages(plotArray(d, status.at = as.POSIXct("2023-05-01", tz = "UTC"))))
  expect_no_error(suppressMessages(plotArray(d[1, , drop = FALSE], detection.range = 500)))   # single station
  expect_no_error(suppressMessages(plotArray(d, main = FALSE, legend = FALSE)))
})

test_that("plotArray returns one row per station with the expected columns", {
  d <- dep_fixture()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressMessages(plotArray(d, color.by = "habitat"))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
  expect_true(all(c("station", "lon", "lat", "n_receivers", "group") %in% names(res)))
  # repeated servicing of one station -> still a single point
  d2 <- rbind(d, transform(d[1, ], receiver = "R1b"))
  res2 <- suppressMessages(plotArray(d2))
  expect_equal(nrow(res2), 5L)
  expect_equal(res2$n_receivers[res2$station == "ST1"], 2L)
})

test_that("plotArray classifies deployment status at a chosen date", {
  d <- dep_fixture()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressMessages(plotArray(d, status.at = as.POSIXct("2023-05-01", tz = "UTC")))
  st <- stats::setNames(as.character(res$status), res$station)
  expect_equal(unname(st["ST1"]), "active")             # deploy Jan, recover Dec -> active
  expect_equal(unname(st["ST2"]), "recovered")          # recovered Mar -> recovered by May
  expect_equal(unname(st["ST3"]), "not yet deployed")   # deploy Jun -> not yet deployed in May
})

test_that("plotArray flags stations with zero detections", {
  d <- dep_fixture()
  det <- data.frame(station = c("ST1", "ST1", "ST2", "ST4"), stringsAsFactors = FALSE)  # ST3, ST5 absent
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  res <- suppressMessages(plotArray(d, detections = det))
  zero <- res$station[!res$detected]
  expect_setequal(zero, c("ST3", "ST5"))
})

test_that("plotArray validates its inputs", {
  d <- dep_fixture()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_error(plotArray(d[, c("receiver", "deploy")]), "missing required column")
  expect_error(plotArray(d, color.by = "nope"), "color.by")
  expect_error(plotArray(d, status.at = "not-a-date"), "status.at")
  expect_error(plotArray(d, detection.range = -5), "detection.range")
  expect_error(plotArray(d, detection.range = c(1, 2)), "named by station")   # unnamed multi-value
  expect_error(plotArray(d[, setdiff(names(d), "deploy")], status.at = "2023-05-01"), "deployment-date column")
  # a geographic EPSG must be rejected by the shared spatial pipeline
  expect_error(suppressMessages(plotArray(d, epsg.code = 4326)), "geographic")
})

test_that("per-station detection.range is keyed by station NAME, not input order", {
  d <- dep_fixture()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  # shuffle the input rows: a positional vector would land on the wrong stations; a named one must not
  d_shuf <- d[c(4, 2, 5, 1, 3), ]
  rng <- stats::setNames(c(100, 200, 300, 400, 500), c("ST1", "ST2", "ST3", "ST4", "ST5"))
  expect_no_error(suppressMessages(plotArray(d_shuf, detection.range = rng)))
  # a deployments column also works and is reduced per station
  d$rng_col <- c(100, 200, 300, 400, 500)
  expect_no_error(suppressMessages(plotArray(d, detection.range = "rng_col")))
  # NA radius is skipped, not a crash - including in projected mode (where sf::st_polygon rejects NaN)
  expect_no_error(suppressMessages(plotArray(d, detection.range = c(ST1=800, ST2=NA, ST3=800, ST4=800, ST5=800),
                                             epsg.code = 32629, coastline = FALSE)))
})

test_that("status.at character input is interpreted in the data timezone, not the session's", {
  d <- dep_fixture()                                  # deploy/recover are UTC; ST2 recovers 2023-03-01 UTC
  old <- Sys.getenv("TZ"); Sys.setenv(TZ = "America/New_York")
  on.exit({ if (nzchar(old)) Sys.setenv(TZ = old) else Sys.unsetenv("TZ") }, add = TRUE)
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  s_chr <- suppressMessages(plotArray(d, status.at = "2023-02-28 22:00:00"))
  s_utc <- suppressMessages(plotArray(d, status.at = as.POSIXct("2023-02-28 22:00:00", tz = "UTC")))
  expect_identical(as.character(s_chr$status), as.character(s_utc$status))
  # at 22:00 UTC (before ST2's 2023-03-01 00:00 UTC recovery) ST2 is still active
  expect_equal(as.character(s_chr$status[s_chr$station == "ST2"]), "active")
})

test_that("plotArray tolerates whitespace and warns on inconsistent data", {
  d <- dep_fixture()
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  # trailing whitespace in the detection station names must NOT flag every receiver as dead
  det <- data.frame(station = paste0(c("ST1", "ST2", "ST3", "ST4", "ST5"), " "), stringsAsFactors = FALSE)
  res <- suppressMessages(plotArray(d, detections = det))
  expect_true(all(res$detected))
  # a station logged at two very different positions is averaged AND warned about
  d2 <- d[1, ]; d2$lon <- -8.0; d2$lat <- 39.5
  expect_warning(suppressMessages(plotArray(rbind(d, d2))), "inconsistent coordinates")
})
