vue_detections <- function() {
  data.frame(check.names = FALSE,
    "Date and Time (UTC)" = c("2023-06-01 00:00:00", "2023-06-01 01:00:00", "2023-06-01 02:00:00"),
    "Receiver" = c("VR2W-1001", "VR2W-1001", "VR2W-1002"),
    "Transmitter" = c("A69-1602-111", "A69-1602-111", "A69-1602-222"),
    "Transmitter Name" = c("shark1", "shark1", "shark2"),
    "Station Name" = c("st A", "st A", "st B"),
    "Latitude" = c(37.0, 37.0, 37.1), "Longitude" = c(-8.0, -8.0, -8.1))
}

test_that("importDetections harmonises a VUE export into a mobyData", {
  f <- tempfile(fileext = ".csv"); on.exit(unlink(f))
  write.csv(vue_detections(), f, row.names = FALSE)
  d <- importDetections(f, source = "vue")
  expect_true(is_moby(d))
  expect_true(all(c("ID", "datetime", "transmitter", "receiver", "station", "lon", "lat") %in% colnames(d)))
  expect_s3_class(d$datetime, "POSIXct")
  expect_true(is.numeric(d$lon))
  # no animal id in VUE -> initialised from transmitter
  expect_true(all(as.character(d$ID) %in% c("A69-1602-111", "A69-1602-222")))
})

test_that("importDetections handles GLATOS and ETN data frames", {
  gl <- data.frame(animal_id = c("fish1", "fish1"),
                   detection_timestamp_utc = c("2023-06-01 00:00:00", "2023-06-01 03:00:00"),
                   transmitter_codespace = c("A69-1602", "A69-1602"), transmitter_id = c("111", "111"),
                   receiver_sn = c("1001", "1001"), station = c("stA", "stA"),
                   deploy_lat = c(37, 37), deploy_long = c(-8, -8))
  d <- importDetections(gl, source = "glatos")
  expect_true("fish1" %in% as.character(d$ID))
  expect_true(all(d$transmitter == "A69-1602-111"))  # codespace + id combined

  etn <- data.frame(animal_id = "eel1", date_time = "2023-06-01 00:00:00",
                    acoustic_tag_id = "A69-1303-9", receiver_id = "R1", station_name = "S1",
                    deploy_latitude = 51, deploy_longitude = 3)
  de <- importDetections(etn, source = "etn")
  expect_equal(as.character(de$ID), "eel1")
  expect_equal(de$station, "S1")
})

test_that("importDetections supports generic col.map and missing optional fields", {
  raw <- data.frame(when = "2023-06-01 00:00:00", tag = "T1", rec = "R1")
  d <- importDetections(raw, source = "generic",
                        col.map = list(datetime = "when", transmitter = "tag", receiver = "rec"))
  expect_true(is_moby(d))
  expect_equal(as.character(d$ID), "T1")
  expect_false("lon" %in% colnames(d))  # absent optional field degrades gracefully
})

test_that("importDeployments produces the canonical deployment schema", {
  raw <- data.frame(check.names = FALSE,
    Receiver = "VR2W-1001", Station = "st A", Latitude = 37, Longitude = -8,
    Deploymentdate = "2023-01-01 00:00:00", Dateout = "2023-04-01 00:00:00")
  f <- tempfile(fileext = ".csv"); on.exit(unlink(f))
  write.csv(raw, f, row.names = FALSE)
  dep <- importDeployments(f, source = "vue")
  expect_true(all(c("receiver", "station", "lon", "lat", "deploy", "recover") %in% colnames(dep)))
  expect_s3_class(dep$deploy, "POSIXct")
})

test_that("checkDeployments flags overlapping deployments and bad coordinates/ranges", {
  dep <- data.frame(
    receiver = c("R1", "R1", "R2", "R3"),
    station  = c("A", "A", "B", "C"),
    lon = c(-8, -8, -8.1, 0), lat = c(37, 37, 37.1, 0),
    deploy = as.POSIXct(c("2023-01-01", "2023-03-01", "2023-01-01", "2023-05-01"), tz = "UTC"),
    recover = as.POSIXct(c("2023-04-01", "2023-06-01", "2023-12-01", "2023-04-01"), tz = "UTC"))
  qc <- checkDeployments(dep, verbose = FALSE)
  expect_s3_class(qc, "mobyQC")
  expect_true("Overlapping deployments" %in% names(qc$counts))   # R1 redeployed before recovery
  expect_true("Implausible coordinates" %in% names(qc$counts))   # R3 at 0,0
  expect_true("Invalid date range" %in% names(qc$counts))        # R3 recover < deploy
})

test_that("checkDeployments cross-checks detections against metadata", {
  f <- tempfile(fileext = ".csv"); on.exit(unlink(f))
  d <- vue_detections()
  d[4, ] <- list("2024-01-01 00:00:00", "VR2W-9999", "A69-1602-222", "shark2", "st X", 37.2, -8.2)
  write.csv(d, f, row.names = FALSE)
  det <- importDetections(f, source = "vue")
  dep <- data.frame(receiver = c("VR2W-1001", "VR2W-1002"), station = c("st A", "st B"),
                    lon = c(-8, -8.1), lat = c(37, 37.1),
                    deploy = as.POSIXct(c("2023-01-01", "2023-01-01"), tz = "UTC"),
                    recover = as.POSIXct(c("2023-12-01", "2023-12-01"), tz = "UTC"))
  qc <- checkDeployments(dep, detections = det, verbose = FALSE)
  expect_true("Receiver missing from metadata" %in% names(qc$counts))  # VR2W-9999
})

test_that("checkDeployments 'checks' selector runs only the requested groups", {
  dep <- data.frame(
    receiver = c("R1", "R1", "R3"), station = c("A", "A", "C"),
    lon = c(-8, -8, 0), lat = c(37, 37, 0),
    deploy = as.POSIXct(c("2023-01-01", "2023-03-01", "2023-05-01"), tz = "UTC"),
    recover = as.POSIXct(c("2023-04-01", "2023-06-01", "2023-04-01"), tz = "UTC"))
  types <- function(...) sort(unique(checkDeployments(dep, ..., verbose = FALSE)$report$type))
  expect_setequal(types(checks = "dates"), "Invalid date range")             # R3 recover<deploy only
  expect_setequal(types(checks = "coordinates"), "Implausible coordinates")  # R3 0,0 only
  expect_true("Overlapping deployments" %in% types(checks = "overlaps"))
  expect_false("Overlapping deployments" %in% types(checks = "coordinates")) # subset excludes it
  # 'all' (default) is the union of the groups
  expect_true(all(c("Invalid date range", "Implausible coordinates") %in% types()))
})

# a big land square (~55 km across) in WGS84, for the on-land coordinate check
qc_land_square <- function()
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(-9.0, 38.0), c(-8.5, 38.0), c(-8.5, 38.5), c(-9.0, 38.5), c(-9.0, 38.0)))), crs = 4326))

test_that("on-land check flags receivers on land, spares near-shore, and is opt-in", {
  skip_if_not_installed("sf")
  sq <- qc_land_square()
  dep <- data.frame(receiver = c("R1", "R2", "R3"),
                    station = c("deep_inland", "at_sea", "near_edge"),
                    lon = c(-8.75, -9.30, -8.503), lat = c(38.25, 38.25, 38.25),
                    deploy = as.POSIXct("2023-01-01", tz = "UTC"),
                    recover = as.POSIXct("2023-06-01", tz = "UTC"))

  # no land.shape -> the check simply does not run
  expect_false("Coordinates on land" %in% checkDeployments(dep, verbose = FALSE)$report$type)

  onland <- function(tol) {
    r <- checkDeployments(dep, land.shape = sq, land.tolerance = tol, verbose = FALSE)$report
    r$station[r$type == "Coordinates on land"]
  }
  expect_setequal(onland(500), "deep_inland")               # near_edge (262 m in) is spared at 500 m
  expect_setequal(onland(0), c("deep_inland", "near_edge")) # a bare intersect catches both
  expect_false("at_sea" %in% onland(0))                     # a point in water is never flagged

  # a deep-inland error is labelled as a likely metadata error; details name the deployment records
  det <- checkDeployments(dep, land.shape = sq, verbose = FALSE)$report
  msg <- det$details[det$type == "Coordinates on land"]
  expect_match(msg, "likely a metadata error")
  expect_match(msg, "deployment records")
})

test_that("on-land check dedups by (station, position) so a moved receiver is still caught", {
  skip_if_not_installed("sf")
  sq <- qc_land_square()
  # one station name, two deployments: first at sea, later moved onto land
  dep <- data.frame(receiver = "R1", station = c("S", "S"),
                    lon = c(-9.30, -8.75), lat = c(38.25, 38.25),
                    deploy = as.POSIXct(c("2023-01-01", "2023-04-01"), tz = "UTC"),
                    recover = as.POSIXct(c("2023-03-01", "2023-06-01"), tz = "UTC"))
  r <- checkDeployments(dep, land.shape = sq, verbose = FALSE)$report
  expect_equal(sum(r$type == "Coordinates on land"), 1L)    # the on-land position, not masked by the at-sea one
})

test_that("min.active.days flags short-effort stations (gap-aware, report-only, opt-in)", {
  # station A: two 15-day deployments 5 months apart -> 30 OPERATIONAL days (not the ~200-day span);
  # station B: one 300-day deployment.
  dep <- data.frame(
    receiver = c("R1", "R1", "R2"), station = c("A", "A", "B"),
    lon = c(-9, -9, -8.9), lat = c(38, 38, 38.1),
    deploy  = as.POSIXct(c("2023-01-01", "2023-06-20", "2023-01-01"), tz = "UTC"),
    recover = as.POSIXct(c("2023-01-16", "2023-07-05", "2023-11-01"), tz = "UTC"))
  det <- data.frame(receiver = c("R1", "R1", "R2"), ID = c("f1", "f2", "f1"),
                    datetime = as.POSIXct(c("2023-01-05", "2023-01-06", "2023-02-01"), tz = "UTC"),
                    station = c("A", "A", "B"))

  r <- checkDeployments(dep, detections = det, min.active.days = 180, verbose = FALSE)$report
  sd <- r[r$type == "Short monitoring duration", ]
  expect_setequal(sd$station, "A")               # gap-aware: A=30 operational days < 180; B=304 not flagged
  expect_match(sd$details, "30 day")             # summed windows, not the calendar span
  expect_equal(sd$n_detections, 2L)              # impact reported when detections supplied
  expect_equal(sd$n_individuals, 2L)

  # opt-in: no threshold -> no such rows; and nothing is ever removed (report-only)
  expect_false("Short monitoring duration" %in% checkDeployments(dep, verbose = FALSE)$report$type)

  # threshold below A's duration -> A no longer flagged
  expect_false("A" %in% checkDeployments(dep, min.active.days = 20, verbose = FALSE
                 )$report$station[checkDeployments(dep, min.active.days = 20, verbose = FALSE)$report$type == "Short monitoring duration"])

  expect_error(checkDeployments(dep, min.active.days = -5, verbose = FALSE), "positive number")
  expect_error(checkDeployments(dep, min.active.days = c(1, 2), verbose = FALSE), "positive number")
})

test_that("on-land check auto-resolves land.shape from the detections mobyData, and fails soft without a CRS", {
  skip_if_not_installed("sf")
  sq <- qc_land_square()
  dep <- data.frame(receiver = "R1", station = "deep_inland", lon = -8.75, lat = 38.25,
                    deploy = as.POSIXct("2023-01-01", tz = "UTC"),
                    recover = as.POSIXct("2023-06-01", tz = "UTC"))
  det <- suppressWarnings(as_moby(
    data.frame(receiver = "R1", ID = "f", datetime = as.POSIXct("2023-02-01", tz = "UTC"),
               station = "deep_inland", lon = -8.75, lat = 38.25),
    land.shape = sq, tagging.dates = as.POSIXct("2023-01-01", tz = "UTC")))
  # land.shape not passed explicitly -> borrowed from the detections metadata
  r <- checkDeployments(dep, detections = det, verbose = FALSE)$report
  expect_true("Coordinates on land" %in% r$type)

  # a land layer with no CRS -> skip with a message, never error, never a false row
  no_crs <- sf::st_set_crs(sq, NA)
  expect_message(checkDeployments(dep, land.shape = no_crs, verbose = FALSE), "could not run the on-land check")
  r2 <- suppressMessages(checkDeployments(dep, land.shape = no_crs, verbose = FALSE))$report
  expect_false("Coordinates on land" %in% r2$type)
})

test_that("deployment.lon.col/lat.col let a non-canonical deployment log be checked without renaming", {
  # same station, coordinates ~1.9 km apart, but columns named Longitude/Latitude/Deploy_Date
  raw <- data.frame(receiver = c("R1", "R2"), station = c("A", "A"),
                    Longitude = c(-8.00, -8.00), Latitude = c(37.000, 37.017),
                    Deploy_Date = as.POSIXct(c("2023-01-01", "2023-01-01"), tz = "UTC"))
  # canonical-named twin, to prove the two resolve identically
  canon <- data.frame(receiver = raw$receiver, station = raw$station,
                      lon = raw$Longitude, lat = raw$Latitude, deploy = raw$Deploy_Date)

  a <- checkDeployments(raw, deployment.deploy.col = "Deploy_Date", deployment.lon.col = "Longitude",
                        deployment.lat.col = "Latitude", verbose = FALSE)$report
  b <- checkDeployments(canon, verbose = FALSE)$report
  expect_equal(a, b)                                                        # identical to the canonical run
  expect_true("Inconsistent station coordinates" %in% a$type)              # coord check actually ran
  expect_match(a$details[a$type == "Inconsistent station coordinates"], "km")

  # an explicitly named coordinate column that is absent is reported, not silently ignored
  expect_message(
    checkDeployments(raw, deployment.deploy.col = "Deploy_Date", deployment.lon.col = "lon_wgs84",
                     deployment.lat.col = "Latitude", verbose = FALSE),
    "lon_wgs84.*not found")
})

test_that("scope = 'detected' restricts metadata checks to receivers present in the detections", {
  # R1 detected a tag; R2 never did. Both have a coverage gap in the log.
  dep <- data.frame(
    receiver = c("R1", "R1", "R2", "R2"),
    station  = c("A", "A", "B", "B"),
    lon = c(-8, -8, -8.1, -8.1), lat = c(37, 37, 37.1, 37.1),
    deploy  = as.POSIXct(c("2023-01-01", "2023-06-01", "2023-01-01", "2023-06-01"), tz = "UTC"),
    recover = as.POSIXct(c("2023-02-01", "2023-07-01", "2023-02-01", "2023-07-01"), tz = "UTC"))
  det <- data.frame(receiver = "R1", ID = "t1",
                    datetime = as.POSIXct("2023-01-10", tz = "UTC"), station = "A")

  full <- checkDeployments(dep, detections = det, verbose = FALSE)              # scope = "all" (default)
  scoped <- checkDeployments(dep, detections = det, scope = "detected", verbose = FALSE)

  gaps_full   <- full$report[full$report$type == "Coverage gap", ]
  gaps_scoped <- scoped$report[scoped$report$type == "Coverage gap", ]
  expect_true("R2" %in% gaps_full$receiver)         # full audit flags the undetected receiver
  expect_false("R2" %in% gaps_scoped$receiver)      # scoped audit drops it
  expect_true("R1" %in% gaps_scoped$receiver)       # ...but keeps the detected one
  expect_lt(nrow(scoped$report), nrow(full$report)) # scoping reduces, never adds
  # $deployments still holds the complete log regardless of scope
  expect_equal(nrow(scoped$deployments), nrow(dep))
})

test_that("scope = 'detected' without detections warns and audits everything", {
  dep <- data.frame(receiver = c("R1", "R1"), station = c("A", "A"),
                    lon = c(-8, -8), lat = c(37, 37),
                    deploy  = as.POSIXct(c("2023-01-01", "2023-06-01"), tz = "UTC"),
                    recover = as.POSIXct(c("2023-02-01", "2023-07-01"), tz = "UTC"))
  expect_message(checkDeployments(dep, scope = "detected", verbose = FALSE), "needs a 'detections'")
  a <- suppressMessages(checkDeployments(dep, scope = "detected", verbose = FALSE))
  b <- checkDeployments(dep, verbose = FALSE)
  expect_equal(a$report, b$report)                  # falls back to a full audit
})

test_that("inconsistent-coordinate messages report metres/kilometres, not 'units'", {
  # same station name, two coordinates ~1.9 km apart (geographic)
  dep <- data.frame(receiver = c("R1", "R2"), station = c("A", "A"),
                    lon = c(-8.00, -8.00), lat = c(37.000, 37.017),
                    deploy  = as.POSIXct(c("2023-01-01", "2023-01-01"), tz = "UTC"),
                    recover = as.POSIXct(c("2023-06-01", "2023-06-01"), tz = "UTC"))
  msg <- checkDeployments(dep, verbose = FALSE)$report
  incoord <- msg$details[msg$type == "Inconsistent station coordinates"]
  expect_length(incoord, 1)
  expect_match(incoord, "\\bkm\\b")                 # ~1.9 km, expressed in km
  expect_false(grepl("units", incoord))             # never the bare "units"
  # the .formatDistance helper picks m below 1 km and km above
  expect_equal(moby:::.formatDistance(320), "320 m")
  expect_equal(moby:::.formatDistance(15765), "15.8 km")
  expect_equal(moby:::.formatDistance(NA_real_), "an undetermined distance")
})

test_that("checkDeployments handles a clean log and a missing 'detections' request without error", {
  clean <- data.frame(receiver = "R1", station = "A", lon = -8, lat = 37,
                      deploy = as.POSIXct("2023-01-01", tz = "UTC"),
                      recover = as.POSIXct("2023-06-01", tz = "UTC"))
  qc <- checkDeployments(clean, verbose = FALSE)                 # empty report must not crash
  expect_s3_class(qc, "mobyQC")
  expect_equal(nrow(qc$report), 0)
  expect_no_error(print(qc))                                     # print handles the empty case
  expect_message(checkDeployments(clean, checks = "detections", verbose = FALSE), "no 'detections' supplied")
  expect_error(checkDeployments(clean, checks = "bogus", verbose = FALSE), "should be one of")
})

test_that("checkDeployments accepts non-canonical deploy/recover column names (deployment.deploy.col/recover.col)", {
  dep <- data.frame(receiver = c("R1", "R2"), station = c("A", "B"), lon = c(-9, -9), lat = c(38, 38),
                    deploy_date = as.POSIXct(c("2023-01-01", "2023-01-01"), tz = "UTC"),
                    recover_date = as.POSIXct(c("2023-06-01", "2023-06-01"), tz = "UTC"))
  expect_error(checkDeployments(dep, verbose = FALSE), "missing required")     # canonical 'deploy' absent by default
  qc <- suppressMessages(checkDeployments(dep, deployment.deploy.col = "deploy_date",
                                          deployment.recover.col = "recover_date", verbose = FALSE))
  expect_s3_class(qc, "mobyQC")
})

test_that("importTags harmonises the transmitter nominal delay (incl. min/max midpoint)", {
  # an explicit nominal delay is carried through as numeric
  tg <- rays_tags; tg$nominal_delay <- "90"
  out <- importTags(tg, source = "generic",
                    col.map = list(ID = "ID", transmitter = "transmitter",
                                   tagging_date = "tagging_date", nominal_delay = "nominal_delay"))
  expect_true(is.numeric(out$nominal_delay))
  expect_true(all(out$nominal_delay == 90))
  # many tag exports give a delay RANGE: the nominal delay is its midpoint
  tg2 <- rays_tags; tg2$min_delay <- 60; tg2$max_delay <- 180
  out2 <- suppressMessages(importTags(tg2, source = "generic",
                                      col.map = list(ID = "ID", transmitter = "transmitter",
                                                     min_delay = "min_delay", max_delay = "max_delay")))
  expect_true(all(out2$nominal_delay == 120))
})

test_that("assignAnimalIDs populates meta$nominal.delay, which filterDetections uses automatically", {
  tg <- rays_tags; tg$nominal_delay <- 120
  tags <- importTags(tg, source = "generic",
                     col.map = list(ID = "ID", transmitter = "transmitter",
                                    tagging_date = "tagging_date", nominal_delay = "nominal_delay"))
  det <- rays_detections[rays_detections$transmitter %in% rays_tags$transmitter, ]
  md <- suppressWarnings(suppressMessages(assignAnimalIDs(det, tags)))

  nd <- mobyMeta(md)$nominal.delay
  expect_false(is.null(nd))
  expect_true(all(nd == 120))
  expect_true(all(names(nd) %in% as.character(unique(md$ID))))

  # the whole point: filterDetections enables its min_lag filter with NO argument supplied,
  # at the documented 30 x nominal.delay threshold (30 * 120 = 3600 s)
  res <- suppressWarnings(suppressMessages(filterDetections(md)))
  expect_true(any(grepl("min_lag > 3600 s", res$data_discarded$reason)))

  # and it can be opted out of
  md_off <- suppressWarnings(suppressMessages(assignAnimalIDs(det, tags, set.nominal.delay = FALSE)))
  expect_null(mobyMeta(md_off)$nominal.delay)
})
