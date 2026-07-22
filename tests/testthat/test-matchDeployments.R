test_that("matchDeployments matches windows, back-fills metadata and flags mismatches", {
  det <- as_moby(data.frame(
    ID = factor(c("A", "A", "A", "A")),
    datetime = as.POSIXct(c("2023-02-01", "2023-05-01", "2023-08-01", "2024-06-01"), tz = "UTC"),
    receiver = c("R1", "R1", "R1", "R9"),
    station = c(NA, NA, NA, NA),
    lon = c(NA, -8.0, -7.0, NA),   # row 3 (-7) disagrees with metadata (-8.1)
    lat = c(NA, 37.0, 37.0, NA)),
    epsg.code = 32629)
  # R1 redeployed: stA (Jan-Jun) then stB (Jul-Dec)
  dep <- data.frame(receiver = c("R1", "R1"), station = c("stA", "stB"),
                    lon = c(-8.0, -8.1), lat = c(37.0, 37.1),
                    deploy = as.POSIXct(c("2023-01-01", "2023-07-01"), tz = "UTC"),
                    recover = as.POSIXct(c("2023-06-30", "2023-12-31"), tz = "UTC"))

  res <- matchDeployments(det, dep, coord.tolerance = 500, verbose = FALSE)
  expect_true(is_moby(res))
  expect_true(all(c("deployment_matched", "coord_mismatch") %in% colnames(res)))

  # R9 detection in 2024 falls outside every window
  expect_false(res$deployment_matched[res$receiver == "R9"])
  expect_true(all(res$deployment_matched[res$receiver == "R1"]))

  feb <- res$datetime == as.POSIXct("2023-02-01", tz = "UTC")
  aug <- res$datetime == as.POSIXct("2023-08-01", tz = "UTC")
  # station back-filled from the correct (time-resolved) deployment window
  expect_equal(res$station[feb], "stA")
  expect_equal(res$station[aug], "stB")
  # missing coordinates back-filled from metadata
  expect_equal(res$lon[feb], -8.0)
  # coordinate disagreement flagged
  expect_true(res$coord_mismatch[aug])
  # original metadata preserved
  expect_equal(mobyMeta(res)$epsg.code, 32629)
})

test_that("matchDeployments can drop unmatched detections", {
  det <- as_moby(data.frame(ID = factor(c("A", "A")),
                            datetime = as.POSIXct(c("2023-02-01", "2024-06-01"), tz = "UTC"),
                            receiver = c("R1", "R9"), station = NA_character_,
                            lon = NA_real_, lat = NA_real_))
  dep <- data.frame(receiver = "R1", station = "stA", lon = -8, lat = 37,
                    deploy = as.POSIXct("2023-01-01", tz = "UTC"),
                    recover = as.POSIXct("2023-12-31", tz = "UTC"))
  res <- matchDeployments(det, dep, drop.unmatched = TRUE, verbose = FALSE)
  expect_equal(nrow(res), 1L)
  expect_true(all(res$deployment_matched))
})

test_that("detection and deployment columns disambiguate: same concept, different names, both honoured", {
  # detection coords live in x/y; deployment coords live in Longitude/Latitude - the namespaced
  # arguments must resolve each against its own dataset without collision.
  dep <- data.frame(receiver = "R1", Site = "A", Longitude = -9, Latitude = 38,
                    deploy_date = as.POSIXct("2023-01-01", tz = "UTC"),
                    recover_date = as.POSIXct("2023-06-01", tz = "UTC"))
  det <- data.frame(receiver = "R1", site = "A", x = NA_real_, y = NA_real_,
                    when = as.POSIXct("2023-02-01", tz = "UTC"), ID = "f1")
  res <- suppressWarnings(suppressMessages(matchDeployments(
    det, dep,
    datetime.col = "when", station.col = "site", lon.col = "x", lat.col = "y",   # detection columns
    deployment.station.col = "Site", deployment.lon.col = "Longitude", deployment.lat.col = "Latitude",
    deployment.deploy.col = "deploy_date", deployment.recover.col = "recover_date", verbose = FALSE)))
  expect_true(res$deployment_matched[1])
  expect_equal(res$x[1], -9)     # detection lon (x) back-filled from deployment lon (Longitude)
  expect_equal(res$y[1], 38)
})

test_that("matchDeployments accepts non-canonical deploy/recover column names (deployment.deploy.col/recover.col)", {
  dep <- data.frame(receiver = "R1", station = "A", lon = -9, lat = 38,
                    deploy_date = as.POSIXct("2023-01-01", tz = "UTC"),
                    recover_date = as.POSIXct("2023-06-01", tz = "UTC"))
  det <- data.frame(receiver = "R1", station = "A", datetime = as.POSIXct("2023-02-01", tz = "UTC"), ID = "f1")
  expect_error(matchDeployments(det, dep, datetime.col = "datetime", verbose = FALSE), "missing required")
  res <- suppressWarnings(suppressMessages(matchDeployments(
    det, dep, deployment.deploy.col = "deploy_date", deployment.recover.col = "recover_date",
    datetime.col = "datetime", verbose = FALSE)))
  expect_equal(nrow(res), 1L)
})
