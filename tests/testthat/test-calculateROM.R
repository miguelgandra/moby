# Build a small tracks-style dataset (consistent coordinates and step distances).
# A moves north then returns to the origin (net displacement = 0); B makes a single
# straight-line step (net displacement = total distance).
mov_dataset <- function() {
  d0 <- data.frame(
    ID = c("A", "A", "A", "A", "B", "B"),
    timebin = as.POSIXct(c("2023-01-01 10:00", "2023-01-01 11:00", "2023-01-01 12:00",
                           "2023-01-01 13:00", "2023-01-01 10:00", "2023-01-01 11:00"), tz = "UTC"),
    lon = c(0, 0, 0, 0, 0, 0),
    lat = c(0, 0.001, 0.002, 0, 0, 0.003),
    stringsAsFactors = FALSE)
  # great-circle distances on geographic coordinates, no epsg.code required
  suppressWarnings(suppressMessages(
    calculateStepDistances(d0, id.col = "ID", lon.col = "lon", lat.col = "lat", verbose = FALSE)))
}

test_that("calculateStepDistances returns a distance-enriched data frame with trajectories as an attribute", {
  d <- mov_dataset()
  expect_s3_class(d, "data.frame")
  expect_true("dist_m" %in% colnames(d))
  # trajectories are exposed via the accessor, not a list element
  expect_null(d$trajectories)
  tr <- getTrajectories(d)
  expect_false(is.null(tr))
  expect_true(all(unique(as.character(d$ID)) %in% names(tr)))
})

test_that("calculateROM aggregates step distances into per-individual metrics", {
  d <- mov_dataset()
  m <- suppressWarnings(suppressMessages(
    calculateROM(d, id.col = "ID", timebin.col = "timebin", dist.col = "dist_m")))

  expect_true(is.data.frame(m))
  expect_true(all(c("ID", "n_steps", "total_distance_m", "mean_rom", "max_rom") %in% colnames(m)))
  expect_equal(attr(m, "interval"), 60)

  distA <- d$dist_m[d$ID == "A"]
  expect_equal(m$n_steps[m$ID == "A"], sum(!is.na(distA)))
  expect_equal(m$total_distance_m[m$ID == "A"], sum(distA, na.rm = TRUE))
  # interval is 60 min, so ROM (per hour) equals the per-step distance unchanged
  expect_equal(m$mean_rom[m$ID == "A"], mean(distA, na.rm = TRUE))
  expect_equal(m$max_rom[m$ID == "A"], max(distA, na.rm = TRUE))

  # B makes a single step
  expect_equal(m$n_steps[m$ID == "B"], 1L)
})

test_that("calculateROM errors when the distance column is missing", {
  d <- mov_dataset()
  d$dist_m <- NULL
  expect_error(
    calculateROM(d, id.col = "ID", timebin.col = "timebin", dist.col = "dist_m"),
    "distance column")
})

test_that("calculateLinearityIndex measures movement directness", {
  d <- mov_dataset()
  li <- suppressWarnings(suppressMessages(
    calculateLinearityIndex(d, id.col = "ID", timebin.col = "timebin",
                            lon.col = "lon", lat.col = "lat", dist.col = "dist_m")))

  expect_true(is.data.frame(li))
  expect_true(all(c("ID", "net_distance_m", "total_distance_m", "linearity_index") %in% colnames(li)))
  # A returns to the origin -> net displacement 0 -> linearity index 0
  expect_equal(li$linearity_index[li$ID == "A"], 0)
  # B moves in a straight line (single step) -> net == total -> linearity index 1
  expect_equal(li$linearity_index[li$ID == "B"], 1, tolerance = 1e-6)
})

test_that("movementTable ROM values match calculateROM (single source of truth)", {
  d <- mov_dataset()
  m <- suppressWarnings(suppressMessages(
    calculateROM(d, id.col = "ID", timebin.col = "timebin", dist.col = "dist_m")))

  kud_mock <- list(summary_table = data.frame(ID = c("A", "B"),
                                              "UD 95 (km2)" = c(10, 20), check.names = FALSE))
  attr(kud_mock, "bandwidth") <- 1

  mt <- suppressWarnings(suppressMessages(
    movementTable(d, ud.results = kud_mock, id.col = "ID", timebin.col = "timebin",
                  lon.col = "lon", lat.col = "lat", dist.col = "dist_m")))

  rom_col <- grep("^ROM", colnames(mt), value = TRUE)
  expect_length(rom_col, 1)
  expect_equal(mt[mt$ID == "A", rom_col], sprintf("%.1f", m$mean_rom[m$ID == "A"]))
  expect_equal(mt[mt$ID == "B", rom_col], sprintf("%.1f", m$mean_rom[m$ID == "B"]))
})
