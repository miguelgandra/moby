# Deterministic synthetic telemetry fixtures used across tests.

# Disable the automatic coastline fallback in plotMaps/plotMovements for the whole suite, so map tests
# stay fast and deterministic and do not depend on rnaturalearth/maps. The coastline path is exercised
# explicitly (and guarded) in test-plotMaps.R / test-plotMovements.R.
options(moby.coastline = FALSE)

# Detections stored in a NON-UTC timezone (UTC+10), for timezone-regression tests. The bundled `rays`
# data is UTC-only, so tz bugs (day/hour derived in UTC instead of the data's zone) are invisible on
# it; this fixture makes them catchable. Animal "A"'s three detections all fall on ONE LOCAL calendar
# day (2023-06-02) but span TWO UTC days; animal "B" spans two local days. So a min.days = 2 filter
# drops A and keeps B iff day-bucketing uses the DATA timezone (correct); UTC bucketing wrongly keeps A.
nonutc_detections <- function(tz = "Etc/GMT-10") {
  dt <- as.POSIXct(c("2023-06-02 08:00", "2023-06-02 09:00", "2023-06-02 11:00",
                     "2023-06-02 10:00", "2023-06-04 10:00", "2023-06-04 12:00"), tz = tz)
  data.frame(ID = factor(c("A", "A", "A", "B", "B", "B")), datetime = dt,
             lon = -8.0, lat = 37.0, stringsAsFactors = FALSE)
}

# Build a simple long-format detection data frame.
make_detections <- function(ids = c("A", "B"),
                            n_per_id = 24,
                            start = "2023-06-01 00:00:00",
                            tz = "UTC",
                            by = "1 hour",
                            stations = c("R1", "R2", "R3"),
                            seed = 1) {
  set.seed(seed)
  times <- seq(as.POSIXct(start, tz = tz), by = by, length.out = n_per_id)
  do.call(rbind, lapply(ids, function(id) {
    data.frame(
      ID = id,
      datetime = times,
      lon = -8.0,
      lat = 37.0,
      station = sample(stations, n_per_id, replace = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}
