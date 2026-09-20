# The package's class convention, pinned in one place:
#
#   as_moby()  CONSTRUCTS  - the only function that turns a plain data frame into a mobyData
#   operations PRESERVE    - mobyData in, mobyData out; plain in, plain out
#   demotion   STRIPS BOTH - class and metadata are two halves of one fact
#
# The third rule is the one that used to be broken: as.data.frame() drops oldClass() but keeps other
# attributes, and every metadata read in moby goes to attr(x, "moby") without checking the class, so
# a "plain" data frame could still be silently driving column resolution.

test_that("demotion removes the metadata, not just the class", {
  for (demote in list(as.data.frame, function(x) head(x, 3), function(x) tail(x, 3))) {
    out <- demote(rays)
    expect_false(is_moby(out))
    expect_null(mobyMeta(out))
    expect_null(attr(out, "moby.land.name"))
    expect_identical(class(out), "data.frame")
  }
})

test_that("a demoted frame no longer drives column resolution", {
  d <- as.data.frame(rays)
  # rays maps id.col to "ID"; once demoted, resolution must fall back to the canonical defaults
  # rather than reading metadata off an object that claims to be plain
  expect_null(attr(d, "moby"))
  expect_identical(.resolveArgs(d, list(id.col = NULL))$id.col, unname(.mobyDefaults[["id.col"]]))
})

test_that("'[' preserves both class and metadata", {
  sub <- rays[1:50, ]
  expect_true(is_moby(sub))
  expect_identical(mobyMeta(sub)$epsg.code, mobyMeta(rays)$epsg.code)
  expect_length(mobyMeta(sub)$tagging.dates, length(mobyMeta(rays)$tagging.dates))
})

test_that(".restoreClass is a no-op when the input carried no metadata", {
  d <- as.data.frame(rays)
  expect_identical(.restoreClass(d, NULL), d)
  expect_false(is_moby(.restoreClass(d, NULL)))
})

test_that("operations preserve the input's class in both directions", {
  plain <- as.data.frame(rays)
  tag_dates <- mobyMeta(rays)$tagging.dates

  ops <- list(
    matchDeployments = list(
      moby  = function() matchDeployments(rays, rays_deployments, station.col = "station", verbose = FALSE),
      plain = function() matchDeployments(plain, rays_deployments, station.col = "station", verbose = FALSE)),
    calculateCOAs = list(
      moby  = function() calculateCOAs(rays, verbose = FALSE),
      plain = function() calculateCOAs(plain, verbose = FALSE)),
    filterDetections = list(
      moby  = function() filterDetections(rays, verbose = FALSE)$data,
      plain = function() filterDetections(plain, tagging.dates = tag_dates, verbose = FALSE)$data)
  )

  for (nm in names(ops)) {
    expect_true(is_moby(ops[[nm]]$moby()), info = paste(nm, "- mobyData in should give mobyData out"))
    out <- ops[[nm]]$plain()
    expect_false(is_moby(out), info = paste(nm, "- plain in should give plain out"))
    expect_null(mobyMeta(out), info = paste(nm, "- plain out must carry no metadata"))
  }
})

test_that("correctPositions preserves the dataset class and metadata in both return paths", {
  fx <- readRDS(test_path("_spatial", "fixtures.rds"))
  coast <- fx$land_sf
  tags <- setNames(as.POSIXct("2023-04-01", tz = "UTC"), "A01")
  make_moby <- function(rows) {
    d <- fx$onland[rows, , drop = FALSE]
    names(d)[names(d) == "ID"] <- "transmitter"
    names(d)[names(d) == "datetime"] <- "when"
    names(d)[names(d) == "lon"] <- "longitude"
    names(d)[names(d) == "lat"] <- "latitude"
    d$timebin <- d$when
    suppressWarnings(as_moby(d, id.col = "transmitter", datetime.col = "when",
                             timebin.col = "timebin", lon.col = "longitude",
                             lat.col = "latitude", tagging.dates = tags,
                             epsg.code = fx$epsg, land.shape = coast, verbose = FALSE))
  }

  for (rows in list(seq_len(nrow(fx$onland)), c(2L, 4L, 6L))) {
    md <- make_moby(rows)
    out <- suppressWarnings(correctPositions(md, spatial.layer = coast, verbose = FALSE))

    expect_identical(attr(md, "moby.land.name"), "coast")
    expect_identical(class(out), "list")
    expect_true(is_moby(out$data))
    expect_identical(mobyMeta(out$data), mobyMeta(md))
    expect_identical(attr(out$data, "moby.land.name"), attr(md, "moby.land.name"))
    expect_identical(names(out$data), names(md))
    expect_identical(attr(out, "points.relocated"), if (length(rows) == 6L) 3L else 0L)
    expect_true(inherits(attr(out, "processing.date"), "POSIXct"))
  }

  plain <- as.data.frame(make_moby(c(2L, 4L, 6L)))
  out_plain <- correctPositions(plain, spatial.layer = coast, lon.col = "longitude",
                                lat.col = "latitude", epsg.code = fx$epsg, verbose = FALSE)
  expect_identical(class(out_plain$data), "data.frame")
  expect_null(mobyMeta(out_plain$data))
})

test_that("matchDeployments keeps the caller's column map instead of the canonical defaults", {
  d <- as.data.frame(rays_detections)
  names(d)[names(d) == "datetime"] <- "dt"
  names(d)[names(d) == "station"]  <- "site"
  md <- as_moby(d, datetime.col = "dt", station.col = "site", verbose = FALSE)

  res <- matchDeployments(md, rays_deployments, datetime.col = "dt", station.col = "site",
                          verbose = FALSE)
  meta <- mobyMeta(res)
  expect_identical(meta$datetime.col, "dt")
  expect_identical(meta$station.col, "site")
  # rebuilding through as_moby() used to overwrite these with "datetime"/"station", leaving the
  # metadata pointing at columns that do not exist
  expect_true(all(c(meta$datetime.col, meta$station.col) %in% names(res)))
})

test_that("matchDeployments no longer demands an ID column it never uses", {
  d <- as.data.frame(rays_detections)
  d <- d[, setdiff(names(d), "ID")]
  expect_no_error(res <- matchDeployments(d, rays_deployments, station.col = "station",
                                          verbose = FALSE))
  expect_false(is_moby(res))
  expect_false("ID" %in% names(res))
})

test_that("as_moby(tags=) keys derived metadata to the declared roster, not the observed IDs", {
  # an animal tagged but never detected stays a declared level; filterDetections() validates
  # tagging.dates against those levels, so the derivation has to cover them too
  tg <- data.frame(transmitter = c("T1", "T2", "T3"), ID = c("a", "b", "ghost"),
                   tagging_date = as.POSIXct(c("2023-01-01", "2023-01-02", "2023-01-03"), tz = "UTC"),
                   nominal_delay = c(90, 90, 90), stringsAsFactors = FALSE)
  d <- data.frame(ID = factor(c("a", "a", "b"), levels = c("a", "b", "ghost")),
                  transmitter = c("T1", "T1", "T2"),
                  datetime = as.POSIXct(c("2023-06-01", "2023-06-02", "2023-06-03"), tz = "UTC"),
                  station = "S1", lon = -8, lat = 37, stringsAsFactors = FALSE)

  md <- as_moby(d, tags = tg, verbose = FALSE)
  expect_setequal(names(mobyMeta(md)$tagging.dates), c("a", "b", "ghost"))
  expect_setequal(names(mobyMeta(md)$nominal.delay), c("a", "b", "ghost"))
  # and the dataset it produces is one filterDetections() will accept
  expect_no_error(filterDetections(md, verbose = FALSE))

  # a plain character ID column declares nothing, so the observed values are the right key
  d2 <- d; d2$ID <- as.character(d2$ID)
  expect_setequal(names(mobyMeta(as_moby(d2, tags = tg, verbose = FALSE))$tagging.dates), c("a", "b"))
})
