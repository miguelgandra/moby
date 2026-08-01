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
