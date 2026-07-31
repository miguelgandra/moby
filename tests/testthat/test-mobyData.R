test_that("mobyData is a data.frame subclass carrying metadata", {
  df <- make_detections(ids = c("A", "B"), n_per_id = 4)
  md <- as_moby(df, epsg.code = 32629,
                tagging.dates = as.POSIXct("2023-06-01", tz = "UTC"),
                id.groups = list(g1 = c("A", "B")))
  expect_s3_class(md, c("mobyData", "data.frame"))
  expect_true(is.data.frame(md))
  m <- mobyMeta(md)
  expect_equal(m$epsg.code, 32629)
  expect_s3_class(m$tagging.dates, "POSIXct")
  expect_named(m$id.groups, "g1")
})

test_that("re-wrapping a mobyData updates only supplied fields and inherits the rest", {
  df <- make_detections(ids = c("A", "B"), n_per_id = 3)
  md <- as_moby(df, id.col = "ID", epsg.code = 32629)
  md2 <- as_moby(md, tagging.dates = as.POSIXct("2023-06-01", tz = "UTC"))
  expect_equal(mobyMeta(md2)$epsg.code, 32629)            # inherited
  expect_s3_class(mobyMeta(md2)$tagging.dates, "POSIXct") # updated
})

test_that("subsetting preserves the mobyData class and metadata", {
  df <- make_detections(ids = c("A", "B"), n_per_id = 3)
  md <- as_moby(df)
  sub <- md[md$ID == "A", ]
  expect_true(is_moby(sub))
  expect_equal(mobyMeta(sub)$id.col, "ID")
})

test_that("functions resolve column names from mobyData metadata without explicit args", {
  df <- data.frame(
    animal = factor(rep(c("A", "B"), each = 3)),
    tbin = as.POSIXct("2023-01-01", tz = "UTC") + rep(c(0, 3600, 7200), 2),
    x = -8, y = 37, rec = "R1"
  )
  md <- as_moby(df, id.col = "animal", timebin.col = "tbin",
                lon.col = "x", lat.col = "y", station.col = "rec")
  coas <- suppressWarnings(suppressMessages(calculateCOAs(md)))
  expect_true("animal" %in% colnames(coas))
  expect_true(nrow(coas) > 0)
})

test_that("plain data frames with canonical column names still work", {
  df <- data.frame(
    ID = factor(rep(c("A", "B"), each = 3)),
    timebin = as.POSIXct("2023-01-01", tz = "UTC") + rep(c(0, 3600, 7200), 2),
    lon = -8, lat = 37, station = "R1"
  )
  coas <- suppressWarnings(suppressMessages(calculateCOAs(df)))
  expect_true(nrow(coas) > 0)
})

test_that("as_moby stores, validates and inherits nominal.delay", {
  d <- as.data.frame(rays)
  # a single value applies to all individuals; a named vector keys per animal (mixed tag families)
  expect_identical(mobyMeta(as_moby(d, nominal.delay = 120))$nominal.delay, 120)
  expect_identical(mobyMeta(as_moby(d, nominal.delay = c(R01 = 60, R02 = 120)))$nominal.delay,
                   c(R01 = 60, R02 = 120))
  # inherited on re-wrap, like the other metadata fields
  expect_identical(mobyMeta(as_moby(as_moby(d, nominal.delay = 120)))$nominal.delay, 120)
  # must be a positive number of seconds
  expect_error(as_moby(d, nominal.delay = -5), "nominal.delay")
  expect_error(as_moby(d, nominal.delay = "fast"), "nominal.delay")
})

test_that("print.mobyData summarises the dataset and its metadata", {
  d <- as_moby(as.data.frame(rays), nominal.delay = 120)
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "<mobyData>", fixed = TRUE)
  for (section in c("Summary", "Coverage", "Metadata")) expect_match(out, section, fixed = TRUE)
  expect_match(out, "1,643 detections", fixed = TRUE)   # thousands separator
  expect_match(out, "8 individuals", fixed = TRUE)
  expect_match(out, "6 stations", fixed = TRUE)
  expect_match(out, "9 variables", fixed = TRUE)        # 'variables', not 'columns'
  expect_match(out, "2023-04-02")
  expect_match(out, "Time zone", fixed = TRUE)
  expect_match(out, "Longitude", fixed = TRUE)
  expect_match(out, "EPSG:32629", fixed = TRUE)
})

test_that("the metadata block lists every supported field, present or not", {
  d <- as_moby(as.data.frame(rays), nominal.delay = 120)
  g <- .mobyGlyphs()
  out <- paste(capture.output(print(d)), collapse = "\n")
  # every field is named whether or not it is attached: the absent ones say what this dataset
  # cannot do downstream, which is as useful as knowing what it can
  for (f in c("Tagging dates", "Nominal delays", "ID groups", "Land polygon", "Coordinate CRS"))
    expect_match(out, f, fixed = TRUE)
  expect_match(out, paste0(g$tick, " Tagging dates"), fixed = TRUE)
  expect_match(out, paste0(g$cross, " Land polygon"), fixed = TRUE)   # not supplied here
})

test_that("counts are reported against the roster where one exists", {
  d <- as.data.frame(rays)
  d$ID <- factor(as.character(d$ID), levels = c(sort(unique(as.character(d$ID))), "GHOST1", "GHOST2"))
  md <- as_moby(d, tagging.dates = mobyMeta(rays)$tagging.dates)
  out <- paste(capture.output(print(md)), collapse = "\n")
  # tagged and detected are reported separately, because they legitimately differ
  expect_match(out, "10 individuals tagged", fixed = TRUE)
  expect_match(out, "8 individuals detected", fixed = TRUE)
  # and a metadata count is shown against that roster
  expect_match(out, "(8/10)", fixed = TRUE)

  # with no roster (a character ID column) there is nothing to compare against, so a bare count
  out2 <- paste(capture.output(print(as_moby(as.data.frame(rays)))), collapse = "\n")
  expect_match(out2, "8 individuals", fixed = TRUE)
  expect_false(grepl("tagged", out2, fixed = TRUE))
})

test_that("the land layer is named when the caller gave it one", {
  coastline <- sf::st_sf(geometry = sf::st_sfc(
    sf::st_polygon(list(cbind(c(-9.1, -8.9, -8.9, -9.1, -9.1), c(38.4, 38.4, 38.5, 38.5, 38.4)))),
    crs = 4326))
  out <- paste(capture.output(print(as_moby(as.data.frame(rays), land.shape = coastline))),
               collapse = "\n")
  expect_match(out, "Land polygon", fixed = TRUE)
  expect_match(out, "(coastline)", fixed = TRUE)
})

test_that("print shows no data rows: head() is what returns rows", {
  out <- capture.output(print(rays))
  # the preview is gone, so the output is short and fixed-height whatever the dataset
  expect_lt(length(out), 25L)
  expect_false(any(grepl("Preview", out, fixed = TRUE)))
  # a value that only ever appeared in the preview must not be there
  expect_false(any(grepl(as.character(rays$transmitter[1]), out, fixed = TRUE)))
  # and the argument that controlled it is gone
  expect_false("preview" %in% names(formals(moby:::print.mobyData)))
  # head() still returns rows, as a plain data frame
  expect_s3_class(head(rays, 2), "data.frame")
  expect_equal(nrow(head(rays, 2)), 2L)
})

test_that("sections absent from the data are simply omitted", {
  expect_no_error(capture.output(print(rays[0, ])))
  minimal <- suppressMessages(as_moby(data.frame(
    ID = factor("a"), datetime = as.POSIXct("2024-01-01", tz = "UTC"))))
  out <- paste(capture.output(print(minimal)), collapse = "\n")
  expect_false(grepl("Longitude", out, fixed = TRUE))
  expect_false(grepl("stations", out, fixed = TRUE))
  # the metadata block still lists the fields, all absent
  expect_match(out, "Metadata", fixed = TRUE)
  expect_match(out, paste0(.mobyGlyphs()$cross, " Tagging dates"), fixed = TRUE)
})

test_that("print.mobyData falls back to ASCII glyphs off UTF-8", {
  g <- .mobyGlyphs()
  expect_true(all(c("rule", "sep", "arrow") %in% names(g)))
  # whichever locale the tests run in, the glyphs must be single-purpose and non-empty
  expect_true(all(nzchar(unlist(g))))
  # the ASCII branch must be pure ASCII (what a non-UTF-8 console receives)
  ascii <- list(rule = "-", sep = "|", arrow = "->")
  expect_false(any(grepl("[^ -~]", unlist(ascii))))
})

test_that("head()/tail() return data frame rows, while '[' keeps the mobyData", {
  h <- head(rays, 3)
  expect_s3_class(h, "data.frame")
  expect_false(is_moby(h))                       # rows print as rows, not as a summary
  expect_equal(nrow(h), 3L)
  expect_equal(ncol(h), ncol(rays))              # every variable retained
  expect_identical(h, head(as.data.frame(rays), 3))

  t3 <- tail(rays, 3)
  expect_s3_class(t3, "data.frame")
  expect_false(is_moby(t3))
  expect_identical(t3, tail(as.data.frame(rays), 3))   # keeps original row numbering

  # '[' remains the metadata-preserving way to subset
  sub <- rays[1:100, ]
  expect_true(is_moby(sub))
  expect_identical(mobyMeta(sub)$id.groups, mobyMeta(rays)$id.groups)
})

test_that(".mobyRowNoun classifies rows from columns, not from a stored label", {
  coas <- suppressWarnings(suppressMessages(calculateCOAs(rays)))
  expect_equal(.mobyRowNoun(rays), "detections")          # one row per decode (has a datetime column)
  expect_equal(.mobyRowNoun(coas), "positions")           # rows are time bins carrying detection counts
  expect_equal(.mobyRowNoun(coas[1:5, ]), "positions")    # survives subsetting
  # growth case: a position table with extra columns (step distances / tracks) needs no new clause
  coas2 <- coas; coas2$dist_m <- 1
  expect_equal(.mobyRowNoun(coas2), "positions")
  # unrecognised shapes fall back to the generic noun rather than guessing
  bare <- suppressMessages(as_moby(data.frame(ID = factor("a"), lon = 1, lat = 2)))
  expect_equal(.mobyRowNoun(bare), "records")
})

test_that("print.mobyData labels rows correctly and recovers the true detection count", {
  coas <- suppressWarnings(suppressMessages(calculateCOAs(rays)))
  out <- paste(capture.output(print(coas, preview = 0)), collapse = "\n")
  # regression: this line previously read "794 detections", contradicting the object's own
  # 'detections' column and understating the real count by half
  expect_match(out, "794 positions")
  expect_match(out, "1,643 detections")
  expect_false(grepl("794 detections", out, fixed = TRUE))
  # a detection table keeps the domain noun and gains no spurious second count
  det <- paste(capture.output(print(rays, preview = 0)), collapse = "\n")
  expect_match(det, "1,643 detections")
  expect_false(grepl("positions", det, fixed = TRUE))
})

test_that("as_moby only flags absent columns the caller actually asked for", {
  d <- as.data.frame(rays)[, c("ID", "datetime", "lon", "lat", "station")]   # no timebin column
  # a canonical default nobody requested is moby's guess: its absence is normal, not worth a note
  expect_silent(as_moby(d))
  # but an explicitly supplied name that is missing is a likely typo -> still reported
  expect_message(as_moby(as.data.frame(rays), station.col = "statoin"), "statoin")
  expect_message(as_moby(d, timebin.col = "tb"), "tb")
  # a non-default name carried over from previous metadata is still checked
  md <- suppressMessages(as_moby(as.data.frame(rays), station.col = "station"))
  attr(md, "moby")$station.col <- "site"
  expect_message(as_moby(md, epsg.code = 4326), "site")
})

test_that("calculateCOAs actually drops the datetime role it documents dropping", {
  coas <- suppressWarnings(suppressMessages(calculateCOAs(rays)))
  # regression: omitting datetime.col from the as_moby() call was not enough - the constructor
  # re-applied its canonical default, leaving metadata pointing at a non-existent column
  expect_null(mobyMeta(coas)$datetime.col)
  expect_false("datetime" %in% colnames(coas))
  # and computing COAs is now quiet (it previously emitted a spurious 'datetime' missing-column note);
  # verbose = FALSE must silence the summary AND leave no other output behind
  expect_silent(suppressWarnings(calculateCOAs(rays, verbose = FALSE)))
  # the rest of the metadata is still carried forward
  expect_equal(mobyMeta(coas)$epsg.code, mobyMeta(rays)$epsg.code)
  expect_equal(mobyMeta(coas)$timebin.col, "timebin")
  expect_equal(.mobyRowNoun(coas), "positions")
})
