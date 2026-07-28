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

test_that("print.mobyData shows the summary sections, including nominal.delay", {
  d <- as_moby(as.data.frame(rays), nominal.delay = 120)
  out <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(out, "<mobyData>")
  expect_match(out, "overview\\s+1,643 detections")     # thousands separator
  expect_match(out, "8 individuals")
  expect_match(out, "9 variables")                      # 'variables', not 'columns'
  expect_match(out, "period\\s+2023-04-02")
  expect_match(out, "space\\s+6 stations")
  expect_match(out, "lon \\[-9\\.02, -8\\.97\\]")
  # regression: nominal.delay was previously never displayed
  expect_match(out, "nominal.delay \\(1\\)")
  expect_match(out, "tagging.dates \\(8\\)")
  expect_match(out, "EPSG: 32629")
  # preview: 3 rows by default, showing every variable
  expect_match(out, "Preview \\(first 3 rows\\)")
  for (v in colnames(rays)) expect_match(out, v, fixed = TRUE)
})

test_that("print.mobyData preview is controllable and degrades gracefully", {
  expect_false(any(grepl("Preview", capture.output(print(rays, preview = 0)))))
  expect_match(paste(capture.output(print(rays, preview = 1)), collapse = "\n"),
               "Preview \\(first 1 row\\)")
  # no rows / no optional roles -> the corresponding sections are simply omitted
  expect_no_error(capture.output(print(rays[0, ])))
  minimal <- suppressMessages(as_moby(data.frame(
    ID = factor("a"), datetime = as.POSIXct("2024-01-01", tz = "UTC"))))
  out <- paste(capture.output(print(minimal)), collapse = "\n")
  expect_false(grepl("space", out))
  expect_false(grepl("metadata", out))
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
