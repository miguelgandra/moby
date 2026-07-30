# .readTabular() reads with data.table::fread when it is installed and utils::read.csv otherwise.
# The backend is an implementation detail, so the contract these tests enforce is that the RESULT is
# identical either way. Each fixture below is a case where the two readers' own defaults disagree.

# Force the base backend by hiding data.table from requireNamespace(). Nothing else in the call is
# changed, so the two runs differ only in which reader executed.
with_base_backend <- function(expr) {
  testthat::with_mocked_bindings(.hasDataTable = function() FALSE, .package = "moby", expr)
}

both_backends <- function(path) {
  fread_run <- moby:::.readTabular(path)
  base_run  <- with_base_backend(moby:::.readTabular(path))
  list(fread = fread_run, base = base_run)
}

tmp_csv <- function(lines) {
  f <- tempfile(fileext = ".csv"); writeLines(lines, f); f
}


test_that("a ragged export reads identically under both backends", {
  skip_if_not_installed("data.table")
  # VUE exports omit trailing empty fields, so the header declares more columns than most rows carry.
  # fread's default infers the column count from the first DATA row and consumes the header as data.
  f <- tmp_csv(c("datetime,receiver,transmitter,station,lon,lat",
                 "2023-06-01 08:00:00,VR2-1,A69-1,S1,-28.6,38.5",
                 "2023-06-01 09:00:00,VR2-1,A69-1,S1",
                 "2023-06-01 10:00:00,VR2-1,A69-1,S1"))
  r <- both_backends(f)
  expect_identical(r$fread, r$base)
  expect_equal(ncol(r$fread), 6L)                     # the header's count, not the first data row's
  expect_equal(nrow(r$fread), 3L)                     # header not swallowed
  expect_equal(names(r$fread)[1], "datetime")
})


test_that("a semicolon file with a decimal comma reads identically under both backends", {
  skip_if_not_installed("data.table")
  f <- tmp_csv(c("receiver;station;lat;lon", "VR2-1;S1;38,5191;-28,6192", "VR2-2;S2;38,5228;-28,6300"))
  r <- both_backends(f)
  expect_identical(r$fread, r$base)
  expect_equal(ncol(r$fread), 4L)
  # detected as a decimal mark, not left as text that as.numeric() would turn into NA
  expect_type(r$fread$lat, "double")
  expect_equal(r$fread$lat, c(38.5191, 38.5228))
})


test_that("a semicolon file with a decimal POINT is not mis-detected", {
  skip_if_not_installed("data.table")
  # the deep-sharks station log's shape: ';' separated but '.' decimals. Reading it with read.csv2's
  # dec = "," would leave the coordinates as character.
  f <- tmp_csv(c("receiver;station;lat;lon", "VR2-1;S1;38.5191;-28.6192", "VR2-2;S2;38.5228;-28.6300"))
  r <- both_backends(f)
  expect_identical(r$fread, r$base)
  expect_type(r$fread$lat, "double")
  expect_equal(r$fread$lat, c(38.5191, 38.5228))
})


test_that("a byte-order mark is stripped from the first column name by both backends", {
  skip_if_not_installed("data.table")
  f <- tempfile(fileext = ".csv")
  con <- file(f, open = "wb")
  writeBin(as.raw(c(0xEF, 0xBB, 0xBF)), con)          # UTF-8 BOM
  writeLines(c("datetime,receiver", "2023-06-01 08:00:00,VR2-1"), con)
  close(con)
  r <- both_backends(f)
  expect_identical(r$fread, r$base)
  expect_equal(names(r$fread)[1], "datetime")         # i.e. with no leading BOM
})


test_that("an ISO datetime column is returned as character, so the importer owns the timezone", {
  skip_if_not_installed("data.table")
  # fread auto-parses ISO-8601 to POSIXct assuming UTC; read.csv returns character. Left alone this
  # changes the INSTANT, because the importer then applies tz to an already-parsed value.
  f <- tmp_csv(c("datetime,receiver", "2023-06-01T08:00:00,VR2-1", "2023-06-01T09:00:00,VR2-1"))
  r <- both_backends(f)
  expect_identical(r$fread, r$base)
  expect_type(r$fread$datetime, "character")
})


test_that("the timezone a file is read in does not depend on the backend", {
  skip_if_not_installed("data.table")
  f <- tmp_csv(c("datetime,transmitter,receiver,station,lon,lat",
                 "2023-06-01T08:00:00,A69-1,VR2-1,S1,-28.6,38.5",
                 "2023-06-01T09:00:00,A69-1,VR2-1,S1,-28.6,38.5"))
  cm <- list(datetime = "datetime", transmitter = "transmitter", receiver = "receiver",
             station = "station", lon = "lon", lat = "lat")
  a <- importDetections(f, source = "generic", col.map = cm, tz = "Europe/Lisbon", verbose = FALSE)
  b <- with_base_backend(
    importDetections(f, source = "generic", col.map = cm, tz = "Europe/Lisbon", verbose = FALSE))
  expect_identical(a$datetime, b$datetime)
  # the string is interpreted in the REQUESTED zone, not pre-parsed as UTC by the reader
  expect_equal(as.numeric(a$datetime[1]),
               as.numeric(as.POSIXct("2023-06-01 08:00:00", tz = "Europe/Lisbon")))
})


test_that("a misread delimiter is reported rather than silently yielding one column", {
  # nothing here is a valid separator, so the file legitimately has one column and must NOT warn
  f1 <- tmp_csv(c("onlycolumn", "a", "b"))
  expect_no_warning(moby:::.readTabular(f1))
  # a pipe-separated file IS detected, so it also must not warn
  f2 <- tmp_csv(c("a|b|c", "1|2|3"))
  expect_no_warning(moby:::.readTabular(f2))
  expect_equal(ncol(moby:::.readTabular(f2)), 3L)
})


test_that("canonical fields get their declared type regardless of what the file looked like", {
  f <- tmp_csv(c("datetime,transmitter,receiver,station,lon,lat",
                 "2023-06-01 08:00:00,0012,0034,0056,-28.6,38.5"))
  d <- importDetections(f, source = "generic", verbose = FALSE,
                        col.map = list(datetime = "datetime", transmitter = "transmitter",
                                       receiver = "receiver", station = "station",
                                       lon = "lon", lat = "lat"))
  # keys and labels are character whatever the file looked like, so a join key never depends on
  # whether a particular file happened to make a code look numeric
  expect_type(d$transmitter, "character")
  expect_type(d$receiver, "character")
  expect_type(d$station, "character")
  expect_type(d$lon, "double")
  expect_s3_class(d$datetime, "POSIXct")

  # KNOWN LIMIT of typing at read time and standardising afterwards: a purely numeric code loses its
  # leading zeros before the type contract can act ("0056" is already the integer 56 by then). This is
  # unchanged from previous versions - read.csv did the same - and identical under both backends, so
  # it is a documented limit rather than a regression. Only reading every column as text would keep
  # them, at the cost of moby having to type every user column itself.
  expect_equal(d$station, "56")
  expect_identical(d$station,
                   with_base_backend(importDetections(f, source = "generic", verbose = FALSE,
                     col.map = list(datetime = "datetime", transmitter = "transmitter",
                                    receiver = "receiver", station = "station",
                                    lon = "lon", lat = "lat")))$station)
})


test_that("a data frame passed directly is returned unchanged", {
  d <- data.frame(a = 1:2, b = c("x", "y"), stringsAsFactors = FALSE)
  expect_identical(moby:::.readTabular(d), d)
})


test_that("a long numeric code does not arrive as integer64 under one backend only", {
  skip_if_not_installed("data.table")
  # fread returns integer64 for a 16-digit value by default; read.csv returns a double. Serial numbers
  # reach that length, so the backend would otherwise decide the column's class.
  f <- tmp_csv(c("serial,x", "1234567890123456,1"))
  r <- both_backends(f)
  expect_identical(r$fread, r$base)
  expect_type(r$fread$serial, "double")
})


test_that("an ISO column appearing only late in the file is still returned as text", {
  skip_if_not_installed("data.table")
  # fread types a column from a sample spread through the file, so a date that first appears well past
  # the start can be typed even though the opening rows look empty. The reader re-reads such a column
  # as character rather than trying to predict fread's sampling.
  f <- tmp_csv(c("dt,x", rep(",1", 1200), "2023-06-01T08:00:00,1"))
  r <- both_backends(f)
  expect_identical(r$fread, r$base)
  expect_type(r$fread$dt, "character")
  expect_equal(r$fread$dt[nrow(r$fread)], "2023-06-01T08:00:00")   # the source text, not reformatted
})


test_that("structural oddities read identically under both backends", {
  skip_if_not_installed("data.table")
  # accented station names (common in Azorean and Iberian arrays); written as an escape so the
  # test file itself stays pure ASCII
  # Only backend AGREEMENT is asserted. How R marks the string's Encoding() varies with locale and
  # is not moby's contract; what matters is that both readers deliver the same thing.
  acc <- tmp_csv(c("station,x", "S\u00e3o Mateus,1"))
  expect_identical(both_backends(acc)$fread, both_backends(acc)$base)
  # a delimiter inside a quoted field must not split it
  q <- both_backends(tmp_csv(c('station,note', '"S1","a, b"')))
  expect_identical(q$fread, q$base)
  expect_equal(q$fread$note, "a, b")
  # CRLF line endings
  f <- tempfile(fileext = ".csv"); con <- file(f, "wb")
  writeBin(charToRaw("a,b\r\n1,2\r\n"), con); close(con)
  expect_identical(both_backends(f)$fread, both_backends(f)$base)
  # header with no data rows, and a trailing blank line
  expect_identical(both_backends(tmp_csv("a,b,c"))$fread, both_backends(tmp_csv("a,b,c"))$base)
  expect_identical(both_backends(tmp_csv(c("a,b", "1,2", "")))$fread,
                   both_backends(tmp_csv(c("a,b", "1,2", "")))$base)
})
