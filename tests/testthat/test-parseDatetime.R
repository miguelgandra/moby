# The datetime parser's contract:
#   - datetime.format supplied  -> used exactly
#   - not supplied              -> inferred, but ONLY when the interpretation is unique
#   - several readings possible -> refuse, and ask for datetime.format
#   - never a silently wrong or silently all-NA column
#   - character input: tz is applied AS the zone; POSIXct input: the instant is preserved

P <- function(...) moby:::.parseDatetime(...)


test_that("unambiguous layouts are inferred without a format argument", {
  expect_equal(format(P(c("2023-06-01 08:00:00", "2023-06-02 09:00:00"), tz = "UTC")),
               c("2023-06-01 08:00:00", "2023-06-02 09:00:00"))
  expect_equal(format(P(c("2023-06-01T08:00:00", "2023-06-02T09:00:00"), tz = "UTC")),
               c("2023-06-01 08:00:00", "2023-06-02 09:00:00"))
  # the layout whose absence silently emptied all 847 rows of a real deployment log
  expect_equal(format(P(c("15/07/2003 00:00", "29/07/2003 12:00"), tz = "UTC")),
               c("2003-07-15 00:00:00", "2003-07-29 12:00:00"))
  # month names, which strptime otherwise reads only in an English locale
  expect_equal(format(P(c("06-Jun-2013 08:00:00", "11-Jul-2013 09:00:00"), tz = "UTC")),
               c("2013-06-06 08:00:00", "2013-07-11 09:00:00"))
})


test_that("a day above 12 settles the day-first / month-first question by itself", {
  # no ambiguity here: "%m/%d/%Y" cannot read "25/06/2013", so only one reading survives
  expect_equal(format(P(c("06/06/2013", "25/06/2013"), tz = "UTC"), "%Y-%m-%d"),
               c("2013-06-06", "2013-06-25"))
})


test_that("a genuinely ambiguous column is refused rather than guessed", {
  # every day <= 12, so day-first and month-first both fit and disagree. The old parser silently
  # preferred month-first because it came first in the list.
  expect_error(P(c("06/06/2013", "07/06/2013", "05/12/2013"), tz = "UTC"),
               "ambiguous", fixed = TRUE)
  expect_error(P(c("06/06/2013", "07/06/2013"), tz = "UTC"), "datetime.format", fixed = TRUE)
  # the message must name both readings, so the user can tell which they meant
  err <- tryCatch(P(c("06/06/2013", "07/06/2013"), tz = "UTC"), error = conditionMessage)
  expect_match(err, "%d/%m/%Y", fixed = TRUE)
  expect_match(err, "%m/%d/%Y", fixed = TRUE)
})


test_that("an explicit format is used exactly, and resolves what inference refuses", {
  v <- c("06/06/2013", "07/06/2013")
  expect_equal(format(P(v, tz = "UTC", format = "%d/%m/%Y"), "%Y-%m-%d"),
               c("2013-06-06", "2013-06-07"))
  expect_equal(format(P(v, tz = "UTC", format = "%m/%d/%Y"), "%Y-%m-%d"),
               c("2013-06-06", "2013-07-06"))
})


test_that("failure is explicit rather than an all-NA column", {
  expect_error(P(c("not a date", "nope"), tz = "UTC"), "could not be parsed")
  # an explicit format that fits nothing is a mistake worth reporting, not silent NAs
  expect_error(P("2023-06-01", tz = "UTC", format = "%d/%m/%Y %H:%M"), "could not be parsed")
  # the failing column is named, since an import maps several of them
  err <- tryCatch(P("xx", tz = "UTC", field = "recover"), error = conditionMessage)
  expect_match(err, "'recover'", fixed = TRUE)
})


test_that("blank and placeholder values become NA without derailing the column", {
  # a column that is entirely empty is a legitimately absent optional field
  expect_true(all(is.na(P(c(NA, "", NA), tz = "UTC"))))
  # station logs use placeholders for "never recovered"; those rows are NA, the rest still parse.
  # A value NO layout can read is junk, not evidence against the layout the other rows fit.
  r <- P(c("15/07/2003 00:00", "/  /", "29/07/2003 12:00"), tz = "UTC")
  expect_equal(sum(is.na(r)), 1L)
  expect_equal(format(r[c(1, 3)]), c("2003-07-15 00:00:00", "2003-07-29 12:00:00"))
})


test_that("a layout that only fits by misreading the string is rejected", {
  # "%Y/%m/%d" will happily read "25/06/2013" as the year 25, because strptime ignores what it did
  # not consume. Telemetry post-dates 1970, so an implausible year proves the layout wrong.
  r <- P(c("25/06/2013", "26/06/2013"), tz = "UTC")
  expect_equal(format(r, "%Y"), c("2013", "2013"))
})


test_that("a date-only layout does not win over one that also reads the time", {
  # "%Y-%m-%d" matches a full timestamp too - it just drops the time - so specificity has to decide
  r <- P(c("2023-06-01 08:30:00", "2023-06-02 09:45:00"), tz = "UTC")
  expect_equal(format(r, "%H:%M"), c("08:30", "09:45"))
})


test_that("the timezone contract differs by input type, on purpose", {
  # character carries no zone, so tz is applied AS the zone
  ch <- P("2023-06-01 08:00:00", tz = "Europe/Lisbon")
  expect_equal(format(ch, "%H:%M"), "08:00")
  expect_equal(as.numeric(ch), as.numeric(as.POSIXct("2023-06-01 08:00:00", tz = "Europe/Lisbon")))

  # POSIXct is already an instant: tz changes only the display, and the instant is preserved.
  # Reinterpreting it would silently move every timestamp in the dataset.
  utc <- as.POSIXct("2023-06-01 08:00:00", tz = "UTC")
  po <- P(utc, tz = "Europe/Lisbon")
  expect_equal(as.numeric(po), as.numeric(utc))
  expect_equal(attr(po, "tzone"), "Europe/Lisbon")
})


test_that("month names parse the same way whatever the session locale", {
  old <- Sys.getlocale("LC_TIME")
  on.exit(Sys.setlocale("LC_TIME", old), add = TRUE)
  skip_if(!nzchar(suppressWarnings(Sys.setlocale("LC_TIME", "fr_FR.UTF-8"))),
          "fr_FR locale unavailable")
  expect_equal(format(P("06-Jun-2013 08:00:00", tz = "UTC")), "2013-06-06 08:00:00")
  # and the session's locale is left as it was found
  expect_equal(Sys.getlocale("LC_TIME"), "fr_FR.UTF-8")
})
