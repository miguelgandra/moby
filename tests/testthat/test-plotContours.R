contour_dataset <- function() {
  set.seed(1)
  n <- 1500
  t0 <- as.POSIXct("2022-11-01", tz = "UTC")
  dt <- t0 + sort(runif(n, 0, 300 * 86400))   # ~10 months, spanning two calendar years
  hr <- as.numeric(format(dt, "%H")) + as.numeric(format(dt, "%M")) / 60
  doy <- as.numeric(format(dt, "%j"))
  data.frame(
    ID = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
    datetime = dt,
    depth = 15 - 8 * cos(2 * pi * hr / 24) + 5 * sin(2 * pi * doy / 365) + rnorm(n, 0, 2),
    temperature = 18 + 4 * sin(2 * pi * (doy - 80) / 365) + rnorm(n, 0, 0.5),
    species = factor(sample(c("sp1", "sp2"), n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
}

coords <- c(-8.9, 37.0)

test_that("plotContours runs across modes, intervals and overlays", {
  d <- contour_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  cc <- function(...) plotContours(d, variables = "depth", id.col = "ID", datetime.col = "datetime",
                                   coords = coords, ...)
  expect_no_error(suppressWarnings(suppressMessages(cc())))                              # true extent
  expect_no_error(suppressWarnings(suppressMessages(cc(annual.cycle = TRUE))))           # annual cycle
  expect_no_error(suppressWarnings(suppressMessages(cc(date.interval = "1 week"))))      # weekly bins
  expect_no_error(suppressWarnings(suppressMessages(cc(time.interval = "30 mins"))))     # sub-hourly
  expect_no_error(suppressWarnings(suppressMessages(cc(diel.lines = 2))))                # 2 diel lines
  expect_no_error(suppressWarnings(suppressMessages(cc(shared.scale = TRUE, split.by = "species"))))
  expect_no_error(suppressWarnings(suppressMessages(cc(legend = FALSE, grid = FALSE))))
  expect_no_error(suppressWarnings(suppressMessages(cc(zlab = "Depth (m)"))))             # colour-bar label
  expect_no_error(suppressWarnings(suppressMessages(cc(zlab = FALSE))))                   # no colour-bar label
  expect_no_error(suppressWarnings(suppressMessages(cc(zlim = c(5, 25)))))                # clamped colour range
  expect_no_error(suppressWarnings(suppressMessages(cc(reverse.scale = TRUE))))           # reversed legend axis
  expect_no_error(suppressWarnings(suppressMessages(cc(sample.size = FALSE))))            # no (n=) annotation
  expect_no_error(suppressWarnings(suppressMessages(cc(sample.size.color = "red"))))      # custom (n=) colour
  expect_no_error(suppressWarnings(suppressMessages(
    plotContours(d, variables = c("depth", "temperature"), var.titles = c("Depth (m)", "Temp"),
                 id.col = "ID", datetime.col = "datetime", coords = coords, ncol = 2))))
})

test_that("plotContours errors clearly when the date.interval yields a single bin", {
  d <- contour_dataset()
  d <- d[d$datetime < min(d$datetime) + 15 * 86400, ]   # ~2 weeks -> one monthly bin
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_error(
    plotContours(d, variables = "depth", id.col = "ID", datetime.col = "datetime", diel.lines = 0),
    "single date bin")
})

test_that("plotContours does not require coords when diel lines are off", {
  d <- contour_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(
    plotContours(d, variables = "depth", id.col = "ID", datetime.col = "datetime", diel.lines = 0))))
  expect_error(
    plotContours(d, variables = "depth", id.col = "ID", datetime.col = "datetime", diel.lines = 4),
    "coords")
})

test_that("plotContours validates inputs and returns invisibly", {
  d <- contour_dataset()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  expect_error(plotContours(d, variables = "bogus", id.col = "ID", datetime.col = "datetime",
                            diel.lines = 0), "variable")
  expect_error(plotContours(d, variables = "depth", var.titles = c("a", "b"), id.col = "ID",
                            datetime.col = "datetime", diel.lines = 0), "var.titles")
  ret <- suppressWarnings(suppressMessages(
    plotContours(d, variables = "depth", id.col = "ID", datetime.col = "datetime", diel.lines = 0)))
  expect_null(ret)
})
