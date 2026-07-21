deploy_log <- function() {
  # R1: two deployments 29-day gap apart (-> broken bar); R2: still active (NA recover)
  data.frame(
    receiver = c("R1", "R1", "R2", "R3"),
    station  = c("A", "A", "B", "C"),
    region   = c("north", "north", "north", "south"),
    deploy   = as.POSIXct(c("2020-01-01", "2020-03-01", "2020-01-01", "2020-02-01"), tz = "UTC"),
    recover  = as.POSIXct(c("2020-02-01", "2020-04-01", NA, "2020-05-01"), tz = "UTC"),
    stringsAsFactors = FALSE)
}

test_that("plotDeployments returns a correct per-row coverage table", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  cov <- suppressMessages(plotDeployments(deploy_log(), end = as.POSIXct("2020-06-01", tz = "UTC")))
  expect_s3_class(cov, "data.frame")
  expect_true(all(c("row", "group", "n_deployments", "first", "last", "active_days", "gap_days", "coverage") %in% names(cov)))
  expect_setequal(cov$row, c("R1", "R2", "R3"))                 # one row per receiver
  r1 <- cov[cov$row == "R1", ]
  expect_equal(r1$n_deployments, 2)                             # 29-day gap > merge.gaps -> 2 segments
  expect_equal(r1$active_days, 62)                              # 31 + 31 days listening
  expect_equal(r1$gap_days, 29)                                 # the servicing gap
  expect_equal(r1$coverage, round(62 / 91, 3))
})

test_that("merge.gaps fuses short gaps into a single bar", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  cov <- suppressMessages(plotDeployments(deploy_log(), merge.gaps = 40))  # 29-day gap now merged
  expect_equal(cov$n_deployments[cov$row == "R1"], 1)
})

test_that("plotDeployments aggregates by row.by and colours/legends by group.by", {
  pdf(tempfile(fileext = ".pdf"), width = 8, height = 5); on.exit(dev.off(), add = TRUE)
  # row.by='station' aggregates R1's two deployments (both at station A) into one row
  cov <- suppressMessages(plotDeployments(deploy_log(), row.by = "station", group.by = "region",
                                          end = as.POSIXct("2020-06-01", tz = "UTC")))
  expect_setequal(cov$row, c("A", "B", "C"))
  expect_setequal(cov$group, c("north", "south"))
  # per-row distinct colours idiom (example-2 style): group.by = row.by
  expect_no_error(suppressMessages(plotDeployments(deploy_log(), row.by = "station", group.by = "station")))
})

test_that("plotDeployments overlays events and honours sort.by", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  ev <- data.frame(receiver = c("R2", "R3"), time = as.POSIXct(c("2020-02-15", "2020-03-01"), tz = "UTC"))
  expect_no_error(suppressMessages(plotDeployments(deploy_log(), events = ev, sort.by = "start")))
  expect_no_error(suppressMessages(plotDeployments(deploy_log(), events = ev, sort.by = "name")))
  # an events frame missing the row key errors clearly
  expect_error(plotDeployments(deploy_log(), events = data.frame(time = as.POSIXct("2020-02-15", tz = "UTC"))), "row key")
})

test_that("plotDeployments supports group dividers and a manual date.interval", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- deploy_log()
  expect_no_error(suppressMessages(plotDeployments(d, group.by = "region", dividers = TRUE)))       # separators
  expect_no_error(suppressMessages(plotDeployments(d, date.interval = 3, date.format = "%b")))       # every 3rd month
  expect_no_error(suppressMessages(plotDeployments(d, date.interval = 1, date.format = "%b %Y")))
  expect_error(plotDeployments(d, date.interval = -1), "positive integer")                           # invalid manual interval
})

test_that("plotDeployments validates inputs and writes files without leaking a device", {
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  d <- deploy_log()
  expect_error(plotDeployments(d, row.by = "nope"), "not found")
  expect_error(plotDeployments(d, deploy.col = "nope"), "not found")
  d2 <- d; d2$deploy <- as.character(d2$deploy)
  expect_error(plotDeployments(d2), "POSIXct")
  before <- length(dev.list()); f <- tempfile(fileext = ".png")
  expect_no_error(suppressMessages(plotDeployments(d, group.by = "region", file = f)))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)                     # no device leaked
})
