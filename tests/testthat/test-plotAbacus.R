aba_dataset <- function() {
  set.seed(1)
  ids <- c("A1", "A2", "B1")
  d <- do.call(rbind, lapply(ids, function(id) {
    n <- sample(20:40, 1)
    data.frame(ID = id,
               datetime = as.POSIXct("2023-04-01", tz = "UTC") + sort(runif(n, 0, 80 * 86400)),
               station = sample(c("S1", "S2", "S3"), n, replace = TRUE),
               stringsAsFactors = FALSE)
  }))
  d$ID <- factor(d$ID)
  d$station <- factor(d$station)
  d
}

tag_dates <- function() as.POSIXct(c(A1 = "2023-04-01", A2 = "2023-04-03", B1 = "2023-04-02"), tz = "UTC")

test_that("plotAbacus runs across core option combinations", {
  d <- aba_dataset(); tg <- tag_dates()
  grp <- list("Group A" = c("A1", "A2"), "Group B" = "B1")
  pal <- c("orange", "red2", "turquoise3")
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressMessages(plotAbacus(d, tagging.dates = tg, id.groups = grp)))
  expect_no_error(suppressMessages(plotAbacus(d, tagging.dates = tg, color.by = "station", color.pal = pal, legend = FALSE)))
  expect_no_error(suppressMessages(plotAbacus(d, tagging.dates = tg, shade = FALSE, main = "title")))
})

test_that("plotAbacus aggregation (bin) and count scaling run", {
  d <- aba_dataset(); tg <- tag_dates()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressMessages(plotAbacus(d, tagging.dates = tg, color.by = "station", bin = "day")))
  expect_no_error(suppressMessages(plotAbacus(d, tagging.dates = tg, color.by = "station", bin = "auto", scale.by.count = TRUE)))
  expect_no_error(suppressMessages(plotAbacus(d, tagging.dates = tg, bin = 720)))  # 12 h
})

test_that("plotAbacus custom shading periods run and validate", {
  d <- aba_dataset(); tg <- tag_dates()
  periods <- data.frame(start = as.POSIXct(c("2023-04-10", "2023-05-20"), tz = "UTC"),
                        end   = as.POSIXct(c("2023-04-30", "2023-06-10"), tz = "UTC"),
                        label = c("phase A", "phase B"))
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressMessages(plotAbacus(d, tagging.dates = tg, shade = periods)))
  # missing required columns should error clearly
  expect_error(suppressMessages(plotAbacus(d, tagging.dates = tg, shade = data.frame(begin = 1, finish = 2))),
               "start.*end")
})

test_that("plotAbacus returns invisibly (no stray autoprint)", {
  d <- aba_dataset(); tg <- tag_dates()
  pdf(tempfile(fileext = ".pdf"), width = 9, height = 6)
  on.exit(dev.off(), add = TRUE)
  ret <- suppressMessages(plotAbacus(d, tagging.dates = tg, color.by = "station"))
  expect_null(ret)
  expect_invisible(suppressMessages(plotAbacus(d, tagging.dates = tg, color.by = "station")))
})
