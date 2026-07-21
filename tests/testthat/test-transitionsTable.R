tt_dataset <- function() {
  d <- data.frame(
    ID = factor(c(rep("A", 6), rep("B", 4))),
    datetime = as.POSIXct("2023-01-01", tz = "UTC") + c(0, 1, 2, 3, 4, 5, 0, 1, 2, 3) * 3600,
    site = c("S1", "S1", "S2", "S2", "S3", "S1", "S2", "S3", "S3", "S2"))
  as_moby(d, station.col = "site",
          tagging.dates = as.POSIXct(c(A = "2023-01-01", B = "2023-01-01"), tz = "UTC"))
}

test_that("transitionsTable formats a movement network into a clean table", {
  net <- calculateTransitions(tt_dataset(), spatial.col = "site")
  tt <- transitionsTable(net)
  expect_s3_class(tt, "data.frame")
  expect_true(all(c("Type", "Movements", "Individuals") %in% colnames(tt)))
  expect_true(any(grepl("Mean duration", colnames(tt))))
  expect_true(all(grepl("-->", tt$Type)))           # transitions only (no stationary rows)
  expect_true(any(grepl("%", tt$Individuals)))       # n (pct%)
})

test_that("transitionsTable summarises id.metadata per transition type", {
  net <- calculateTransitions(tt_dataset(), spatial.col = "site")
  meta <- data.frame(ID = c("A", "B"), sex = c("F", "M"), length = c(120, 135))
  tt <- transitionsTable(net, id.metadata = meta)
  expect_true(any(grepl("Mean Length", colnames(tt))))
  expect_true("Sex" %in% colnames(tt))
})

test_that("transitionsTable adds group rows for id.groups and rejects non-movement input", {
  net <- calculateTransitions(tt_dataset(), spatial.col = "site",
                              id.groups = list(grp1 = "A", grp2 = "B"))
  tt <- transitionsTable(net)
  expect_true(any(tt$Type %in% c("grp1", "grp2")))
  expect_error(transitionsTable(data.frame(a = 1)), "movement")
})
