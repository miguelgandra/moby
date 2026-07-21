# Regression: the parallel (cores > 1) branch of randomizeAssociations() must randomise EACH subset,
# exactly like the serial branch. Before the fix it randomised the full input table, so per-subset
# null distributions were computed from the wrong data (silently wrong p-values).

test_that("randomizeAssociations parallel branch randomises the subset, not the full table", {
  skip_on_cran()
  skip_if_not_installed("doRNG")
  skip_if_not_installed("doSNOW")
  skip_if_not_installed("foreach")
  skip_if_not_installed("parallel")
  skip_if(parallel::detectCores() < 2, "needs >= 2 cores")

  # P1 (bins 1-30): only A,B detected. P2 (bins 31-60): only C,D detected. Hence in period P1 any dyad
  # involving C or D CANNOT co-occur, so its null overlap must be NA. The bug randomised the full table
  # (where C,D exist in P2), which let them "co-occur" in P1 and inflated those nulls to ~100.
  tb <- as.POSIXct("2023-01-01", tz = "UTC") + (0:59) * 3600
  mk <- function(id, idx, st) data.frame(ID = id, timebin = tb[idx], station = st, stringsAsFactors = FALSE)
  df <- rbind(mk("A", 1:30, "R1"), mk("B", 1:30, "R1"),
              mk("C", 31:60, "R2"), mk("D", 31:60, "R2"))
  df$ID <- factor(df$ID, levels = c("A", "B", "C", "D"))
  md <- as_moby(df, timebin.col = "timebin", station.col = "station",
                tagging.dates = as.POSIXct("2023-01-01", tz = "UTC"))
  wt <- suppressWarnings(suppressMessages(createWideTable(md, value.col = "station", verbose = FALSE)))
  wt$period <- ifelse(wt$timebin < tb[31], "P1", "P2")
  ov <- suppressWarnings(suppressMessages(calculateAssociations(wt, subset = "period")))

  run <- function(cores) suppressWarnings(suppressMessages(
    randomizeAssociations(wt, ov, iterations = 30, cores = cores, random.seed = 1, verbose = FALSE)))
  serial   <- run(1)
  parallel <- run(2)

  # the NA-pattern of the null distribution must be identical between serial and parallel runs
  expect_identical(is.na(serial$pairwise_results$mean_null),
                   is.na(parallel$pairwise_results$mean_null))

  # specifically, the C-D dyad in period P1 is structurally impossible and must stay NA in parallel
  pr <- parallel$pairwise_results
  cd_p1 <- pr$subset == "P1" &
    ((pr$id1 == "C" & pr$id2 == "D") | (pr$id1 == "D" & pr$id2 == "C"))
  expect_true(all(is.na(pr$mean_null[cd_p1])))
})
