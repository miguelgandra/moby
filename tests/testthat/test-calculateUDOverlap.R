# Pairwise UD overlap. KDE (adehabitatHR) is fast + reliable and drives most checks;
# AKDE (ctmm) is guarded and kept minimal (model fitting is slow).

make_track <- function(id, n = 80, seed = 1) {
  set.seed(seed)
  t <- as.POSIXct("2023-01-01", tz = "UTC") + (1:n) * 3600
  data.frame(ID = id, timebin = t,
             lon = 5e5 + cumsum(rnorm(n, 0, 80)),
             lat = 4.1e6 + cumsum(rnorm(n, 0, 80)), stringsAsFactors = FALSE)
}

moby_track <- function(ids, seeds = seq_along(ids)) {
  d <- do.call(rbind, Map(function(i, s) make_track(i, seed = s), ids, seeds))
  d$ID <- factor(d$ID)
  as_moby(d, timebin.col = "timebin", lon.col = "lon", lat.col = "lat", epsg.code = 32629)
}

kde_ud <- function(ids = c("A", "B", "C")) {
  suppressWarnings(calculateUDs(moby_track(ids), method = "kde", bandwidth = 200, verbose = FALSE))
}

# ---- KDE ---------------------------------------------------------------------

test_that("KDE overlap returns tidy unordered pairs bounded in [0, 1]", {
  skip_if_not_installed("adehabitatHR")
  ud  <- kde_ud()
  res <- suppressWarnings(calculateUDOverlap(ud, index = "BA", verbose = FALSE))

  expect_s3_class(res, "data.frame")
  expect_true(all(c("id1", "id2", "BA") %in% names(res)))
  expect_equal(nrow(res), choose(3, 2))               # 3 animals -> 3 pairs
  expect_true(all(res$id1 != res$id2))                # no self-pairs
  # each unordered pair appears exactly once
  keys <- apply(res[, c("id1", "id2")], 1, function(p) paste(sort(p), collapse = "-"))
  expect_false(any(duplicated(keys)))
  expect_true(all(res$BA >= 0 & res$BA <= 1))
  expect_equal(attr(res, "method"), "kde")
  expect_equal(attr(res, "index"), "BA")
  expect_equal(attr(res, "contour"), 95)

  # BA is symmetric; the attached "matrix" is the estimator's own output verbatim (not a rebuild),
  # so its diagonal is the conditional self-overlap (~contour), not a forced 1.
  coll <- if (inherits(ud$ud, "estUDm")) ud$ud else ud$ud[[1]]
  m <- attr(res, "matrix")
  expect_true(isSymmetric(unname(m)))
  expect_equal(m, adehabitatHR::kerneloverlaphr(coll, method = "BA", percent = 95, conditional = TRUE))
})

test_that("KDE overlap honours alternative indices and rejects invalid ones", {
  skip_if_not_installed("adehabitatHR")
  ud <- kde_ud()
  res <- suppressWarnings(calculateUDOverlap(ud, index = "udoi", verbose = FALSE))  # case-insensitive
  expect_equal(attr(res, "index"), "UDOI")
  expect_true("UDOI" %in% names(res))
  expect_error(calculateUDOverlap(ud, index = "NOPE"), "Invalid 'index'")
})

test_that("directional KDE indices (HR/PHR) return both directions and stay asymmetric", {
  skip_if_not_installed("adehabitatHR")
  ud   <- kde_ud()                                   # A, B, C in one unit
  res  <- suppressWarnings(calculateUDOverlap(ud, index = "PHR", verbose = FALSE))

  expect_equal(nrow(res), 3 * 2)                     # one row per ORDERED pair: n*(n-1) = 6
  keys <- paste(res$id1, res$id2)
  expect_true(all(c("A B", "B A") %in% keys))        # both directions present
  expect_false(any(res$id1 == res$id2))              # still no self-pairs

  # ground truth: the estimator's own (asymmetric) matrix on the same estUDm
  coll <- if (inherits(ud$ud, "estUDm")) ud$ud else ud$ud[[1]]
  m <- adehabitatHR::kerneloverlaphr(coll, method = "PHR", percent = 95, conditional = TRUE)
  expect_false(isSymmetric(unname(m)))               # PHR is genuinely directional
  for (r in seq_len(nrow(res))) expect_equal(res$PHR[r], m[res$id1[r], res$id2[r]])

  # attached matrix is the raw asymmetric estimator output, not a symmetrised rebuild
  am <- attr(res, "matrix")
  expect_false(isSymmetric(unname(am)))
  expect_equal(am, m)
})

test_that("an all-singleton AKDE result keeps the CI columns (stable schema)", {
  skip_if_not_installed("ctmm")
  # two AKDE UDs in SEPARATE units -> every unit is a singleton -> 0 pairs, but the AKDE 5-col schema
  # (incl. BA_lower/BA_upper) must be preserved so downstream rbind()s stay consistent.
  fake_ud <- structure(list(), class = "UD")
  ud <- structure(list(ud = stats::setNames(list(fake_ud, fake_ud), c("U1::A", "U2::B"))),
                  method = "akde")
  res <- suppressWarnings(calculateUDOverlap(ud, index = "BA", verbose = FALSE))
  expect_equal(nrow(res), 0L)
  expect_true(all(c("id1", "id2", "BA", "BA_lower", "BA_upper") %in% names(res)))
})

test_that("a single individual yields zero pairs without error", {
  skip_if_not_installed("adehabitatHR")
  ud  <- kde_ud(ids = "A")
  res <- suppressWarnings(calculateUDOverlap(ud, verbose = FALSE))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0L)
})

test_that("id.groups annotates pairs as within / between", {
  skip_if_not_installed("adehabitatHR")
  res <- suppressWarnings(calculateUDOverlap(kde_ud(), index = "BA",
           id.groups = list(g1 = c("A", "B"), g2 = "C"), verbose = FALSE))
  expect_true(all(c("group1", "group2", "pair_type") %in% names(res)))
  expect_true(all(res$pair_type %in% c("within", "between")))
  ab <- res[res$id1 == "A" & res$id2 == "B", ]
  expect_equal(ab$pair_type, "within")
  ac <- res[res$id1 == "A" & res$id2 == "C", ]
  expect_equal(ac$pair_type, "between")
})

# ---- AKDE (ctmm) -------------------------------------------------------------

test_that("AKDE overlap reports BA with a bracketing confidence interval", {
  skip_if_not_installed("ctmm")
  ud  <- suppressWarnings(calculateUDs(moby_track(c("A", "B")), verbose = FALSE))  # akde default
  res <- suppressWarnings(calculateUDOverlap(ud, index = "BA", level = 0.95, verbose = FALSE))

  expect_equal(nrow(res), 1L)
  expect_true(all(c("BA", "BA_lower", "BA_upper") %in% names(res)))
  expect_true(all(res$BA_lower <= res$BA & res$BA <= res$BA_upper))
  expect_true(all(res$BA >= 0 & res$BA <= 1))
  expect_equal(attr(res, "method"), "akde")
  expect_equal(attr(res, "level"), 0.95)
  # AKDE (ctmm::overlap) only offers the Bhattacharyya coefficient
  expect_error(calculateUDOverlap(ud, index = "UDOI"), "AKDE")
})

# ---- input validation --------------------------------------------------------

test_that("bad input is rejected clearly", {
  # argument validation fires before the UD is inspected, so a stub with a $ud element suffices
  stub <- list(ud = 1)
  expect_error(calculateUDOverlap(list()), "\\$ud")
  expect_error(calculateUDOverlap(stub, level = 1.5), "level")
  expect_error(calculateUDOverlap(stub, contour = 0), "contour")
})
