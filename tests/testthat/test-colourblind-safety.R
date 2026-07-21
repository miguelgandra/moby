# Turns the package's "colourblind-safe" claim from an assertion into a verified property: simulate
# dichromacy (deuteranopia, protanopia) and check the palette colours stay perceptually distinguishable.
# Distance = CIE76 deltaE (Euclidean distance in CIE Lab), computed with base grDevices.
#
# NOTE: categorical colourblind-safety is only achievable up to ~8 colours -- a limit of dichromatic
# vision, not of moby. We therefore verify Okabe-Ito across its guaranteed 2..8 range, the 4-season
# shade palette, and viridis's sequential lightness ordering; large-n categorical safety is *not*
# claimed or tested (no palette can deliver it).

skip_if_not_installed("colorspace")

.lab        <- function(hex) grDevices::convertColor(t(grDevices::col2rgb(hex)) / 255, "sRGB", "Lab")
.min_pair_dE <- function(hex) { d <- as.matrix(stats::dist(.lab(hex))); min(d[upper.tri(d)]) }
.cvd        <- function(hex, f) colorspace::hex(f(colorspace::hex2RGB(substr(hex, 1, 7))))
.views      <- list(normal = identity,
                    deutan = function(h) .cvd(h, colorspace::deutan),
                    protan = function(h) .cvd(h, colorspace::protan))

test_that("Okabe-Ito categorical palette stays distinguishable under dichromacy (n = 2..8)", {
  for (n in 2:8) {
    cols <- moby:::.okabe_ito_pal(n)
    for (nm in names(.views)) {
      dE <- .min_pair_dE(.views[[nm]](cols))
      expect_gt(dE, 10)   # >10 = a clearly perceptible categorical difference
    }
  }
})

test_that("season-shade palette (economist, n = 4) stays distinguishable under dichromacy", {
  cols <- moby:::.economist_pal(4)
  for (nm in names(.views)) expect_gt(.min_pair_dE(.views[[nm]](cols)), 12)
})

test_that("viridis sequential ramp preserves monotonic lightness under dichromacy", {
  cols <- moby:::.viridis_pal(12)
  mono <- function(hex) { L <- .lab(hex)[, 1]; all(diff(L) > 0) || all(diff(L) < 0) }
  expect_true(mono(cols))                               # perceptually ordered for normal vision
  expect_true(mono(.cvd(cols, colorspace::deutan)))     # ...and still ordered under deuteranopia
  expect_true(mono(.cvd(cols, colorspace::protan)))     # ...and protanopia
})

test_that("the >8 categorical fallback beats the old interpolation for normal vision", {
  # Beyond 8 categories nothing is CVD-safe; the fallback ('Dark 3') at least keeps colours well
  # separated for normal vision, and must not regress below the interpolation it replaced.
  fallback <- moby:::.okabe_ito_pal(10)
  interp   <- grDevices::colorRampPalette(moby:::.okabe_ito_pal(8))(10)
  expect_gt(.min_pair_dE(fallback), .min_pair_dE(interp))
})
