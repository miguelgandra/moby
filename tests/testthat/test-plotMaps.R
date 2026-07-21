# UD maps rely on adehabitatHR (kde) for the home-range estimation; keep small and guarded.
kud_dataset <- function(ids = c("A", "B"), n = 120, seed = 1, offset = 0) {
  set.seed(seed)
  do.call(rbind, lapply(seq_along(ids), function(k) {
    t <- as.POSIXct("2023-01-01", tz = "UTC") + (1:n) * 3600
    data.frame(ID = ids[k], timebin = t,
               lon = 5e5 + offset + cumsum(rnorm(n, 0, 90)),
               lat = 4.1e6 + offset + cumsum(rnorm(n, 0, 90)), stringsAsFactors = FALSE)
  }))
}

make_kud <- function(d) {
  d$ID <- factor(d$ID)
  md <- as_moby(d, timebin.col = "timebin", lon.col = "lon", lat.col = "lat", epsg.code = 32629)
  list(md = md, kud = suppressWarnings(calculateUDs(md, method = "kde", bandwidth = 250, verbose = FALSE)))
}

test_that("plotMaps draws UD isopleths and returns per-individual home-range areas", {
  skip_if_not_installed("adehabitatHR")
  x <- make_kud(kud_dataset())
  pdf(tempfile(fileext = ".pdf"), width = 8, height = 4); on.exit(dev.off(), add = TRUE)
  res <- suppressWarnings(suppressMessages(plotMaps(x$md, uds = x$kud, epsg.code = 32629)))
  expect_s3_class(res, "data.frame")
  expect_true(all(c("id", "isopleth", "area_km2", "n_detections") %in% names(res)))
  expect_equal(res$isopleth, c(0.95, 0.95))
  expect_true(all(res$area_km2 > 0))                           # real (metric) home-range areas
})

test_that("plotMaps restores par/layout when file = NULL (device-hygiene fix)", {
  skip_if_not_installed("adehabitatHR")
  x <- make_kud(kud_dataset())
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 1)); mar0 <- par("mar")
  suppressWarnings(suppressMessages(plotMaps(x$md, uds = x$kud, epsg.code = 32629, ncol = 2)))
  expect_equal(par("mar"), mar0)                               # margins restored
  expect_equal(par("mfrow"), c(1L, 1L))                        # layout restored
})

test_that("plotMaps survives local coordinates where easting < northing (axis-swap regression)", {
  skip_if_not_installed("adehabitatHR")
  # small local coordinates: raster::extent(bbox) previously scrambled the axes -> crash/mis-crop
  d <- kud_dataset(); d$lon <- d$lon - 5e5 + 15; d$lat <- d$lat - 4.1e6 + 2
  x <- make_kud(d)
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(plotMaps(x$md, uds = x$kud, epsg.code = 32629))))
})

test_that("plotMaps validates ud.contour and writes files without leaking a device", {
  skip_if_not_installed("adehabitatHR")
  x <- make_kud(kud_dataset())
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_error(plotMaps(x$md, uds = x$kud, epsg.code = 32629, ud.contour = 95), "proportion")
  before <- length(dev.list()); f <- tempfile(fileext = ".png")
  expect_no_error(suppressWarnings(suppressMessages(
    plotMaps(x$md, uds = x$kud, epsg.code = 32629, file = f))))
  expect_true(file.exists(f) && file.info(f)$size > 0)
  expect_equal(length(dev.list()), before)
})

# --- AKDE home ranges (ctmm) with confidence envelope -------------------------------------------
make_akde <- function(d) {
  d$ID <- factor(d$ID)
  md <- as_moby(d, timebin.col = "timebin", lon.col = "lon", lat.col = "lat", epsg.code = 32629)
  list(md = md, ud = suppressWarnings(suppressMessages(calculateUDs(md, verbose = FALSE))))  # akde default
}

test_that("plotMaps draws AKDE home ranges with a confidence envelope (band / outline / none)", {
  skip_if_not_installed("ctmm")
  x <- make_akde(kud_dataset(ids = "A"))
  expect_equal(attr(x$ud, "method"), "akde")               # exercising the AKDE (sf-contour) path
  for (style in c("band", "outline", "none")) {
    pdf(tempfile(fileext = ".pdf"))
    res <- suppressWarnings(suppressMessages(
      plotMaps(x$md, uds = x$ud, ud.uncertainty = style, epsg.code = 32629)))
    dev.off()
    expect_s3_class(res, "data.frame")
    expect_true(all(c("id", "isopleth", "area_km2", "n_detections") %in% names(res)))
    expect_gt(res$area_km2[1], 0)                            # real AKDE home-range area drawn
  }
  expect_error(plotMaps(x$md, uds = x$ud, ud.uncertainty = "nope"), "should be one of")
})

# --- coastline fallback --------------------------------------------------------------------------
test_that(".defaultCoastline fetches land for the extent and returns it in the map CRS", {
  skip_if_not_installed("maps")
  # a projected window over the SW Iberian coast -> land found, returned in the map CRS
  land <- moby:::.defaultCoastline(c(xmin = 480000, ymin = 4080000, xmax = 560000, ymax = 4160000),
                                   32629, coastline = TRUE, verbose = FALSE)
  expect_s3_class(land, "sfc")
  expect_gt(length(land), 0)
  expect_equal(sf::st_crs(land)$epsg, 32629L)
  # a window over open ocean -> no land
  expect_null(moby:::.defaultCoastline(c(xmin = 1e5, ymin = 4.0e6, xmax = 1.2e5, ymax = 4.02e6),
                                       32629, coastline = TRUE, verbose = FALSE))
})

test_that("plotMaps draws a default coastline when land.shape is NULL and coastline = TRUE", {
  skip_if_not_installed("maps"); skip_if_not_installed("adehabitatHR")
  x <- make_kud(kud_dataset())
  pdf(tempfile(fileext = ".pdf")); on.exit(dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(suppressMessages(
    plotMaps(x$md, uds = x$kud, epsg.code = 32629, coastline = TRUE))))
})
