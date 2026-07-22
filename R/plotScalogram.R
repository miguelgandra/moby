##################################################################################################
## Scalogram: time-resolved rhythms via the continuous wavelet transform ##########################
##################################################################################################

#' Plot wavelet scalograms (time-resolved rhythms)
#'
#' @description Computes and plots a continuous wavelet transform (CWT) scalogram per individual - a
#' time (x) by period (y) heat map of spectral power - showing how dominant rhythms in a detection or
#' environmental time series strengthen, weaken or shift through the study. It is the time-resolved
#' complement of \code{\link{plotPeriodogram}} (which gives one global spectrum per animal). This is a
#' wrapper around \code{\link[wavScalogram]{cwt_wst}} (Bolos & Benitez, 2022).
#'
#' @details Each series is first laid on a complete, regularly-sampled time grid (via the shared
#' signal builder), so gaps are represented rather than deleted, then gap-filled, detrended and
#' (optionally) standardised before the CWT. The cone of influence (COI) - the region where edge
#' effects distort the estimate - is *always* excluded when computing the colour scale, so edge
#' artefacts never dominate the colours; by default it is then drawn as a translucent veil
#' (`mask.coi = FALSE`) or, if preferred, hard-masked in solid grey (`mask.coi = TRUE`).
#'
#' Power is shown on a log scale by default (`power.scaling = "log"`): raw wavelet power is extremely
#' peaked, so a linear scale hides everything but the single strongest feature.
#'
#' @inheritParams as_moby
#' @param data A `mobyData` or data frame of binned detections / time-based measurements. Gaps are
#' allowed (see `gap.handling`).
#' @param variable Name of the numeric column to analyse.
#' @param id.groups Optional named list of ID groups (one block of panels each).
#' @param wavelet.type Wavelet passed to `cwt_wst`'s `wname`: "MORLET" (default), "DOG", "PAUL",
#' "HAAR" or "HAAR2".
#' @param gap.handling How internal gaps are filled: "zero" (absence; default), "mean", "locf"
#' (last observation carried forward) or "interpolate" (linear interpolation).
#' @param detrend One of "none" (default; demean), "linear" (remove OLS trend), "diff" (first
#' difference) or "loess" (residuals of a LOESS fit). `diff`/`loess` are high-pass - use with care.
#' @param loess.span Span for `detrend = "loess"`. Defaults to 0.75.
#' @param standardize Logical; z-score each series after detrending. Defaults to FALSE.
#' @param power.scaling Colour scaling of wavelet power: "log" (default), "linear", "sqrt" or
#' "quantile" (exploratory - maximises contrast but can look noisy).
#' @param upper.value,upper.quant Optional winsorising thresholds on raw power (an absolute value, or
#' a quantile in (0, 1]) applied before scaling, to tame extreme outliers.
#' @param shared.scale Logical; share the colour scale across individuals for density comparison.
#' Defaults to FALSE.
#' @param mask.coi Logical; if TRUE the COI is hard-masked in solid grey, otherwise drawn as a
#' translucent veil. The COI is excluded from the colour scale either way. Defaults to FALSE.
#' @param period.range Length-2 numeric period range (y-axis), in `time.unit`. Defaults to c(3, 48).
#' @param axis.periods Periods highlighted on the y-axis, in `time.unit`. Defaults to c(6,12,16,24,48).
#' @param time.unit Unit for `period.range`/`axis.periods`: "hours" (default), "mins" or "days".
#' @param color.pal Colour palette for power. If NULL, the perceptually-uniform viridis palette.
#' @param min.days Discard individuals detected on fewer than this many distinct days.
#' @param date.format Optional \code{\link[base]{strptime}} format for x-axis labels. If NULL
#' (default), calendar-aware labels are chosen automatically.
#' @param main Overall plot title. If NULL, generated from `variable`.
#' @param cex Global expansion factor for all plot text. Defaults to 1.
#' @param ncol Number of panel columns. Defaults to 1.
#' @param cores Number of CPU cores for the CWT computation (needs \pkg{parallel}, \pkg{doSNOW},
#' \pkg{foreach} when > 1). Defaults to 1.
#' @template deviceArgs
#' @param ... Further arguments passed to \code{\link[wavScalogram]{cwt_wst}}.
#'
#' @return Invisibly, a tidy data frame with one row per individual: the dominant (global) period
#' (in `time.unit`) and its power. The per-individual global wavelet spectra and their period axis
#' are attached as attributes `"spectra"` and `"periods"`.
#' @references Bolos, V. J., & Benitez, R. (2022). wavScalogram: an R package with scalogram wavelet
#' tools for time series analysis. The R Journal, 14(2), 164-185.
#' @seealso \code{\link{plotPeriodogram}}, \code{\link{plotActograms}}
#' @examples
#' \donttest{
#' # Hourly detection counts per individual, then a time-resolved wavelet scalogram
#' binned <- aggregate(list(detections = rep(1L, nrow(rays))),
#'                     by = list(ID = rays$ID, timebin = rays$timebin), FUN = sum)
#' one <- binned[binned$ID %in% rays_tags$ID[rays_tags$species == "Raja clavata"], ]
#' if (requireNamespace("wavScalogram", quietly = TRUE)) {
#'   plotScalogram(one, variable = "detections", id.col = "ID",
#'                 timebin.col = "timebin", ncol = 2)
#' }
#' }
#' @export

plotScalogram <- function(data,
                          variable,
                          id.col = NULL,
                          timebin.col = NULL,
                          id.groups = NULL,
                          wavelet.type = "MORLET",
                          gap.handling = c("zero", "mean", "locf", "interpolate"),
                          detrend = c("none", "linear", "diff", "loess"),
                          loess.span = 0.75,
                          standardize = FALSE,
                          power.scaling = c("log", "linear", "sqrt", "quantile"),
                          upper.value = NULL,
                          upper.quant = NULL,
                          shared.scale = FALSE,
                          mask.coi = FALSE,
                          period.range = c(3, 48),
                          axis.periods = c(6, 12, 16, 24, 48),
                          time.unit = c("hours", "mins", "days"),
                          color.pal = NULL,
                          min.days = NULL,
                          date.format = NULL,
                          main = NULL,
                          ncol = 1,
                          cores = 1,
                          cex = 1,
                          file = NULL,
                          width = NULL,
                          height = NULL,
                          res = 300,
                          ...) {

  ##############################################################################
  # Checks #####################################################################
  ##############################################################################

  reviewed <- .validateArguments()
  data <- as.data.frame(reviewed$data)
  gap.handling <- match.arg(gap.handling); detrend <- match.arg(detrend)
  power.scaling <- match.arg(power.scaling); time.unit <- match.arg(time.unit)

  errors <- c()
  if(!requireNamespace("wavScalogram", quietly = TRUE))
    errors <- c(errors, "The 'wavScalogram' package is required (install.packages('wavScalogram')).")
  if(!inherits(data[[variable]], c("numeric", "integer")))
    errors <- c(errors, "The 'variable' column must be numeric or integer.")
  if(!is.numeric(period.range) || length(period.range) != 2)
    errors <- c(errors, "'period.range' must be a length-2 numeric (min, max).")
  if(detrend == "loess" && (!is.numeric(loess.span) || loess.span <= 0))
    errors <- c(errors, "'loess.span' must be a positive number.")
  if(!is.null(upper.value) && (!is.numeric(upper.value) || upper.value < 0))
    errors <- c(errors, "'upper.value' must be a non-negative number.")
  if(!is.null(upper.quant) && (!is.numeric(upper.quant) || upper.quant <= 0 || upper.quant > 1))
    errors <- c(errors, "'upper.quant' must be in (0, 1].")
  if(cores > 1 && !all(vapply(c("parallel", "doSNOW", "foreach"), requireNamespace, logical(1), quietly = TRUE)))
    errors <- c(errors, "cores > 1 requires the 'parallel', 'doSNOW' and 'foreach' packages.")
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)

  if(is.null(main)) main <- paste("Wavelet power spectrum -", tools::toTitleCase(variable))
  if(is.null(id.groups)) id.groups <- list(levels(factor(data[[id.col]])))
  if(is.null(color.pal)) color.pal <- .viridis_pal(100)

  cex_title <- 1.2 * cex; cex_lab <- 1.1 * cex; cex_axis <- 1.0 * cex; cex_legend <- 0.9 * cex

  time_factor <- switch(time.unit, mins = 1, hours = 60, days = 1440)
  unit_abbrev <- switch(time.unit, mins = "min", hours = "h", days = "days")


  ##############################################################################
  # Build series ###############################################################
  ##############################################################################

  built <- .buildRhythmSeries(data, id.col, timebin.col, value.col = variable, binary = FALSE,
                              gap.handling = gap.handling, detrend.method = detrend,
                              loess.span = loess.span, standardize = standardize, min.days = min.days)
  if(length(built$ids) == 0) stop("No individuals passed the 'min.days' filter.", call. = FALSE)

  # panel order: grouped, retaining only individuals that survived the filter
  id.groups <- lapply(id.groups, function(x) x[x %in% built$ids]); id.groups <- id.groups[lengths(id.groups) > 0]
  ord_ids <- unlist(id.groups, use.names = FALSE)
  n_ind <- length(ord_ids)

  # low-side scale padding so the shortest period is not truncated at the axis edge
  scales_arg <- c(period.range[1] * time_factor * 0.95, period.range[2] * time_factor, 20)


  ##############################################################################
  # CWT computation ############################################################
  ##############################################################################

  .printScalogramSummary(n_ids = n_ind, n_total = built$n_total, dt = built$dt, wavelet = wavelet.type,
                         gap.handling = gap.handling, detrend = detrend, power.scaling = power.scaling,
                         period.range = period.range, unit = unit_abbrev)

  data_ts <- lapply(ord_ids, function(id) stats::ts(built$series[[id]]$values))
  cwt_fun <- function(x) wavScalogram::cwt_wst(x, dt = built$dt, scales = scales_arg, powerscales = TRUE,
                                               wname = wavelet.type, border_effects = "BE",
                                               makefigure = FALSE, energy_density = TRUE, figureperiod = TRUE, ...)
  pb <- utils::txtProgressBar(min = 0, max = n_ind, style = 3)
  if(cores > 1){
    cl <- parallel::makeCluster(cores); doSNOW::registerDoSNOW(cl); on.exit(parallel::stopCluster(cl), add = TRUE)
    `%dopar%` <- foreach::`%dopar%`
    opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
    cwts <- foreach::foreach(i = seq_along(data_ts), .options.snow = opts, .packages = "wavScalogram") %dopar% cwt_fun(data_ts[[i]])
  }else{
    cwts <- vector("list", n_ind)
    for(i in seq_along(data_ts)){ cwts[[i]] <- cwt_fun(data_ts[[i]]); utils::setTxtProgressBar(pb, i) }
  }
  close(pb)


  ##############################################################################
  # Power, scale limits, global spectra ########################################
  ##############################################################################

  power_of <- function(cwt) t(t(abs(cwt$coefs)^2) / cwt$scales)          # power density: |W|^2 / scale
  coi_mask_of <- function(cwt, Z){                                        # TRUE inside COI (valid), per cell
    m <- matrix(TRUE, nrow = nrow(Z), ncol = ncol(Z))
    ff <- cwt$fourierfactor
    for(j in seq_len(nrow(Z))) m[j, cwt$scales * ff > cwt$coi_maxscale[j] * ff] <- FALSE
    m
  }
  winsorize <- function(Z){
    if(!is.null(upper.quant)) Z[Z > stats::quantile(Z, upper.quant, na.rm = TRUE)] <- stats::quantile(Z, upper.quant, na.rm = TRUE)
    if(!is.null(upper.value)) Z[Z > upper.value] <- upper.value
    Z
  }
  rescale_power <- function(Z) switch(power.scaling, log = log10(Z + 1e-12), sqrt = sqrt(Z), Z)

  masks <- lapply(cwts, function(cwt){ Z <- power_of(cwt); coi_mask_of(cwt, Z) })
  # scale limits are ALWAYS computed inside the COI (edge effects never drive the colours)
  valid_scaled <- lapply(seq_along(cwts), function(i){ Z <- rescale_power(winsorize(power_of(cwts[[i]]))); Z[masks[[i]]] })
  density_range <- range(unlist(valid_scaled), na.rm = TRUE)

  # per-individual global spectrum (time-average of power within the COI) and dominant period
  spectra <- list(); periods_list <- list(); dom <- data.frame()
  for(i in seq_along(cwts)){
    cwt <- cwts[[i]]; Zlin <- winsorize(power_of(cwt)); m <- masks[[i]]
    gs <- colSums(Zlin * m, na.rm = TRUE) / pmax(colSums(m), 1)
    periods <- cwt$scales * cwt$fourierfactor / time_factor
    spectra[[ord_ids[i]]] <- gs; periods_list[[ord_ids[i]]] <- periods
    k <- if(all(!is.finite(gs))) NA_integer_ else which.max(gs)
    dom <- rbind(dom, data.frame(id = ord_ids[i], dominant_period = if(is.na(k)) NA_real_ else periods[k],
                                 power = if(is.na(k)) NA_real_ else gs[k], stringsAsFactors = FALSE))
  }


  ##############################################################################
  # Layout #####################################################################
  ##############################################################################

  layout_params <- .setLayout(ncol, id.groups, plots.height = 6, dividers.height = 2, legend = FALSE)
  nplots <- max(layout_params$matrix, na.rm = TRUE)
  layout_params$matrix[is.na(layout_params$matrix)] <- nplots + 1

  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 5.5, slope = 4, n = max(0, ncol - 1), lo = 5.5, hi = 30),
                     h.rule = list(base = 1.5, slope = 2.6, n = ceiling(n_ind / ncol), lo = 4, hi = 30),
                     crowd.unit = "individuals")
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  layout(mat = layout_params$matrix, heights = layout_params$heights)
  oma_left <- if(length(id.groups) > 1) 4 else 2
  par(mar = c(4, 5, 3, 6), oma = c(2, oma_left, 3, 2), mgp = c(3, 0.8, 0))


  ##############################################################################
  # Draw panels ################################################################
  ##############################################################################

  for(i in seq_len(nplots)){
    if(i > n_ind){ plot.new(); next }
    cwt <- cwts[[i]]; id_name <- ord_ids[i]
    times <- built$series[[id_name]]$times
    Z <- rescale_power(winsorize(power_of(cwt)))

    zlim <- if(shared.scale) density_range else range(Z[masks[[i]]], na.rm = TRUE)
    breaks <- .powerBreaks(Z, zlim, masks[[i]], power.scaling, color.pal)

    graphics::image(x = seq_len(nrow(cwt$coefs)), y = seq_len(ncol(cwt$coefs)), z = Z,
                    zlim = zlim, breaks = breaks, col = color.pal, axes = FALSE, useRaster = TRUE,
                    main = "", xlab = "Date", ylab = paste0("Period (", unit_abbrev, ")"),
                    frame.plot = TRUE, cex.lab = cex_lab, xaxs = "i", yaxs = "i")
    title(main = id_name, cex.main = cex_title, line = 1)

    # x-axis: calendar-aware dates mapped onto the column-index axis
    ax <- .prettyDateAxis(min(times), max(times), n = 7, format = date.format)
    at_major <- vapply(as.numeric(ax$at), function(t) which.min(abs(as.numeric(times) - t)), integer(1))
    at_minor <- vapply(as.numeric(ax$minor), function(t) which.min(abs(as.numeric(times) - t)), integer(1))
    keep <- !duplicated(at_major)
    axis(1, at = at_major[keep], labels = ax$labels[keep], cex.axis = cex_axis)
    if(length(at_minor)) axis(1, at = unique(at_minor), labels = FALSE, tck = -0.02, lwd.ticks = 0.5)

    # y-axis: reference periods on the log-scale wavelet axis
    period_idx <- .rescale(log2(axis.periods * time_factor),
                           from = log2(range(cwt$scales * cwt$fourierfactor)), to = c(1, length(cwt$scales)))
    axis(2, at = period_idx, labels = axis.periods, cex.axis = cex_axis, las = 1)
    abline(h = period_idx, lty = 5, col = "grey70", lwd = 0.7)

    .drawCOI(cwt, mask.coi)

    if(!shared.scale || i == n_ind)
      .drawPowerLegend(color.pal, Z, zlim, masks[[i]], breaks, power.scaling, cex_legend)
    box()
  }

  mtext(main, side = 3, line = 1.6, outer = TRUE, cex = cex_title, font = 2)
  if(length(id.groups) > 1){
    label_pos <- grconvertY(1 - (layout_params$group_positions / sum(layout_params$heights)), "ndc", "user")
    text(x = grconvertX(0.01, "ndc", "user"), y = label_pos, labels = names(id.groups),
         srt = 90, cex = cex_title + 0.2, font = 2, xpd = NA, adj = c(0.5, 0.5))
  }

  attr(dom, "spectra") <- spectra; attr(dom, "periods") <- periods_list
  invisible(dom)
}


##################################################################################################
## Internal helpers ##############################################################################

#' Colour breaks for one scalogram panel (quantile scaling needs bespoke breaks)
#' @keywords internal
#' @noRd
.powerBreaks <- function(Z, zlim, mask, power.scaling, color.pal){
  n <- length(color.pal) + 1
  if(power.scaling == "quantile"){
    v <- Z[mask]; v <- v[is.finite(v)]; pos <- v[v > 0]
    if(length(pos) > 0){
      br <- unique(stats::quantile(pos, probs = seq(0, 1, length.out = n)))
      br[1] <- min(zlim); br <- sort(unique(br))
      if(length(br) == n) return(br)
    }
  }
  seq(min(zlim), max(zlim), length.out = n)
}

#' Draw the cone of influence (translucent veil or hard grey mask)
#' @keywords internal
#' @noRd
.drawCOI <- function(cwt, mask.coi){
  ff <- cwt$fourierfactor
  x <- seq_len(nrow(cwt$coefs))
  yr <- log2(range(cwt$scales * ff))
  coi_y <- .rescale(log2(cwt$coi_maxscale * ff + 1e-20), from = yr, to = c(1, length(cwt$scales)))
  y_top <- par("usr")[4]; y_bot <- par("usr")[3]
  coi_y <- pmin(pmax(coi_y, y_bot), y_top)
  col <- if(mask.coi) "grey20" else adjustcolor("white", alpha.f = 0.5)
  polygon(c(x, rev(x)), c(coi_y, rep(y_top, length(x))), col = col, border = TRUE)
  min_y <- .rescale(log2(cwt$coi_minscale * ff + 1e-20), from = yr, to = c(1, length(cwt$scales)))
  if(length(min_y) == 1 && is.finite(min_y) && min_y > y_bot)
    rect(par("usr")[1], y_bot, par("usr")[2], min_y, col = col, border = TRUE)
}

#' Draw the power colour bar for one panel (log/quantile ticks anti-transformed to raw power)
#' @keywords internal
#' @noRd
.drawPowerLegend <- function(color.pal, Z, zlim, mask, breaks, power.scaling, cex_legend){
  lab_sci <- FALSE; zlab <- NULL; digit <- 3
  fmt <- function(v){ out <- v
    sci <- which(abs(v) > 0 & (abs(v) < 0.1 | abs(v) >= 1000)); norm <- setdiff(seq_along(v), sci)
    out[sci] <- formatC(v[sci], format = "e", digits = 2); out[norm] <- formatC(v[norm], format = "f", digits = 2); out }
  if(power.scaling == "quantile"){
    zlim2 <- zlim; zval <- breaks[round(seq(1, length(breaks), length.out = 6))]; zlab <- fmt(zval)
  }else if(power.scaling == "log"){
    zlim2 <- zlim; zval <- pretty(zlim); zval <- zval[zval >= min(zlim) & zval <= max(zlim)]; zlab <- fmt(10^zval)
  }else{
    zlim2 <- zlim; zval <- pretty(zlim); zval <- zval[zval >= min(zlim) & zval <= max(zlim)]
    d <- .decimalPlaces(zval); digit <- if(all(is.na(d))) 2 else max(d, na.rm = TRUE); if(digit > 3) lab_sci <- TRUE
  }
  .colorlegend(col = color.pal, zlim = zlim2, zval = zval, zlab = zlab, posx = c(0.9, 0.915),
               posy = c(0.15, 0.85), main = "", digit = digit, lab.scientific = lab_sci,
               cex = cex_legend, xpd = NA)
}

#' @keywords internal
#' @noRd
.printScalogramSummary <- function(n_ids, n_total, dt, wavelet, gap.handling, detrend, power.scaling,
                                   period.range, unit){
  kv <- .kv
  .summaryOpen("Wavelet scalogram")
  kv("Individuals", sprintf("%d of %d (min.days filter)", n_ids, n_total))
  kv("Sampling", sprintf("dt = %g min", dt))
  kv("Wavelet", sprintf("%s; power: %s", wavelet, power.scaling))
  kv("Pre-processing", sprintf("gaps: %s; detrend: %s", gap.handling, detrend))
  kv("Period range", sprintf("%g-%g %s", period.range[1], period.range[2], unit))
  .summaryClose()
}

##################################################################################################
##################################################################################################
##################################################################################################
