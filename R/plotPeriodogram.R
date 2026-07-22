##################################################################################################
## Periodogram: dominant rhythms in per-individual detection series ###############################
##################################################################################################

#' Plot detection periodograms (dominant rhythms)
#'
#' @description Computes and plots a periodogram per individual - spectral power against period - to
#' reveal dominant cyclic patterns in a detection time series, such as tidal (~6-12 h), diel (24 h)
#' or lunar rhythms in habitat use. One small panel is drawn per animal, with the tidal and diel
#' bands shaded and the dominant period marked. It is the global-spectrum complement of
#' \code{\link{plotScalogram}} (which resolves how those rhythms change through time).
#'
#' @details The detection series is first laid on a complete, regularly-sampled time grid (via the
#' shared signal builder), so temporal gaps are *represented* rather than deleted - the ordinary FFT
#' assumes even sampling, and dropping absent bins corrupts the period axis. By default `method =
#' "fft"` computes a raw periodogram (\code{\link[stats]{spec.pgram}}) on that grid, treating an
#' unobserved bin as absence (0). For data where absence-as-zero is not appropriate, `method =
#' "lomb"` runs a Lomb-Scargle periodogram on the observed (time, value) pairs (requires the
#' \pkg{lomb} package) and reports an analytic false-alarm probability for the dominant peak.
#'
#' Detrending defaults to `"none"` (mean removal only): the common first-difference detrend is a
#' high-pass filter that *suppresses* exactly the low-frequency diel/tidal peaks this tool looks for.
#'
#' @inheritParams as_moby
#' @param data A `mobyData` or data frame of binned detections (must contain a `detections` column).
#' @param type Response analysed: `"detections"` (counts) or `"presences"` (binary per bin).
#' @param method `"fft"` (default; periodogram on the gap-filled regular grid) or `"lomb"`
#' (Lomb-Scargle on the observed samples; needs the \pkg{lomb} package).
#' @param detrend One of `"none"` (default; demean), `"linear"` (remove OLS trend), or `"diff"`
#' (first difference - high-pass, use with care).
#' @param gap.handling How internal gaps in the regular grid are filled for `method = "fft"`: `"zero"`
#' (absence; default), `"mean"`, `"locf"` or `"interpolate"`.
#' @param period.range Optional length-2 numeric (hours) giving the period range to display and search
#' for the dominant peak. If NULL (default), spans ~2 sampling intervals to twice the largest
#' `axis.periods` value.
#' @param axis.periods Reference periods (hours) marked on the x-axis. Defaults to c(48, 24, 12, 8, 6).
#' @param highlight.bands Logical; shade the tidal (~6-12.4 h) and diel (~24 h) bands. Defaults to TRUE.
#' @param min.days Discard individuals detected on fewer than this many distinct days. Defaults to 10.
#' @param id.groups Optional named list of ID groups (one block of panels each).
#' @param shared.scale Logical; share the y-axis (spectral power) across panels. Defaults to TRUE.
#' @param color.pal Colour for the spectrum. If NULL, a colourblind-safe colour is used.
#' @param background.color Panel background colour. Defaults to "grey96".
#' @param cex Global expansion factor for all plot text. Defaults to 1.
#' @param ncol Number of panel columns. If NULL, set from the number of individuals.
#' @template deviceArgs
#'
#' @return Invisibly, a tidy data frame with one row per individual: the dominant period (h), its
#' power, the false-alarm probability (`method = "lomb"` only), and whether it falls in the tidal or
#' diel band.
#' @seealso \code{\link{plotScalogram}}, \code{\link{plotActograms}}, \code{\link{plotChronogram}}
#' @examples
#' # Collapse detections into hourly counts per individual, then find dominant rhythms
#' binned <- aggregate(list(detections = rep(1L, nrow(rays))),
#'                     by = list(ID = rays$ID, timebin = rays$timebin), FUN = sum)
#' plotPeriodogram(binned, id.col = "ID", timebin.col = "timebin")
#' @export


plotPeriodogram <- function(data,
                            id.col = NULL,
                            timebin.col = NULL,
                            id.groups = NULL,
                            type = "detections",
                            method = c("fft", "lomb"),
                            detrend = c("none", "linear", "diff"),
                            gap.handling = c("zero", "mean", "locf", "interpolate"),
                            period.range = NULL,
                            axis.periods = c(48, 24, 12, 8, 6),
                            highlight.bands = TRUE,
                            min.days = 10,
                            shared.scale = TRUE,
                            color.pal = NULL,
                            background.color = "grey96",
                            ncol = NULL,
                            cex = 1,
                            file = NULL,
                            width = NULL,
                            height = NULL,
                            res = 300) {

  ##############################################################################
  # Checks #####################################################################
  ##############################################################################

  reviewed <- .validateArguments()
  data <- as.data.frame(reviewed$data)
  method <- match.arg(method); detrend <- match.arg(detrend); gap.handling <- match.arg(gap.handling)

  errors <- c()
  if(!type %in% c("detections", "presences")) errors <- c(errors, "'type' must be 'detections' or 'presences'.")
  if(!"detections" %in% colnames(data)) errors <- c(errors, "A 'detections' column is required (binned detection counts).")
  if(method == "lomb" && !requireNamespace("lomb", quietly = TRUE))
    errors <- c(errors, "method = 'lomb' requires the 'lomb' package (install.packages('lomb')).")
  if(!is.null(period.range) && (!is.numeric(period.range) || length(period.range) != 2))
    errors <- c(errors, "'period.range' must be a length-2 numeric (hours).")
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)

  if(is.null(id.groups)) id.groups <- list(levels(factor(data[, id.col])))

  cex_title <- 2.0 * cex; cex_lab <- 1.5 * cex; cex_axis <- 1.1 * cex
  if(is.null(color.pal)) color.pal <- .okabe_ito_pal(1)


  ##############################################################################
  # Build series + compute periodograms ########################################
  ##############################################################################

  built <- .buildRhythmSeries(data, id.col, timebin.col, value.col = "detections",
                              binary = (type == "presences"), gap.handling = gap.handling,
                              detrend.method = detrend, min.days = min.days)
  ids <- built$ids
  if(length(ids) == 0) stop("No individuals passed the 'min.days' filter.", call. = FALSE)
  id.groups <- lapply(id.groups, function(x) x[x %in% ids]); id.groups <- id.groups[lengths(id.groups) > 0]

  dt_h <- built$dt / 60
  if(is.null(period.range)) period.range <- c(2 * dt_h, 2 * max(axis.periods))

  pgrams <- lapply(ids, function(id) .computePeriodogram(built$series[[id]], built$dt, method, period.range))
  names(pgrams) <- ids
  ymax <- max(vapply(pgrams, function(p) if(length(p$power)) max(p$power) else 0, numeric(1)), na.rm = TRUE)


  ##############################################################################
  # Layout #####################################################################
  ##############################################################################

  n_ind <- length(ids)
  if(is.null(ncol)) ncol <- if(n_ind == 1) 1 else 2
  layout_params <- .setLayout(ncol, id.groups, plots.height = 6, dividers.height = 2, legend = FALSE)
  bottom_indices <- sort(unlist(apply(layout_params$matrix, 2, function(x){
    r <- rle(!is.na(x)); x[cumsum(r$lengths)[r$values]] }, simplify = FALSE)))
  nplots <- max(layout_params$matrix, na.rm = TRUE)
  layout_params$matrix[is.na(layout_params$matrix)] <- nplots + 1
  first_col <- layout_params$matrix[, 1]

  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 4.5, slope = 3.2, n = max(0, ncol - 1), lo = 4.5, hi = 30),
                     h.rule = list(base = 2, slope = 2.4, n = ceiling(n_ind / ncol), lo = 4.5, hi = 30),
                     crowd.unit = "individuals")
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  layout(mat = layout_params$matrix, heights = layout_params$heights)
  par(mar = c(1.6, if(length(id.groups) > 1) 6 else 5.5, 0.4, 1), oma = c(4, if(length(id.groups) > 1) 4 else 1, 1, 1), mgp = c(3, 0.8, 0))


  ##############################################################################
  # Summary + draw #############################################################
  ##############################################################################

  band_of <- function(p) if(is.na(p)) NA_character_ else if(p >= 6 & p <= 12.4) "tidal" else if(p >= 22 & p <= 26) "diel" else "other"
  results <- do.call(rbind, lapply(ids, function(id){
    p <- pgrams[[id]]
    data.frame(id = id, dominant_period = p$dom_period, power = p$dom_power, fap = p$fap,
               band = band_of(p$dom_period), stringsAsFactors = FALSE)
  }))

  .printPeriodogramSummary(n_ids = n_ind, n_total = built$n_total, dt = built$dt, method = method,
                           detrend = detrend, period.range = period.range, results = results)

  tidal <- c(6, 12.4); diel <- c(22, 26)
  for(i in seq_len(nplots)){
    id <- ids[i]; p <- pgrams[[id]]
    ylim <- if(shared.scale) c(0, ymax) else c(0, if(length(p$power)) max(p$power) else 1)

    plot(NA, xlim = period.range, ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = background.color, border = NA)
    if(highlight.bands){
      rect(tidal[1], par("usr")[3], tidal[2], par("usr")[4], col = adjustcolor("#0072B2", 0.12), border = NA)
      rect(diel[1],  par("usr")[3], diel[2],  par("usr")[4], col = adjustcolor("#D55E00", 0.12), border = NA)
    }
    if(length(p$power)){
      graphics::segments(p$period, 0, p$period, p$power, col = color.pal[1], lwd = 1)
      if(!is.na(p$dom_period)) points(p$dom_period, p$dom_power, pch = 25, bg = "grey15", col = "grey15", cex = 0.9 * cex, xpd = NA)
    }

    if(!shared.scale || i %in% first_col){
      title(ylab = "Spectral power", line = 4, cex.lab = cex_lab, xpd = NA)
      axis(2, at = pretty(ylim), las = 1, cex.axis = cex_axis, col = NA, col.ticks = 1)
    }
    if(i %in% bottom_indices){
      title(xlab = "Period (h)", line = 2.6, cex.lab = cex_lab, xpd = NA)
      axis(1, at = axis.periods, labels = axis.periods, cex.axis = cex_axis, col = NA, col.ticks = 1)
    }
    legend("topright", legend = id, text.font = 2, cex = cex_title, bty = "n")
    box()
  }

  if(length(id.groups) > 1){
    label_pos <- grconvertY(1 - (layout_params$group_positions / sum(layout_params$heights)), "ndc", "user")
    text(x = grconvertX(0.01, "ndc", "user"), y = label_pos,
         labels = names(id.groups), srt = 90, cex = cex_title + 0.2, font = 2, xpd = NA, adj = c(0.5, 0.5))
  }

  invisible(results)
}


##################################################################################################
## Internal helpers ##############################################################################

#' Compute a periodogram (fft or Lomb-Scargle) restricted to a period range (hours)
#' @keywords internal
#' @noRd
.computePeriodogram <- function(series, dt, method, period.range){
  x <- series$values
  if(length(x) < 4 || stats::sd(x, na.rm = TRUE) == 0)
    return(list(period = numeric(0), power = numeric(0), fap = NA_real_, dom_period = NA_real_, dom_power = NA_real_))

  if(method == "lomb"){
    tnum <- as.numeric(series$times) / 60                     # minutes
    ls <- lomb::lsp(x = x, times = tnum, type = "period",
                    from = period.range[1] * 60, to = period.range[2] * 60, plot = FALSE)
    period_h <- ls$scanned / 60; power <- ls$power; fap <- ls$p.value
  }else{
    sp <- stats::spec.pgram(stats::ts(x), taper = 0.05, detrend = FALSE, fast = TRUE, plot = FALSE)
    period_h <- (1 / sp$freq) * dt / 60; power <- sp$spec; fap <- NA_real_
  }
  keep <- is.finite(period_h) & period_h >= period.range[1] & period_h <= period.range[2]
  period_h <- period_h[keep]; power <- power[keep]
  ord <- order(period_h); period_h <- period_h[ord]; power <- power[ord]
  dom <- if(length(power)) which.max(power) else NA_integer_
  list(period = period_h, power = power, fap = fap,
       dom_period = if(length(power)) period_h[dom] else NA_real_,
       dom_power  = if(length(power)) power[dom] else NA_real_)
}

#' @keywords internal
#' @noRd
.printPeriodogramSummary <- function(n_ids, n_total, dt, method, detrend, period.range, results){
  kv <- .kv
  .summaryOpen("Detection periodogram")
  kv("Individuals", sprintf("%d of %d (min.days filter)", n_ids, n_total))
  kv("Sampling", sprintf("dt = %g min", dt))
  kv("Method", sprintf("%s; detrend: %s", if(method == "lomb") "Lomb-Scargle" else "FFT periodogram", detrend))
  kv("Period range", sprintf("%.1f-%.1f h", period.range[1], period.range[2]))
  nt <- sum(results$band == "tidal", na.rm = TRUE); nd <- sum(results$band == "diel", na.rm = TRUE)
  kv("Dominant peak", sprintf("%d tidal, %d diel, %d other", nt, nd, n_ids - nt - nd))
  .summaryClose()
}

##################################################################################################
##################################################################################################
##################################################################################################
