##################################################################################################
## Internal helpers: signal / time-series construction ###########################################
## Shared by the spectral plots (plotPeriodogram, plotScalogram): turn a gappy per-individual      ##
## detection log into an analysis-ready, regularly-sampled series.                                 ##
##################################################################################################
## The single source of truth for how a binned detection series becomes a signal. Both spectral
## functions call it, so gap handling, detrending and the sampling-interval inference live in ONE
## place (previously duplicated and inconsistent - the periodogram's copy dropped NA bins, which
## collapses temporal gaps and corrupts the frequency axis).


#' Build regularly-sampled per-individual signal series from binned detections
#'
#' @description Infers the sampling interval, lays every individual on one complete regular time grid
#' (so gaps are *represented*, not deleted), optionally binarises to presence, fills internal gaps,
#' detrends and standardises. Pre/first- and post/last-detection regions are trimmed (they are not
#' internal gaps). Returns one ready-to-analyse series per retained individual, plus the sampling
#' interval, timezone and per-individual diagnostics.
#'
#' @param data A data frame of binned detections.
#' @param id.col,timebin.col Column names for the individual ID and the (POSIXct) time bin.
#' @param value.col Name of the numeric column to analyse (e.g. "detections", or any metric).
#' @param binary Logical; if TRUE, observed values are reduced to presence (>=1 -> 1) before gap fill.
#' @param gap.handling One of "zero" (absence = 0), "mean", "locf", "interpolate".
#' @param detrend.method One of "none" (demean only), "linear" (remove OLS trend), "diff" (first
#' difference) or "loess". `diff` and `loess` are high-pass and can suppress the very low-frequency
#' (diel/tidal) peaks these tools look for - hence "none" is the default upstream.
#' @param loess.span Span for `detrend.method = "loess"`.
#' @param standardize Logical; z-score each series after detrending.
#' @param min.days Optional; drop individuals detected on fewer than this many distinct days.
#' @return A list: `series` (named list per retained id of `list(values, times)`), `dt` (sampling
#' interval, minutes), `tz`, `ids` (retained), `n_total`, `diagnostics` (per-id data frame with the
#' number of bins, internal-gap fraction and days detected).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.buildRhythmSeries <- function(data, id.col, timebin.col, value.col, binary = FALSE,
                               gap.handling = "zero", detrend.method = "none", loess.span = 0.75,
                               standardize = FALSE, min.days = NULL) {

  tz <- .dataTZ(data[[timebin.col]])
  data[[id.col]] <- droplevels(as.factor(data[[id.col]]))

  # sampling interval (minutes) from the smallest positive gap between distinct time bins
  sorted <- sort(unique(data[[timebin.col]]))
  if(length(sorted) < 2) stop("At least two unique time bins are required to infer the sampling interval.", call. = FALSE)
  d <- as.numeric(diff(sorted), units = "mins")
  dt <- min(d[d > 0])
  if(!is.finite(dt) || dt <= 0) stop("Could not determine a valid sampling interval (dt > 0). Check 'timebin.col'.", call. = FALSE)

  # complete regular grid across the study span (so gaps are represented, not deleted)
  grid <- seq(min(sorted), max(sorted), by = dt * 60)

  n_total <- nlevels(data[[id.col]])
  ids_all <- levels(data[[id.col]])

  # per-individual diagnostics (days detected, for the min.days filter and reporting)
  day <- strftime(data[[timebin.col]], "%Y-%m-%d", tz = tz)
  days_detected <- tapply(day, data[[id.col]], function(x) length(unique(x)))

  keep <- if(is.null(min.days)) ids_all else ids_all[which(as.numeric(days_detected) >= min.days)]

  series <- list(); diag_rows <- list()
  for(id in keep){
    sub <- data[data[[id.col]] == id, c(timebin.col, value.col)]
    if(any(duplicated(sub[[timebin.col]])))
      sub <- stats::aggregate(sub[[value.col]], by = list(sub[[timebin.col]]), FUN = mean, na.rm = TRUE)
    colnames(sub) <- c("t", "v")

    v <- rep(NA_real_, length(grid))
    v[match(sub$t, grid)] <- sub$v
    if(binary){ v[!is.na(v) & v >= 1] <- 1; v[!is.na(v) & v < 1] <- 0 }

    # trim pre/post-monitoring region (outer NAs) - not internal gaps
    obs <- which(!is.na(v))
    if(length(obs) == 0) next
    rng <- obs[1]:obs[length(obs)]
    v <- v[rng]; times <- grid[rng]
    n_gap <- sum(is.na(v)); gap_frac <- n_gap / length(v)

    # fill internal gaps
    v <- switch(gap.handling,
      zero        = { v[is.na(v)] <- 0; v },
      mean        = { v[is.na(v)] <- mean(v, na.rm = TRUE); v },
      locf        = .naLocf(v),
      interpolate = .naApprox(v),
      v)

    # detrend (demean, or remove a trend) then optionally standardise
    ti <- seq_along(v)
    v <- switch(detrend.method,
      none   = v - mean(v, na.rm = TRUE),
      linear = stats::residuals(stats::lm(v ~ ti)),
      diff   = diff(v),
      loess  = { fit <- try(stats::loess(v ~ ti, span = loess.span), silent = TRUE)
                 if(inherits(fit, "try-error")) v - mean(v, na.rm = TRUE) else stats::residuals(fit) },
      v - mean(v, na.rm = TRUE))
    if(standardize) v <- as.numeric(scale(v))

    series[[id]] <- list(values = v, times = if(detrend.method == "diff") times[-1] else times)
    diag_rows[[id]] <- data.frame(id = id, n = length(v), gap_frac = gap_frac,
                                  days = as.numeric(days_detected[id]), stringsAsFactors = FALSE)
  }

  list(series = series, dt = dt, tz = tz, ids = names(series), n_total = n_total,
       diagnostics = if(length(diag_rows)) do.call(rbind, diag_rows) else NULL)
}

##################################################################################################
##################################################################################################
##################################################################################################
