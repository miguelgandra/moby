#######################################################################################################
# Diagnostic for the min_lag false-detection threshold ################################################
#######################################################################################################

#' Diagnose the min_lag false-detection threshold
#'
#' @description
#' Plots the empirical distribution of `min_lag` - each detection's time gap to the nearest other
#' detection of the same transmitter on the same receiver - normalised by the transmitter nominal
#' delay, and reports the proportion of detections that the short-interval false-detection filter of
#' [filterDetections()] would flag at several candidate thresholds. It is a quick check on whether the
#' default `min.lag.factor = 30` is appropriate for a given dataset, rather than assuming it.
#'
#' @details
#' A genuinely present tag produces bursts of closely spaced detections (small `min_lag`), whereas an
#' isolated spurious decode has no nearby companion (large `min_lag`). The two typically appear as
#' separate modes on the log axis, and a defensible threshold sits in the valley between them. Read the
#' plot rather than a single number: if the proportion flagged changes little across the reference
#' multipliers, the choice is robust; if it changes sharply, the threshold matters and should be
#' justified for that dataset.
#'
#' `min_lag` is computed with the same internal helper used by [filterDetections()], so - on the same
#' input - the flagged proportions reported here match what the filter would remove. The diagnostic
#' applies only the optional duplicate-removal step (not the pre-tagging / cut-off steps), so run it on
#' the detections you intend to filter.
#'
#' @param data A `mobyData` object or a data frame of detections.
#' @param nominal.delay Transmitter nominal (mean) delay, in seconds. A single value applied to all
#'   individuals, or a vector named by `id.col` for mixed tag families. Read from the `mobyData`
#'   metadata (`nominal.delay`) when not supplied. Required: the diagnostic normalises `min_lag` by it.
#' @param factors Reference threshold multipliers (of the nominal delay) to overlay on the plot and
#'   tabulate. Defaults to `c(20, 30, 50)`; the middle value is the [filterDetections()] default.
#' @param id.col Name of the column with animal IDs. Resolved from the metadata / canonical default
#'   ("ID") when not supplied.
#' @param datetime.col Name of the column with detection timestamps (`POSIXct`). Resolved from the
#'   metadata / canonical default ("datetime") when not supplied.
#' @param station.col Name of the receiver/station column. Resolved from the metadata / canonical
#'   default ("station") when not supplied.
#' @param remove.duplicates Logical; drop exact-duplicate records (same ID, timestamp and station)
#'   before computing `min_lag`, matching the stage-0 behaviour of [filterDetections()]. Defaults to TRUE.
#' @param bar.color Fill colour for the histogram. Defaults to a colourblind-safe blue.
#' @param background.color Panel background colour. Defaults to "grey96"; `NA` draws none.
#' @param breaks Number of histogram bins on the log10 axis. Defaults to 40.
#' @param cex Global expansion factor for labels. Defaults to 1.
#' @param main Optional plot title.
#' @template deviceArgs
#'
#' @return Invisibly, a data frame with one row per `factors` value: the multiplier (`factor`), the
#'   threshold in seconds (`threshold_s`; `NA` when `nominal.delay` varies across individuals), and the
#'   number (`n_flagged`) and percentage (`pct_flagged`) of detections flagged - i.e. with `min_lag`
#'   above the threshold, or a lone decode - at that multiplier.
#'
#' @seealso [filterDetections()]
#'
#' @examples
#' # assess whether the default 30x min_lag threshold suits the data (120 s nominal delay)
#' flagged <- plotMinLag(rays, nominal.delay = 120)
#' flagged
#'
#' @export
plotMinLag <- function(data,
                       nominal.delay = NULL,
                       factors = c(20, 30, 50),
                       id.col = NULL,
                       datetime.col = NULL,
                       station.col = NULL,
                       remove.duplicates = TRUE,
                       bar.color = "#0072B2",
                       background.color = "grey96",
                       breaks = 40,
                       cex = 1,
                       main = NULL,
                       file = NULL,
                       width = NULL,
                       height = NULL,
                       res = 300) {

  ##############################################################################
  # Resolve arguments ##########################################################
  ##############################################################################

  prev_meta <- attr(data, "moby")
  if (is.null(nominal.delay) && !is.null(prev_meta) && !is.null(prev_meta$nominal.delay))
    nominal.delay <- prev_meta$nominal.delay

  # resolves id.col / datetime.col / station.col from metadata or canonical defaults, validates the
  # columns exist, and coerces id.col to a factor
  reviewed <- .validateArguments()
  data <- as.data.frame(reviewed$data)

  if (is.null(nominal.delay))
    .mobyAbort("'nominal.delay' is required (a single value, a vector named by ID, or stored in the ",
               "mobyData metadata): the diagnostic normalises min_lag by the nominal delay.")
  if (!is.numeric(factors) || length(factors) < 1 || any(!is.finite(factors)) || any(factors <= 0))
    .mobyAbort("'factors' must be one or more positive numbers (threshold multipliers of the nominal delay).")
  if (!inherits(data[[datetime.col]], "POSIXct"))
    .mobyAbort("'", datetime.col, "' must be in POSIXct format.")
  if (anyNA(data[[datetime.col]]))
    .mobyAbort("'", datetime.col, "' contains missing value(s); please remove or fix them first.")

  ##############################################################################
  # Compute min_lag ############################################################
  ##############################################################################

  # stage 0: exact-duplicate removal (matches filterDetections)
  if (isTRUE(remove.duplicates)) {
    key_cols <- c(id.col, datetime.col, station.col)
    data <- data[!duplicated(data[, key_cols, drop = FALSE]), , drop = FALSE]
  }

  # per-ID nominal delay (tolerant of unknown tags -> NA, those detections are excluded)
  cp <- .checkAnimalParams(nominal.delay, "Nominal delay", expected_class = "numeric", data, id.col,
                           allow.missing = TRUE)
  if (!is.null(cp$errors)) .mobyAbort(paste(cp$errors, collapse = "\n"))
  nd_by_level <- cp$vector
  names(nd_by_level) <- levels(data[, id.col])

  # order by (ID, time) so min_lag sees each tag's detections in ascending order (as the filter does)
  data <- data[order(data[, id.col], data[, datetime.col]), , drop = FALSE]
  nd_row <- nd_by_level[as.character(data[, id.col])]

  keep <- !is.na(nd_row)
  if (!any(keep))
    .mobyAbort("No detections have a known nominal delay; cannot build the diagnostic.")
  if (!all(keep))
    message("moby: ", sum(!keep), " detection(s) from tag(s) with no nominal delay were excluded.")
  data <- data[keep, , drop = FALSE]; nd_row <- nd_row[keep]

  # nearest same-receiver gap per detection, computed within each individual (shared helper)
  minlag <- rep(NA_real_, nrow(data))
  ids_chr <- as.character(data[, id.col])
  for (lv in unique(ids_chr)) {
    idx <- which(ids_chr == lv)
    minlag[idx] <- .minLagValues(as.numeric(data[idx, datetime.col]), as.character(data[idx, station.col]))
  }
  norm <- minlag / nd_row                               # min_lag in units of the nominal delay

  ##############################################################################
  # Flagged proportion at each reference multiplier ############################
  ##############################################################################

  n_total <- nrow(data)
  n_lone  <- sum(is.na(minlag))                          # lone decodes: flagged at every threshold
  single_nominal <- length(unique(nd_row)) == 1
  results <- data.frame(
    factor      = factors,
    threshold_s = if (single_nominal) factors * nd_row[1] else NA_real_,
    n_flagged   = vapply(factors, function(f) sum(norm > f, na.rm = TRUE) + n_lone, numeric(1)),
    row.names   = NULL)
  results$pct_flagged <- round(100 * results$n_flagged / n_total, 2)

  ##############################################################################
  # Draw #######################################################################
  ##############################################################################

  lognl <- log10(norm[is.finite(norm) & norm > 0])
  if (length(lognl) < 2)
    .mobyAbort("Too few finite min_lag values to plot (need at least 2 detections that share a receiver).")

  if (!is.null(file)) {
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 7, slope = 0, n = 0, lo = 5, hi = 12),
                     h.rule = list(base = 5, slope = 0, n = 0, lo = 4, hi = 10))
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  par(mar = c(4.4, 4.4, if (is.null(main)) 2.4 else 3.2, 1.2), mgp = c(2.6, 0.7, 0))

  rng <- range(lognl); if (diff(rng) < .Machine$double.eps^0.5) rng <- rng + c(-0.5, 0.5)
  br <- seq(rng[1], rng[2], length.out = max(2L, breaks) + 1)
  h  <- hist(lognl, breaks = br, plot = FALSE)

  plot(NA, xlim = range(br), ylim = c(0, max(h$counts) * 1.12), xaxt = "n",
       xlab = "min_lag / nominal delay", ylab = "detections",
       main = if (is.null(main)) "min_lag distribution" else main,
       cex.lab = 1.1 * cex, cex.axis = cex, cex.main = 1.2 * cex)
  if (!is.na(background.color))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = background.color, border = NA)
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$counts, col = adjustcolor(bar.color, 0.7), border = NA)

  # log x-axis at human-readable multipliers of the nominal delay
  ticks <- c(0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000)
  ticks <- ticks[log10(ticks) >= min(br) & log10(ticks) <= max(br)]
  axis(1, at = log10(ticks), labels = ticks, cex.axis = cex)

  # reference threshold lines (the filterDetections default, 30x, is emphasised)
  is_default <- abs(factors - 30) < 1e-9
  for (k in seq_along(factors)) {
    xf <- log10(factors[k])
    if (xf < min(br) || xf > max(br)) next
    abline(v = xf, col = "grey20", lwd = if (is_default[k]) 2.2 else 1.6, lty = if (is_default[k]) 1 else 2)
  }
  box()

  # legend reports the proportion flagged at each threshold (robust to close/overlapping factors,
  # unlike inline labels)
  lab <- paste0(formatC(factors, format = "g"), "x")
  lab[is_default] <- paste0(lab[is_default], " (default)")
  legend("topright", bty = "n", cex = 0.9 * cex, title = "detections flagged", title.adj = 0,
         legend = paste0(lab, ":  ", sprintf("%.1f%%", results$pct_flagged)),
         lty = ifelse(is_default, 1, 2), lwd = ifelse(is_default, 2.2, 1.6), col = "grey20",
         seg.len = 1.6, text.col = "grey15")

  invisible(results)
}
