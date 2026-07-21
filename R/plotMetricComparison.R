##################################################################################################
## Compare behavioural metrics across the levels of a grouping factor #############################
##################################################################################################

#' Compare metrics across the levels of a grouping factor
#'
#' @description Draws a multi-panel comparison of one or more per-individual behavioural metrics (e.g.
#' activity rate, distance moved, home-range area, association index) across the levels of a single
#' grouping factor - typically a within-individual time-frame such as diel phase or season. One panel
#' per metric shows the distribution of the per-individual values at each level (box, violin or
#' points), and, optionally, annotates a design-appropriate group comparison test.
#'
#' @details **Statistical testing.** Comparisons of a per-individual metric across time-frames are
#' *repeated-measures by design* (the same animals are measured under each level), so the default
#' tests respect within-individual pairing: **Wilcoxon signed-rank** for 2 levels and the
#' **Friedman test** for >2 levels (both on the individuals present in *all* compared levels;
#' Skillings-Mack is used for incomplete blocks if that package is installed). Independent-groups
#' tests (Mann-Whitney / Kruskal-Wallis) are used only with `paired = FALSE`. An **effect size**
#' (rank-biserial for signed-rank, Kendall's W for Friedman) is always computed, and p-values are
#' corrected across the metric panels (`p.adjust.method`). Tests are suppressed below `min.n`
#' complete pairs/blocks.
#'
#' This annotation is a **descriptive/exploratory aid, not confirmatory inference**: it is computed on
#' the complete-case subset, and telemetry non-detection is often *informative* (an animal absent at
#' night because it left the array), so the comparison can be biased in ways no in-plot test corrects.
#' For confirmatory inference, fit a mixed-effects model (individual random effect; a distribution
#' suited to the metric) and check its diagnostics. Significance flagging is therefore **off by
#' default**; the full test results are returned invisibly.
#'
#' The metrics must already be columns in `data`; compute them upstream (e.g. with the `calculate*`
#' family) so this function stays a pure visualiser. Multiple rows per individual x level are
#' aggregated by `agg.fun` (mean by default).
#'
#' @note The layout is sized in inch units, so the figure stays consistent across devices.
#'
#' @inheritParams as_moby
#' @param data A `mobyData` or data frame containing the metric column(s), an individual column and
#' the grouping column.
#' @param metrics Character vector of numeric column names to compare (one panel each).
#' @param split.by Name of the grouping column (the within-individual factor: diel phase, season, ...).
#' @param metric.labels Optional display names for the panels (defaults to the column names).
#' @param agg.fun Function aggregating multiple rows per individual x level to one value. Defaults to
#' the mean.
#' @param plot.type One of "box" (default), "violin" or "points".
#' @param paired Logical; treat the design as repeated-measures (within-individual). Defaults to TRUE.
#' Set FALSE only for genuine independent-groups factors (e.g. sex, population).
#' @param test One of "auto" (default; run the design-appropriate test), or "none".
#' @param p.adjust.method Multiple-comparison correction across the metric panels, passed to
#' \code{\link[stats]{p.adjust}}. Defaults to "holm".
#' @param min.n Minimum complete pairs/blocks required to run a test. Defaults to 6.
#' @param flag.significant Logical; if TRUE, outline panels whose adjusted p is below `alpha`.
#' Defaults to FALSE (avoids dichotomous "significant/not" reading).
#' @param alpha Significance threshold for `flag.significant`. Defaults to 0.05.
#' @param discard.incomplete Logical; if TRUE (default), drop individuals missing from any level for
#' the *plotted distributions* too (so the boxes match the tested set). If FALSE, all available
#' values are plotted while the test still uses complete blocks.
#' @param display.n Logical; append the per-level sample size to the x-axis labels. Defaults to FALSE.
#' @param color.pal Fill colours, one per level. If NULL, a colourblind-safe palette is used.
#' @param background.color Panel background colour. Defaults to "grey96".
#' @param outliers Logical; show boxplot outliers. Defaults to FALSE.
#' @param main Optional overall title.
#' @param cex Global expansion factor for all plot text. Defaults to 1.
#' @param ncol Number of columns in the panel layout. Defaults to 2.
#' @template deviceArgs
#'
#' @return Invisibly, a tidy per-metric data frame of the test results (test, paired, n's, statistic,
#' df, raw and adjusted p, effect size, method), with the per-(individual, level, metric) values
#' attached as `attr(., "values")`.
#' @seealso \code{\link{plotStationStats}}, \code{\link{plotGroupSizeDistribution}}
#' @examples
#' # Compare per-individual residency indices between the two species
#' res <- calculateResidency(rays, last.monitoring.date = max(rays$datetime))
#' res$species <- rays_tags$species[match(res$ID, rays_tags$ID)]
#' plotMetricComparison(res, metrics = c("IR1", "IR2"), split.by = "species",
#'                      paired = FALSE)
#' @export


plotMetricComparison <- function(data,
                                 metrics,
                                 split.by,
                                 id.col = NULL,
                                 metric.labels = NULL,
                                 agg.fun = function(x) mean(x, na.rm = TRUE),
                                 plot.type = c("box", "violin", "points"),
                                 paired = TRUE,
                                 test = c("auto", "none"),
                                 p.adjust.method = "holm",
                                 min.n = 6,
                                 flag.significant = FALSE,
                                 alpha = 0.05,
                                 discard.incomplete = TRUE,
                                 display.n = FALSE,
                                 color.pal = NULL,
                                 background.color = "grey96",
                                 outliers = FALSE,
                                 main = NULL,
                                 cex = 1,
                                 ncol = 2,
                                 file = NULL,
                                 width = NULL,
                                 height = NULL,
                                 res = 300) {

  ######################################################################################
  # Initial checks #####################################################################
  ######################################################################################

  data <- as.data.frame(data)
  id.col <- .resolveArgs(data, list(id.col = id.col))$id.col
  plot.type <- match.arg(plot.type)
  test <- match.arg(test)

  errors <- c()
  if(!id.col %in% colnames(data)) errors <- c(errors, sprintf("ID column '%s' not found; set 'id.col'.", id.col))
  if(missing(metrics) || length(metrics) == 0) errors <- c(errors, "'metrics' must name at least one column.")
  if(missing(split.by) || is.null(split.by)) errors <- c(errors, "'split.by' (the grouping column) is required.")
  miss <- setdiff(metrics, colnames(data))
  if(length(miss) > 0) errors <- c(errors, sprintf("metric(s) not found in data: %s", paste(miss, collapse = ", ")))
  if(!is.null(split.by) && !split.by %in% colnames(data)) errors <- c(errors, sprintf("'split.by' column '%s' not found.", split.by))
  if(!is.null(metric.labels) && length(metric.labels) != length(metrics)) errors <- c(errors, "'metric.labels' must match 'metrics' in length.")
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)

  if(!is.factor(data[[split.by]])){
    warning("Converting 'split.by' to a factor.", call. = FALSE)
    data[[split.by]] <- as.factor(data[[split.by]])
  }
  levels_ <- levels(droplevels(data[[split.by]]))
  if(length(levels_) < 2) stop("'split.by' must have at least 2 levels.", call. = FALSE)
  labels_ <- if(is.null(metric.labels)) metrics else metric.labels


  ######################################################################################
  # Aggregate to per-(individual x level) values + run tests ###########################
  ######################################################################################

  ids <- levels(factor(data[[id.col]]))
  mats <- lapply(metrics, function(m)
    tapply(data[[m]], list(factor(data[[id.col]], levels = ids), factor(data[[split.by]], levels = levels_)),
           agg.fun))
  names(mats) <- metrics

  # design-appropriate test per metric, then correct across panels
  stats_tab <- do.call(rbind, lapply(metrics, function(m){
    r <- .metricTest(mats[[m]], levels_, paired = paired, min.n = min.n)
    r <- cbind(metric = m, r, stringsAsFactors = FALSE); r
  }))
  if(test == "none"){ stats_tab$test <- "none"; stats_tab$p_raw <- NA_real_ }
  stats_tab$p_adj <- stats::p.adjust(stats_tab$p_raw, method = p.adjust.method)
  rownames(stats_tab) <- NULL

  # tidy per-(id, level, metric) values (for the return + plotting)
  values <- do.call(rbind, lapply(metrics, function(m){
    mm <- mats[[m]]
    data.frame(id = rep(rownames(mm), ncol(mm)), level = factor(rep(colnames(mm), each = nrow(mm)), levels = levels_),
               metric = m, value = as.vector(mm), stringsAsFactors = FALSE)
  }))
  values <- values[!is.na(values$value), ]


  ######################################################################################
  # Appearance & layout ################################################################
  ######################################################################################

  cex_title <- 1.15 * cex; cex_lab <- 1.0 * cex; cex_axis <- 0.85 * cex; cex_annot <- 0.75 * cex
  n_lev <- length(levels_)
  if(is.null(color.pal)) color.pal <- if(n_lev <= 8) .okabe_ito_pal(n_lev) else grDevices::hcl.colors(n_lev, "Dark 3")

  n_panels <- length(metrics)
  rows <- ceiling(n_panels / ncol)

  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 2.2, slope = 2.8, n = ncol, lo = 4.5, hi = 28),
                     h.rule = list(base = 1.8, slope = 2.8, n = rows, lo = 3.5, hi = 28),
                     crowd.unit = "panels")
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  par(mfrow = c(rows, ncol),
      mai = c(0.55, 0.75, 0.5, 0.2),
      omi = c(0.1, 0.1, if(!is.null(main)) 0.34 else 0.1, 0.1),
      mgp = c(2.4, 0.6, 0), lwd = 0.7)


  ######################################################################################
  # Console summary ####################################################################
  ######################################################################################

  .printMetricComparisonSummary(metrics = metrics, split.by = split.by, levels = levels_,
                                paired = paired, test = test, stats = stats_tab,
                                p.adjust.method = p.adjust.method)


  ######################################################################################
  # Draw panels ########################################################################
  ######################################################################################

  fmt_p <- function(p) if(is.na(p)) "NA" else if(p < 0.001) "<0.001" else sprintf("%.3f", p)

  for(i in seq_len(n_panels)){
    m <- metrics[i]
    st <- stats_tab[stats_tab$metric == m, ]
    # per-level value vectors (optionally restricted to complete blocks for the drawn distribution)
    mm <- mats[[m]]
    keep_rows <- if(discard.incomplete) rowSums(is.na(mm)) == 0 else rep(TRUE, nrow(mm))
    vals <- lapply(seq_len(n_lev), function(j){ v <- mm[keep_rows, j]; v[!is.na(v)] })
    ns <- vapply(vals, length, integer(1))

    finite <- unlist(vals); finite <- finite[is.finite(finite)]
    ylim <- if(length(finite)) range(finite) else c(0, 1)
    ypad <- diff(ylim) * 0.08; ylim <- c(ylim[1] - ypad, ylim[2] + ypad * 2.4)  # headroom for annotation

    plot(NA, xlim = c(0.5, n_lev + 0.5), ylim = ylim, axes = FALSE, xlab = "", ylab = "")
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = background.color, border = NA)

    for(j in seq_len(n_lev)){
      v <- vals[[j]]; col <- color.pal[((j - 1) %% length(color.pal)) + 1]
      if(length(v) == 0) next
      if(plot.type == "box"){
        graphics::boxplot(v, at = j, add = TRUE, axes = FALSE, col = col, border = "grey25",
                          boxwex = 0.5, outline = outliers, lwd = 0.8)
      }else if(plot.type == "violin"){
        .drawViolin(v, j, col)
        graphics::segments(j - 0.02, stats::median(v), j + 0.02, stats::median(v), lwd = 2)
      }else{
        graphics::points(jitter(rep(j, length(v)), amount = 0.12), v, pch = 16,
                         col = adjustcolor(col, 0.7), cex = 0.8 * cex)
        graphics::segments(j - 0.2, stats::median(v), j + 0.2, stats::median(v), lwd = 2, col = "grey25")
      }
    }

    # axes
    lab_x <- if(display.n) paste0(levels_, "\n(", ns, ")") else levels_
    axis(1, at = seq_len(n_lev), labels = lab_x, cex.axis = cex_axis, tick = FALSE, line = -0.2)
    yat <- pretty(ylim); yat <- yat[yat >= ylim[1] & yat <= ylim[2]]
    axis(2, at = yat, las = 1, cex.axis = cex_axis)
    title(main = labels_[i], cex.main = cex_title, font.main = 2, line = 1.1)

    # neutral test annotation (effect size + p at least as prominent as p; no stars)
    if(test != "none"){
      note <- if(!is.na(st$note)) st$note
              else sprintf("%s: %s=%.2f, p=%s (n=%d)", st$test, .effLabel(st$effect_type),
                           st$effect, fmt_p(st$p_adj), st$n_complete)
      mtext(note, side = 3, line = 0.1, cex = cex_annot, col = "grey20")
      if(flag.significant && !is.na(st$p_adj) && st$p_adj < alpha) box(col = "grey20", lwd = 2)
      else box(col = "grey55")
    }else box(col = "grey55")
  }

  if(!is.null(main)) mtext(main, side = 3, outer = TRUE, font = 2, cex = cex_title * 1.1, line = 0.2)

  attr(stats_tab, "values") <- values
  invisible(stats_tab)
}


##################################################################################################
## Internal helpers ##############################################################################

#' Short effect-size symbol for the annotation
#' @keywords internal
#' @noRd
.effLabel <- function(type) switch(as.character(type),
  "rank-biserial" = "r", "Kendall W" = "W", "Kendall W (approx)" = "W", "eta2_H" = "eta2", "e")

#' Draw a symmetric (kernel-density) violin for one level
#' @keywords internal
#' @noRd
.drawViolin <- function(v, at, col, width = 0.42){
  if(length(v) < 2 || diff(range(v)) == 0){
    graphics::segments(at - width * 0.6, stats::median(v), at + width * 0.6, stats::median(v), col = col, lwd = 3)
    return(invisible(NULL))
  }
  d <- stats::density(v)
  w <- d$y / max(d$y) * width
  graphics::polygon(c(at + w, rev(at - w)), c(d$x, rev(d$x)), col = adjustcolor(col, 0.85), border = "grey25", lwd = 0.7)
}


##################################################################################################
## Console summary (internal) ####################################################################

#' @keywords internal
#' @noRd
.printMetricComparisonSummary <- function(metrics, split.by, levels, paired, test, stats, p.adjust.method){
  kv <- .kv
  .summaryOpen("Metric comparison")
  kv("Metrics", paste(metrics, collapse = ", "))
  kv("Grouping", sprintf("%s (%s): %s", split.by, if(paired) "repeated-measures" else "independent groups",
                         paste(levels, collapse = ", ")))
  if(test != "none"){
    tests_used <- unique(stats$test[!is.na(stats$test)])
    kv("Test", sprintf("%s; correction: %s", paste(tests_used, collapse = "/"), p.adjust.method))
    dropped <- stats$metric[stats$n_dropped > 0 & !is.na(stats$n_dropped)]
    if(length(dropped) > 0){
      frac <- max(stats$n_dropped / pmax(stats$n_total, 1), na.rm = TRUE)
      kv("Incomplete", sprintf("%d/%d individuals dropped (max %.0f%%); complete-case only",
                               max(stats$n_dropped, na.rm = TRUE), max(stats$n_total, na.rm = TRUE), 100 * frac))
      if(frac > 0.3) cat("  ! large data loss - non-detection may be informative; consider a mixed model.\n")
    }
  } else kv("Test", "none")
  .summaryClose()
}

##################################################################################################
##################################################################################################
##################################################################################################
