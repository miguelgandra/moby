##################################################################################################
## Internal helpers: statistics ##################################################################
## Design-appropriate group-comparison tests + effect sizes for plotMetricComparison();          ##
## permutation-based empirical p-values for the network randomisers.                             ##
##################################################################################################


#' Empirical (permutation) p-values with multiple-comparison adjustment
#'
#' @description For each observed statistic, computes a one- or two-sided empirical p-value against
#' its null distribution (a row of `null_mat`), applies a multiple-comparison adjustment across the
#' set, and labels the direction of the deviation. Shared by \code{randomizeTransitions()} and
#' \code{randomizeAssociations()} so both network randomisers use one correct inference routine.
#'
#' @details The one-sided p-values use the standard permutation-test tail counts:
#' `"greater"` = (\#{null >= obs} + 1) / (n + 1) (right tail; observed significantly greater than the
#' null), `"less"` = (\#{null <= obs} + 1) / (n + 1). The two-sided p-value is `min(2 * min(left,
#' right), 1)` (capped at 1). NA null values are dropped per row, and rows with no valid data (or an
#' NA observed value) return NA. Set `p.adjust.method = "none"` to keep raw p-values.
#'
#' @param observed Numeric vector of observed statistics.
#' @param null_mat Numeric matrix of null draws, one row per observed statistic (statistics x iterations).
#' @param alternative One of "two.sided" (default), "greater", "less".
#' @param p.adjust.method Passed to \code{\link[stats]{p.adjust}} ("fdr" default; "none" to skip).
#' @param alpha Significance threshold applied to the ADJUSTED p-value for the direction label.
#' @param labels Named length-2 character vector giving the labels for the two deviation directions,
#' as `c(more = ..., less = ...)`. Defaults to `c(more = "more", less = "less")`; the association
#' randomiser passes `c(more = "positive", less = "negative")`.
#' @return A data frame with `mean_null`, `p_value`, `p_adjusted` and `association` (the direction
#' label where the adjusted p-value is below `alpha`, else "non-significant"; NA where undefined).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.empiricalPvalues <- function(observed, null_mat, alternative = "two.sided",
                              p.adjust.method = "fdr", alpha = 0.05,
                              labels = c(more = "more", less = "less")) {
  n <- length(observed)
  pvals <- numeric(n); direction <- character(n)
  for (i in seq_len(n)) {
    nd <- null_mat[i, ]; nd <- nd[!is.na(nd)]; obs <- observed[i]
    if (is.na(obs) || length(nd) == 0) { pvals[i] <- NA_real_; direction[i] <- NA_character_; next }
    m <- length(nd)
    if (alternative == "greater")      pvals[i] <- (sum(nd >= obs) + 1) / (m + 1)
    else if (alternative == "less")    pvals[i] <- (sum(nd <= obs) + 1) / (m + 1)
    else {
      pl <- (sum(nd >= obs) + 1) / (m + 1); pr <- (sum(nd <= obs) + 1) / (m + 1)
      pvals[i] <- min(2 * min(pl, pr), 1)
    }
    direction[i] <- if (obs > mean(nd)) labels[["more"]] else labels[["less"]]
  }
  padj <- stats::p.adjust(pvals, method = p.adjust.method)
  assoc <- ifelse(is.na(padj), NA_character_, ifelse(padj < alpha, direction, "non-significant"))
  data.frame(mean_null = rowMeans(null_mat, na.rm = TRUE), p_value = pvals, p_adjusted = padj,
             association = assoc, stringsAsFactors = FALSE)
}


##################################################################################################
## Design-appropriate group-comparison tests + effect sizes for plotMetricComparison().          ##
##################################################################################################
## Telemetry metric comparisons across diel/seasonal levels are REPEATED-MEASURES by design (the
## same individuals measured under each level). The default tests therefore respect within-individual
## pairing: Wilcoxon signed-rank (2 levels) / Friedman (>2 levels), both requiring complete blocks;
## Skillings-Mack is used for incomplete blocks when the (Suggested) package is available; small
## samples fall below a guard. Independent-groups tests (Mann-Whitney / Kruskal-Wallis) are used only
## when paired = FALSE. Effect sizes are always computed. Nothing here fits a model - that is left to
## the user's own confirmatory analysis.


##################################################################################################
## Effect sizes ##################################################################################

#' Matched-pairs rank-biserial correlation (paired signed-rank effect size)
#'
#' @description Effect size for a Wilcoxon signed-rank comparison of `y` vs `x` (paired). Positive
#' values mean `y` tends to exceed `x`. Range -1..+1. Zeros (ties) are dropped, matching the test.
#' The sign convention is pinned by tests (a hand-rolled formula is a known footgun).
#' @keywords internal
#' @noRd
.rankBiserialPaired <- function(x, y){
  d <- y - x
  d <- d[!is.na(d) & d != 0]
  n <- length(d)
  if(n == 0) return(NA_real_)
  r <- rank(abs(d))
  2 * sum(r[d > 0]) / (n * (n + 1) / 2) - 1
}

#' Rank-biserial correlation for an independent-groups (Mann-Whitney) comparison of `y` vs `x`
#' @keywords internal
#' @noRd
.rankBiserialUnpaired <- function(x, y){
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if(length(x) == 0 || length(y) == 0) return(NA_real_)
  u <- suppressWarnings(stats::wilcox.test(y, x))$statistic
  2 * unname(u) / (length(x) * length(y)) - 1
}


##################################################################################################
## One-metric comparison test ####################################################################

#' Run the design-appropriate comparison test for one metric
#'
#' @description Given an `individual x level` matrix of values (`NA` = missing observation), runs the
#' appropriate test and computes an effect size. `paired = TRUE` (repeated measures) uses signed-rank
#' (2 levels) / Friedman (>2 levels) on complete blocks - or Skillings-Mack on all data when blocks
#' are incomplete and the `Skillings.Mack` package is installed. `paired = FALSE` uses
#' Mann-Whitney / Kruskal-Wallis on all available data. Returns a one-row summary; never throws
#' (degenerate cases yield `p = NA` with a note).
#' @param mat Numeric matrix, rows = individuals, cols = levels (in order), `NA` for missing.
#' @param levels Character level labels (same order as columns).
#' @param paired Logical; repeated-measures design.
#' @param min.n Minimum complete pairs/blocks below which the test is suppressed.
#' @return A one-row `data.frame` of the test result and effect size.
#' @keywords internal
#' @noRd
.metricTest <- function(mat, levels, paired = TRUE, min.n = 6){

  k <- ncol(mat)
  complete <- rowSums(is.na(mat)) == 0
  n_total <- nrow(mat); n_complete <- sum(complete)

  res <- data.frame(test = NA_character_, paired = paired, n_total = n_total,
                    n_complete = n_complete, n_dropped = n_total - n_complete,
                    statistic = NA_real_, df = NA_real_, p_raw = NA_real_,
                    effect = NA_real_, effect_type = NA_character_,
                    p_method = NA_character_, note = NA_character_, stringsAsFactors = FALSE)

  tc_ <- function(expr) tryCatch(expr, error = function(e) NULL)

  if(paired){
    incomplete <- n_complete < n_total
    use_sm <- incomplete && k > 2 && requireNamespace("Skillings.Mack", quietly = TRUE)

    if(use_sm){
      # Skillings-Mack: repeated-measures test that tolerates incomplete blocks (reduces to Friedman
      # when complete), so partially-observed individuals are used rather than discarded.
      long <- data.frame(y = as.vector(mat),
                         block = factor(rep(seq_len(n_total), k)),
                         grp   = factor(rep(levels, each = n_total), levels = levels))
      sm <- tc_(Skillings.Mack::Ski.Mack(y = long$y, groups = long$grp, blocks = long$block))
      if(!is.null(sm)){
        res$test <- "Skillings-Mack"; res$statistic <- unname(sm$statistic %||% sm$SM.stat)
        res$df <- k - 1; res$p_raw <- unname(sm$p.value %||% sm$p.value.simulated)
        res$p_method <- "asymptotic/simulated"; res$n_complete <- n_total
        res$effect <- if(!is.na(res$statistic)) res$statistic / (n_total * (k - 1)) else NA_real_
        res$effect_type <- "Kendall W (approx)"
        return(res)
      }
    }

    if(n_complete < min.n){ res$note <- sprintf("n=%d complete < min.n=%d", n_complete, min.n); return(res) }
    m <- mat[complete, , drop = FALSE]

    if(k == 2){
      x <- m[, 1]; y <- m[, 2]
      tt <- tc_(suppressWarnings(stats::wilcox.test(y, x, paired = TRUE, exact = n_complete < 50)))
      res$test <- "Wilcoxon signed-rank"
      res$p_method <- if(n_complete < 50) "exact/normal" else "normal"
      if(!is.null(tt)){ res$statistic <- unname(tt$statistic); res$p_raw <- tt$p.value }
      res$effect <- .rankBiserialPaired(x, y); res$effect_type <- "rank-biserial"
    }else{
      long <- data.frame(y = as.vector(m), grp = factor(rep(levels, each = nrow(m)), levels = levels),
                         block = factor(rep(seq_len(nrow(m)), k)))
      tt <- tc_(stats::friedman.test(y ~ grp | block, data = long))
      res$test <- "Friedman"; res$p_method <- "asymptotic"
      if(!is.null(tt)){ res$statistic <- unname(tt$statistic); res$df <- unname(tt$parameter); res$p_raw <- tt$p.value
                        res$effect <- unname(tt$statistic) / (nrow(m) * (k - 1)); res$effect_type <- "Kendall W" }
    }

  }else{
    # independent-groups design (explicit paired = FALSE): use all available observations
    vals <- lapply(seq_len(k), function(j) mat[!is.na(mat[, j]), j])
    ns <- vapply(vals, length, integer(1))
    if(any(ns < 2)){ res$note <- "a level has <2 observations"; return(res) }
    res$n_complete <- sum(ns)
    if(k == 2){
      tt <- tc_(suppressWarnings(stats::wilcox.test(vals[[2]], vals[[1]])))
      res$test <- "Mann-Whitney"; res$p_method <- "exact/normal"
      if(!is.null(tt)){ res$statistic <- unname(tt$statistic); res$p_raw <- tt$p.value }
      res$effect <- .rankBiserialUnpaired(vals[[1]], vals[[2]]); res$effect_type <- "rank-biserial"
    }else{
      long <- data.frame(y = unlist(vals), grp = factor(rep(levels, ns), levels = levels))
      tt <- tc_(stats::kruskal.test(y ~ grp, data = long))
      res$test <- "Kruskal-Wallis"; res$p_method <- "asymptotic"
      if(!is.null(tt)){ H <- unname(tt$statistic); N <- sum(ns)
                        res$statistic <- H; res$df <- unname(tt$parameter); res$p_raw <- tt$p.value
                        res$effect <- max(0, (H - k + 1) / (N - k)); res$effect_type <- "eta2_H" }
    }
  }
  res
}

`%||%` <- function(a, b) if(is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

##################################################################################################
##################################################################################################
##################################################################################################
