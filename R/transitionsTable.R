#######################################################################################################
## Movement transitions summary table #################################################################
#######################################################################################################

#' Summarise a movement network as a transitions table
#'
#' @description Produces a publication-ready summary table from a movement network (the output of
#' \code{\link{calculateTransitions}}). Each row is a directed transition between two locations,
#' with the number of movements, the number (and percentage) of distinct individuals performing
#' it, and the mean transit duration. When per-animal metadata is supplied, numeric variables are
#' summarised as mean +/- error and categorical variables as level counts, per transition type.
#'
#' This is the formatting counterpart to \code{\link{calculateTransitions}} (which holds the
#' numeric network), mirroring the \code{\link{calculateResidency}} / \code{\link{summaryTable}}
#' split. It is purely a table: network visualisation is handled by \code{plot()} on the network
#' object, and temporal distributions of transition timing are available from the network's
#' `transition_records` attribute.
#'
#' @param network A `mobyNetwork` object of type `"movement"`, from \code{\link{calculateTransitions}}.
#' @param id.metadata Optional data frame of per-animal metadata. Must contain an animal-ID column
#' matching the network's `id.col`. Numeric columns are summarised as mean +/- error and
#' categorical columns as counts, per transition type.
#' @param error.stat Error statistic for numeric metadata summaries: `"se"` (standard error,
#' default) or `"sd"` (standard deviation).
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A \code{mobyTable}: a TYPED data frame (counts stay integer, indices numeric, dates
#' POSIXct), one row per directed transition, so the result can be computed on directly. Presentation - fixed
#' precision, the display-only `mean +/- error` row, group headings - is applied by
#' \code{\link[=format.mobyTable]{format}} and \code{\link[=print.mobyTable]{print}}. Export the
#' rendered version with \code{write.csv(format(x), file, row.names = FALSE)}.
#'
#' Columns use stable snake_case names, which is what `format(decimals=)` and `format(group.by=)`
#' are keyed on; the publication headers live in `format(style = "report")`.
#' \itemize{
#'   \item `transition` (`"A --> B"`), and a `group` factor when the network was built with `id.groups`
#'   \item `n_movements`, `n_individuals`, `pct_individuals`
#'   \item `mean_duration`, `error_duration` (hours, or days when transits are long - the unit is
#'     recorded on the table and named in the header)
#'   \item `mean_<var>` for each numeric `id.metadata` column, and the raw column for each
#'     categorical one
#' }
#'
#' The count and the share of individuals are separate columns, as are the mean duration and its
#' error: written as one string (`"3 (75%)"`, `"30.7 +/- 12.6"`) neither could be computed on.
#'
#' @seealso \code{\link{calculateTransitions}}, \code{\link{summaryTable}}
#'
#' @examples
#' data(rays)
#' trans <- calculateTransitions(rays, spatial.col = "station")
#' # publication-ready summary of the directed transitions
#' transitionsTable(trans)
#'
#' @export

transitionsTable <- function(network, id.metadata = NULL, error.stat = "se",
                             verbose = getOption("moby.verbose", TRUE)) {

  if (!inherits(network, "mobyNetwork") || !identical(attr(network, "network.type"), "movement")) {
    stop("'network' must be a movement 'mobyNetwork' object (see calculateTransitions()).", call. = FALSE)
  }
  error.stat <- match.arg(tolower(error.stat), c("se", "sd"))
  errFun <- if (error.stat == "se") function(x) .stdError(x) else function(x) stats::sd(x, na.rm = TRUE)

  edges <- networkEdges(network)
  records <- attr(network, "transition_records")
  group_sizes <- attr(network, "group.sizes")
  id.groups <- attr(network, "id.groups")
  ordered_sites <- attr(network, "ordered.sites")
  id.col <- if (!is.null(id.metadata)) {
    # the id column name used by the network (from its construction metadata)
    nm <- attr(network, "id.col")
    if (is.null(nm)) "ID" else nm
  } else NULL

  group_levels <- if (is.null(id.groups)) "all" else names(id.groups)

  # ---- header ---------------------------------------------------------------------------------
  # Placed after the argument checks, so a call that is going to error never prints a banner first.
  # This function prints a header and nothing else: the returned table IS the summary, so a closing
  # line would only restate what the user is about to read. The one criterion is the statistic behind
  # every "+/-" the table reports (worded exactly as in summaryTable(), so the two never disagree);
  # the duration unit is derived from the data and already spelled out in the column name.
  # Count the sites the network actually HAS. attr(, "ordered.sites") is the factor levels of the
  # spatial column, which includes never-visited stations, so using it here reported more sites than
  # calculateTransitions() and networkMetrics() do for the same object.
  n_sites <- length(unique(as.character(networkNodes(network)$site)))
  .mobyHeader("transitionsTable()", "Summarising directed transitions between sites",
              input = paste0(.fmtCount(n_sites, "site"), " ", .mobyGlyph("mid"), " ",
                             .fmtCount(nrow(edges), "transition")),
              criteria = c("error" = if (error.stat == "se") "standard error (se)"
                                     else "standard deviation (sd)"),
              verbose = verbose)

  # overall duration unit (days if mean transit time is long). A display choice, so it travels as
  # table metadata and the column keeps one stable name whatever the data does.
  all_dur <- edges$mean_duration_h
  use_days <- length(all_dur) > 0 && mean(all_dur, na.rm = TRUE) > 72
  dur_units <- if (use_days) "d" else "h"

  # metadata column types
  if (!is.null(id.metadata)) {
    if (!id.col %in% colnames(id.metadata)) {
      # fall back to a sensible guess
      id.col <- colnames(id.metadata)[1]
    }
    meta_cols <- setdiff(colnames(id.metadata), id.col)
    numeric_cols <- meta_cols[vapply(id.metadata[meta_cols], function(x) is.numeric(x), logical(1))]
    character_cols <- setdiff(meta_cols, numeric_cols)
  }

  tables <- list()
  for (g in group_levels) {
    e_g <- edges[edges$group == g, , drop = FALSE]
    rec_g <- records[[g]]
    gsize <- if (!is.null(group_sizes)) group_sizes[[g]] else NA

    if (nrow(e_g) == 0) {
      tab <- data.frame(transition = character(0), n_movements = integer(0),
                        n_individuals = integer(0), pct_individuals = numeric(0),
                        mean_duration = numeric(0), error_duration = numeric(0),
                        check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      type <- paste(e_g$from, "-->", e_g$to)
      # count and share are two facts, so they are two columns. As one string ("3 (75%)") neither
      # could be computed on, and the share was unrecoverable once the group size was gone.
      pct <- if (!is.na(gsize) && gsize > 0) e_g$n_individuals / gsize * 100 else NA_real_

      # per-transition duration (mean +/- error), from the records
      dur_mean <- dur_err <- rep(NA_real_, nrow(e_g))
      for (k in seq_len(nrow(e_g))) {
        d <- rec_g$duration_h[rec_g$from == e_g$from[k] & rec_g$to == e_g$to[k]]
        d <- d[is.finite(d)]
        if (length(d) > 0) { dur_mean[k] <- mean(d); dur_err[k] <- if (length(d) > 1) errFun(d) else NA }
      }
      if (use_days) { dur_mean <- dur_mean / 24; dur_err <- dur_err / 24 }

      tab <- data.frame(transition = type, n_movements = e_g$n_movements,
                        n_individuals = e_g$n_individuals,
                        pct_individuals = as.numeric(rep(pct, length.out = nrow(e_g))),
                        mean_duration = dur_mean, error_duration = dur_err,
                        check.names = FALSE, stringsAsFactors = FALSE)

      # per-transition metadata summaries
      if (!is.null(id.metadata)) {
        tx_ids <- lapply(seq_len(nrow(e_g)), function(k)
          unique(rec_g$id[rec_g$from == e_g$from[k] & rec_g$to == e_g$to[k]]))
        for (nc in numeric_cols) {
          vals <- vapply(tx_ids, function(ids) {
            v <- id.metadata[[nc]][as.character(id.metadata[[id.col]]) %in% ids]
            v <- v[!is.na(v)]
            if (length(v) == 0) return(NA_real_)
            mean(v)
          }, numeric(1))
          tab[[paste0("mean_", nc)]] <- vals
        }
        for (cc in character_cols) {
          vals <- vapply(tx_ids, function(ids) {
            v <- id.metadata[[cc]][as.character(id.metadata[[id.col]]) %in% ids]
            v <- v[!is.na(v)]
            if (length(v) == 0) return(NA_character_)
            tb <- table(v)
            paste(paste0(as.integer(tb), " ", names(tb)), collapse = " | ")
          }, character(1))
          tab[[cc]] <- vals
        }
      }

      # order by site sequence
      ord <- order(factor(e_g$from, levels = ordered_sites), factor(e_g$to, levels = ordered_sites))
      tab <- tab[ord, , drop = FALSE]
    }

    if (length(group_levels) > 1 && nrow(tab) > 0) tab$group <- g
    tables[[g]] <- tab
  }

  out <- do.call(.rbindFill, tables)
  if (length(group_levels) > 1 && "group" %in% names(out))
    out$group <- factor(out$group, levels = group_levels)
  rownames(out) <- NULL
  .newMobyTable(out, kind = "transitions", error.stat = error.stat, label.col = "transition",
                units = list(duration = dur_units),
                extra = list(id.groups = id.groups, processing.date = Sys.time()))
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
