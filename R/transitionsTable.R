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
#'
#' @return A data frame with one row per directed transition (and group-label rows when the
#' network was built with `id.groups`).
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

transitionsTable <- function(network, id.metadata = NULL, error.stat = "se") {

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

  # overall duration unit (days if mean transit time is long)
  all_dur <- edges$mean_duration_h
  use_days <- length(all_dur) > 0 && mean(all_dur, na.rm = TRUE) > 72
  dur_label <- if (use_days) "Mean duration (d)" else "Mean duration (h)"

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
      tab <- data.frame(Type = character(0), Movements = integer(0), Individuals = character(0),
                        check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      type <- paste(e_g$from, "-->", e_g$to)
      pct <- if (!is.na(gsize) && gsize > 0) round(e_g$n_individuals / gsize * 100) else NA
      individuals <- if (all(is.na(pct))) as.character(e_g$n_individuals) else
        paste0(e_g$n_individuals, " (", pct, "%)")

      # per-transition duration (mean +/- error), from the records
      dur_mean <- dur_err <- rep(NA_real_, nrow(e_g))
      for (k in seq_len(nrow(e_g))) {
        d <- rec_g$duration_h[rec_g$from == e_g$from[k] & rec_g$to == e_g$to[k]]
        d <- d[is.finite(d)]
        if (length(d) > 0) { dur_mean[k] <- mean(d); dur_err[k] <- if (length(d) > 1) errFun(d) else NA }
      }
      if (use_days) { dur_mean <- dur_mean / 24; dur_err <- dur_err / 24 }
      duration <- paste0(sprintf("%.1f", dur_mean),
                         ifelse(is.na(dur_err), "", paste0(" \u00b1 ", sprintf("%.1f", dur_err))))

      tab <- data.frame(Type = type, Movements = e_g$n_movements, Individuals = individuals,
                        check.names = FALSE, stringsAsFactors = FALSE)
      tab[[dur_label]] <- duration

      # per-transition metadata summaries
      if (!is.null(id.metadata)) {
        tx_ids <- lapply(seq_len(nrow(e_g)), function(k)
          unique(rec_g$id[rec_g$from == e_g$from[k] & rec_g$to == e_g$to[k]]))
        for (nc in numeric_cols) {
          vals <- vapply(tx_ids, function(ids) {
            v <- id.metadata[[nc]][as.character(id.metadata[[id.col]]) %in% ids]
            v <- v[!is.na(v)]
            if (length(v) == 0) return(NA_character_)
            m <- mean(v); e <- if (length(v) > 1) errFun(v) else NA
            digits <- max(.decimalPlaces(id.metadata[[nc]]), na.rm = TRUE) + 1
            paste0(sprintf(paste0("%.", digits, "f"), m),
                   ifelse(is.na(e), "", paste0(" \u00b1 ", sprintf(paste0("%.", digits, "f"), e))))
          }, character(1))
          tab[[paste("Mean", tools::toTitleCase(nc))]] <- vals
        }
        for (cc in character_cols) {
          vals <- vapply(tx_ids, function(ids) {
            v <- id.metadata[[cc]][as.character(id.metadata[[id.col]]) %in% ids]
            v <- v[!is.na(v)]
            if (length(v) == 0) return(NA_character_)
            tb <- table(v)
            paste(paste0(as.integer(tb), " ", names(tb)), collapse = " | ")
          }, character(1))
          tab[[tools::toTitleCase(cc)]] <- vals
        }
      }

      # order by site sequence
      ord <- order(factor(e_g$from, levels = ordered_sites), factor(e_g$to, levels = ordered_sites))
      tab <- tab[ord, , drop = FALSE]
    }

    # group title row when several groups
    if (length(group_levels) > 1) {
      title_row <- tab[0, , drop = FALSE]
      title_row[1, ] <- ""
      title_row$Type <- g
      tab <- rbind(title_row, tab)
    }
    tables[[g]] <- tab
  }

  out <- do.call(.rbindFill, tables)
  out[is.na(out)] <- "-"
  rownames(out) <- NULL
  out
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
