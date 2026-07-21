#######################################################################################################
## Plot association significance / overlap matrix ####################################################
#######################################################################################################

#' Plot an association significance / overlap matrix
#'
#' @description Draws an id x id matrix of pairwise association results from a
#' \code{\link{randomizeAssociations}} null model - coloured either by permutation significance
#' (positive / negative / non-significant) or by mean overlap percentage. Individuals can be sorted by
#' a supplied metric (e.g. body size).
#'
#' @param random.results The null-model output of \code{\link{randomizeAssociations}}.
#' @param type One of `"significance"` (default; +/-/ns from the FDR-adjusted permutation test) or
#' `"mean overlap"` (the mean pairwise overlap percentage).
#' @param color.pal Colour palette. For `"significance"`, a length-3 vector for
#' positive / negative / non-significant; for `"mean overlap"`, a continuous gradient. If NULL,
#' colourblind-safe defaults are used.
#' @param sort.by Optional numeric/integer/factor/POSIXct vector (length = number of IDs) used to
#' reorder rows and columns (e.g. by size). Values are shown in a header row/column: numbers to a
#' common number of decimal places, factors as their labels (ordered by level) and dates formatted
#' (ordered chronologically).
#' @param sort.by.title Header label for the `sort.by` metric. If omitted, it is inferred from a
#' `data$column` expression.
#' @param full.scale Logical; for `"mean overlap"`, map the colour ramp over the fixed 0-100% domain
#' instead of the observed range. Defaults to FALSE.
#' @param discard.missing Logical; drop rows/columns with no valid pairwise interactions. Defaults to TRUE.
#' @param main Optional overall title.
#' @param cex Global expansion factor for all text (cell labels and legend). Defaults to 1.
#' @template deviceArgs
#'
#' @return Invisibly, the id x id matrix of plotted values (significance labels or overlap percentages).
#' @seealso \code{\link{calculateAssociations}}, \code{\link{randomizeAssociations}}, \code{\link{plotAssociations}}
#' @examples
#' \donttest{
#' # Significance / overlap matrix from a permutation null model
#' wide  <- createWideTable(rays, value.col = "station")
#' assoc <- calculateAssociations(wide)
#' rand  <- randomizeAssociations(wide, assoc, iterations = 100, random.seed = 1)
#' plotAssociationMatrix(rand)
#' plotAssociationMatrix(rand, type = "mean overlap")
#' }
#' @export


plotAssociationMatrix <- function(random.results,
                                  type = c("significance", "mean overlap"),
                                  color.pal = NULL,
                                  sort.by = NULL,
                                  sort.by.title = NULL,
                                  full.scale = FALSE,
                                  discard.missing = TRUE,
                                  main = NULL,
                                  cex = 1,
                                  file = NULL,
                                  width = NULL,
                                  height = NULL,
                                  res = 300) {

  ##############################################################################
  # Checks #####################################################################
  ##############################################################################

  type <- match.arg(type)
  errors <- .validateRandomResults(random.results)
  if(!is.null(sort.by)){
    if(length(sort.by) != length(attributes(random.results)$ids))
      errors <- c(errors, "'sort.by' must have the same length as the number of IDs in random.results.")
    else if(!inherits(sort.by, c("numeric", "integer", "factor", "POSIXct")))
      errors <- c(errors, "'sort.by' must contain sortable values (numeric, integer, factor or POSIXct).")
  }
  if(!is.null(color.pal) && type == "significance" && length(color.pal) != 3)
    errors <- c(errors, "'color.pal' must contain 3 colours when type = 'significance'.")
  if(!is.null(color.pal) && type == "mean overlap" && length(color.pal) < 10)
    warning("- Consider supplying more colours in 'color.pal' for a smoother gradient.", call. = FALSE)
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)

  # infer sort.by.title from a `data$column` expression when not supplied
  if(is.null(sort.by.title) && !is.null(sort.by)){
    e <- substitute(sort.by)
    if(is.call(e) && as.character(e[[1]]) == "$") sort.by.title <- as.character(e[[3]])
    else stop("Please provide a 'sort.by.title'.", call. = FALSE)
  }

  cex_table <- 0.8 * cex; cex_legend <- 0.7 * cex

  ##############################################################################
  # Prepare data ###############################################################
  ##############################################################################

  pairwise_results <- random.results$pairwise_results
  id.groups <- attributes(random.results)$id.groups
  ids <- as.character(attributes(random.results)$ids)
  group_comparisons <- attributes(random.results)$group.comparisons

  # default palettes: colourblind-safe triple for significance, viridis for overlap
  if(is.null(color.pal)){
    if(type == "significance") color.pal <- adjustcolor(c("#0072B2", "#D55E00", "grey85"), alpha.f = 0.35)
    else color.pal <- c(adjustcolor(.viridis_pal(100)[1], alpha.f = 0.65), adjustcolor(.viridis_pal(100)[2:100], alpha.f = 0.7))
  }

  # filter to within-/between-group comparisons when requested
  if(!is.null(id.groups)){
    id_lookup <- .meltList(id.groups); colnames(id_lookup) <- c("id", "group")
    grp <- function(x) as.character(.mapValues(x, from = id_lookup$id, to = as.character(id_lookup$group), warn_missing = FALSE))
    g1 <- grp(pairwise_results$id1); g2 <- grp(pairwise_results$id2)
    if(identical(group_comparisons, "within"))  pairwise_results <- pairwise_results[g1 == g2, , drop = FALSE]
    if(identical(group_comparisons, "between")) pairwise_results <- pairwise_results[g1 != g2, , drop = FALSE]
  }

  # id x id matrix (shared reshaper; drop=FALSE keeps labels even with a single surviving row/col)
  if(type == "mean overlap"){
    mat <- .associationMatrix(pairwise_results, "association", discard.missing)
    overlap_range <- range(pairwise_results$association, na.rm = TRUE)     # (was range() on the whole df -> crash)
  }else{
    mat <- .associationMatrix(pairwise_results, "significance", discard.missing)
    mat[mat == "positive"] <- "+"; mat[mat == "negative"] <- "-"; mat[mat == "non-significant"] <- "ns"
  }
  result_matrix <- mat

  ids1 <- rownames(mat); ids2 <- colnames(mat)
  n_ids <- length(unique(c(ids1, ids2)))

  # assemble the display table: IDs as the first row and column
  contingency_table <- cbind(ids1, mat)
  contingency_table <- rbind(c(NA, ids2), contingency_table)
  rownames(contingency_table) <- NULL; colnames(contingency_table) <- NULL
  contingency_table[1, 1] <- "IDs"

  # optional sort of rows/columns by a supplied metric
  if(!is.null(sort.by)){
    # subset with [ ] to keep sort.by's class: factors and POSIXct must reach .sortByLabels() intact
    # (unlist() drops them to codes / epoch seconds)
    vars1 <- sort.by[match(ids1, ids)]
    vars2 <- sort.by[match(ids2, ids)]
    # sorting below uses these underlying values; only the printed labels are class-specific
    labels <- .sortByLabels(vars1, vars2)
    contingency_table <- cbind(contingency_table, c(NA, labels$labs1))
    contingency_table <- rbind(contingency_table, c(NA, labels$labs2, NA))
    cols <- ncol(contingency_table); rows <- nrow(contingency_table)
    contingency_table <- contingency_table[, c(1, cols, 2:(cols - 1))]
    contingency_table <- contingency_table[c(1, rows, 2:(rows - 1)), ]
    contingency_table[-c(1:2), ] <- contingency_table[-c(1:2), ][order(vars1), ]
    contingency_table[, -c(1:2)] <- contingency_table[, -c(1:2)][, order(vars2)]
    contingency_table[2, 2] <- sort.by.title
    rownames(contingency_table) <- NULL; colnames(contingency_table) <- NULL
  }

  n_cols <- ncol(contingency_table); n_rows <- nrow(contingency_table)

  ##############################################################################
  # Device + summary ###########################################################
  ##############################################################################

  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 5, slope = 0.35, n = n_cols, lo = 4, hi = 36),
                     h.rule = list(base = 3, slope = 0.30, n = n_rows, lo = 4, hi = 36),
                     crowd.unit = "individuals")
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  .printAssociationMatrixSummary(n_ids = n_ids, type = type, mat = result_matrix,
                                 metric = attributes(random.results)$metric)

  ##############################################################################
  # Cell colours ###############################################################
  ##############################################################################

  values <- expand.grid(row = seq_len(n_rows), col = seq_len(n_cols))
  values$value <- unlist(apply(values, 1, function(x) contingency_table[x[1], x[2]]))
  values$font <- 1; values$font[values$row == 1 | values$col == 1] <- 2
  values$color <- NA

  if(type == "significance"){
    values$color[values$value == "+"]  <- color.pal[1]
    values$color[values$value == "-"]  <- color.pal[2]
    values$color[values$value == "ns"] <- color.pal[3]
    values$color[is.na(values$value) & values$row > 1 & values$col > 1] <- color.pal[3]
  }else{
    if(!is.null(sort.by)) overlap_cells <- which(values$row > 2 & values$col > 2)
    else overlap_cells <- which(values$row > 1 & values$col > 1)
    overlaps <- as.numeric(values$value[overlap_cells])
    if(full.scale) idx <- round(.rescale(overlaps, from = c(0, 100), to = c(1, length(color.pal) - 1)))
    else idx <- round(.rescale(overlaps, from = range(overlaps, na.rm = TRUE), to = c(1, length(color.pal) - 1)))
    idx[is.na(idx)] <- 0
    values$color[overlap_cells] <- color.pal[idx + 1]
    values$value[overlap_cells] <- sprintf("%.1f", as.numeric(values$value[overlap_cells]))
    values$value[values$value == "NA"] <- ""
  }

  ##############################################################################
  # Draw #######################################################################
  ##############################################################################

  par(mar = c(1, 1, if(is.null(main)) 1 else 3, 12))
  plot(0, 0, xlim = c(0.5, n_cols + 0.5), ylim = c(n_rows + 0.5, 0.5), xaxs = "i", yaxs = "i",
       type = "n", axes = FALSE, main = "", ylab = "", xlab = "")
  rect(values$col - 0.5, values$row - 0.5, values$col + 0.5, values$row + 0.5, col = values$color, border = NA)
  text(values$col, values$row, labels = values$value, font = values$font, cex = cex_table, xpd = TRUE, adj = 0.5)

  # separator gridlines + bounding box
  sep <- if(!is.null(sort.by)) 2.5 else 1.5
  segments(x0 = sep, x1 = n_cols + 0.5, y0 = sep); segments(y0 = sep, y1 = n_rows + 0.5, x0 = sep)
  segments(x0 = 0.5, x1 = n_cols + 0.5, y0 = 0.5, xpd = TRUE)
  segments(x0 = 0.5, x1 = n_cols + 0.5, y0 = n_rows + 0.5, xpd = TRUE)
  segments(y0 = 0.5, y1 = n_rows + 0.5, x0 = 0.5, xpd = TRUE)
  segments(y0 = 0.5, y1 = n_rows + 0.5, x0 = n_cols + 0.5, xpd = TRUE)

  # legend
  if(type == "significance"){
    labs <- c("overlap > random\n(joint resource utilization)", "overlap < random\n(spatiotemporal segregation)", "non-significant overlap")
    .legend(x = par("usr")[2] + diff(par("usr")[1:2]) * 0.02, y = mean(par("usr")[3:4]),
            legend = labs, box.lwd = 0.1, fill = color.pal, box.cex = c(1.5, 1), bty = "n",
            cex = cex_legend, y.intersp = 1.6, xpd = TRUE)
  }else{
    if(full.scale){ overlap_range <- c(0, 100); overlap_labs <- pretty(0:100, min.n = 4) }
    else { overlap_labs <- pretty(overlap_range, min.n = 4); overlap_labs <- overlap_labs[overlap_labs >= min(overlap_range) & overlap_labs <= max(overlap_range)] }
    .colorlegend(col = color.pal, zlim = overlap_range, zval = overlap_labs, posx = c(0.88, 0.895),
                 posy = c(0.2, 0.65), main = "Overlap (%)", main.cex = cex_legend + 0.1,
                 digit = max(0, max(.decimalPlaces(overlap_labs), na.rm = TRUE)), main.adj = 0, cex = cex_legend)
  }
  if(!is.null(main)) mtext(main, side = 3, outer = FALSE, font = 2, cex = cex_table + 0.3, line = 1)

  invisible(result_matrix)
}


#######################################################################################################
## Internal helpers ###################################################################################

# Render the 'sort.by' values shown in the header row/column. Numeric metrics share a common number of
# decimal places (so the column reads as a scale); factors and dates print as their labels instead of
# their underlying codes/epoch seconds.
#' @keywords internal
#' @noRd
.sortByLabels <- function(vars1, vars2){
  if(is.factor(vars1)){
    return(list(labs1 = as.character(vars1), labs2 = as.character(vars2)))
  }
  if(inherits(vars1, "POSIXct")){
    return(list(labs1 = format(vars1), labs2 = format(vars2)))
  }
  decimals <- .decimalPlaces(c(vars1, vars2))
  digits <- if(all(is.na(decimals))) 0 else max(0, max(decimals, na.rm = TRUE))
  fmt <- paste0("%.", digits, "f")
  list(labs1 = sprintf(fmt, vars1), labs2 = sprintf(fmt, vars2))
}

#' @keywords internal
#' @noRd
.printAssociationMatrixSummary <- function(n_ids, type, mat, metric){
  kv <- .kv
  .summaryOpen("Association matrix")
  kv("Individuals", n_ids)
  kv("Metric", if(is.null(metric)) "association index" else metric)
  kv("Display", type)
  if(type == "significance"){
    v <- unlist(mat)
    kv("Significant", sprintf("%d positive, %d negative, %d ns",
                              sum(v == "+", na.rm = TRUE), sum(v == "-", na.rm = TRUE), sum(v == "ns", na.rm = TRUE)))
  }
  .summaryClose()
}

#######################################################################################################
#######################################################################################################
