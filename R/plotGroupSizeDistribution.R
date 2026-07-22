#######################################################################################################
## Plot the distribution of co-occurring group sizes ##################################################
#######################################################################################################

#' Plot the distribution of co-occurring group sizes
#'
#' @description Computes and plots the frequency distribution of **co-occurring group sizes** - how
#' often 2, 3, 4, ... individuals were detected together (at the same station within a time bin) over
#' the study. It is the aggregation-size view of co-occurrence, complementing the pairwise view
#' (\code{\link{plotAssociations}}) and the per-location view
#' (\code{\link{plotStationStats}} with `type = "co-occurrences"`). The distribution can be split by a
#' time-bin variable (e.g. diel phase) into a stacked bar chart, and computed within/between animal
#' classes via `id.groups` and `group.comparisons`.
#'
#' @details The input is standard detection data (a `mobyData` or a data frame); the wide
#' `time-bin x individual` table needed for co-occurrence is built internally via
#' \code{\link{createWideTable}}. For each time bin, clusters of individuals sharing a station are
#' found and their sizes tallied. When `id.groups` is set, a `"within"` comparison counts clusters of
#' same-class individuals and a `"between"` comparison counts clusters spanning 2+ classes (checked
#' per cluster). Bar height is the overall frequency of each group size (as a percentage of all
#' co-occurrence events in that panel); `split.by` partitions each bar by a time-bin category.
#'
#' The computed distribution is returned invisibly as a tidy data frame.
#'
#' @note The layout and legend are sized in inch units, so the figure stays consistent across
#' datasets and graphics devices.
#'
#' @inheritParams as_moby
#' @param data A `mobyData` object or a data frame of binned detections.
#' @param id.groups Optional named list of ID groups (e.g. species), each drawn in its own panel.
#' @param group.comparisons For `id.groups`: `"within"` (intra-class), `"between"` (inter-class pairs)
#' or `"all"` (both). Defaults to `"all"`.
#' @param split.by Optional name of a time-bin variable (e.g. diel phase, season) by which to split
#' each distribution into a stacked bar. If NULL (default), the distribution is computed over the
#' whole study.
#' @param levels.order Optional integer vector giving the preferred stacking order of the `split.by`
#' levels.
#' @param color.pal Fill colour(s): one per `split.by` level (or a single colour when `split.by` is
#' NULL). If NULL, a colourblind-safe palette is used.
#' @param background.color Panel background colour. Defaults to "grey96".
#' @param annotate Logical; annotate each bar with its total count (`n=`). Defaults to TRUE.
#' @param main Optional overall title.
#' @param legend Logical; draw the `split.by` legend. Defaults to TRUE.
#' @param legend.pos Keyword position for the legend. Defaults to "topright".
#' @param legend.horiz Logical; draw the legend horizontally. Defaults to FALSE.
#' @param cex Global expansion factor for all plot text. Defaults to 1.
#' @param ncol Number of columns in the panel layout. Defaults to 1.
#' @template deviceArgs
#'
#' @return Invisibly, a tidy data frame of the plotted distribution (columns `series`, `group_size`,
#' `split_level`, `count`, `freq`).
#' @examples
#' # Frequency of co-occurring group sizes (one panel per ID group in 'rays')
#' plotGroupSizeDistribution(rays)
#' @export


plotGroupSizeDistribution <- function(data,
                                      id.col = NULL,
                                      timebin.col = NULL,
                                      station.col = NULL,
                                      id.groups = NULL,
                                      split.by = NULL,
                                      group.comparisons = "all",
                                      levels.order = NULL,
                                      color.pal = NULL,
                                      background.color = "grey96",
                                      annotate = TRUE,
                                      main = NULL,
                                      legend = TRUE,
                                      legend.pos = "topright",
                                      legend.horiz = FALSE,
                                      ncol = 1,
                                      cex = 1,
                                      file = NULL,
                                      width = NULL,
                                      height = NULL,
                                      res = 300) {

  ######################################################################################
  # Initial checks #####################################################################
  ######################################################################################

  reviewed <- .validateArguments()
  data <- as.data.frame(reviewed$data)

  errors <- c()
  if(length(group.comparisons) != 1 || !group.comparisons %in% c("within", "between", "all"))
    errors <- c(errors, "'group.comparisons' must be one of 'within', 'between' or 'all'.")
  if(!is.null(split.by) && !split.by %in% colnames(data))
    errors <- c(errors, paste0("The 'split.by' variable '", split.by, "' was not found in the data."))
  if(!is.null(id.groups) && any(duplicated(unlist(id.groups))))
    errors <- c(errors, "Repeated ID(s) in id.groups.")
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse="\n"), call.=FALSE)


  ######################################################################################
  # Build the wide co-occurrence table #################################################
  ######################################################################################

  wide <- suppressWarnings(createWideTable(data, id.col=id.col, timebin.col=timebin.col,
                                           value.col=station.col, verbose=FALSE))
  ids <- attr(wide, "ids")
  timebins <- wide[[timebin.col]]
  det_mat <- as.matrix(wide[, ids, drop=FALSE])           # time-bin x individual, values = station

  # per-time-bin split value (assumed constant within a bin, e.g. diel phase)
  split_vec <- NULL
  if(!is.null(split.by)){
    sm <- data[!duplicated(data[[timebin.col]]), c(timebin.col, split.by)]
    sv <- sm[[split.by]][match(timebins, sm[[timebin.col]])]
    lv <- if(is.factor(data[[split.by]])) levels(data[[split.by]]) else sort(unique(stats::na.omit(sv)))
    split_vec <- factor(as.character(sv), levels=as.character(lv))
  }


  ######################################################################################
  # Compute the group-size distribution ################################################
  ######################################################################################

  dist <- .computeGroupSizeDistribution(det_mat, split_vec, ids, id.groups, group.comparisons)
  if(is.null(dist))
    stop("No co-occurrences were found (no time bin had 2+ individuals sharing a station).", call.=FALSE)

  series_levels <- levels(dist$series)
  split_levels  <- levels(dist$split_level)
  all_sizes     <- sort(unique(dist$group_size))
  n_series      <- length(series_levels)

  # optional user reordering of split levels
  if(!is.null(levels.order)){
    if(is.null(split.by)) warning("'levels.order' is ignored when 'split.by' is not set.", call.=FALSE)
    else split_levels <- split_levels[levels.order]
  }


  ######################################################################################
  # Appearance & layout (device-stable) ################################################
  ######################################################################################

  cex_title <- 1.1*cex; cex_lab <- 1.0*cex; cex_axis <- 0.8*cex; cex_legend <- 0.7*cex; cex_annot <- 0.7*cex

  if(is.null(color.pal)){
    color.pal <- if(is.null(split.by)) "grey45"
                 else if(length(split_levels) <= 8) .okabe_ito_pal(length(split_levels))
                 else grDevices::hcl.colors(length(split_levels), "Dark 3")
  }

  rows <- ceiling(n_series / ncol)
  cw <- par("cin")[1]
  show_leg <- legend && !is.null(split.by)
  # right margin for an in-panel legend; a little extra top for the overall title
  right_in <- if(show_leg) max(nchar(as.character(split_levels))) * cw * cex_legend + 0.45 else 0.15

  # optional file output (before layout / device-metric reads)
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule=list(base=3, slope=2.6, n=ncol, lo=4.5, hi=28),
                     h.rule=list(base=1.6, slope=2.6, n=rows, lo=3.5, hi=28),
                     crowd.unit="panels")
    on.exit(grDevices::dev.off(), add=TRUE, after=FALSE)
  }

  original_par <- .savePar(); on.exit(.restorePar(original_par), add=TRUE)
  par(mfrow=c(rows, ncol),
      mai=c(0.62, 0.72, 0.30, right_in),
      omi=c(0.30, 0.20, if(!is.null(main)) 0.34 else 0.10, 0.10),
      mgp=c(2.4, 0.6, 0), lwd=0.7)


  ######################################################################################
  # Console summary ####################################################################
  ######################################################################################

  .printGroupSizeDistributionSummary(
    n_ids = length(ids), n_bins = length(timebins),
    split.by = split.by, n_split = length(split_levels), id.groups = id.groups,
    group.comparisons = group.comparisons, n_series = n_series,
    size_range = range(all_sizes), n_events = sum(dist$count))


  ######################################################################################
  # Draw panels (one per series) #######################################################
  ######################################################################################

  faceted <- n_series > 1
  for(i in seq_len(n_series)){
    sub <- dist[dist$series == series_levels[i], ]

    # freq / count matrices: rows = split levels (stack segments), cols = group sizes (bars)
    fm <- matrix(0, nrow=length(split_levels), ncol=length(all_sizes),
                 dimnames=list(split_levels, as.character(all_sizes)))
    cm <- fm
    if(nrow(sub) > 0){
      idx <- cbind(match(as.character(sub$split_level), split_levels), match(sub$group_size, all_sizes))
      keep <- !is.na(idx[,1])
      fm[idx[keep,, drop=FALSE]] <- sub$freq[keep]
      cm[idx[keep,, drop=FALSE]] <- sub$count[keep]
    }

    ylim <- c(0, 110)
    bp <- barplot(fm, col=NA, border=NA, axes=FALSE, axisnames=FALSE, ylim=ylim, space=0.4, xpd=FALSE)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.color, border=NA)
    barplot(fm, col=color.pal, border="grey30", axes=FALSE, axisnames=FALSE, ylim=ylim, space=0.4, add=TRUE, xpd=FALSE)

    # y-axis (%) on the leftmost column
    if((i - 1) %% ncol == 0){
      axis(2, at=seq(0, 100, by=20), labels=paste0(seq(0, 100, by=20), "%"), las=1, cex.axis=cex_axis)
      title(ylab="Frequency", cex.lab=cex_lab, line=2.6, xpd=NA)
    }
    # x-axis (group sizes) on the bottom row
    if(i > n_series - ncol){
      axis(1, at=bp, labels=all_sizes, cex.axis=cex_axis, tick=FALSE, line=-0.4)
    }
    # panel title (series) when faceting
    if(faceted){
      g <- series_levels[i]; g <- paste0(toupper(substring(g,1,1)), substring(g,2))
      title(main=g, font.main=2, cex.main=cex_title, line=0.4)
    }
    # count annotation above each bar
    if(annotate){
      totals <- colSums(cm); heights <- colSums(fm)
      lab <- ifelse(totals > 0, paste0("n=", totals), "")
      text(x=bp, y=heights + 4, labels=lab, cex=cex_annot, xpd=FALSE)
    }
    # split legend (once, in the first panel)
    if(show_leg && i == 1){
      .legend(legend.pos, legend=split_levels, fill=color.pal, horiz=legend.horiz,
              box.cex=c(1.4, 1.1), y.intersp=1.3, bty="n", cex=cex_legend)
    }
    box(col="grey40")
  }

  # shared x-axis title and optional overall title
  mtext("Co-occurring group size", side=1, outer=TRUE, line=0.4, cex=cex_lab)
  if(!is.null(main)) mtext(main, side=3, outer=TRUE, font=2, cex=cex_title*1.1, line=0.3)

  invisible(dist)
}


#######################################################################################################
## Distribution computation (internal) ################################################################
#######################################################################################################

#' Compute the co-occurring group-size distribution as a tidy long table
#' @keywords internal
#' @noRd
.computeGroupSizeDistribution <- function(det_mat, split_vec, ids, id.groups, group.comparisons){

  if(is.null(id.groups)){
    series <- list("All" = ids)
    ids_table <- data.frame(ID=ids, group="All", stringsAsFactors=FALSE)
  }else{
    series <- .buildGroupComparisonList(id.groups, group.comparisons)
    ids_table <- .meltList(id.groups); colnames(ids_table) <- c("ID", "group")
  }
  split_levels <- if(is.null(split_vec)) "All" else levels(split_vec)

  rows <- list()
  for(s in names(series)){
    cols <- colnames(det_mat) %in% series[[s]]
    if(!any(cols)) next
    sub <- det_mat[, cols, drop=FALSE]
    grp <- .mapValues(colnames(sub), ids_table$ID, ids_table$group, warn_missing=FALSE)
    comparison <- if(grepl("<->", s, fixed=TRUE)) "between" else "within"
    per_bin <- apply(sub, 1, .countCooccurrences, comparison=comparison, groups=grp, return="sizes")
    if(!is.list(per_bin)) per_bin <- as.list(per_bin)
    for(lev in split_levels){
      keep <- if(is.null(split_vec)) rep(TRUE, length(per_bin)) else (!is.na(split_vec) & split_vec == lev)
      sizes <- unlist(per_bin[keep]); sizes <- sizes[!is.na(sizes) & sizes > 0]
      if(length(sizes) == 0) next
      tab <- table(sizes)
      rows[[length(rows) + 1]] <- data.frame(series=s, group_size=as.integer(names(tab)),
                                             split_level=lev, count=as.numeric(tab), stringsAsFactors=FALSE)
    }
  }
  if(length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  out$freq <- stats::ave(out$count, out$series, FUN=function(x) x / sum(x) * 100)
  out$series      <- factor(out$series, levels=names(series))
  out$split_level <- factor(out$split_level, levels=split_levels)
  rownames(out) <- NULL
  out
}


#######################################################################################################
## Console summary (internal) #########################################################################
#######################################################################################################

#' @keywords internal
#' @noRd
.printGroupSizeDistributionSummary <- function(n_ids, n_bins, split.by, n_split, id.groups,
                                               group.comparisons, n_series, size_range, n_events){
  kv <- .kv
  .summaryOpen("Co-occurring group-size distribution")
  kv("Individuals", format(n_ids, big.mark=","))
  kv("Time bins", format(n_bins, big.mark=","))
  if(!is.null(id.groups)) kv("Groups", sprintf("%d (%s -> %d series)", length(id.groups), group.comparisons, n_series))
  if(!is.null(split.by)) kv("Split by", sprintf("%s (%d levels)", split.by, n_split))
  kv("Group sizes", sprintf("%d to %d", size_range[1], size_range[2]))
  kv("Co-occurrences", format(n_events, big.mark=","))
  .summaryClose()
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
