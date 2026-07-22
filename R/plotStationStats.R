#######################################################################################################
## Plot station (receiver) statistics #################################################################
#######################################################################################################

#' Plot receiver-based statistics
#'
#' @description Computes and plots summary statistics across receivers (or any coarser spatial unit
#' via `aggregate.by`) from passive acoustic telemetry data, as bar charts. Available statistics are
#' the number of `"detections"`, the mean within-individual detection share (`"average detections"`),
#' the number of unique `"individuals"`, and the number of `"co-occurrences"` (time-bins where two or
#' more animals share a location). Data can be split by animal class through `id.groups`.
#'
#' @details Each requested `type` is drawn in its own row of panels with its own axis and unit, so
#' quantities with different denominators are never conflated on a shared scale. When `id.groups` is
#' supplied, each group (and, for co-occurrences, each between-group pair when `group.comparisons` is
#' `"between"`/`"all"`) becomes a column, sharing the per-type y-scale so groups are directly
#' comparable. Bar height encodes counts for count-based statistics and a proportion for
#' `"average detections"` by default (`value.scale` overrides this); the complementary quantity is shown as
#' a bar annotation.
#'
#' The computed statistics are returned invisibly as a tidy data frame, so the plotted numbers are
#' available programmatically.
#'
#' @note The layout (margins, the optional numbered-location key) is sized in inch units, so it stays
#' consistent across datasets and graphics devices.
#'
#' @inheritParams as_moby
#' @param data A data frame containing animal detections with corresponding time-bins.
#' @param type One or more statistics to plot (each drawn in its own panel row): `"detections"`,
#' `"average detections"`, `"individuals"`, `"co-occurrences"`. Defaults to `"detections"`.
#' @param value.scale What the bar height encodes: `"natural"` (default; counts for count-based statistics,
#' a proportion for `"average detections"`), `"count"`, or `"proportion"`.
#' @param aggregate.by Optional column name to summarise by (e.g. habitat or region); defaults to
#' `station.col`. Each station must map to a single `aggregate.by` value.
#' @param id.groups Optional named list of ID groups, used to split animals into classes (e.g.
#' species); each becomes a panel column.
#' @param group.comparisons For co-occurrences with `id.groups`: `"within"` (intra-group), `"between"`
#' (inter-group pairs) or `"all"` (both). Defaults to `"all"`.
#' @param color.pal Fill colours, one per `type`. If NULL, a colourblind-safe palette is used.
#' @param background.color Panel background colour. Defaults to "grey96".
#' @param station.labels How to label locations on the x-axis: `"names"` (default), `"rotated"`
#' (vertical names) or `"numbered"` (numbers, with a name key drawn in the right margin).
#' @param annotate Bar annotation: `"auto"` (default; the quantity the axis does not show),
#' `"count"`, `"proportion"`, `"both"` or `"none"`.
#' @param annot.min Optional numeric; hide annotations for bars below this count (reduces clutter).
#' @param main Optional overall title.
#' @param xlab Optional x-axis title (defaults to the title-cased `aggregate.by` name).
#' @param legend Logical; when `station.labels = "numbered"`, draw the number-to-name key. Defaults
#' to TRUE.
#' @param cex Global expansion factor for all plot text. Defaults to 1.
#' @template deviceArgs
#' @param ... Further arguments passed to \code{\link[graphics]{barplot}}.
#'
#' @return Invisibly, a tidy data frame of the plotted statistics (columns `series`, `type`,
#' `location`, `count`, `proportion`).
#' @examples
#' # number of detections and of unique individuals per receiver
#' plotStationStats(rays, type = c("detections", "individuals"))
#' @export


plotStationStats <- function(data,
                             id.col = NULL,
                             timebin.col = NULL,
                             station.col = NULL,
                             id.groups = NULL,
                             group.comparisons = "all",
                             aggregate.by = NULL,
                             type = "detections",
                             value.scale = c("natural", "count", "proportion"),
                             color.pal = NULL,
                             background.color = "grey96",
                             station.labels = c("names", "rotated", "numbered"),
                             annotate = c("auto", "count", "proportion", "both", "none"),
                             annot.min = NULL,
                             main = NULL,
                             xlab = NULL,
                             legend = TRUE,
                             cex = 1,
                             file = NULL,
                             width = NULL,
                             height = NULL,
                             res = 300,
                             ...) {

  ######################################################################################
  # Initial checks #####################################################################
  ######################################################################################

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  value.scale    <- match.arg(value.scale)
  station.labels <- match.arg(station.labels)
  annotate       <- match.arg(annotate)

  valid_types <- c("detections", "average detections", "individuals", "co-occurrences")
  errors <- c()
  if(any(!type %in% valid_types)) errors <- c(errors, paste0("Invalid 'type'. Choose one or more of: ", paste(shQuote(valid_types), collapse=", "), "."))
  if(length(group.comparisons) != 1 || !group.comparisons %in% c("within", "between", "all"))
    errors <- c(errors, "'group.comparisons' must be one of 'within', 'between' or 'all'.")
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse="\n"), call.=FALSE)

  original_par <- .savePar()
  on.exit(.restorePar(original_par))

  # coerce station column to factor
  if(!is.factor(data[, station.col])){
    warning("- Converting 'station.col' to factor.", call.=FALSE)
    data[, station.col] <- as.factor(data[, station.col])
  }

  # resolve the aggregation column (defaults to the station column)
  if(is.null(aggregate.by)){
    aggregate.by <- station.col
  }else{
    counts <- table(data[[station.col]], data[[aggregate.by]])
    multi  <- rownames(counts)[rowSums(counts > 0) > 1]
    if(length(multi) > 0)
      stop("Each station must map to a single 'aggregate.by' group. Offending stations: ",
           paste(multi, collapse=", "), call.=FALSE)
  }
  if(!is.factor(data[, aggregate.by])){
    warning("- Converting 'aggregate.by' to factor.", call.=FALSE)
    data[, aggregate.by] <- as.factor(data[, aggregate.by])
  }


  ######################################################################################
  # Compute statistics (tidy long table) ###############################################
  ######################################################################################

  stats <- .computeStationStats(data, type, id.groups, group.comparisons,
                                id.col, timebin.col, station.col, aggregate.by)

  type          <- levels(stats$type)          # requested order, validated
  series_levels <- levels(stats$series)
  n_types       <- length(type)
  n_series      <- length(series_levels)
  loc_levels    <- levels(stats$location)

  # optional file output: size the device to the panel grid (rows = types, columns = series)
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule=list(base=3,   slope=3.2, n=n_series, lo=5, hi=30),
                     h.rule=list(base=1.2, slope=2.7, n=n_types,  lo=4, hi=30),
                     crowd.unit="panels")
    on.exit(grDevices::dev.off(), add=TRUE, after=FALSE)
  }


  ######################################################################################
  # Layout & appearance ################################################################
  ######################################################################################

  cex_title <- 1.1 * cex; cex_lab <- 1.0 * cex; cex_axis <- 0.8 * cex; cex_annot <- 0.7 * cex

  if(is.null(color.pal)) color.pal <- if(n_types <= 8) .okabe_ito_pal(n_types) else grDevices::hcl.colors(n_types, "Dark 3")

  disp_names <- c("detections"="Detections", "average detections"="Mean detection share",
                  "individuals"="Individuals", "co-occurrences"="Co-occurrences")

  # per-type scale ("average detections" is intrinsically a proportion)
  type_scale <- function(tp) if(tp == "average detections") "proportion" else if(value.scale == "natural") "count" else value.scale
  # per-type y-max (shared across the row's series), with headroom for annotations
  ymax_by_type <- vapply(type, function(tp){
    col <- if(type_scale(tp) == "count") "count" else "proportion"
    m <- suppressWarnings(max(stats[stats$type == tp, col], na.rm=TRUE))
    if(!is.finite(m) || m <= 0) 1 else m
  }, numeric(1))

  # x-axis labels
  labs_x <- if(station.labels == "numbered") as.character(seq_along(loc_levels)) else loc_levels
  cw <- par("cin")[1]

  # device-stable margins (inches). Panels keep compact, UNIFORM margins; the shared decorations
  # (column headers, x-labels, axis titles, key) live in the outer margins and are drawn with
  # xpd=NA, so all panels stay the same height regardless of label length.
  csi <- par("csi")
  xlab_in    <- if(station.labels == "rotated") max(nchar(loc_levels)) * cw * cex_axis + 0.10 else 0.24
  show_key   <- legend && station.labels == "numbered"
  key_labs   <- paste0(seq_along(loc_levels), ". ", loc_levels)
  omi_bottom <- xlab_in + 0.34                                   # x-labels + x-title
  omi_top    <- (if(n_series > 1) 0.28 else 0.06) + (if(!is.null(main)) 0.34 else 0)
  omi_right  <- if(show_key) max(nchar(key_labs)) * cw * cex_axis + 0.25 else 0.12
  par(mfrow = c(n_types, n_series),
      mai = c(0.14, 0.85, 0.14, 0.12),
      omi = c(omi_bottom, 0.20, omi_top, omi_right),
      mgp = c(2.6, 0.6, 0), lwd = 0.7)


  ######################################################################################
  # Console summary ####################################################################
  ######################################################################################

  .printStationStatsSummary(type=type, aggregate.by=aggregate.by, n_loc=length(loc_levels),
                            n_ids=nlevels(data[, id.col]), id.groups=id.groups, n_series=n_series,
                            group.comparisons=group.comparisons, scale=value.scale)


  ######################################################################################
  # Draw panels (rows = type, columns = series) ########################################
  ######################################################################################

  for(r in seq_len(n_types)){
    sc <- type_scale(type[r]); is_count <- sc == "count"
    ymax <- ymax_by_type[r]; ylim <- c(0, ymax * 1.18)
    ylab <- paste0(disp_names[type[r]], if(is_count) " (n)" else " (%)")
    fill <- color.pal[((r - 1) %% length(color.pal)) + 1]

    for(cc in seq_len(n_series)){
      sub <- stats[stats$type == type[r] & stats$series == series_levels[cc], ]
      sub <- sub[order(sub$location), ]
      height <- if(is_count) sub$count else sub$proportion
      height[is.na(height)] <- 0

      # blank barplot to fix coordinates, background, then the bars
      bp <- barplot(height, ylim=ylim, col=NA, border=NA, axes=FALSE, names.arg=NA, space=0.35, xpd=FALSE)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.color, border=NA)
      barplot(height, ylim=ylim, col=fill, border="grey30", axes=FALSE, names.arg=NA, space=0.35, add=TRUE, xpd=FALSE, ...)

      # y-axis + title on the leftmost column (shared per-type scale)
      if(cc == 1){
        if(is_count){
          at <- pretty(c(0, ymax)); at <- at[at <= ymax * 1.05]
          axis(2, at=at, las=1, cex.axis=cex_axis)
        }else{
          at <- seq(0, 1, by=0.2); at <- at[at <= ymax * 1.05 + 1e-9]
          axis(2, at=at, labels=paste0(at * 100, "%"), las=1, cex.axis=cex_axis)
        }
        title(ylab=ylab, cex.lab=cex_lab, line=2.9, xpd=NA)
      }
      usr <- par("usr"); yspan <- usr[4] - usr[3]
      # column header (series name) above the top row, drawn into the outer top margin
      if(r == 1 && n_series > 1)
        text(mean(bp), usr[4] + yspan * 0.07, series_levels[cc], font=2, cex=cex_title, xpd=NA)
      # x-axis (locations) below the bottom row, drawn into the outer bottom margin
      if(r == n_types){
        if(station.labels == "rotated")
          text(bp, usr[3] - yspan * 0.03, loc_levels, srt=90, adj=c(1, 0.5), cex=cex_axis, xpd=NA)
        else
          text(bp, usr[3] - yspan * 0.05, labs_x, adj=c(0.5, 1), cex=cex_axis, xpd=NA)
      }

      # bar annotation (the quantity the axis does not show, by default)
      if(annotate != "none"){
        cnt <- sub$count; prp <- sub$proportion
        lab <- switch(annotate,
          "auto"       = if(is_count) paste0(round(prp * 100), "%")
                         else ifelse(is.na(cnt), paste0(round(prp * 100, 1), "%"), paste0("n=", cnt)),
          "count"      = ifelse(is.na(cnt), "", paste0("n=", cnt)),
          "proportion" = paste0(round(prp * 100), "%"),
          "both"       = ifelse(is.na(cnt), paste0(round(prp * 100, 1), "%"),
                                paste0("n=", cnt, " (", round(prp * 100), "%)")))
        if(!is.null(annot.min)){
          hide <- ifelse(is.na(cnt), prp * 100 < annot.min, cnt < annot.min)
          lab[hide] <- ""
        }
        text(x=bp, y=height + ymax * 0.045, labels=lab, cex=cex_annot, xpd=FALSE)
      }
      box(col="grey40")
    }
  }

  # shared x-axis title and optional overall title (placed at the edges of the outer margins)
  mtext(if(!is.null(xlab)) xlab else tools::toTitleCase(aggregate.by), side=1, outer=TRUE,
        line=xlab_in / csi + 0.3, cex=cex_lab)
  if(!is.null(main)) mtext(main, side=3, outer=TRUE, font=2, cex=cex_title * 1.1, line=omi_top / csi - 0.9)

  # numbered-location key (right outer margin)
  if(show_key){
    x_key <- grconvertX(1 - omi_right / par("din")[1] + 0.04, "ndc", "user")
    y_key <- grconvertY(0.97, "ndc", "user")
    legend(x=x_key, y=y_key, legend=key_labs, bty="n", cex=cex_axis, xpd=NA, y.intersp=1.15)
  }

  invisible(stats)
}


#######################################################################################################
## Statistic computation (internal) ###################################################################
#######################################################################################################

#' Compute per-location statistics as a tidy long table
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.computeStationStats <- function(data, type, id.groups, group.comparisons,
                                 id.col, timebin.col, station.col, aggregate.by){

  all_levels <- levels(data[, aggregate.by])
  station_groups <- if(aggregate.by != station.col) unique(data[, c(station.col, aggregate.by)]) else NULL

  # ID series (groups / comparisons) or a single "All" series
  if(is.null(id.groups)){
    series <- list("All" = levels(data[, id.col]))
    ids_table <- data.frame(ID=levels(data[, id.col]), group="All", stringsAsFactors=FALSE)
  }else{
    series <- .buildGroupComparisonList(id.groups, group.comparisons)
    ids_table <- .meltList(id.groups); colnames(ids_table) <- c("ID", "group")
  }

  rows <- list()
  for(s in names(series)){
    sd <- data[data[, id.col] %in% series[[s]], , drop=FALSE]
    sd[, id.col] <- droplevels(sd[, id.col])
    is_between <- grepl("<->", s, fixed=TRUE)
    for(tp in type){
      v <- .stationMetric(tp, sd, id.col, timebin.col, station.col, aggregate.by,
                          all_levels, ids_table, is_between, station_groups)
      rows[[length(rows) + 1]] <- data.frame(series=s, type=tp, location=all_levels,
                                             count=v$count, proportion=v$proportion,
                                             stringsAsFactors=FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out$series   <- factor(out$series, levels=names(series))
  out$type     <- factor(out$type, levels=type)
  out$location <- factor(out$location, levels=all_levels)
  rownames(out) <- NULL
  out
}


#' Compute a single statistic for one ID series, indexed by location level
#' @keywords internal
#' @noRd

.stationMetric <- function(tp, sd, id.col, timebin.col, station.col, aggregate.by,
                           all_levels, ids_table, is_between, station_groups){

  n_loc <- length(all_levels)
  if(nrow(sd) == 0) return(list(count=rep(0, n_loc), proportion=rep(0, n_loc)))
  loc <- factor(sd[, aggregate.by], levels=all_levels)

  if(tp == "detections"){
    count <- as.numeric(table(loc))
    prop  <- if(sum(count) > 0) count / sum(count) else count

  }else if(tp == "individuals"){
    ids   <- sd[, id.col]
    count <- as.numeric(tapply(ids, loc, function(x) length(unique(x)))); count[is.na(count)] <- 0
    n_ind <- length(unique(ids))
    prop  <- if(n_ind > 0) count / n_ind else count

  }else if(tp == "average detections"){
    tab <- table(sd[, id.col], loc)                    # individual x location detection counts
    tab <- tab[rowSums(tab) > 0, , drop=FALSE]          # individuals actually detected
    val <- if(nrow(tab) > 0) colMeans(tab / rowSums(tab)) else rep(0, n_loc)
    count <- rep(NA_real_, n_loc)                       # no meaningful raw count
    prop  <- as.numeric(val)

  }else{  # co-occurrences
    count <- .stationCooccurrences(sd, ids_table, if(is_between) "between" else "within",
                                   id.col, timebin.col, station.col, aggregate.by, station_groups, all_levels)
    prop  <- if(sum(count) > 0) count / sum(count) else count
  }
  list(count=count, proportion=prop)
}


#' Count co-occurrence events per location for one ID series
#' @keywords internal
#' @noRd

.stationCooccurrences <- function(sd, ids_table, comparison, id.col, timebin.col, station.col,
                                  aggregate.by, station_groups, all_levels){
  zero <- rep(0, length(all_levels))
  if(length(unique(sd[, id.col])) < 2) return(zero)
  wide <- suppressWarnings(createWideTable(sd, id.col=id.col, timebin.col=timebin.col,
                                           value.col=station.col, verbose=FALSE))
  wide <- wide[, -1, drop=FALSE]                        # drop timebin column
  colnames(wide) <- .mapValues(colnames(wide), ids_table$ID, ids_table$group, warn_missing=FALSE)
  events <- unlist(apply(wide, 1, .countCooccurrences, comparison=comparison, groups=colnames(wide)))
  if(is.null(events) || length(events) == 0) return(zero)
  if(aggregate.by != station.col)
    events <- as.character(station_groups[match(events, station_groups[, station.col]), aggregate.by])
  as.numeric(table(factor(events, levels=all_levels)))
}


#######################################################################################################
## Console summary (internal) #########################################################################
#######################################################################################################

#' @keywords internal
#' @noRd
.printStationStatsSummary <- function(type, aggregate.by, n_loc, n_ids, id.groups, n_series,
                                      group.comparisons, scale){
  kv <- .kv
  .summaryOpen("Station statistics")
  kv("Statistics", paste(type, collapse=", "))
  kv("Aggregate", sprintf("%s (%d levels)", aggregate.by, n_loc))
  kv("Individuals", format(n_ids, big.mark=","))
  if(!is.null(id.groups)) kv("Groups", sprintf("%d (%s -> %d series)", length(id.groups), group.comparisons, n_series))
  kv("Bar height", if(scale == "natural") "counts (share for avg. detections)" else scale)
  .summaryClose()
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
