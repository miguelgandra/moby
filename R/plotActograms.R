##################################################################################################
## Function to generate actogram plots ###########################################################
##################################################################################################

#' Actogram plot (per-individual diel activity)
#'
#' @description Produces an actogram: for each tagged individual, detections are plotted by
#' time-of-day (y) across the monitoring period (x), so diel patterns and their seasonal shift are
#' visible per animal. It is the raw, per-individual counterpart of \code{\link{plotChronogram}}
#' (which aggregates the population into a single panel). Detections can be coloured by any variable
#' (e.g. receiver or habitat) via `color.by`, animals can be grouped with `id.groups`, the diel
#' cycle (sunrise/sunset, optionally dawn/dusk) is overlaid, and each panel marks the release date
#' and, when `tag.durations` are supplied, the estimated tag-expiry date.
#'
#' @details Panels are laid out in the order of `id.groups` (and, within a group, in the factor-level
#' order of `id.col`); to reorder, reorder the `id.col` factor levels or the IDs within `id.groups`.
#'
#' Date-axis labels are, by default, chosen automatically to suit the temporal span of the data (from
#' a few days to several years), with calendar-aligned breaks. Supplying a numeric `date.interval`
#' switches to manual mode, where every n-th formatted date (per `date.format`, starting at
#' `date.start`) is labelled. The timezone for the time-of-day axis is taken from the data.
#'
#' @inheritParams as_moby
#' @param data A data frame containing animal detections.
#' @param tagging.dates A POSIXct vector of tagging/release dates (single value or named by ID).
#' Inherited from the `mobyData` metadata when available.
#' @param tag.durations Optional numeric vector of estimated tag battery durations (in days), used to
#' mark each animal's estimated tag-expiry date. Either a single value (applied to all IDs) or a
#' named vector matching the IDs in `id.col`.
#' @param id.groups Optional named list of ID groups used to visually aggregate animals belonging to
#' the same class (e.g. species, sex or age). Each element is a vector of IDs in that group.
#' @param discard.missing Logical. If TRUE (default), animals without detections are omitted; if
#' FALSE, an empty panel is drawn for each.
#' @param color.by Name of the column used to colour-code detections (e.g. station or habitat). If
#' NULL, all detections are drawn in a single colour and no colour legend is shown.
#' @param color.pal Colours used to plot detections, one per `color.by` level. May be a plain vector
#' (matched to the levels by position) or a named vector (matched by name). If NULL, a colourblind-
#' safe palette is used (Okabe-Ito for up to 8 levels, HCL "Dark 3" beyond).
#' @param sunriset.coords Longitude/latitude (length-2 numeric, matrix or `SpatialPoints`) used for
#' the diel-phase times. Required only when `diel.lines > 0`.
#' @param diel.lines Number of diel-boundary lines to overlay: 0 (none), 2 (sunrise/sunset) or 4
#' (dawn, sunrise, sunset, dusk). Defaults to 2.
#' @param solar.depth Sun angle below the horizon (degrees) defining twilight. Defaults to 18.
#' @param date.interval Controls the x-axis date labels. Either `"auto"` (default), which generates
#' calendar-aligned "pretty" breaks appropriate to the temporal span, or a numeric value selecting
#' every n-th formatted date (manual mode; used with `date.start` and `date.format`).
#' @param date.format A date format (\code{\link[base]{strftime}}) for the x-axis labels. If NULL
#' (default), it is chosen automatically in `"auto"` mode and falls back to `"%b/%y"` in manual mode.
#' @param date.start Integer controlling label placement. In manual mode, the index of the first
#' displayed date. In `"auto"` mode, a phase/offset that re-anchors the automatically chosen labels
#' within each period (e.g. for yearly labels it selects the month). Defaults to 1.
#' @param pch Plotting symbol for detections. Defaults to 16 (filled circle).
#' @param pt.cex Expansion factor for the detection points. Defaults to 1.6.
#' @param alpha Opacity of the detection points, from 0 (transparent) to 1 (opaque). Defaults to 1.
#' @param highlight.isolated Logical. If TRUE (default) and `color.by` is set, isolated detections are
#' brought forward (plotted on top) so they are not hidden behind denser point clusters.
#' @param background.color Background colour of each panel. Defaults to "grey96".
#' @param grid Logical; draw faint hour/date guide lines. Defaults to FALSE.
#' @param grid.color Colour of the grid lines. Defaults to "white".
#' @param main Optional overall title above the whole panel.
#' @param legend Logical. If TRUE (default) and `color.by` is set, a colour legend is drawn.
#' @param legend.cols Integer number of columns for the `color.by` legend. If NULL (default), chosen
#' automatically.
#' @param ncol Number of columns in the panel grid. If NULL (default), set automatically (1 for a
#' single individual, 2 otherwise).
#' @param cex Global expansion factor applied to all plot text. Point size is controlled separately
#' via `pt.cex`. Defaults to 1.
#' @template deviceArgs
#' @param ... Further graphical parameters passed to \code{\link[graphics]{points}}.
#'
#' @return Called for its side effect: it generates actogram panel(s).
#' @examples
#' # per-individual actograms with the diel cycle overlaid
#' # (coordinates are required to compute sunrise/sunset times)
#' plotActograms(rays, sunriset.coords = c(-9, 38.4))
#' @export


plotActograms <- function(data,
                          id.col = NULL,
                          datetime.col = NULL,
                          tagging.dates = NULL,
                          tag.durations = NULL,
                          id.groups = NULL,
                          discard.missing = TRUE,
                          color.by = NULL,
                          color.pal = NULL,
                          sunriset.coords = NULL,
                          diel.lines = 2,
                          solar.depth = 18,
                          date.interval = "auto",
                          date.format = NULL,
                          date.start = 1,
                          pch = 16,
                          pt.cex = 1.6,
                          alpha = 1,
                          highlight.isolated = TRUE,
                          background.color = "grey96",
                          grid = FALSE,
                          grid.color = "white",
                          main = NULL,
                          legend = TRUE,
                          legend.cols = NULL,
                          ncol = NULL,
                          cex = 1,
                          file = NULL,
                          width = NULL,
                          height = NULL,
                          res = 300,
                          ...) {

  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  tag.durations <- reviewed_params$tag.durations

  if(diel.lines > 0 && is.null(sunriset.coords))
    stop("'sunriset.coords' must be provided when diel.lines > 0.", call.=FALSE)
  if(!(identical(date.interval, "auto") || (is.numeric(date.interval) && length(date.interval)==1)))
    stop("'date.interval' must be either \"auto\" or a single numeric value.", call.=FALSE)
  if(!is.null(color.by) && !is.null(color.pal)){
    n_lev <- nlevels(data[, color.by])
    if(length(color.pal) < n_lev) stop("The number of supplied colors needs to be greater than or equal to the number of color.by levels.", call.=FALSE)
    if(length(color.pal) > n_lev && is.null(names(color.pal))) warning("- The number of specified colors exceeds the number of levels in color.by.", call.=FALSE)
  }

  # per-element text sizes derive from the global 'cex' (keeps the historical proportions)
  cex_title  <- 2.2 * cex
  cex_lab    <- 1.8 * cex
  cex_axis   <- 1.4 * cex
  cex_legend <- 1.6 * cex

  tz <- .dataTZ(data[, datetime.col])
  single_color <- is.null(color.by)
  show_legend  <- legend && !single_color

  original_par <- .savePar()
  on.exit(.restorePar(original_par))


  ##############################################################################
  # Prepare data ###############################################################
  ##############################################################################

  if(is.null(id.groups)) id.groups <- list(levels(data[, id.col]))

  # name the per-animal vectors so every later lookup is by ID name (robust to ordering)
  names(tagging.dates) <- levels(data[, id.col])
  if(!is.null(tag.durations)) names(tag.durations) <- levels(data[, id.col])
  end.dates <- if(!is.null(tag.durations)) tagging.dates + tag.durations*86400 else NULL

  # handle missing animals (no detections)
  orig_levels <- levels(data[, id.col])
  missing_IDs <- names(which(table(data[, id.col])==0))
  n_total_ids <- length(orig_levels); n_missing <- length(missing_IDs)
  if(discard.missing && length(missing_IDs) > 0){
    data[, id.col] <- droplevels(data[, id.col])
    keep <- levels(data[, id.col])
    tagging.dates <- tagging.dates[keep]
    if(!is.null(end.dates)) end.dates <- end.dates[keep]
    id.groups <- lapply(id.groups, function(x) x[!x %in% missing_IDs])
    id.groups <- id.groups[lengths(id.groups) > 0]
  }
  # ordered IDs to plot (id.groups order, restricted to existing levels)
  plot_ids <- unlist(id.groups, use.names=FALSE)
  plot_ids <- plot_ids[plot_ids %in% levels(data[, id.col])]

  # time-of-day (hours) and day, in the data's timezone
  data <- data[order(data[, datetime.col]), ]
  data$day  <- lubridate::floor_date(data[, datetime.col], "day")
  data$hour <- as.numeric(difftime(data[, datetime.col], data$day, units="hours"))

  # colour data by variable (name-aware palette matching), as in plotAbacus
  if(!single_color){
    groups  <- levels(data[, color.by]); ngroups <- nlevels(data[, color.by])
    if(is.null(color.pal)) color.pal <- if(ngroups <= 8) .okabe_ito_pal(ngroups) else grDevices::hcl.colors(ngroups, "Dark 3")
    if(!is.null(names(color.pal))){
      codes <- match(as.character(data[, color.by]), names(color.pal))
      group_colors <- color.pal[match(groups, names(color.pal))]
    }else{
      codes <- as.integer(data[, color.by]); group_colors <- color.pal[seq_along(groups)]
    }
    data$plot_color <- color.pal[codes]
  }else{
    groups <- NULL; ngroups <- NA
    data$plot_color <- if(is.null(color.pal)) "black" else color.pal[1]
  }
  data_individual <- split(data, f=data[, id.col])

  n_det <- sum(!is.na(data[, datetime.col]))


  ##############################################################################
  # Axis & diel preparation ####################################################
  ##############################################################################

  xmin <- min(c(tagging.dates, data[, datetime.col]), na.rm=TRUE)   # include release dates so markers show
  xmax <- max(data[, datetime.col], na.rm=TRUE)
  date_lim1 <- lubridate::floor_date(xmin, "day")
  date_lim2 <- lubridate::ceiling_date(xmax, "day")
  complete_dates <- seq.POSIXt(date_lim1, date_lim2, by="day", tz=tz)

  # diel-phase boundary times (per day), only when requested
  daytimes_table <- NULL
  if(diel.lines > 0){
    daytimes_table <- getSunTimes(sunriset.coords, date_lim1, date_lim2, by="%Y-%m-%d", solar.depth=solar.depth)
    daytimes_table$interval <- as.POSIXct(daytimes_table$interval, "%Y-%m-%d", tz=tz)
  }

  # x-axis breaks/labels: calendar-aligned pretty breaks (auto) or every n-th formatted date (manual)
  ax <- if(identical(date.interval, "auto")) .prettyDateAxis(date_lim1, date_lim2, format=date.format, phase=date.start) else NULL
  if(identical(date.interval, "auto")){
    ax_at <- ax$at; ax_lab <- ax$labels; ax_minor <- ax$minor
  }else{
    fmt <- if(is.null(date.format)) "%b/%y" else date.format
    all_dates <- strftime(complete_dates, fmt)
    consec <- rle(all_dates)
    consec_dates <- paste0(all_dates, "_", rep(seq_along(consec$lengths), consec$lengths))
    unique_dates <- unique(consec_dates)
    ax_minor <- complete_dates[vapply(unique_dates, function(x) min(which(consec_dates==x)), integer(1))]
    disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
    ax_at <- complete_dates[vapply(disp_dates, function(x) min(which(consec_dates==x)), integer(1))]
    ax_lab <- sub("\\_.*", "", disp_dates)
  }

  hour_labels <- sprintf("%02dh", seq(0, 24, by=2))


  ##############################################################################
  # Layout #####################################################################
  ##############################################################################

  if(is.null(ncol)) ncol <- if(length(plot_ids) == 1) 1 else 2

  # optional file output: size the device to the panel grid (before .setLayout/layout below)
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule=list(base=5, slope=3.25, n=max(0, ncol-1), lo=5, hi=30),
                     h.rule=list(base=2, slope=2.75, n=ceiling(length(plot_ids)/ncol), lo=4, hi=30),
                     crowd.unit="individuals")
    on.exit(grDevices::dev.off(), add=TRUE, after=FALSE)
  }
  if(show_legend && is.null(legend.cols)) legend.cols <- min(3, max(1, ngroups))

  layout_params <- .setLayout(ncol, id.groups, plots.height=6, dividers.height=1, legend=show_legend)
  nplots <- max(layout_params$matrix, na.rm=TRUE)
  layout_params$matrix[is.na(layout_params$matrix)] <- nplots + 1   # blank cells -> empty regions
  layout(mat=layout_params$matrix, heights=layout_params$heights)

  oma <- if(length(id.groups) > 1) c(1, 8, 1, 2) else c(1, 4, 1, 2)
  if(!is.null(main)) oma[3] <- oma[3] + 2
  par(mar=c(3, 2, 3, 0), oma=oma)

  n_panels <- if(show_legend) nplots - 1 else nplots
  first_col <- layout_params$matrix[, 1]


  ##############################################################################
  # Console summary ############################################################
  ##############################################################################

  .printActogramSummary(
    n_ids = length(plot_ids), n_groups = length(id.groups), n_total = n_total_ids,
    n_missing = n_missing, discard = discard.missing, n_det = n_det, xmin = xmin, xmax = xmax,
    color.by = color.by, ngroups = ngroups, diel.lines = diel.lines,
    date.interval = date.interval, date.format = date.format, ax = ax, legend = show_legend)


  ##############################################################################
  # Draw panels ################################################################
  ##############################################################################

  for(i in seq_len(n_panels)){
    selected_id <- plot_ids[i]
    data_plot   <- data_individual[[selected_id]]

    # empty panel canvas (fixed hour range; common date range across panels)
    plot(x=date_lim1, y=0, type="n", axes=FALSE, xaxs="i",
         xlim=c(date_lim1, date_lim2), ylim=c(0, 24), xlab="", ylab="",
         main=selected_id, cex.main=cex_title)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.color, border=NA)

    if(grid) abline(h=0:24, v=ax_minor, lwd=0.05, col=grid.color)

    # x-axis (date)
    axis(side=1, at=ax_at, labels=ax_lab, cex.axis=cex_axis, tck=-0.04)
    axis(side=1, at=ax_minor, labels=FALSE, tck=-0.02, lwd.ticks=0.5)
    # y-axis (hour-of-day) on first-column panels only
    if(i %in% first_col){
      title(ylab="Hour", cex.lab=cex_lab, line=4, xpd=NA)
      axis(2, at=seq(0, 24, by=2), labels=hour_labels, tck=-0.04, cex.axis=cex_axis, las=1)
      axis(2, at=0:24, labels=FALSE, tck=-0.02, lwd.ticks=0.5)
    }

    # detections
    if(nrow(data_plot) > 0){
      det <- data_plot
      if(highlight.isolated && !single_color){
        det <- det[!is.na(det[, color.by]), , drop=FALSE]
        det <- det[order(det[, datetime.col]), ]
        runs <- rle(det$plot_color)
        det <- det[order(rep(runs$lengths, runs$lengths), decreasing=TRUE), ]
      }
      points(x=det$day, y=det$hour, col=adjustcolor(det$plot_color, alpha.f=alpha), pch=pch, cex=pt.cex, ...)
    }

    # release marker line, and estimated tag-expiry line when available
    abline(v=tagging.dates[selected_id], lwd=1.4)
    if(!is.null(end.dates) && !is.na(end.dates[selected_id])) abline(v=end.dates[selected_id], lwd=1.4)

    # diel-phase lines
    if(diel.lines > 0){
      lines(daytimes_table$interval, daytimes_table$sunrises, lty=2)
      lines(daytimes_table$interval, daytimes_table$sunsets,  lty=2)
      if(diel.lines == 4){
        lines(daytimes_table$interval, daytimes_table$dawns, lty=2)
        lines(daytimes_table$interval, daytimes_table$dusks, lty=2)
      }
    }
    box()
  }

  ##############################################################################
  # Colour legend ##############################################################
  ##############################################################################

  if(show_legend){
    par(mar=c(1, 1, 3, 1))
    plot.new()
    legend("center", legend=groups, pch=pch, col=adjustcolor(group_colors, alpha.f=alpha),
           ncol=legend.cols, cex=cex_legend, bg=background.color, pt.cex=cex_legend+0.6, y.intersp=1.4)
  }

  # id.group labels (rotated, in the left outer margin)
  if(length(id.groups) > 1){
    label_pos <- layout_params$group_positions
    layout_height <- sum(layout_params$heights)
    label_pos <- grconvertY(1 - (label_pos / layout_height), "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels=names(id.groups),
         srt=90, cex=cex_title + 0.2, font=2, xpd=NA, adj=c(0.5, 0.5))
  }

  # overall title
  if(!is.null(main)) mtext(main, side=3, outer=TRUE, font=2, cex=cex_lab*1.1, line=0.2)

  invisible(NULL)
}


##################################################################################################
## Console summary (internal) ####################################################################
##################################################################################################

#' Print a concise diagnostic summary for an actogram plot
#'
#' @description Reports the input data and plotting decisions made by \code{\link{plotActograms}} as
#' a compact, aligned console block.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.printActogramSummary <- function(n_ids, n_groups, n_total, n_missing, discard, n_det, xmin, xmax,
                                  color.by, ngroups, diel.lines, date.interval, date.format, ax, legend){
  dur <- as.integer(difftime(xmax, xmin, units="days"))

  ids_desc <- format(n_ids, big.mark=",")
  if(n_groups > 1) ids_desc <- paste0(ids_desc, sprintf(" (%d groups)", n_groups))
  if(n_missing > 0) ids_desc <- paste0(ids_desc, sprintf("; %d missing (%s)", n_missing, if(discard) "removed" else "shown"))

  date_desc <- if(identical(date.interval, "auto")){
    unit_lab <- c(hour="hourly", day="daily", week="weekly", month="monthly", year="yearly")[ax$unit]
    sprintf("auto -> %s (\"%s\")", unname(unit_lab), ax$format)
  } else sprintf("every %g (\"%s\")", date.interval, if(is.null(date.format)) "%b/%y" else date.format)

  kv <- .kv
  .summaryOpen("Actogram plot")
  kv("Individuals", ids_desc)
  kv("Detections", format(n_det, big.mark=","))
  kv("Period", sprintf("%s to %s (%d d)", format(xmin, "%Y-%m-%d"), format(xmax, "%Y-%m-%d"), dur))
  if(!is.null(color.by)) kv("Colour", sprintf("%s (%d levels)", color.by, ngroups))
  kv("Diel", if(diel.lines > 0) sprintf("%d lines", diel.lines) else "off")
  kv("Date axis", date_desc)
  kv("Legend", if(legend) "shown" else "none")
  .summaryClose()
}

##################################################################################################
##################################################################################################
##################################################################################################
