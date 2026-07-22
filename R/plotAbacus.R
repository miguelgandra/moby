#######################################################################################################
## Function to generate abacus plots  #################################################################
#######################################################################################################

#' Abacus plot
#'
#' @description Produces an abacus plot: a timeline of colour-coded detections for each tagged
#' individual over the monitoring period. Animals can be grouped (e.g. by species, sex or age) via
#' `id.groups`, detections can be coloured by any variable (e.g. receiver/site) through `color.by`,
#' and per-animal release and estimated tag-expiry dates are marked with dedicated symbols. Optional
#' background shading (annual seasons or custom periods) and a labelled top time-band provide
#' temporal context.
#'
#' @details Rows are laid out in the order of `id.groups` (and, within a single group, in the
#' factor-level order of `id.col`); to reorder animals, reorder the `id.col` factor levels or the
#' IDs within `id.groups`.
#'
#' Date-axis labels are, by default, chosen automatically to suit the temporal span of the data
#' (from a few days to several years). Supplying a numeric `date.interval` switches to manual mode,
#' where every n-th formatted date (per `date.format`, starting at `date.start`) is labelled. When
#' `highlight.isolated = TRUE` and a `color.by` variable is set, isolated detections are drawn on
#' top of denser point clusters so they remain visible.
#'
#' @note The layout is device-stable: decorative elements (legend, background shading, the top band and
#' the inter-legend spacing) are sized in physical (inch) units and positioned with
#' \code{\link[graphics]{grconvertX}} / \code{\link[graphics]{grconvertY}}, while the right margin is
#' sized from the measured legend width. Only the data panel stretches with the device, so the
#' appearance stays consistent across datasets and graphics devices of differing dimensions.
#'
#' @inheritParams as_moby
#' @param data A data frame containing animal detections.
#' @param tagging.dates A POSIXct vector of tagging/release dates (single value or named by ID).
#' Inherited from the `mobyData` metadata when available.
#' @param tag.durations Optional numeric vector of estimated tag battery durations (in days), used
#' to mark the estimated tag-expiry date of each animal. Either a single value (applied to all IDs)
#' or a named vector matching the IDs in `id.col`.
#' @param id.groups Optional named list of ID groups used to visually aggregate animals belonging to
#' the same class (e.g. species, sex or age). Each element is a vector of IDs in that group.
#' @param color.by Name of the column used to colour-code detections (e.g. station, habitat or a
#' temporal class). If NULL, all detections are drawn in a single colour.
#' @param color.pal Colours used to plot detections, one per `color.by` level. May be a plain vector
#' (matched to the levels by position) or a **named** vector (matched by name, which is safer when
#' the palette order may not follow the factor-level order).
#' @param discard.missing Logical. If TRUE, animals without detections are omitted. Defaults to FALSE.
#' @param shade Controls the background shading. One of: `TRUE` (default) for annual-season shading
#' (see \code{\link{shadeSeasons}}); `FALSE` for a plain `background.color`; or a data frame of custom
#' periods with columns `start` and `end` (POSIXct or coercible) and, optionally, `label` (used for
#' the legend and to assign a colour per category) and `color`.
#' @param top.band A date format (\code{\link[base]{strftime}}) for the labelled time band drawn
#' above the panel (e.g. `"%Y"` for years). Set to FALSE to omit it. Defaults to `"%Y"`.
#' @param legend Logical. If TRUE (default), a legend (shading periods, release/expiry symbols and
#' `color.by` levels) is drawn in the right margin.
#' @param main Optional plot title.
#' @param bin Optional aggregation of detections before plotting, useful for large datasets (reduces
#' overplotting, file size and render time). One of: `NULL` (default, one point per detection); a
#' unit string (`"hour"`, `"day"`, `"week"`, `"month"`); a number of minutes; or `"auto"` (a bin
#' width chosen from the temporal span). Detections are aggregated by `id.col`, `color.by` (if set)
#' and time bin, so one point represents all detections of that group within the interval.
#' @param scale.by.count Logical. When `bin` is set, if TRUE the point *area* scales with the number
#' of detections in each bin (`pt.cex` is the reference size); if FALSE (default) all points share
#' a constant size. Ignored when `bin = NULL`.
#' @param date.interval Controls the x-axis date labels. Either `"auto"` (default), which generates
#' "pretty" breaks appropriate to the temporal span, or a numeric value selecting every n-th
#' formatted date (manual mode; used together with `date.start` and `date.format`).
#' @param date.format A date format (\code{\link[base]{strftime}}) for the x-axis labels. If NULL
#' (default), it is chosen automatically in `"auto"` mode and falls back to `"%b"` in manual mode.
#' @param date.start Integer controlling label placement. In manual mode (numeric `date.interval`),
#' the index of the first displayed date. In `"auto"` mode, a phase/offset that re-anchors the
#' automatically chosen labels within each period without changing the interval - e.g. for yearly
#' labels it selects the month (1 = January ... 12 = December); for monthly labels, the day of month.
#' Defaults to 1.
#' @param pch Plotting symbol for detections. Defaults to 16 (filled circle).
#' @param pt.cex Expansion factor for the detection points, or `"auto"` to scale the points to the
#' available row height (keeping density consistent across datasets and devices). Defaults to 1.
#' @param alpha Opacity of the detection points, from 0 (transparent) to 1 (opaque). Defaults to 1.
#' @param event.pch Length-2 vector of symbols for the release and tag-expiry markers, ideally named
#' `c(release = , expiry = )`. Defaults to `c(release = 8, expiry = 4)`.
#' @param highlight.isolated Logical. If TRUE and `color.by` is set, isolated detections are brought
#' forward (plotted on top) so they are not hidden behind denser point clusters. Defaults to TRUE.
#' @param background.color Background colour used when `shade = FALSE` (or for uncovered gaps with custom periods). Defaults to "grey96".
#' @param legend.cols Integer number of columns for the `color.by` legend. If NULL (default), it is
#' chosen automatically from the estimated legend height versus the device height.
#' @param cex Global expansion factor applied to all plot text (axis labels, titles, legend and the
#' top band). Detection point size is controlled separately via `pt.cex`. Defaults to 1.
#' @template deviceArgs
#' @param ... Further graphical parameters passed to \code{\link[graphics]{points}} when drawing the
#' detections (e.g. `lwd`, `bg`).
#'
#' @return Called for its side effect: it generates an abacus plot.
#' @examples
#' # abacus plot of the bundled ray detections
#' # (release dates are read from the mobyData metadata)
#' plotAbacus(rays)
#'
#' # colour the detections by receiver station
#' plotAbacus(rays, color.by = "station")
#' @export


plotAbacus <- function(data,
                       id.col = NULL,
                       datetime.col = NULL,
                       id.groups = NULL,
                       tagging.dates = NULL,
                       tag.durations = NULL,
                       color.by = NULL,
                       color.pal = NULL,
                       discard.missing = FALSE,
                       shade = TRUE,
                       top.band = "%Y",
                       legend = TRUE,
                       main = NULL,
                       bin = NULL,
                       scale.by.count = FALSE,
                       date.interval = "auto",
                       date.format = NULL,
                       date.start = 1,
                       pch = 16,
                       pt.cex = 1,
                       alpha = 1,
                       event.pch = c(release = 8, expiry = 4),
                       highlight.isolated = TRUE,
                       background.color = "grey96",
                       legend.cols = NULL,
                       cex = 1,
                       file = NULL,
                       width = NULL,
                       height = NULL,
                       res = 300,
                       ...) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  tag.durations <- reviewed_params$tag.durations

  # check color.by / color.pal consistency
  if(!is.null(color.by) && !is.null(color.pal)) {
    n_lev <- nlevels(data[, color.by])
    if(length(color.pal) < n_lev) stop("The number of supplied colors needs to be greater than or equal to the number of color.by levels", call.=FALSE)
    if(length(color.pal) > n_lev && is.null(names(color.pal))) warning("- The number of specified colors exceeds the number of levels in color.by", call.=FALSE)
  }
  # validate date.interval
  if(!(identical(date.interval, "auto") || (is.numeric(date.interval) && length(date.interval)==1))){
    stop("'date.interval' must be either \"auto\" or a single numeric value.", call.=FALSE)
  }

  # per-element text sizes derive from the global 'cex' (keeps the historical proportions)
  cex_lab    <- 0.8 * cex
  cex_axis   <- 0.7 * cex
  cex_legend <- 0.7 * cex
  cex_mural  <- 0.7 * cex

  # release / expiry symbols (named or positional)
  has_nm      <- !is.null(names(event.pch))
  release_pch <- if(has_nm && "release" %in% names(event.pch)) unname(event.pch[["release"]]) else event.pch[1]
  expiry_pch  <- if(has_nm && "expiry"  %in% names(event.pch)) unname(event.pch[["expiry"]])  else event.pch[2]

  # save the current par settings and ensure they are restored upon function exit
  original_par <- .savePar()
  on.exit(.restorePar(original_par))


  ##############################################################################
  ## Prepare data ##############################################################
  ##############################################################################

  # define single id.group if needed
  if(is.null(id.groups)){
    id.groups <- list(levels(data[,id.col]))
  }

  # name the per-animal vectors so every later lookup is by ID *name* (robust to ordering)
  names(tagging.dates) <- levels(data[,id.col])
  if(!is.null(tag.durations)) names(tag.durations) <- levels(data[,id.col])

  # estimate tag-lifetime end dates (vectorised, class-safe, named)
  end.dates <- if(!is.null(tag.durations)) tagging.dates + tag.durations*86400 else NULL

  # facts for the console summary (captured before any subsetting/aggregation)
  n_det_raw   <- sum(!is.na(data[,datetime.col]))

  # handle missing animals (no detections)
  orig_levels <- levels(data[,id.col])
  missing_IDs <- names(which(table(data[,id.col])==0))
  n_total_ids <- length(orig_levels)
  n_missing   <- length(missing_IDs)
  if(discard.missing){
    if(length(missing_IDs)>0){
      data[,id.col] <- droplevels(data[,id.col])
      keep <- levels(data[,id.col])
      tagging.dates <- tagging.dates[keep]
      if(!is.null(end.dates)) end.dates <- end.dates[keep]
      id.groups <- lapply(id.groups, function(x) x[!x %in% missing_IDs])
      id.groups <- id.groups[lengths(id.groups) > 0]                 # drop now-empty groups
    }
  } else if(length(missing_IDs)>0){
    # add empty placeholder rows so missing animals still appear as empty rows
    dummy_data <- data.frame(x=missing_IDs); colnames(dummy_data) <- id.col
    data <- .rbindFill(data, dummy_data)
    data[,id.col] <- factor(as.character(data[,id.col]), levels=orig_levels)  # preserve original level order
  }

  # row layout: animals grouped by id.groups, with a blank spacer row between groups
  row.groups <- id.groups
  for(i in seq_along(id.groups)){
    if(i==length(id.groups)){break}
    row.groups[[i]] <- c(id.groups[[i]], paste0("blank",i))
  }
  data$row_index <- factor(as.character(data[,id.col]), levels=unlist(row.groups))
  total_rows <- nlevels(data$row_index)
  data$row_index <- as.numeric(data$row_index)

  # optional file output: open a device sized to the data (height grows with the number of rows).
  # Done before the inch-based layout below, so its device-metric queries read the file device.
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base=9,   slope=0,    n=0,          lo=6, hi=12),
                     h.rule = list(base=1.5, slope=0.22, n=total_rows, lo=3, hi=22),
                     crowd.unit="individuals")
    on.exit(grDevices::dev.off(), add=TRUE, after=FALSE)
  }

  # color data by variable (name-aware palette matching)
  if(!is.null(color.by)) {
    groups <- levels(data[,color.by])
    ngroups <- nlevels(data[,color.by])
    if(is.null(color.pal)) {
      # colourblind-safe Okabe-Ito for up to 8 levels; even, distinguishable HCL hues beyond
      color.pal <- if(ngroups <= 8) .okabe_ito_pal(ngroups) else grDevices::hcl.colors(ngroups, "Dark 3")
    }
    # match by name when the palette is named, otherwise by factor-level position
    if(!is.null(names(color.pal))){
      codes <- match(as.character(data[,color.by]), names(color.pal))
      group_colors <- color.pal[match(groups, names(color.pal))]
    }else{
      codes <- as.integer(data[,color.by])
      group_colors <- color.pal[seq_along(groups)]
    }
    data$plot_color <- color.pal[codes]
  }else{
    data$plot_color <- if(is.null(color.pal)) "black" else color.pal[1]
  }

  # optional aggregation: collapse detections into time bins (one point per id x color.by x bin),
  # carrying the per-bin detection count in '.count' (used for optional count-scaled sizing)
  data$.count <- 1L
  if(!is.null(bin)){
    tz <- .dataTZ(data[,datetime.col])
    bin_sec <- .binSeconds(bin, data[,datetime.col])
    centers <- as.POSIXct((floor(as.numeric(data[,datetime.col]) / bin_sec) + 0.5) * bin_sec,
                          origin="1970-01-01", tz=tz)
    cby <- if(!is.null(color.by)) as.character(data[,color.by]) else ""
    key <- paste(data$row_index, cby, as.numeric(centers), sep="\r")   # NA-datetime rows stay unique (distinct row_index)
    first <- !duplicated(key)
    counts <- as.integer(table(key)[key[first]])
    data <- data[first, , drop=FALSE]
    data[,datetime.col] <- centers[first]                              # bin-centre time (NA preserved for placeholders)
    data$.count <- counts
  }
  n_points <- sum(!is.na(data[,datetime.col]))                         # plotted points (= detections if no binning)


  ##############################################################################
  ## Layout: device-stable margins (inches) ####################################
  ##############################################################################

  # character width/height in inches at cex = 1 (device property, available before plotting)
  char_w_in <- par("cin")[1]
  char_h_in <- par("cin")[2]

  # fixed inch geometry for the optional top band (constant height + gap above the panel)
  band_gap_in    <- 0.10
  band_height_in <- 0.17

  # margins (inches): bottom (axis + title), left (ID + group labels), top (band + title), right (legend)
  bottom_in <- 0.80
  left_in   <- if(length(id.groups) > 1) 0.95 else 0.85
  top_in    <- 0.30 + (if(!isFALSE(top.band)) band_gap_in + band_height_in else 0) +
                      (if(!is.null(main)) 0.32 else 0)

  # shading legend labels (seasons, custom-period categories, or none)
  shade_labels <- if(isTRUE(shade)) c("spring","summer","autumn","winter")
                  else if(is.data.frame(shade) && "label" %in% colnames(shade)) unique(as.character(shade$label))
                  else character(0)

  # right margin sized from the (measured) legend content, or minimal when no legend
  if(legend){
    leg_labels <- c(shade_labels,
                    "release date", if(!is.null(end.dates)) "tag lifetime",
                    if(!is.null(color.by)) groups)
    # choose number of legend columns from the estimated legend height vs. available height
    if(!is.null(color.by) && is.null(legend.cols)){
      avail_h_in <- par("din")[2] - top_in - bottom_in
      est_h_in   <- ngroups * cex_legend * char_h_in * 1.2
      legend.cols <- if(est_h_in > 0.9*avail_h_in) 2L else 1L
    }
    if(is.null(legend.cols)) legend.cols <- 1L
    lab_w_in <- max(nchar(leg_labels)) * char_w_in * cex_legend
    right_in <- lab_w_in * legend.cols + 0.22 + 0.18
  }else{
    right_in <- 0.30
  }

  par(mai=c(bottom_in, left_in, top_in, right_in), mgp=c(2.5,0.6,0), xpd=TRUE)


  ##############################################################################
  ## Finalise adaptive decisions & report summary ##############################
  ##############################################################################

  xmin <- min(tagging.dates, na.rm=TRUE)
  xmax <- max(data[,datetime.col], na.rm=TRUE)

  # auto date-axis breaks/labels (computed once; reused for both the summary and drawing).
  # date.start acts as a phase/offset, re-anchoring the auto-chosen labels within each period.
  ax <- if(identical(date.interval, "auto"))
          .prettyDateAxis(xmin, xmax, format=date.format, band.format=top.band, phase=date.start) else NULL

  # resolve density-adaptive point size now that the plot region (par('pin')) is known
  pt_cex_auto <- identical(pt.cex, "auto")
  if(pt_cex_auto){
    row_h_in <- par("pin")[2] / (total_rows + 1)
    pt.cex <- max(0.3, min(1.6, 0.9 * row_h_in / char_h_in))
  }

  .printAbacusSummary(
    n_ids = nlevels(data[,id.col]), n_groups = length(id.groups),
    n_total = n_total_ids, n_missing = n_missing, discard = discard.missing,
    n_det = n_det_raw, n_points = n_points, xmin = xmin, xmax = xmax,
    color.by = color.by, ngroups = if(!is.null(color.by)) ngroups else NA,
    bin = bin, bin_sec = if(!is.null(bin)) bin_sec else NA, scale.by.count = scale.by.count,
    pt.cex = pt.cex, pt_cex_auto = pt_cex_auto, date.interval = date.interval,
    date.format = date.format, ax = ax, top.band = top.band, shade = shade,
    legend = legend, legend.cols = legend.cols)


  ##############################################################################
  ## Generate plot #############################################################
  ##############################################################################

  # empty plot using the default axis expansion (xaxs/yaxs = "r"): R proportionally pads the
  # data range, so edge elements (release/expiry markers, points at the limits) sit comfortably
  # inside the panel without device-specific manual insets. All decorations (box, shading,
  # guides, band) are then drawn to the expanded panel extent given by par("usr").
  plot(x=as.numeric(xmin), y=0, type="n",
       xlim=c(as.numeric(xmin), as.numeric(xmax)), ylim=c(total_rows+1, 0),
       axes=FALSE, xlab="", ylab="")
  usr <- par("usr")
  x0 <- usr[1]; x1 <- usr[2]                       # panel x-extent (data + default expansion)

  date_lim1 <- lubridate::ceiling_date(xmin, "day")
  date_lim2 <- lubridate::floor_date(xmax, "day")
  complete_dates <- seq.POSIXt(date_lim1, date_lim2, by="day", tz=.dataTZ(date_lim1))

  # axis titles (x-title placed clear of the date labels)
  mtext("Date", side=1, line=2.0, cex=cex_lab)
  if(length(id.groups)==1){
    title(ylab="Animal ID", cex.lab=cex_lab, line=2.4)
  }else{
    label_pos <- sapply(id.groups, function(x) mean(unique(data$row_index[as.character(data[,id.col]) %in% x])))
    ok <- is.finite(label_pos)
    if(any(ok)) mtext(names(id.groups)[ok], side=2, at=label_pos[ok], line=2.4, cex=cex_lab)
  }

  # background shading: custom periods (data frame), annual seasons (TRUE), or plain (FALSE)
  shade_legend <- NULL
  if(is.data.frame(shade)){
    sp <- .validateShade(shade, .dataTZ(data[,datetime.col]))
    sp <- sp[as.numeric(sp$end) >= x0 & as.numeric(sp$start) <= x1, , drop=FALSE]   # keep overlapping periods
    rect(xleft=x0, xright=x1, ybottom=0, ytop=total_rows+1, col=background.color, border=NA)   # base layer
    if(nrow(sp) > 0){
      rect(xleft=pmax(as.numeric(sp$start), x0), xright=pmin(as.numeric(sp$end), x1),
           ybottom=0, ytop=total_rows+1, col=sp$color, border=NA)
      if("label" %in% colnames(sp)){
        shade_legend <- unique(sp[, c("label","color")])
        colnames(shade_legend) <- c("label","color")
      }
    }
  } else if(isTRUE(shade)){
    seasons_table <- shadeSeasons(date_lim1, date_lim2, interval=60*24)
    seasons_table$start <- as.numeric(seasons_table$start)
    seasons_table$end <- as.numeric(seasons_table$end)
    seasons_table$start[1] <- x0
    seasons_table$end[nrow(seasons_table)] <- x1
    rect(xleft=seasons_table$start, xright=seasons_table$end, ybottom=0, ytop=total_rows+1, col=seasons_table$color, border=NA)
    shade_legend <- seasons_table[!duplicated(seasons_table$season), c("season","color")]
    colnames(shade_legend) <- c("label","color")
    shade_legend <- shade_legend[order(match(shade_legend$label, c("spring", "summer", "autumn", "winter"))),]
  } else {
    rect(xleft=x0, xright=x1, ybottom=0, ytop=total_rows+1, col=background.color, border=NA)
  }

  # reference for optional count-scaled point sizes (binned data only)
  max_count <- suppressWarnings(max(data$.count, na.rm=TRUE))
  if(!is.finite(max_count) || max_count < 1) max_count <- 1

  # plot detections of each individual
  for(i in seq_len(total_rows)) {
    pts <- data[data$row_index==i, , drop=FALSE]
    if(nrow(pts)==0){next}                       # blank spacer row between groups
    id_name <- as.character(pts[,id.col][1])
    # horizontal guide
    segments(x0=x0, x1=x1, y0=i, lty=2, lwd=0.5)
    # detections (drop placeholder NA rows)
    det <- pts[!is.na(pts[,datetime.col]), , drop=FALSE]
    if(nrow(det)>0){
      if(highlight.isolated){
        # reorder so isolated detections (short colour runs) are plotted last / in front
        det <- det[order(det[,datetime.col]),]
        runs <- rle(det$plot_color)
        runs_length <- rep(runs$lengths, runs$lengths)
        det <- det[order(runs_length, decreasing=TRUE),]
      }
      # point size: optionally scale area with the per-bin detection count
      cex_pts <- if(scale.by.count && !is.null(bin)) pt.cex * (0.8 + 1.7 * sqrt(det$.count / max_count)) else pt.cex
      points(x=det[,datetime.col], y=det$row_index, pch=pch,
             col=adjustcolor(det$plot_color, alpha.f=alpha), cex=cex_pts, ...)
    }
    # release marker (looked up by ID name)
    points(x=tagging.dates[id_name], y=i, pch=release_pch, cex=1.2)
    # estimated tag-expiry marker, if available and within the plotted span
    if(!is.null(end.dates) && !is.na(end.dates[id_name]) && end.dates[id_name]<=date_lim2){
      points(x=end.dates[id_name], y=i, pch=expiry_pch, cex=1.2)
    }
  }

  ##############################################################################
  ## X-axis date labels ########################################################
  ##############################################################################

  if(identical(date.interval, "auto")){
    # pretty, span- and layout-aware breaks (kept finer than, and non-redundant with, top.band)
    axis(side=1, at=ax$at, labels=ax$labels, cex.axis=cex_axis, pos=total_rows+1, tcl=-0.45)
    axis(side=1, at=ax$minor, labels=FALSE, pos=total_rows+1, tcl=-0.24, lwd.ticks=0.75)
  }else{
    # manual mode: every n-th formatted date, beginning at date.start
    fmt <- if(is.null(date.format)) "%b" else date.format
    all_dates <- strftime(complete_dates, fmt)
    consec_dates <- rle(all_dates)
    consec_dates <- paste0(all_dates, "_", rep(seq_along(consec_dates$lengths), consec_dates$lengths))
    unique_dates <- unique(consec_dates)
    indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
    detec_dates <- strftime(seq.POSIXt(xmin, xmax, "day"), fmt)
    start <- min(which(sub("\\_.*", "", unique_dates)==detec_dates[1]))
    end <- max(which(sub("\\_.*", "", unique_dates)==detec_dates[length(detec_dates)]))
    indexes <- indexes[start:end]
    unique_dates <- unique_dates[start:end]
    disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]
    disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
    disp_dates <- sub("\\_.*", "", disp_dates)
    axis(side=1, labels=disp_dates, at=complete_dates[disp_indexes], cex.axis=cex_axis, pos=total_rows+1)
    axis(side=1, labels=FALSE, at=complete_dates[indexes], pos=total_rows+1, tck=-0.016, lwd.ticks=0.8)
  }

  # y-axis (animal IDs) and plot box. pos=x0 keeps the axis flush with the panel's left
  # border (the data extent), matching the x-axis (pos=total_rows+1); otherwise the axis would
  # sit at the padded plot-region edge, detached from the box.
  id_rows <- stats::aggregate(data$row_index, by=list(data[,id.col]), function(z) z[1])$x
  axis(side=2, labels=levels(data[,id.col]), at=id_rows, las=1, cex.axis=cex_axis, pos=x0)
  rect(xleft=x0, xright=x1, ybottom=0, ytop=total_rows+1, col=NA, border="black", xpd=TRUE)

  # panel top in inches (used by both the top band and the title)
  box_top_in <- grconvertY(0, "user", "inches")


  ##############################################################################
  ## Top band (fixed inch thickness, drawn above the panel) ####################
  ##############################################################################

  if(!isFALSE(top.band)){
    # convert a constant inch band above the panel top (user y = 0) into user coordinates.
    # On the inverted y-axis the panel top maps to a *high* inch value, so we ADD inches to
    # move further up (into the top margin) -- never into the data panel.
    y0m <- grconvertY(box_top_in + band_gap_in, "inches", "user")
    y1m <- grconvertY(box_top_in + band_gap_in + band_height_in, "inches", "user")

    # group consecutive equal band labels (e.g. years)
    band_vals <- strftime(complete_dates, top.band)
    r <- rle(band_vals)
    ends    <- cumsum(r$lengths)
    starts  <- c(1, utils::head(ends, -1) + 1)
    centers <- (starts + ends) / 2
    labels  <- r$values
    width_frac <- r$lengths / length(complete_dates)
    divs <- ends[-length(ends)]                      # interior separators only

    # which labels to display: drop too-narrow edge groups; keep only those present in the data
    keep <- labels %in% unique(strftime(data[,datetime.col], top.band, tz=.dataTZ(data[,datetime.col])))
    if(length(keep) >= 1 && width_frac[1] < 0.10) keep[1] <- FALSE
    if(length(keep) >= 1 && width_frac[length(labels)] < 0.10) keep[length(labels)] <- FALSE

    rect(xleft=x0, xright=x1, ybottom=y0m, ytop=y1m, col="black", border="black", xpd=NA)
    if(length(divs) > 0) segments(x0=complete_dates[divs], y0=y0m, y1=y1m, col="white", lwd=1.5, xpd=NA)
    if(any(keep)) text(x=complete_dates[round(centers[keep])], y=mean(c(y0m, y1m)),
                       labels=labels[keep], col="white", cex=cex_mural, adj=0.5, xpd=NA)
  }

  # plot title (above the top band)
  if(!is.null(main)){
    title_y_in <- box_top_in + (if(!isFALSE(top.band)) band_gap_in + band_height_in else 0) + 0.16
    text(x=mean(c(x0, x1)), y=grconvertY(title_y_in, "inches", "user"),
         labels=main, font=2, cex=cex_lab*1.25, xpd=NA)
  }


  ##############################################################################
  ## Legend (stacked in inch units in the right margin) ########################
  ##############################################################################

  if(legend){
    # x just outside the panel right edge (fixed inch pad); cursor starts aligned with the panel top
    x_leg <- grconvertX(grconvertX(x1, "user", "inches") + 0.12, "inches", "user")
    y_cursor_in <- grconvertY(0, "user", "inches")

    # place a legend at the current inch cursor and return the next cursor (below it).
    # NOTE: the drawing-function argument is named '.draw' (not e.g. 'legend_fun') on purpose,
    # so the caller's 'legend=' argument cannot partial-match and clobber it.
    place_legend <- function(y_in, .draw, ...){
      co <- .draw(x=x_leg, y=grconvertY(y_in, "inches", "user"), ...)
      h_in <- grconvertY(co$rect$top, "user", "inches") -
              grconvertY(co$rect$top - co$rect$h, "user", "inches")
      y_in - h_in - 0.10
    }

    if(!is.null(shade_legend)){
      y_cursor_in <- place_legend(y_cursor_in, .legend, legend=shade_legend$label,
                                  fill=shade_legend$color, bty="n", border="black",
                                  box.cex=c(1.6, 1.2), y.intersp=1.4, cex=cex_legend)
    }

    rel_lab <- c("release date", if(!is.null(end.dates)) "tag lifetime")
    rel_pch <- c(release_pch, if(!is.null(end.dates)) expiry_pch)
    y_cursor_in <- place_legend(y_cursor_in, graphics::legend, legend=rel_lab, pch=rel_pch,
                                pt.cex=1.2, bty="n", y.intersp=1.4, cex=cex_legend)

    if(!is.null(color.by)){
      place_legend(y_cursor_in, graphics::legend, legend=groups, pch=pch,
                   col=adjustcolor(group_colors, alpha.f=alpha),
                   bty="n", border=NA, ncol=legend.cols, pt.cex=pt.cex*1.2,
                   y.intersp=1.2, cex=cex_legend)
    }
  }

  invisible(NULL)
}


#######################################################################################################
## Console summary (internal) #########################################################################
#######################################################################################################

#' Print a concise diagnostic summary for an abacus plot
#'
#' @description Reports the input data and the (adaptive) plotting decisions made by
#' \code{\link{plotAbacus}} as a compact, aligned console block.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.printAbacusSummary <- function(n_ids, n_groups, n_total, n_missing, discard, n_det, n_points,
                                xmin, xmax, color.by, ngroups, bin, bin_sec, scale.by.count,
                                pt.cex, pt_cex_auto, date.interval, date.format, ax, top.band,
                                shade, legend, legend.cols) {

  dur <- as.integer(difftime(xmax, xmin, units="days"))

  ids_desc <- format(n_ids, big.mark=",")
  if(n_groups > 1) ids_desc <- paste0(ids_desc, sprintf(" (%d groups)", n_groups))
  if(n_missing > 0) ids_desc <- paste0(ids_desc, sprintf("; %d missing (%s)",
                       n_missing, if(discard) "removed" else "shown"))

  bin_label <- if(is.null(bin)) NULL else {
    lab <- switch(as.character(bin_sec), "3600"="hourly", "86400"="daily",
                  "604800"="weekly", "2629800"="monthly", sprintf("%g-min", bin_sec/60))
    if(identical(bin, "auto")) paste0(lab, " [auto]") else lab
  }

  date_desc <- if(identical(date.interval, "auto")){
    unit_lab <- c(hour="hourly", day="daily", week="weekly", month="monthly", year="yearly")[ax$unit]
    sprintf("auto -> %s (\"%s\")", unname(unit_lab), ax$format)
  } else sprintf("every %g (\"%s\")", date.interval, if(is.null(date.format)) "%b" else date.format)

  shade_desc <- if(is.data.frame(shade)) sprintf("custom (%d periods)", nrow(shade)) else
                if(isTRUE(shade)) "seasonal" else "none"

  # narrow, width-robust layout: short labels, compact values, long info split over two lines
  kv   <- .kv
  cont <- .kvCont
  .summaryOpen("Abacus plot")
  kv("Individuals", ids_desc)
  kv("Detections", format(n_det, big.mark=","))
  kv("Period", sprintf("%s to %s (%d d)", format(xmin, "%Y-%m-%d"), format(xmax, "%Y-%m-%d"), dur))
  if(!is.null(color.by)) kv("Colour", sprintf("%s (%d levels)", color.by, ngroups))
  if(!is.null(bin)){
    kv("Binning", sprintf("%s bins", bin_label))
    cont(sprintf("%s pts; mean %.1f det/pt%s", format(n_points, big.mark=","),
                 n_det/n_points, if(scale.by.count) "; area-scaled" else ""))
  }
  kv("Points", sprintf("cex %.2f%s", pt.cex, if(pt_cex_auto) " [auto]" else ""))
  kv("Date axis", date_desc)
  if(!isFALSE(top.band)) kv("Top band", sprintf("\"%s\"", top.band))
  kv("Shading", shade_desc)
  if(legend) kv("Legend", sprintf("%d col%s", legend.cols, if(legend.cols > 1) "s" else ""))
  .summaryClose()
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
