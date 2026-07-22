#######################################################################################################
# Function to generate chronograms (hour x date heatmaps) #############################################
#######################################################################################################

#' Chronogram (hour x date heatmap)
#'
#' @description Draws a chronogram: a two-dimensional hour-of-day by date plot of a chosen metric
#' (number of detections, individuals or co-occurring animals) per time bin. Detections can be shown
#' as size-coded points (`style = "points"`) or colour-coded raster cells (`style = "raster"`), with
#' the diel cycle (dawn / sunrise / sunset / dusk boundaries, and optional day/night/season shading)
#' overlaid to reveal temporal patterns of habitat use. When `split.by` is supplied, a separate panel
#' is drawn for each group (e.g. species).
#'
#' @details Date-axis labels are chosen automatically to suit the temporal span (`date.interval =
#' "auto"`); a numeric `date.interval` reverts to manual mode, and `date.start` shifts the labels (a
#' phase offset in auto mode). For `style = "points"`, `color.by` colours points by any variable while
#' point size encodes the metric; `shade` shades the diel cycle or seasons in the background.
#'
#' The `"co-occurrences"` metric counts the number of distinct individuals sharing the same station
#' within the same time bin (i.e. the local group size, restricted to cells with two or more animals).
#' It is therefore an exploratory descriptor whose magnitude is sensitive to the time-bin width and
#' the spatial resolution of `station.col` (coarser bins/stations inflate apparent co-occurrence) and
#' to detection effort; it is not corrected for the number of animals available. For quantitative
#' co-occurrence/association analysis (with permutation-based null models) use
#' \code{\link{calculateAssociations}}.
#'
#' @note The right margin is sized from the measured legend width and legends are positioned in inch
#' units (via \code{\link[graphics]{grconvertX}} / \code{\link[graphics]{grconvertY}}), so the layout
#' stays consistent across datasets and graphics devices.
#'
#' @inheritParams as_moby
#' @param data A data frame of detections including a time-bin column (`timebin.col`).
#' @param metric The metric to plot: one of `"detections"` (the default), `"individuals"` or
#' `"co-occurrences"` (see Details for the co-occurrence definition and its caveats).
#' @param split.by Optional column name; when supplied, a separate panel is drawn for each of its
#' levels (e.g. species, life stage, habitat).
#' @param style `"points"` (metric encoded by point size) or `"raster"` (metric encoded by colour).
#' @param color.by For `style = "points"`, a column used to colour the points (factor or numeric).
#' @param color.pal Colours (or a palette function) for the points/raster. If NULL, a colourblind-safe
#' default is used (Okabe-Ito for categorical, viridis for continuous/raster).
#' @param coords Longitude/latitude (a length-2 numeric, matrix or `SpatialPoints`) at which
#' to compute diel-phase times. Required only when `diel.lines > 0` or `shade = "diel"`.
#' @param solar.depth Sun angle below the horizon (degrees) defining twilight. Defaults to 18.
#' @param diel.lines Number of diel-boundary lines to draw: 0, 2 (sunrise/sunset) or 4 (dawn,
#' sunrise, sunset, dusk). Defaults to 4.
#' @param shade Background shading for `style = "points"`. One of `"diel"` (day/night/twilight),
#' `"season"`, FALSE (none, default), or a data frame of custom periods with columns `start` and
#' `end` (POSIXct or coercible) and, optionally, `label` and `color` (same structure as
#' \code{\link{plotAbacus}}).
#' @param date.interval Date-axis labels: `"auto"` (default) for span-appropriate "pretty" breaks, or
#' a numeric value selecting every n-th formatted date (manual mode).
#' @param date.format Date format for the labels (\code{\link[base]{strftime}}). If NULL (default),
#' chosen automatically in `"auto"` mode and defaults to `"%d/%b"` in manual mode.
#' @param date.start In manual mode, the index of the first labelled date; in `"auto"` mode, a phase
#' offset re-anchoring the labels within each period. Defaults to 1.
#' @param pt.cex Length-2 numeric giving the min and max point sizes (`style = "points"`). Defaults
#' to c(0.5, 3).
#' @param alpha Opacity of the points, from 0 (transparent) to 1 (opaque). Defaults to 0.7.
#' @param background.color Background colour (`style = "points"`, when `shade` is not used).
#' @param grid Logical; draw faint hour/date guide lines. Defaults to TRUE.
#' @param grid.color Colour of the grid lines. Defaults to "white".
#' @param highlight.isolated Logical; for `style = "points"` with `color.by`, plot isolated colour
#' levels on top so they are not hidden behind denser clusters. Defaults to FALSE.
#' @param shared.scale Use a common metric (colour/size) scale across panels. Defaults to FALSE.
#' @param shared.dates Use a common date range across panels. Defaults to TRUE.
#' @param main Panel title(s). If NULL (default), titles are generated automatically (the group name
#' when `split.by` is used, otherwise the metric); a character vector overrides them (recycled across
#' panels); FALSE omits titles.
#' @param legend Logical; draw the legend(s). Defaults to TRUE.
#' @param legend.cols Number of columns for the `color.by` legend. If NULL, chosen automatically.
#' @param cex Global expansion factor for all plot text. Defaults to 1.
#' @param ncol Number of columns in the panel layout. If NULL, panels are stacked (one column).
#' @template deviceArgs
#' @param ... Further arguments passed to \code{\link[graphics]{image}} (raster) or
#' \code{\link[graphics]{points}} (points).
#'
#' @return Called for its side effect: it generates chronogram panel(s).
#' @examples
#' # chronogram of detections with the diel cycle overlaid
#' plotChronogram(rays, coords = c(-9, 38.4))
#'
#' # without the diel overlay (no coordinates required)
#' plotChronogram(rays, diel.lines = 0)
#' @export


plotChronogram <- function(data,
                           id.col = NULL,
                           timebin.col = NULL,
                           station.col = NULL,
                           metric = "detections",
                           split.by = NULL,
                           style = "points",
                           color.by = NULL,
                           color.pal = NULL,
                           coords = NULL,
                           solar.depth = 18,
                           diel.lines = 4,
                           shade = FALSE,
                           date.interval = "auto",
                           date.format = NULL,
                           date.start = 1,
                           pt.cex = c(0.5, 3),
                           alpha = 0.7,
                           background.color = "white",
                           grid = TRUE,
                           grid.color = "white",
                           highlight.isolated = FALSE,
                           shared.scale = FALSE,
                           shared.dates = TRUE,
                           main = NULL,
                           legend = TRUE,
                           legend.cols = NULL,
                           cex = 1,
                           ncol = NULL,
                           file = NULL,
                           width = NULL,
                           height = NULL,
                           res = 300,
                           ...) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  if(length(metric) != 1 || !metric %in% c("detections", "individuals", "co-occurrences"))
    stop("'metric' must be one of: 'detections', 'individuals', 'co-occurrences'.", call.=FALSE)
  if(!style %in% c("points", "raster")) stop("'style' must be either 'points' or 'raster'.", call.=FALSE)
  if(length(pt.cex) != 2) stop("'pt.cex' should be a vector of length 2 (min and max).", call.=FALSE)
  if(!(identical(date.interval, "auto") || (is.numeric(date.interval) && length(date.interval) == 1)))
    stop("'date.interval' must be either \"auto\" or a single numeric value.", call.=FALSE)
  if(!isFALSE(shade) && !is.data.frame(shade) && !(is.character(shade) && length(shade)==1 && shade %in% c("diel", "season")))
    stop("'shade' must be \"diel\", \"season\", FALSE, or a data frame of periods (start/end, optional label/color).", call.=FALSE)

  # diel-phase times (and 'coords') are only needed when diel lines or diel shading are used
  need_sun <- (is.numeric(diel.lines) && diel.lines > 0) || identical(shade, "diel")
  if(need_sun && is.null(coords))
    stop("'coords' is required when diel.lines > 0 or shade = 'diel'. Set diel.lines = 0 to omit.", call.=FALSE)

  # per-element text sizes derive from the global 'cex'
  cex_axis   <- 1.0 * cex
  cex_legend <- 0.8 * cex
  cex_main   <- 1.2 * cex

  original_par <- .savePar()
  on.exit(.restorePar(original_par))


  #####################################################################################
  # Resolve colour palette ############################################################
  #####################################################################################

  color_is_numeric <- !is.null(color.by) && is.numeric(data[, color.by])
  color_is_factor  <- !is.null(color.by) && !color_is_numeric

  if(style == "raster"){
    if(!is.null(color.by)) warning("- 'color.by' is ignored for raster style.", call.=FALSE)
    raster_fun <- if(is.null(color.pal)) .viridis_pal
                  else if(inherits(color.pal, "function")) color.pal
                  else grDevices::colorRampPalette(color.pal)
  }else{
    if(color_is_factor){
      data[, color.by] <- as.factor(data[, color.by])
      ngroups <- nlevels(data[, color.by])
      if(is.null(color.pal)){
        color.pal <- if(ngroups <= 8) .okabe_ito_pal(ngroups) else grDevices::hcl.colors(ngroups, "Dark 3")
      }else{
        if(inherits(color.pal, "function")) color.pal <- color.pal(ngroups)
        if(length(color.pal) < ngroups) stop("Fewer colours than 'color.by' levels.", call.=FALSE)
      }
    }else if(color_is_numeric){
      color.pal <- if(is.null(color.pal)) .viridis_pal(100)
                   else if(inherits(color.pal, "function")) color.pal(100)
                   else grDevices::colorRampPalette(color.pal)(100)
    }else{
      color.pal <- if(is.null(color.pal)) "grey25"
                   else if(inherits(color.pal, "function")) color.pal(1) else color.pal[1]
    }
  }

  if(is.null(legend.cols)){
    legend.cols <- if(color_is_factor && nlevels(data[, color.by]) >= 15) 2L else 1L
  }


  #####################################################################################
  # Aggregate data (one panel per group) ##############################################
  #####################################################################################

  # timezone of the data: the hour-of-day axis, day columns and shade mapping are all built in the
  # data's own timezone (falls back to UTC when the column carries no tzone), so a study logged in
  # local time gets a correctly-labelled diel axis rather than a UTC-shifted one
  chrono_tz <- .dataTZ(data[, timebin.col])

  interval <- difftime(data[, timebin.col], .lag(data[, timebin.col]), units="min")
  interval <- as.numeric(min(interval[interval > 0], na.rm=TRUE))

  if(!is.null(split.by)){
    if(!inherits(data[, split.by], "factor")) data[, split.by] <- as.factor(data[, split.by])
    data_list <- split(data, f=data[, split.by])
    data_list <- lapply(data_list, function(x){ x[, id.col] <- droplevels(x[, id.col]); x })
    groups <- levels(data[, split.by])
  }else{
    data_list <- list("complete"=data); groups <- "complete"
  }

  plot_data_list <- list(); plot_template <- list(); n_det_total <- 0L

  for(g in seq_along(groups)){
    data_group <- data_list[[groups[g]]]

    src <- if(shared.dates) data else data_group
    first_dates <- lubridate::floor_date(min(src[, timebin.col], na.rm=TRUE), unit="day")
    last_dates  <- lubridate::ceiling_date(max(src[, timebin.col], na.rm=TRUE), unit="day") - 60*60*(interval/60)

    all_timebins <- data.frame(timebin=seq.POSIXt(first_dates, last_dates, by=paste(interval, "mins")))
    all_timebins$day  <- strftime(all_timebins$timebin, "%Y-%m-%d", tz=chrono_tz)
    all_timebins$time <- strftime(all_timebins$timebin, "%H:%M:%S", tz=chrono_tz)
    days <- unique(all_timebins$day); bins <- unique(all_timebins$time)
    plot_matrix <- matrix(NA, nrow=length(days), ncol=length(bins), dimnames=list(days, bins))

    data_group$.row <- if("detections" %in% colnames(data_group)) data_group$detections else 1L
    det_tab <- stats::aggregate(data_group$.row, by=list(data_group[, timebin.col], data_group[, id.col], data_group[, station.col]), sum)
    colnames(det_tab) <- c("timebin", "id", "station", "detections")
    n_det_total <- n_det_total + sum(det_tab$detections, na.rm=TRUE)

    if(metric == "detections"){
      plot_data <- stats::aggregate(det_tab$detections, by=list(det_tab$timebin, det_tab$station), sum)
    }else{
      plot_data <- stats::aggregate(det_tab$id, by=list(det_tab$timebin, det_tab$station), function(x) length(unique(x)))
    }
    colnames(plot_data) <- c("timebin", "station", "var")
    if(metric == "co-occurrences") plot_data <- plot_data[plot_data$var > 1, ]
    plot_data <- plot_data[order(plot_data$timebin), ]
    plot_data$day  <- strftime(plot_data$timebin, "%Y-%m-%d", tz=chrono_tz)
    plot_data$hour <- strftime(plot_data$timebin, "%H:%M:%S", tz=chrono_tz)

    if(color_is_factor || color_is_numeric){
      cd <- stats::aggregate(data_group[, color.by], by=list(data_group[, timebin.col], data_group[, id.col], data_group[, station.col]), .aggFun)
      colnames(cd) <- c("timebin", "id", "station", "level")
      cd <- stats::aggregate(cd$level, by=list(cd$timebin, cd$station), .aggFun)
      colnames(cd) <- c("timebin", "station", "level")
      if(color_is_numeric){
        cd$color <- round(.rescale(cd$level, to=c(1, 100)))
      }else{
        cd$level <- factor(cd$level, levels=levels(data_group[, color.by]))
        cd$color <- as.integer(cd$level)
      }
      plot_data <- .joinKeep(plot_data, cd[, c("timebin","station","color")], by=c("timebin","station"), type="left")
    }else{
      plot_data$color <- 1L
    }

    plot_template[[g]] <- plot_matrix
    plot_data_list[[g]] <- plot_data
  }

  # extend each template by a 4% date margin (breathing room on the index axis)
  extend_range <- function(template){
    dates <- as.Date(rownames(template))
    ext <- floor(as.numeric(max(dates) - min(dates)) * 0.04)
    new_dates <- seq(min(dates) - ext, max(dates) + ext, by="day")
    out <- matrix(NA, nrow=length(new_dates), ncol=ncol(template), dimnames=list(as.character(new_dates), colnames(template)))
    common <- intersect(rownames(out), rownames(template))
    out[common, ] <- template[common, ]
    out
  }
  plot_template <- lapply(plot_template, extend_range)

  var_range_all <- range(unlist(lapply(plot_data_list, function(x) range(x$var, na.rm=TRUE))), na.rm=TRUE)
  if(!all(is.finite(var_range_all))) var_range_all <- c(0, 1)
  same_dates <- length(unique(vapply(plot_data_list, function(x) paste(range(as.character(x$timebin)), collapse="/"), character(1)))) == 1


  #####################################################################################
  # Layout & device-stable margins ####################################################
  #####################################################################################

  nplots <- length(groups)
  if(is.null(ncol)) ncol <- 1L
  rows <- ceiling(nplots / ncol)

  # optional file output: size the device to the panel grid (before the device-metric reads below)
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule=list(base=2.3, slope=3.2, n=ncol, lo=5, hi=30),
                     h.rule=list(base=1.5, slope=3,   n=rows, lo=3.5, hi=30),
                     crowd.unit="panels")
    on.exit(grDevices::dev.off(), add=TRUE, after=FALSE)
  }

  # legend labels that drive the right margin width
  leg_labels <- c(format(pretty(var_range_all)), if(color_is_factor) levels(data[, color.by]), if(color_is_numeric) color.by)
  char_w_in <- par("cin")[1]
  if(legend){
    base_w <- max(nchar(leg_labels), 4) * char_w_in * cex_legend
    right_in <- base_w * (if(color_is_factor) legend.cols else 1) + 0.45
    if(style == "raster") right_in <- max(right_in, 0.8)
  }else right_in <- 0.30

  plot_grid <- matrix(seq_len(nplots), ncol=ncol, byrow=TRUE)
  graphics::layout(plot_grid)
  bottom_plots <- as.numeric(plot_grid)
  right_plots  <- apply(plot_grid, 1, min, na.rm=TRUE)
  # bottom: x-title close to labels; top roomier so stacked panels get a touch more vertical separation
  mai_panel <- c(0.62, 0.8, 0.62, right_in)


  #####################################################################################
  # Console summary ###################################################################
  #####################################################################################

  .printChronogramSummary(
    n_ids = nlevels(data[, id.col]), n_det = n_det_total,
    xmin = min(unlist(lapply(plot_data_list, function(x) min(x$timebin))), na.rm=TRUE),
    xmax = max(unlist(lapply(plot_data_list, function(x) max(x$timebin))), na.rm=TRUE),
    interval = interval, metric = metric, split.by = split.by, n_groups = length(groups),
    style = style, color.by = color.by, color_is_numeric = color_is_numeric,
    diel.lines = if(need_sun) diel.lines else 0, shade = shade,
    date.interval = date.interval, date.start = date.start)


  #####################################################################################
  # Panels ############################################################################
  #####################################################################################

  for(i in seq_len(nplots)){

    plot_data <- plot_data_list[[i]]
    template  <- plot_template[[i]]
    # panel title: user override (recycled / FALSE), else the group (split.by) or the metric
    plot_title <- if(!is.null(main)){
                    if(isFALSE(main)) NULL else main[((i - 1) %% length(main)) + 1]
                  }else if(!is.null(split.by)){
                    paste0(toupper(substring(groups[i], 1, 1)), substring(groups[i], 2))
                  }else{
                    switch(metric, detections="n\u00ba of detections", individuals="n\u00ba of individuals",
                           "co-occurrences"="n\u00ba of co-occurring animals")
                  }

    days <- as.POSIXct(rownames(template), "%Y-%m-%d", tz=chrono_tz)
    bins <- colnames(template)
    day_keys <- strftime(days, "%Y-%m-%d", tz=chrono_tz)

    # diel-phase boundary times (in bin-row units), only if needed
    daytimes <- NULL
    if(need_sun){
      daytimes <- getSunTimes(coords, min(days), max(days), by="%Y-%m-%d", solar.depth)
      daytimes[,-1] <- daytimes[,-1] * (60/interval) + 1
    }

    # x-axis date breaks (auto, span-aware) mapped onto the day-index axis
    if(identical(date.interval, "auto")){
      ax <- .prettyDateAxis(min(days), max(days), format=date.format, phase=date.start)
      maj_idx <- match(strftime(ax$at, "%Y-%m-%d"), day_keys); maj_lab <- ax$labels[!is.na(maj_idx)]; maj_idx <- maj_idx[!is.na(maj_idx)]
      min_idx <- match(strftime(ax$minor, "%Y-%m-%d"), day_keys); min_idx <- min_idx[!is.na(min_idx)]
    }else{
      fmt <- if(is.null(date.format)) "%d/%b" else date.format
      lab <- strftime(days, fmt); rl <- rle(lab); grp_lab <- paste0(lab, "_", rep(seq_along(rl$lengths), rl$lengths))
      uniq <- unique(grp_lab); first_idx <- vapply(uniq, function(x) min(which(grp_lab==x)), integer(1))
      disp <- uniq[seq(date.start, length(uniq), by=date.interval)]
      maj_idx <- vapply(disp, function(x) min(which(grp_lab==x)), integer(1)); maj_lab <- sub("\\_.*", "", disp)
      min_idx <- first_idx[-1]
    }

    hr <- lubridate::hour(lubridate::hms(bins))
    hour_idx <- vapply(0:23, function(x){ w <- which(hr==x); if(length(w)) min(w) else NA_integer_ }, integer(1))

    # vertical grid at every date tick (major + minor), matching the hour grid which spans all hours
    grid_x <- sort(unique(c(maj_idx, min_idx)))

    var_range <- if(shared.scale) var_range_all else range(plot_data$var, na.rm=TRUE)
    if(!all(is.finite(var_range))) var_range <- c(0, 1)

    par(mai=mai_panel, mgp=c(3, 0.6, 0), xpd=FALSE)

    ##########################################################
    if(style == "raster"){
      raster_pal <- raster_fun(256)
      pm <- .castWide(plot_data, "day", "hour", "var", fun.aggregate=max, fill=0)
      miss_h <- base::setdiff(colnames(template), colnames(pm))
      miss_d <- base::setdiff(rownames(template), pm$day)
      if(length(miss_h) > 0){ pm[miss_h] <- NA; pm <- pm[, order(colnames(pm))] }
      if(length(miss_d) > 0) pm <- .rbindFill(pm, data.frame(day=miss_d))
      rownames(pm) <- pm$day; pm <- .dropCols(pm, "day"); pm <- as.matrix(pm[order(rownames(pm)), ])
      plot(1, type="n", xlim=c(1, length(days)), ylim=c(0.5, length(bins)+0.5), xlab="", ylab="", axes=FALSE, xaxs="i", yaxs="i")
      if(grid){
        segments(par("usr")[1], hour_idx, par("usr")[2], hour_idx, lwd=0.3, col=grid.color)
        segments(grid_x, par("usr")[3], grid_x, par("usr")[4], lwd=0.3, col=grid.color)
      }
      image(y=seq_along(bins), x=seq_along(days), z=pm, zlim=var_range, col=raster_pal, add=TRUE, ...)

    ##########################################################
    }else{
      # yaxs="r" adds proportional y-padding so points at hour 0 / 23 are not clipped at the box edge
      plot(0, 0, type="n", xlim=c(1, nrow(template)), ylim=c(1, ncol(template)), xlab="", ylab="", axes=FALSE, xaxs="i", yaxs="r")
      shade_legend <- NULL

      if(is.data.frame(shade)){
        # custom user-defined periods (same input structure as plotAbacus), mapped to day indexes
        sp <- .validateShade(shade, .dataTZ(data[, timebin.col]))
        d_dates <- as.Date(day_keys)
        sp$x0 <- sapply(as.Date(sp$start), function(s){ w <- which(d_dates >= s); if(length(w)) min(w) else NA_integer_ })
        sp$x1 <- sapply(as.Date(sp$end),   function(e){ w <- which(d_dates <= e); if(length(w)) max(w) else NA_integer_ })
        sp <- sp[!is.na(sp$x0) & !is.na(sp$x1) & sp$x1 >= sp$x0, , drop=FALSE]
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.color, border=NA)   # base
        if(nrow(sp) > 0){
          rect(sp$x0, par("usr")[3], sp$x1, par("usr")[4], col=sp$color, border=NA)
          if("label" %in% colnames(sp)){
            shade_legend <- unique(sp[, c("label","color")]); colnames(shade_legend) <- c("label","color")
          }
        }
      }else if(identical(shade, "season")){
        st <- shadeSeasons(min(days), max(days), interval)
        st$start <- vapply(strftime(st$start, "%Y-%m-%d", tz=chrono_tz), function(x){ w<-which(day_keys==x); if(length(w)) min(w) else 1 }, numeric(1))
        st$end   <- vapply(strftime(st$end,   "%Y-%m-%d", tz=chrono_tz), function(x){ w<-which(day_keys==x); if(length(w)) max(w) else length(days) }, numeric(1))
        st$start[1] <- par("usr")[1]; st$end[nrow(st)] <- par("usr")[2]
        rect(st$start, par("usr")[3], st$end, par("usr")[4], col=st$color, border=NA)
        shade_legend <- st[!duplicated(st$season), c("season","color")]
        shade_legend <- shade_legend[order(match(shade_legend$season, c("spring","summer","autumn","winter"))), ]
        colnames(shade_legend) <- c("label","color")
      }else if(identical(shade, "diel")){
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="grey98", border=NA)
        x <- seq_along(daytimes$sunrises)
        nb <- adjustcolor("#CAD9EA", 0.8); dy <- adjustcolor("#FFFBE1", 0.8); cr <- adjustcolor("#EDE8EC", 0.8)
        if(diel.lines == 2){
          polygon(c(x, rev(x)), c(daytimes$sunrises, rep(par("usr")[3], length(x))), col=nb, border=FALSE)
          polygon(c(x, rev(x)), c(daytimes$sunsets,  rep(par("usr")[4], length(x))), col=nb, border=FALSE)
          polygon(c(x, rev(x)), c(daytimes$sunrises, rev(daytimes$sunsets)), col=dy, border=FALSE)
          shade_legend <- data.frame(label=c("day","night"), color=c(dy, nb))
        }else{
          polygon(c(x, rev(x)), c(daytimes$dawns, rep(par("usr")[3], length(x))), col=nb, border=FALSE)
          polygon(c(x, rev(x)), c(daytimes$dusks, rep(par("usr")[4], length(x))), col=nb, border=FALSE)
          polygon(c(x, rev(x)), c(daytimes$sunrises, rev(daytimes$sunsets)), col=dy, border=FALSE)
          polygon(c(x, rev(x)), c(daytimes$dawns, rev(daytimes$sunrises)), col=cr, border=FALSE)
          polygon(c(x, rev(x)), c(daytimes$sunsets, rev(daytimes$dusks)), col=cr, border=FALSE)
          shade_legend <- data.frame(label=c("day","night","twilight"), color=c(dy, nb, cr))
        }
      }else{
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.color, border=NA)
      }

      if(grid){
        segments(par("usr")[1], hour_idx, par("usr")[2], hour_idx, lwd=0.3, col=grid.color)
        segments(grid_x, par("usr")[3], grid_x, par("usr")[4], lwd=0.3, col=grid.color)
      }

      plot_data$cex <- .rescale(plot_data$var, from=var_range, to=pt.cex)
      if(highlight.isolated){
        rl <- rle(plot_data$color); plot_data <- plot_data[order(rep(rl$lengths, rl$lengths), decreasing=TRUE), ]
      }
      xi <- match(plot_data$day, day_keys)
      yi <- match(plot_data$hour, bins)
      pt_col <- if(is.null(color.by)) color.pal[1] else color.pal[plot_data$color]
      points(xi, yi, pch=16, cex=plot_data$cex, col=adjustcolor(pt_col, alpha.f=alpha), ...)
    }

    # diel-boundary lines
    if(need_sun && diel.lines > 0){
      lines(seq_along(daytimes$sunrises), daytimes$sunrises, lty=2)
      lines(seq_along(daytimes$sunsets),  daytimes$sunsets,  lty=2)
      if(diel.lines == 4){
        lines(seq_along(daytimes$dawns), daytimes$dawns, lty=2)
        lines(seq_along(daytimes$dusks), daytimes$dusks, lty=2)
      }
    }

    # title & axes (x- and y-titles drawn with the same cex and offset for visual consistency)
    if(!is.null(plot_title)) title(main=plot_title, cex.main=cex_main, font.main=2, line=1)
    if((!same_dates) || i %in% bottom_plots){
      axis(side=1, at=maj_idx, labels=maj_lab, cex.axis=cex_axis, tcl=-0.45)
      axis(side=1, at=min_idx, labels=FALSE, tcl=-0.24, lwd.ticks=0.6)
      title(xlab="Date", cex.lab=cex_axis, line=1.9)
    }
    if(i %in% right_plots){
      axis(2, at=hour_idx[seq(1, 24, by=2)], labels=paste0(seq(0, 23, by=2), "h"), cex.axis=cex_axis, las=1, tcl=-0.45)
      axis(2, at=hour_idx[seq(2, 24, by=2)], labels=FALSE, tcl=-0.24, lwd.ticks=0.6)
      title(ylab="Hour", cex.lab=cex_axis, line=2.6)
    }
    box()

    #########################################################
    # legends (inch-stacked in the right margin; xpd=NA so margin drawing isn't clipped)
    if(legend){
      par(xpd=NA)
      x_leg <- grconvertX(grconvertX(par("usr")[2], "user", "inches") + 0.12, "inches", "user")
      y_cur <- grconvertY(par("usr")[4], "user", "inches")
      place <- function(y_in, .draw, ...){
        co <- .draw(x=x_leg, y=grconvertY(y_in, "inches", "user"), ...)
        h <- grconvertY(co$rect$top, "user", "inches") - grconvertY(co$rect$top - co$rect$h, "user", "inches")
        y_in - h - 0.10
      }
      if(style == "raster"){
        # colour scale via the package's own colour legend (no 'fields' dependency)
        dx <- par("usr")[2] - par("usr")[1]
        lx <- grconvertX(par("usr")[2] + dx*0.04, "user", "ndc"); lx <- c(lx, lx + 0.018)
        sl <- unique(round(pretty(var_range, min.n=4))); sl <- sl[sl >= min(var_range) & sl <= max(var_range)]
        .colorlegend(col=raster_fun(256), zlim=var_range, zval=sl, digit=0, xpd=TRUE,
                     posx=lx, posy=c(0.3, 0.9), main=.metricLabel(metric), main.cex=cex_legend, cex=cex_legend-0.1)
      }else{
        # all legend blocks share the same left edge (x_leg), base legend geometry and left-aligned
        # titles (title.adj=0), so symbols, points and titles line up consistently
        if(!is.null(shade_legend))
          # bordered squares (pch=22) instead of fill=, so the swatch size is controllable via pt.cex
          y_cur <- place(y_cur, graphics::legend, legend=shade_legend$label, pch=22, pt.bg=shade_legend$color,
                         col="black", pt.cex=1.9, pt.lwd=0.6, bty="n", y.intersp=1.4, cex=cex_legend)
        size_labs <- unique(round(pretty(var_range)))
        size_labs <- size_labs[size_labs >= min(var_range) & size_labs <= max(var_range) & size_labs > 0]
        if(length(size_labs) > 0)
          y_cur <- place(y_cur, graphics::legend, legend=size_labs, pch=16, pt.cex=.rescale(size_labs, from=var_range, to=pt.cex),
                         col=adjustcolor("grey25", alpha.f=alpha), bty="n", y.intersp=1.3, cex=cex_legend,
                         title=.metricLabel(metric), title.adj=0)
        if(color_is_factor){
          place(y_cur, graphics::legend, legend=levels(data[, color.by]), pch=16, pt.cex=1.7,
                col=adjustcolor(color.pal[seq_len(nlevels(data[, color.by]))], alpha.f=alpha),
                bty="n", ncol=legend.cols, y.intersp=1.2, cex=cex_legend, title=color.by, title.adj=0)
        }else if(color_is_numeric){
          dx <- par("usr")[2] - par("usr")[1]
          lx <- grconvertX(par("usr")[2] + dx*0.04, "user", "ndc"); lx <- c(lx, lx + 0.016)
          sl <- pretty(data[, color.by], min.n=4); sl <- sl[sl >= min(data[, color.by]) & sl <= max(data[, color.by])]
          .colorlegend(col=color.pal, zlim=range(data[, color.by], na.rm=TRUE), zval=sl, digit=max(.decimalPlaces(sl)),
                       xpd=TRUE, posx=lx, posy=c(0.4, 0.9), main=color.by, main.cex=cex_legend, cex=cex_legend-0.1)
        }
      }
    }
  }

  invisible(NULL)
}


#######################################################################################################
## Internal helpers ###################################################################################
#######################################################################################################

# most-common (categorical) or mean (numeric) aggregator
.aggFun <- function(x){
  if(is.numeric(x)) mean(x, na.rm=TRUE) else names(which.max(table(x)))
}

# short metric label for the legends
.metricLabel <- function(var) switch(var, detections="detections", individuals="individuals", "co-occurrences"="co-occur.")

#' Print a concise diagnostic summary for a chronogram (internal)
#' @keywords internal
#' @noRd
.printChronogramSummary <- function(n_ids, n_det, xmin, xmax, interval, metric, split.by, n_groups,
                                    style, color.by, color_is_numeric, diel.lines, shade,
                                    date.interval, date.start){
  xmin <- as.POSIXct(xmin, origin="1970-01-01"); xmax <- as.POSIXct(xmax, origin="1970-01-01")
  dur <- as.integer(difftime(xmax, xmin, units="days"))
  kv <- .kv

  diel_desc <- if(diel.lines > 0) sprintf("%d lines", diel.lines) else "off"
  if(!isFALSE(shade)) diel_desc <- paste0(diel_desc, sprintf(" + %s shading", if(is.data.frame(shade)) "custom" else shade))
  date_desc <- if(identical(date.interval, "auto")) "auto" else sprintf("every %g", date.interval)
  if(date.start != 1) date_desc <- paste0(date_desc, sprintf(" (phase %d)", date.start))

  .summaryOpen("Chronogram")
  kv("Individuals", format(n_ids, big.mark=","))
  kv("Detections", format(n_det, big.mark=","))
  kv("Period", sprintf("%s to %s (%d d)", format(xmin, "%Y-%m-%d"), format(xmax, "%Y-%m-%d"), dur))
  kv("Time bin", sprintf("%g min", interval))
  kv("Metric", metric)
  if(!is.null(split.by)) kv("Split by", sprintf("%s (%d groups)", split.by, n_groups))
  kv("Style", style)
  if(!is.null(color.by)) kv("Colour", sprintf("%s (%s)", color.by, if(color_is_numeric) "continuous" else "categorical"))
  kv("Diel", diel_desc)
  kv("Date axis", date_desc)
  .summaryClose()
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
