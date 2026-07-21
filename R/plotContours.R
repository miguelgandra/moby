#######################################################################################################
# Function to generate contour plots  #################################################################
#######################################################################################################

#' Contour plot (hour x date heatmap of a continuous variable)
#'
#' @description Draws a filled-contour plot of one or more continuous variables (e.g. depth,
#' temperature) over hour-of-day (y) and date (x), aggregated into time and date bins. It is the
#' continuous-variable counterpart of \code{\link{plotChronogram}}: the diel cycle (dawn / sunrise /
#' sunset / dusk boundaries) is overlaid to reveal how the variable varies with time of day across
#' the study period. A separate panel is drawn for each variable and, when `split.by` is supplied,
#' for each group.
#'
#' @details The temporal x-axis can be built in two ways, controlled by `annual.cycle`:
#' \itemize{
#'   \item `annual.cycle = FALSE` (default) uses the **true temporal extent** of the data - the
#'   x-axis spans the observed date range (which may be a few months or several years).
#'   \item `annual.cycle = TRUE` **wraps the data into a single standardised annual cycle**
#'   (Jan-Dec), collapsing all years together. This is useful for comparing the typical seasonal
#'   pattern across individuals or datasets on a common framework, regardless of monitoring period.
#' }
#' Both modes share the same binning (`time.interval`, `date.interval`), diel overlay, grid, colour
#' scaling and legend logic; only the date axis (range and labels) differs.
#'
#' @note The colour-scale legend and panel margins are sized in inch units, so the layout stays
#' consistent across datasets and graphics devices.
#'
#' @inheritParams as_moby
#' @param data A data frame containing the variable(s) to plot and a datetime column.
#' @param variables Column name(s) of the continuous variable(s) to plot (one panel per variable).
#' @param var.titles Optional display names for the variables (e.g. "Depth (m)"). If NULL, the column
#' names are used.
#' @param main Panel title(s). If NULL (default), each panel is titled automatically (the group name
#' when `split.by` is set, the variable name, or both when several variables and groups are shown).
#' Pass a character vector to override (recycled across panels), or FALSE to suppress panel titles.
#' @param plot.title Optional overall title above the whole panel.
#' @param split.by Optional column name(s); a separate set of panels is drawn for each level (or
#' combination of levels) - e.g. species, sex, life stage.
#' @param aggregate.fun Function used to aggregate values within each (id, time bin, date bin).
#' Defaults to the mean.
#' @param color.pal Colours (or a palette function) for the contour fill. If NULL, viridis is used.
#' @param time.interval Bin width along the y-axis (hour of day), as a lubridate period string
#' (e.g. "30 mins", "1 hour", "2 hours"). Defaults to "hour".
#' @param date.interval Bin width along the x-axis (date), as a lubridate period string (e.g. "1 day",
#' "1 week", "1 month"). Defaults to "month".
#' @param date.format Optional \code{\link[base]{strftime}} format for the x-axis labels (e.g. "%b",
#' "%b %Y"). If NULL (default), the label dates and format are chosen automatically from the axis
#' span; the tick positions are always calendar-aligned (whole years/months) rather than arbitrary.
#' @param annual.cycle Logical. If FALSE (default), the x-axis spans the true date range of the data;
#' if TRUE, all data are wrapped into a single standardised annual cycle (Jan-Dec). See Details.
#' @param diel.lines Number of diel-boundary lines to draw: 0, 2 (sunrise/sunset) or 4 (dawn, sunrise,
#' sunset, dusk). Defaults to 4.
#' @param diel.lines.color Colour of the diel lines. Defaults to "black".
#' @param sunriset.coords Longitude/latitude (length-2 numeric, matrix or `SpatialPoints`) for the
#' diel-phase times. Required only when `diel.lines > 0`.
#' @param solar.depth Sun angle below the horizon (degrees) defining twilight. Defaults to 18.
#' @param zlim Optional length-2 numeric setting the range of the colour scale (e.g. \code{c(-10, 10)}
#' for a symmetric scale centred on zero), overriding the data-driven range. Note that values lying
#' \emph{outside} this range are not hidden: they are clamped ("squished") to the nearest limit and
#' drawn in the corresponding extreme colour. A narrower \code{zlim} therefore saturates the extremes
#' rather than leaving those cells blank, so the chosen range controls colour mapping only and never
#' removes data. Defaults to NULL (range taken from the data).
#' @param zlab Label for the colour-scale bar (e.g. units such as "Depth (m)"). If NULL (default), the
#' variable's display name (`var.titles`, or the column name) is used; pass a character vector to set
#' it per variable, or FALSE to omit it. To reverse the colour scale, supply a reversed `color.pal`.
#' @param shared.scale Logical; use a common colour scale across all panels. Defaults to FALSE.
#' @param reverse.scale Logical; reverse the orientation of the colour-bar's value axis so the largest
#' values sit at the \emph{bottom} of the legend (akin to reversing an axis, e.g. \code{ylim(max, min)}),
#' which is often more intuitive for variables such as depth. This flips the legend only - each value
#' keeps its colour, and the plotted colours are unchanged. To also change which colour represents high
#' vs low values, supply a reversed `color.pal` (e.g. \code{rev(...)}). Defaults to FALSE.
#' @param grid Logical; draw faint hour/date guide lines. Defaults to TRUE.
#' @param grid.color Colour of the grid lines. Defaults to a translucent black.
#' @param na.color Background colour for cells with no data (e.g. months outside the monitoring
#' period when `annual.cycle = TRUE`). Defaults to a very light grey.
#' @param cex Global expansion factor for all plot text. Defaults to 1.
#' @param legend Logical; draw the colour-scale legend. Defaults to TRUE.
#' @param sample.size Logical; annotate each panel with the number of individuals contributing
#' (`(n=...)`, top-right). Defaults to TRUE.
#' @param sample.size.color Colour of the `sample.size` annotation. "auto" (default) switches between
#' black and white according to the brightness of the fill behind it; any colour can be given instead.
#' @param ncol Number of columns in the panel layout. Defaults to 1.
#' @param disable.par Logical. If TRUE, the function does not set up its own multi-panel layout (so it
#' can be embedded in a user-managed layout); `ncol` then has no effect. Defaults to FALSE.
#' @details The timezone for binning is taken from the `datetime.col` data (falling back to UTC).
#' @template deviceArgs
#' @param ... Further arguments passed to the underlying filled-contour routine.
#'
#' @return Called for its side effect: it generates contour panel(s).
#' @examples
#' # contour plot of a continuous variable (here latitude) over hour-of-day and date,
#' # with the diel cycle overlaid
#' plotContours(rays, variables = "lat", var.titles = "Latitude",
#'              sunriset.coords = c(-9, 38.4))
#' @export


plotContours <- function(data,
                         variables,
                         var.titles = NULL,
                         main = NULL,
                         plot.title = NULL,
                         split.by = NULL,
                         id.col = NULL,
                         datetime.col = NULL,
                         aggregate.fun = function(x) mean(x, na.rm=TRUE),
                         color.pal = NULL,
                         time.interval = "hour",
                         date.interval = "month",
                         date.format = NULL,
                         annual.cycle = FALSE,
                         diel.lines = 4,
                         diel.lines.color = "black",
                         sunriset.coords = NULL,
                         solar.depth = 18,
                         zlim = NULL,
                         zlab = NULL,
                         shared.scale = FALSE,
                         reverse.scale = FALSE,
                         grid = TRUE,
                         grid.color = adjustcolor("black", 0.10),
                         na.color = grey(0.94),
                         cex = 1,
                         legend = TRUE,
                         sample.size = TRUE,
                         sample.size.color = "auto",
                         ncol = 1,
                         disable.par = FALSE,
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
  data <- as.data.frame(data)

  # timezone is taken from the data (most moby functions do this), avoiding a silent hour shift
  tz <- .dataTZ(data[[datetime.col]])

  errors <- c()
  if(!is.null(var.titles) && length(variables) != length(var.titles)) errors <- c(errors, "Number of 'variables' and 'var.titles' do not match.")
  if(diel.lines > 0 && is.null(sunriset.coords)) errors <- c(errors, "'sunriset.coords' must be provided when diel.lines > 0.")
  valid_period <- function(s) tryCatch({ lubridate::as.period(lubridate::period(s)); TRUE }, error=function(e) FALSE)
  if(!valid_period(time.interval)) errors <- c(errors, "'time.interval' must be a valid lubridate period string (e.g. '30 mins', '1 hour').")
  if(!valid_period(date.interval)) errors <- c(errors, "'date.interval' must be a valid lubridate period string (e.g. '1 day', '1 week', '1 month').")
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse="\n"), call.=FALSE)

  cex_main   <- 1.1 * cex
  cex_lab    <- 1.0 * cex
  cex_axis   <- 0.8 * cex
  cex_legend <- 0.7 * cex

  if(is.null(color.pal)) color.pal <- .viridis_pal
  if(!inherits(color.pal, "function")) color.pal <- grDevices::colorRampPalette(color.pal)

  if(!disable.par){
    original_par <- .savePar()
    on.exit(.restorePar(original_par))
  }


  #####################################################################################
  # Grouping ##########################################################################
  #####################################################################################

  if(!is.null(split.by)){
    for(s in split.by) if(!inherits(data[, s], "factor")) data[, s] <- as.factor(data[, s])
    if(length(split.by) > 1){
      data$dummy_group <- as.factor(apply(data[split.by], 1, paste, collapse=" "))
      split.by <- "dummy_group"
    }
  }else{
    data$dummy_group <- factor("1"); split.by <- "dummy_group"
  }


  #####################################################################################
  # Time- and date-binning ############################################################
  #####################################################################################

  # y-axis: time-of-day bins
  start_time <- as.POSIXct("00:00", format="%H:%M", tz=tz); end_time <- as.POSIXct("23:59", format="%H:%M", tz=tz)
  time_breaks <- format(seq(start_time, end_time, by=time.interval), "%H:%M")
  data$time_bin <- factor(format(lubridate::floor_date(data[[datetime.col]], unit=time.interval), "%H:%M"), levels=time_breaks)
  time_bin_hours <- as.numeric(lubridate::duration(time.interval)) / 3600

  # x-axis: date bins. In 'true extent' mode the axis spans the observed date range; in 'annual
  # cycle' mode the data are wrapped into a single reference year (1972, a leap year) and the axis is
  # forced to a full Jan-Dec span, so different individuals/datasets always align on a common frame.
  ref_times <- if(annual.cycle) as.POSIXct(format(data[[datetime.col]], "1972-%m-%d %H:%M:%S"), tz=tz) else data[[datetime.col]]
  floored <- lubridate::floor_date(ref_times, unit=date.interval)
  if(annual.cycle){
    yr_lo <- lubridate::floor_date(as.POSIXct("1972-01-01 00:00:00", tz=tz), unit=date.interval)
    yr_hi <- lubridate::floor_date(as.POSIXct("1972-12-31 23:59:59", tz=tz), unit=date.interval)
    date_break_times <- seq(yr_lo, yr_hi, by=date.interval)
  }else{
    date_break_times <- seq(min(floored, na.rm=TRUE), max(floored, na.rm=TRUE), by=date.interval)
  }
  bin_key <- function(x) format(x, "%Y-%m-%d %H:%M:%S")
  date_breaks <- bin_key(date_break_times)
  data$date_bin <- factor(bin_key(floored), levels=date_breaks)

  # a contour needs at least a 2x2 grid; catch over-coarse intervals with a clear message
  if(length(date_breaks) < 2)
    stop("The chosen 'date.interval' yields a single date bin. Use a finer date.interval or a longer time series.", call.=FALSE)
  if(length(time_breaks) < 2)
    stop("The chosen 'time.interval' yields a single time bin. Use a finer time.interval.", call.=FALSE)

  # x-axis tick positions/labels. The bins live on an index axis (1..n_bk), so labels are placed on
  # calendar-meaningful dates chosen by .prettyDateAxis (e.g. whole years, anchored to January) and
  # then mapped onto the nearest bin index. This keeps ticks both evenly spaced AND on a consistent
  # calendar unit, instead of subsampling bin indices linearly (which makes the labelled month drift
  # over long ranges, because an even index spacing is rarely a whole number of months).
  n_bk <- length(date_breaks)
  month_unit <- grepl("month", date.interval)
  axis_xlab <- if(month_unit) "Month" else "Date"
  # in 'annual cycle' mode the axis is a single normalised year, so label every month without a year
  ax_fmt <- if(!is.null(date.format)) date.format else if(annual.cycle) "%b" else NULL
  ax <- .prettyDateAxis(min(date_break_times), max(date_break_times), n=if(annual.cycle) 12 else 8, format=ax_fmt)
  axis_at <- vapply(as.numeric(ax$at), function(t) which.min(abs(as.numeric(date_break_times) - t)), integer(1))
  keep <- !duplicated(axis_at)
  axis_at <- axis_at[keep]
  axis_lab <- ax$labels[keep]


  #####################################################################################
  # Aggregate #########################################################################
  #####################################################################################

  # stage 1: summarise within each (id, date bin, time bin, group). drop=TRUE keeps only observed
  # combinations (much smaller than the full id x date x time x group crossing); the empty cells are
  # restored downstream because date/time remain factors with their complete level sets.
  agg <- stats::aggregate(data[, variables], by=list(data[, id.col], data$date_bin, data$time_bin, data[, split.by]),
                          aggregate.fun, simplify=TRUE, drop=TRUE)
  colnames(agg)[1:4] <- c("id", "date", "time", "group")
  colnames(agg)[5:ncol(agg)] <- variables

  # per-(group, variable) sample sizes
  nids <- do.call(rbind, lapply(variables, function(v){
    vd <- agg[!is.na(agg[, v]), c("id", "group")]
    n <- stats::aggregate(vd$id, by=list(group=vd$group), function(x) length(unique(x)))
    data.frame(group=n$group, nids=n$x, var=v)
  }))

  # average across individuals
  agg <- stats::aggregate(agg[, variables], by=list(agg$group, agg$date, agg$time), mean, na.rm=TRUE, drop=FALSE)
  colnames(agg)[1:3] <- c("group", "date", "time"); colnames(agg)[4:ncol(agg)] <- variables
  agg <- agg[order(agg$group, agg$date, agg$time), ]
  agg <- split(agg, f=agg$group, drop=TRUE)
  complete_data <- do.call(rbind, agg)
  groups <- names(agg)


  #####################################################################################
  # Diel-phase boundary lines (per date bin) ##########################################
  #####################################################################################

  daytimes <- NULL
  if(diel.lines > 0){
    day_seq <- seq(min(date_break_times), max(date_break_times), by="1 day")
    st <- getSunTimes(sunriset.coords, min(day_seq), max(day_seq), by="%Y-%m-%d", solar.depth=solar.depth)
    st_date <- as.POSIXct(st$interval, format="%Y-%m-%d", tz=tz)
    st$bin <- match(bin_key(lubridate::floor_date(st_date, date.interval)), date_breaks)
    st <- st[!is.na(st$bin), ]
    daytimes <- stats::aggregate(st[, c("dawns","sunrises","sunsets","dusks")], by=list(bin=st$bin), mean, na.rm=TRUE)
    daytimes[, c("dawns","sunrises","sunsets","dusks")] <- daytimes[, c("dawns","sunrises","sunsets","dusks")] / time_bin_hours
  }

  # y-axis (hours) tick positions
  hour_pos    <- (0:23) / time_bin_hours
  even_pos    <- seq(0, 22, by=2) / time_bin_hours
  even_labs   <- sprintf("%02dh", seq(0, 22, by=2))


  #####################################################################################
  # Layout & device-stable margins ####################################################
  #####################################################################################

  nplots <- length(groups) * length(variables)

  # optional file output: size the device to the panel grid (registered independently of disable.par)
  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule=list(base=2,   slope=3.2, n=ncol, lo=4.5, hi=30),
                     h.rule=list(base=1.2, slope=2.8, n=ceiling(nplots/ncol), lo=4, hi=30),
                     crowd.unit="panels")
    on.exit(grDevices::dev.off(), add=TRUE, after=FALSE)
  }

  show_zlab <- legend && !isFALSE(zlab)
  right_in <- if(show_zlab) 1.30 else if(legend) 0.95 else 0.30
  if(!disable.par){
    rows <- ceiling(nplots / ncol)
    par(mfrow=c(rows, ncol), mgp=c(2.2, 0.6, 0),
        mai=c(0.78, 0.7, 0.45, right_in),
        oma=c(0, 0, if(!is.null(plot.title)) 1.4 else 0, 0))
  }


  #####################################################################################
  # Console summary ###################################################################
  #####################################################################################

  .printContoursSummary(
    n_ids = nlevels(data[, id.col]), variables = variables,
    xmin = min(data[[datetime.col]], na.rm=TRUE), xmax = max(data[[datetime.col]], na.rm=TRUE),
    time.interval = time.interval, date.interval = date.interval, annual.cycle = annual.cycle,
    split.by = if(identical(split.by, "dummy_group")) NULL else split.by, n_groups = length(groups),
    diel.lines = diel.lines, shared.scale = shared.scale)


  #####################################################################################
  # Panels ############################################################################
  #####################################################################################

  faceted <- length(groups) > 1
  multivar <- length(variables) > 1
  p <- 0L

  for(i in seq_along(groups)){
    for(v in seq_along(variables)){
      var <- variables[v]
      p <- p + 1L

      contour_matrix <- matrix(agg[[i]][, var], nrow=length(time_breaks), ncol=length(date_breaks),
                               dimnames=list(time_breaks, date_breaks))
      contour_matrix <- rbind(contour_matrix, contour_matrix[1, ])   # wrap midnight for a closed y-range

      # colour scale
      scale <- if(!is.null(zlim)) zlim
               else if(shared.scale) range(complete_data[, var], na.rm=TRUE)
               else range(contour_matrix, na.rm=TRUE)
      if(!all(is.finite(scale))) scale <- c(0, 1)
      # clamp ("squish") values to the scale so out-of-range cells saturate at the extreme colour
      # instead of being dropped (left blank); a no-op when 'scale' already spans the data
      contour_matrix[] <- pmin(pmax(contour_matrix, min(scale)), max(scale))

      # panel title: user override (recycled / FALSE), else group and/or variable name
      vt <- if(!is.null(var.titles)) var.titles[v] else var
      grp <- paste0(toupper(substring(groups[i], 1, 1)), substring(groups[i], 2))
      plot_title <- if(!is.null(main)){
                      if(isFALSE(main)) "" else main[((p - 1) %% length(main)) + 1]
                    }else if(faceted && multivar){
                      paste0(grp, " - ", vt)
                    }else if(faceted){
                      grp
                    }else{
                      vt
                    }

      .filled.contour(x=seq_along(date_breaks), y=0:length(time_breaks), z=t(contour_matrix),
                      main=plot_title, cex.main=cex_main, background=na.color,
                      xlab=axis_xlab, ylab="Hour", nlevels=100, color.palette=color.pal, cex.lab=cex_lab, zlim=scale,
                      plot.axes = {
                        axis(1, at=axis_at, labels=axis_lab, pos=0, cex.axis=cex_axis)
                        axis(2, at=even_pos, labels=even_labs, pos=1, cex.axis=cex_axis, las=1, tcl=-0.45)
                        axis(2, at=hour_pos, labels=FALSE, pos=1, tcl=-0.24, lwd.ticks=0.6)
                        if(diel.lines > 0){
                          lines(daytimes$bin, daytimes$sunrises, lty=2, col=diel.lines.color)
                          lines(daytimes$bin, daytimes$sunsets,  lty=2, col=diel.lines.color)
                        }
                        if(diel.lines == 4){
                          lines(daytimes$bin, daytimes$dawns, lty=2, col=diel.lines.color)
                          lines(daytimes$bin, daytimes$dusks, lty=2, col=diel.lines.color)
                        }
                        if(grid){
                          abline(v=axis_at, lwd=0.4, col=grid.color)
                          abline(h=hour_pos, lwd=0.4, col=grid.color)
                        }}, ...)

      # sample size (top-right); "auto" colour adapts to the fill brightness behind it for readability
      if(sample.size){
        n_v <- nids$nids[nids$group == groups[i] & nids$var == var]
        if(identical(sample.size.color, "auto")){
          nr <- nrow(contour_matrix); nc <- ncol(contour_matrix)
          bg_val <- mean(contour_matrix[max(1, floor(nr * 0.80)):nr, max(1, floor(nc * 0.72)):nc], na.rm=TRUE)
          denom <- diff(range(scale)); if(denom == 0) denom <- 1
          n_col <- .contrastColor(if(is.nan(bg_val)) na.color
                                  else color.pal(100)[min(max(round((bg_val - min(scale)) / denom * 99) + 1, 1), 100)])
        }else n_col <- sample.size.color
        legend("topright", legend=paste0("(n=", n_v, ")"), bty="n", cex=0.85 * cex, text.col=n_col)
      }

      # colour-scale legend (inch-anchored in the right margin)
      if(legend){
        scale_labs <- pretty(scale, min.n=4); scale_labs <- scale_labs[scale_labs >= min(scale) & scale_labs <= max(scale)]
        # reverse.scale flips the legend's value axis (largest values at the bottom): the colours and
        # value range are reversed together so each value still sits beside its own colour; the plot
        # fill is untouched (the variable is colour-encoded, so its axis lives only in the legend)
        leg_col <- color.pal(100); leg_zlim <- scale
        if(reverse.scale){ leg_col <- rev(leg_col); leg_zlim <- rev(scale) }
        fin <- par("fin")
        bx0 <- (fin[1] - right_in) / fin[1] + 0.12 / fin[1]
        bx1 <- bx0 + 0.15 / fin[1]
        bar_top <- 0.85
        # per-variable colour-bar title (units): explicit zlab (recycled / FALSE) else the variable name.
        # Drawn horizontally just above the bar, left-aligned to its left edge (done before .colorlegend,
        # while the main panel's coordinate system is still active). x uses figure fractions, y plot ones.
        zlab_v <- if(isFALSE(zlab)) NULL else if(!is.null(zlab)) zlab[((v - 1) %% length(zlab)) + 1] else vt
        if(!is.null(zlab_v))
          text(grconvertX(bx0, "nfc", "user"), grconvertY(bar_top + 0.035, "npc", "user"),
               labels=zlab_v, adj=c(0, 0), xpd=NA, cex=cex_legend)
        .colorlegend(col=leg_col, zlim=leg_zlim, zval=scale_labs, digit=max(.decimalPlaces(scale_labs)), xpd=TRUE,
                     posx=c(bx0, bx1), posy=c(0.20, bar_top), main="", main.cex=cex_legend, cex=cex_legend)
      }
    }
  }

  if(!is.null(plot.title)) mtext(plot.title, side=3, cex=1.25*cex, font=2, line=-0.5, outer=TRUE)

  invisible(NULL)
}


#######################################################################################################
## Console summary (internal) #########################################################################
#######################################################################################################

#' @keywords internal
#' @noRd
.printContoursSummary <- function(n_ids, variables, xmin, xmax, time.interval, date.interval,
                                  annual.cycle, split.by, n_groups, diel.lines, shared.scale){
  dur <- as.integer(difftime(xmax, xmin, units="days"))
  kv <- .kv
  .summaryOpen("Contour plot")
  kv("Individuals", format(n_ids, big.mark=","))
  kv("Variables", paste(variables, collapse=", "))
  kv("Period", sprintf("%s to %s (%d d)", format(xmin, "%Y-%m-%d"), format(xmax, "%Y-%m-%d"), dur))
  kv("Binning", sprintf("%s x %s", time.interval, date.interval))
  kv("Time axis", if(annual.cycle) "annual cycle (Jan-Dec)" else "true extent")
  if(!is.null(split.by)) kv("Split by", sprintf("%s (%d groups)", paste(split.by, collapse=", "), n_groups))
  kv("Diel", if(diel.lines > 0) sprintf("%d lines", diel.lines) else "off")
  kv("Scale", if(shared.scale) "shared across panels" else "per panel")
  .summaryClose()
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
