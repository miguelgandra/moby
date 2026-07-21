#######################################################################################################
## Internal helpers: graphics ########################################################################
## Low-level plotting primitives (filled contours, colour legends, panel layout, vertex sizing).
#######################################################################################################


##################################################################################################
## Contrasting text colour #######################################################################

#' Pick a legible text colour for a given background
#'
#' @description Returns "black" or "white" depending on the perceived brightness (relative
#' luminance) of a background colour, so overlaid text stays readable on light or dark fills.
#' @param bg A colour (name, hex or any value accepted by [grDevices::col2rgb()]).
#' @return Either "black" or "white".
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.contrastColor <- function(bg){
  if(is.na(bg)) return("black")
  rgb <- grDevices::col2rgb(bg) / 255
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  ifelse(lum < 0.5, "white", "black")
}


##################################################################################################
## Validate custom shading periods ###############################################################

#' Validate and complete a custom shading-periods data frame
#'
#' @description Checks a user-supplied data frame of background-shading periods (used by the
#' time-series plotting functions), coercing `start`/`end` to POSIXct and filling in a `color`
#' column when absent (a soft palette keyed on `label`, if present).
#' @param shade A data frame with at least `start` and `end` columns, optionally `label` and `color`.
#' @param tz Timezone used to coerce the dates.
#' @return The completed data frame (`start`, `end`, `color`, and `label` if supplied).
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.validateShade <- function(shade, tz) {
  if (!all(c("start", "end") %in% colnames(shade))) {
    stop("A custom 'shade' data frame must contain 'start' and 'end' columns.", call. = FALSE)
  }
  shade$start <- as.POSIXct(shade$start, tz = tz)
  shade$end   <- as.POSIXct(shade$end, tz = tz)
  if (any(is.na(shade$start) | is.na(shade$end))) {
    stop("Some 'start'/'end' values in 'shade' could not be parsed as dates.", call. = FALSE)
  }
  if (!"color" %in% colnames(shade)) {
    if ("label" %in% colnames(shade)) {
      labs <- unique(shade$label)
      pal  <- grDevices::adjustcolor(.economist_pal(length(labs)), alpha.f = 0.35)
      shade$color <- pal[match(shade$label, labs)]
    } else {
      shade$color <- "grey85"
    }
  }
  shade
}



##################################################################################################
## Pretty date axis ##############################################################################

#' Pretty, layout-aware date-axis breaks and labels
#'
#' @description Returns sensible major and minor date breaks (and labels) for a time axis spanning
#' an arbitrary range (hours to years). The break *unit* is derived from the supplied `format` when
#' given (so a `"%b"` format yields monthly breaks rather than mismatched ones), otherwise from the
#' span. The result is coordinated with an optional coarser "band" axis (e.g. a year band drawn
#' above the panel): when the band already encodes the year, the x-axis stays at a finer unit and
#' omits the year from its labels, avoiding redundant duplication.
#' @param start.time,end.time POSIXct limits of the axis.
#' @param n Approximate number of major (labelled) breaks. Defaults to 7.
#' @param format Optional date format (\code{\link[base]{strftime}}). If NULL, both the break unit
#' and the label format are chosen automatically.
#' @param band.format Optional date format used by a coarser companion axis/band (e.g. `"%Y"`). Used
#' to keep the x-axis finer than, and non-redundant with, that band.
#' @param phase Integer phase/offset that re-anchors the breaks *within* each period, after the
#' interval has been chosen. For a yearly unit it selects the month (1 = January ... 12 = December);
#' for finer units it shifts the breaks by `phase - 1` sub-units (days for month/week, hours for
#' day). `1` (default) leaves the natural period start unchanged.
#' @return A list with `at` (major breaks, POSIXct), `labels` (formatted), `minor` (minor breaks,
#' POSIXct), and the chosen `unit` and `format`.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.prettyDateAxis <- function(start.time, end.time, n = 7, format = NULL, band.format = FALSE, phase = 1) {
  rng <- c(as.POSIXct(start.time), as.POSIXct(end.time))
  tz <- .dataTZ(rng[1])
  span_days <- as.numeric(difftime(rng[2], rng[1], units = "days"))

  # does the companion band already encode the year (and nothing finer)?
  band_year <- is.character(band.format) &&
               grepl("%Y|%y", band.format) && !grepl("%m|%b|%B|%d|%H|%U|%W|%j", band.format)

  # target temporal unit: inferred from the user 'format' when supplied (so breaks match the
  # label resolution), otherwise from the span. Never coarser than years; capped below the band.
  unit <- if (!is.null(format)) {
    if (grepl("%H|%M|%S", format)) "hour"
    else if (grepl("%d|%j", format)) "day"
    else if (grepl("%U|%W", format)) "week"
    else if (grepl("%m|%b|%B|%q", format)) "month"
    else if (grepl("%Y|%y", format)) "year"
    else "day"
  } else {
    if (span_days <= 2.5) "hour"
    else if (span_days <= 31) "day"
    else if (span_days <= 150) "week"
    else if (span_days <= 365 * 2.5 || band_year) "month"
    else "year"
  }

  # major step (number of 'unit's) sized to give roughly 'n' labels
  by <- switch(unit,
    "hour"  = max(1, round(span_days * 24 / n)),
    "day"   = max(1, round(span_days / n)),
    "week"  = max(1, round(span_days / 7 / n)),
    "month" = { cand <- c(1, 2, 3, 6, 12); if (band_year) cand <- cand[cand < 12]
                cand[which.min(abs(span_days / 30.44 / cand - n))] },
    "year"  = max(1, round(span_days / 365 / n)))
  unit_word <- c(hour = "hours", day = "days", week = "weeks", month = "months", year = "years")[unit]

  # phase offset re-anchors breaks within each period by (phase - 1) sub-steps. The sub-step is one
  # level finer than the labelling step: months when the step spans several months (e.g. a year, or
  # a 6-month label cadence under a year band) -> shifts which month is labelled; otherwise days
  # (within a month/week), hours (within a day) or minutes (within an hour).
  applyPhase <- function(x) {
    if (phase == 1) return(x)
    k <- phase - 1
    if (unit == "year" || (unit == "month" && by > 1)) return(x + lubridate::period(months = k))
    mult <- switch(unit, "month" = 86400, "week" = 86400, "day" = 3600, "hour" = 60, 0)
    x + k * mult
  }

  base <- lubridate::floor_date(rng[1], unit)
  cap  <- lubridate::ceiling_date(rng[2], unit)

  # major (labelled) breaks on the phased grid
  at <- applyPhase(seq(base, cap, by = paste(by, unit_word)))
  at <- at[at >= rng[1] & at <= rng[2]]

  # default label format for the chosen unit (the year is omitted when the band already shows it)
  if (is.null(format)) {
    format <- switch(unit,
      "hour"  = "%d %b %H:%M",
      "day"   = "%d %b",
      "week"  = "%d %b",
      "month" = if (band_year) "%b" else "%b %Y",
      "year"  = "%Y")
  }

  # minor (unlabelled) ticks: a regular subdivision of the major unit at *meaningful* dates
  # (e.g. monthly majors -> weekly minors), on the SAME phased grid; majors are dropped.
  minor_by <- switch(unit,
    "hour"  = "30 mins",
    "day"   = if (by == 1) "12 hours" else "1 day",
    "week"  = if (by == 1) "1 day" else "1 week",
    "month" = if (by == 1) "1 week" else "1 month",
    "year"  = if (by == 1) "3 months" else "1 year")
  minor <- applyPhase(seq(base, cap, by = minor_by))
  minor <- minor[minor >= rng[1] & minor <= rng[2]]
  minor <- minor[!as.numeric(minor) %in% as.numeric(at)]

  list(at = at, labels = strftime(at, format, tz = tz), minor = minor, unit = unit, format = format)
}



##################################################################################################
## Updated filled.contour function  ##############################################################
## Adapted from https://gist.github.com/epijim/6514388 ###########################################

# Modification by Ian Taylor of the filled.contour function to remove the key and
# facilitate overplotting with contour(). Further modified by Carey McGilliard and
# Bridget Ferris to allow multiple plots on one page

#' Filled Contour Plot with Modified Options
#'
#' @description Produces a filled contour plot with additional modifications for flexibility.
#' @note This function is adapted from the original `filled.contour` function, modified by Ian Taylor, Carey McGilliard, and Bridget Ferris.
#' @keywords internal
#' @noRd

.filled.contour <- function(x = seq(0, 1, length.out = nrow(z)),
                            y = seq(0, 1, length.out = ncol(z)),
                            z,
                            xlim = range(x, finite = TRUE),
                            ylim = range(y, finite = TRUE),
                            zlim = range(z, finite = TRUE),
                            levels = pretty(zlim, nlevels),
                            nlevels = 20,
                            color.palette = cm.colors,
                            col = color.palette(length(levels) - 1),
                            plot.title,
                            plot.axes,
                            key.title,
                            key.axes,
                            asp = NA,
                            xaxs = "i",
                            yaxs = "i",
                            las = 1,
                            axes = TRUE,
                            frame.plot = axes,
                            mar,
                            invert.scale=FALSE,
                            background=NA,
                            ...){

  # input validation
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      } else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    } else stop("no 'z' matrix specified")
  } else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }

  # verify x and y values are increasing
  if (any(diff(x) <= 0) || any(diff(y) <= 0))
    stop("increasing 'x' and 'y' values expected")

  # plot initialization
  plot.new()

  # invert scale if needed
  if(invert.scale){
    zscale <- rev(range(levels))
  }else{
    zscale <- range(levels)
  }

  # set up plotting window
  plot.window(xlim = xlim, ylim = ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)

  # optional background fill (drawn behind the contour, so empty/NA cells show this colour)
  if (!is.na(background)) rect(xlim[1], ylim[1], xlim[2], ylim[2], col = background, border = NA)

  # validate matrix structure of z
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) stop("no proper 'z' matrix specified")
  if (!is.double(z)) storage.mode(z) <- "double"

  # draw filled contour
  graphics::.filled.contour(
    as.double(x), as.double(y), z, as.double(levels),
    col = col
  )


  # draw axes and title
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  } else {
    plot.axes
  }

  # draw frame and title
  if (frame.plot) box()
  if (missing(plot.title)) title(...) else plot.title

  invisible()

}


##################################################################################################
## Updated colorlegend function  #################################################################
## Adapted from 'shape' package (https://rdrr.io/cran/shape/src/R/colorlegend.R)
## - added main.adj + main.inset + support for scientific notation
## - added tick.length
## - added horizontal option
## - added zlab option

#' Color Legend
#'
#' @description Creates a color legend for a plot, adapted from https://rdrr.io/cran/shape/src/R/colorlegend.R with additional features such as main.adj, main.inset, and support for scientific notation.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.colorlegend <- function(col = .viridis_pal(100),
                         zlim,
                         zlevels = 5,
                         dz = NULL,
                         zval = NULL,
                         zlab = NULL,
                         log = FALSE,
                         posx = c(0.9, 0.93),
                         posy = c(0.05, 0.9),
                         main = NULL,
                         main.cex = 1.0,
                         main.col = "black",
                         main.adj = 0.5,
                         lab.col = "black",
                         main.inset = 1,
                         digit = 0,
                         left = FALSE,
                         tick.length = 0.3,
                         lab.scientific = FALSE,
                         horizontal = FALSE,
                         ...) {

  ## Set the number of colors
  ncol <- length(col)

  # Set up a new plot layer without erasing the existing plot
  par(new=TRUE)

  ## Save original margin settings and initialize new margins
  omar <- nmar <- par("mar")

  ## Adjust margins based on orientation
  if (horizontal) {
    nmar[1] <- max(nmar[1], 4)  # Increase bottom margin for horizontal legend
  } else {
    nmar[c(2, 4)] <- 0  # Remove left and right margins for vertical legend
  }

  ## Apply updated margin settings
  par(mar=nmar)

  ## Create an empty plot for the legend without axes or frames
  plot(0, type="n", xlab="", ylab="", asp=1, axes=FALSE, frame.plot=FALSE,
       xlim=c(0, 1), ylim=c(0, 1), xaxs="i", yaxs="i")

  ## Get the plotting area coordinates
  pars <- par("usr")

  ## Calculate the size of the plot area
  dx <- pars[2] - pars[1]
  dy <- pars[4] - pars[3]


  ## Draw the color legend (horizontal or vertical)
  if (horizontal) {
    ## ---- Horizontal Legend ----

    ## Set the position of the legend
    ymin <- pars[3] + posy[1] * dy
    ymax <- pars[3] + posy[2] * dy
    xmin <- pars[1] + posx[1] * dx
    xmax <- pars[1] + posx[2] * dx

    ## Create colored rectangles for the horizontal legend
    X <- seq(xmin, xmax, length.out=ncol+1)
    rect(X[-(ncol+1)], ymin, X[-1], ymax, col=col, border=NA)
    rect(xmin, ymin, xmax, ymax, border=lab.col)

    ## Determine tick labels (either provided by zval or calculated)
    if (!is.null(zval)) {
      zz <- zval
    } else if (is.null(dz) & !is.null(zlevels)) {
      zz <- pretty(zlim, n=(zlevels + 1))  # Generate pretty tick labels
    }

    ## Apply logarithmic scaling if specified
    if (log) zz <- log10(zz)

    ## Draw tick marks and labels for the horizontal legend
    if (!is.null(zz)) {
      Xpos <- xmin + (zz - zlim[1]) / (zlim[2] - zlim[1]) * (xmax - xmin)

      ## Draw ticks
      tick.ystart <- ymin - (tick.length * (ymax - ymin))  # Tick length as fraction of height
      tick.yend <- ymin  # Y position of tick end
      segments(Xpos, tick.ystart, Xpos, tick.yend, col = lab.col)

      ## Format labels (scientific or fixed-point)
      if (lab.scientific) {
        labels <- format(zz, scientific = TRUE)
      } else {
        labels <- formatC(zz, digits = digit, format = "f")
      }

      # Replace with custom labels
      if (!is.null(zlab)) {
        labels[match(zz, zval)] <- zlab
      }


      ## Adjust y-position to ensure visibility
      text(Xpos, tick.ystart - 0.02 * dy, labels, col = lab.col, adj=c(0.5, 1), ...)
    }


  } else {
    ## ---- Vertical Legend ----

    ## Set the position of the legend
    ymin <- pars[3] + posy[1] * dy
    ymax <- pars[3] + posy[2] * dy
    xmin <- pars[1] + posx[1] * dx
    xmax <- pars[1] + posx[2] * dx

    ## Create colored rectangles for the vertical legend
    Y <- seq(ymin, ymax, length.out=ncol+1)
    rect(xmin, Y[-(ncol+1)], xmax, Y[-1], col=col, border=NA)
    rect(xmin, ymin, xmax, ymax, border=lab.col)

    ## Determine tick labels (either provided by zval or calculated)
    if (!is.null(zval)) {
      zz <- zval
    } else if (is.null(dz) & !is.null(zlevels)) {
      zz <- pretty(zlim, n=(zlevels + 1))  # Generate pretty tick labels
    }

    ## Apply logarithmic scaling if specified
    if (log) zz <- log10(zz)

    ## Draw tick marks and labels for the vertical legend
    if (!is.null(zz)) {
      Ypos <- ymin + (zz - zlim[1]) / (zlim[2] - zlim[1]) * (ymax - ymin)

      ## Draw ticks
      tick.xstart <- if (left) {xmin - (tick.length * (xmax - xmin))} else {xmax}
      tick.xend <- if (left) {xmin} else {xmax + (tick.length * (xmax - xmin))}
      segments(tick.xstart, Ypos, tick.xend, Ypos, col = lab.col)

      ## Format labels (scientific or fixed-point)
      if (lab.scientific) {
        labels <- format(zz, scientific = TRUE)
      } else {
        labels <- formatC(zz, digits = digit, format = "f")
      }

      # Replace with custom labels
      if (!is.null(zlab)) {
        labels[match(zz, zval)] <- zlab
      }

      # Adjust horizontal alignment based on legend position
      # 1 for right-alignment, 0 for left-alignment
      adj_value <- if (left) 1 else 0

      ## Add labels next to ticks
      #pos <- if (left) 2 else 4
      text(tick.xend + 0.01 * dx, Ypos, labels, col = lab.col, adj = adj_value, ...)
    }
  }

  ## ---- Main Title ----
  if (!is.null(main)) {
    if (horizontal) {
      # Place the title above the horizontal color scale and center it
      text(mean(c(xmin, xmax)), ymax + 0.05 * dy,  # y-position adjusted to ymax
           labels = main, adj = c(0.5, 0.5), cex = main.cex, col = main.col)
    } else {
      # Keep the title below the vertical color scale
      text(mean(c(xmin, xmax)), ymax + 0.05 * dy,
           labels=main, adj=c(main.adj, 0.5), cex=main.cex, col=main.col)
    }
  }

  ## Reset the margin settings
  par(new=FALSE)
  par(mar=omar)
}


##################################################################################################
## Legend function ###############################################################################
## Adapted by Ben Bolker (http://ms.mcmaster.ca/~bolker/R/misc/legendx.R)

#' Legend
#'
#' @description Creates a customized legend for a plot, with options for color scales, logarithmic scaling, and positioning.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.legend <- function(x, y = NULL, legend, fill = NULL, col = par("col"), border="black",
                   lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"),
                   box.lwd = par("lwd"), box.lty = par("lty"), box.col = par("fg"),
                   box.cex = c(0.8,0.5), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd,
                   xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5),
                   text.width = NULL, text.col = par("col"), text.font = NULL,
                   merge = do.lines && has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE, title = NULL,
                   inset = 0, xpd, title.col = text.col, title.adj = 0.5, seg.len = 2){

    ## the 2nd arg may really be `legend'
    if(missing(legend) && !missing(y) &&
       (is.character(y) || is.expression(y))) {
      legend <- y
      y <- NULL
    }
    mfill <- !missing(fill) || !missing(density)

    if(!missing(xpd)) {
      op <- par("xpd")
      on.exit(par(xpd=op))
      par(xpd=xpd)
    }
    title <- as.graphicsAnnot(title)
    if(length(title) > 1) stop("invalid 'title'")
    legend <- as.graphicsAnnot(legend)
    n.leg <- if(is.call(legend)) 1 else length(legend)
    if(n.leg == 0) stop("'legend' is of length 0")
    auto <-
      if (is.character(x))
        match.arg(x, c("bottomright", "bottom", "bottomleft", "left",
                       "topleft", "top", "topright", "right", "center"))
    else NA

    if (is.na(auto)) {
      xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
      nx <- length(x)
      if (nx < 1 || nx > 2) stop("invalid coordinate lengths")
    } else nx <- 0

    xlog <- par("xlog")
    ylog <- par("ylog")

    rect2 <- function(left, top, dx, dy, density = NULL, angle, ...) {
      r <- left + dx; if(xlog) { left <- 10^left; r <- 10^r }
      b <- top  - dy; if(ylog) {  top <- 10^top;  b <- 10^b }
      rect(left, top, r, b, angle = angle, density = density, ...)
    }
    segments2 <- function(x1, y1, dx, dy, ...) {
      x2 <- x1 + dx; if(xlog) { x1 <- 10^x1; x2 <- 10^x2 }
      y2 <- y1 + dy; if(ylog) { y1 <- 10^y1; y2 <- 10^y2 }
      segments(x1, y1, x2, y2, ...)
    }
    points2 <- function(x, y, ...) {
      if(xlog) x <- 10^x
      if(ylog) y <- 10^y
      points(x, y, ...)
    }
    text2 <- function(x, y, ...) {
      ##--- need to adjust  adj == c(xadj, yadj) ?? --
      if(xlog) x <- 10^x
      if(ylog) y <- 10^y
      text(x, y, ...)
    }
    if(trace)
      catn <- function(...)
        do.call("cat", c(lapply(list(...),formatC), list("\n")))

    cin <- par("cin")
    Cex <- cex * par("cex")		# = the `effective' cex for text

    ## at this point we want positive width even for reversed x axis.
    if(is.null(text.width))
      text.width <- max(abs(strwidth(legend, units="user",
                                     cex=cex, font = text.font)))
    else if(!is.numeric(text.width) || text.width < 0)
      stop("'text.width' must be numeric, >= 0")

    xc <- Cex * xinch(cin[1L], warn.log=FALSE) # [uses par("usr") and "pin"]
    yc <- Cex * yinch(cin[2L], warn.log=FALSE)
    if(xc < 0) text.width <- -text.width

    xchar  <- xc
    xextra <- 0
    yextra <- yc * (y.intersp - 1)
    ## watch out for reversed axis here: heights can be negative
    ymax   <- yc * max(1, strheight(legend, units="user", cex=cex)/yc)
    ychar <- yextra + ymax
    if(trace) catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra,ychar))

    if(mfill) {
      ##= sizes of filled boxes.
      xbox <- xc * box.cex[1]
      ybox <- yc * box.cex[2]
      dx.fill <- xbox ## + x.intersp*xchar
    }
    do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 0))
    ) || !missing(lwd)

    ## legends per column:
    n.legpercol <-
      if(horiz) {
        if(ncol != 1)
          warning(gettextf("horizontal specification overrides: Number of columns := %d",
                           n.leg), domain = NA)
        ncol <- n.leg
        1
      } else ceiling(n.leg / ncol)

    has.pch <- !missing(pch) && length(pch) > 0 # -> default 'merge' is available
    if(do.lines) {
      x.off <- if(merge) -0.7 else 0
    } else if(merge)
      warning("'merge = TRUE' has no effect when no line segments are drawn")

    if(has.pch) {
      if(is.character(pch) && !is.na(pch[1L]) &&
         nchar(pch[1L], type="c") > 1) {
        if(length(pch) > 1)
          warning("not using pch[2..] since pch[1L] has multiple chars")
        np <- nchar(pch[1L], type="c")
        pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
      }
      ##D	if(!merge) dx.pch <- x.intersp/2 * xchar
    }

    if (is.na(auto)) {
      ##- Adjust (x,y) :
      if (xlog) x <- log10(x)
      if (ylog) y <- log10(y)
    }
    if(nx == 2) {
      ## (x,y) are specifiying OPPOSITE corners of the box
      x <- sort(x)
      y <- sort(y)
      left <- x[1L]
      top  <- y[2L]
      w <- diff(x)# width
      h <- diff(y)# height
      w0 <- w/ncol # column width

      x <- mean(x)
      y <- mean(y)
      if(missing(xjust)) xjust <- 0.5
      if(missing(yjust)) yjust <- 0.5

    }
    else {## nx == 1  or  auto
      ## -- (w,h) := (width,height) of the box to draw -- computed in steps
      h <- (n.legpercol + !is.null(title)) * ychar + yc
      w0 <- text.width + (x.intersp + 1) * xchar
      if(mfill)	w0 <- w0 + dx.fill
      ##D	if(has.pch && !merge)	w0 <- w0 + dx.pch
      if(do.lines)		w0 <- w0 + (seg.len + x.off)*xchar
      w <- ncol*w0 + .5* xchar
      if (!is.null(title)
          && (abs(tw <- strwidth(title, units="user", cex=cex) + 0.5*xchar)) > abs(w)) {
        xextra <- (tw - w)/2
        w <- tw
      }

      ##-- (w,h) are now the final box width/height.

      if (is.na(auto)) {
        left <- x - xjust * w
        top	 <- y + (1 - yjust) * h
      } else {
        usr <- par("usr")
        inset <- rep_len(inset, 2)
        insetx <- inset[1L]*(usr[2L] - usr[1L])
        left <- switch(auto, "bottomright"=,
                       "topright"=, "right" = usr[2L] - w - insetx,
                       "bottomleft"=, "left"=, "topleft"= usr[1L] + insetx,
                       "bottom"=, "top"=, "center"= (usr[1L] + usr[2L] - w)/2)
        insety <- inset[2L]*(usr[4L] - usr[3L])
        top <- switch(auto, "bottomright"=,
                      "bottom"=, "bottomleft"= usr[3L] + h + insety,
                      "topleft"=, "top"=, "topright" = usr[4L] - insety,
                      "left"=, "right"=, "center" = (usr[3L] + usr[4L] + h)/2)
      }
    }

    if (plot && bty != "n") { ## The legend box :
      if(trace)
        catn("  rect2(",left,",",top,", w=",w,", h=",h,", ...)",sep="")
      rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
            lwd = box.lwd, lty = box.lty, border = box.col)
    }

    ## (xt[],yt[]) := `current' vectors of (x/y) legend text
    xt <- left + xchar + xextra +
      (w0 * rep.int(0:(ncol-1), rep.int(n.legpercol,ncol)))[1L:n.leg]
    yt <- top -	0.5 * yextra - ymax -
      (rep.int(1L:n.legpercol,ncol)[1L:n.leg] - 1 + !is.null(title)) * ychar

    if (mfill) {		#- draw filled boxes -------------
      if(plot) {
        if(!is.null(fill)) fill <- rep_len(fill, n.leg)
        rect2(left = xt, top=yt+ybox/2, dx = xbox, dy = ybox,
              col = fill,
              density = density, angle = angle, border = border)
      }
      xt <- xt + dx.fill
    }
    if(plot && (has.pch || do.lines))
      col <- rep_len(col, n.leg)

    ## NULL is not documented but people use it.
    if(missing(lwd) || is.null(lwd))
      lwd <- par("lwd") # = default for pt.lwd
    if (do.lines) {			#- draw lines ---------------------
      ## NULL is not documented
      if(missing(lty) || is.null(lty)) lty <- 1
      lty <- rep_len(lty, n.leg)
      lwd <- rep_len(lwd, n.leg)
      ok.l <- !is.na(lty) & (is.character(lty) | lty > 0) & !is.na(lwd)
      if(trace)
        catn("  segments2(",xt[ok.l] + x.off*xchar, ",", yt[ok.l],
             ", dx=", seg.len*xchar, ", dy=0, ...)")
      if(plot)
        segments2(xt[ok.l] + x.off*xchar, yt[ok.l], dx= seg.len*xchar, dy=0,
                  lty = lty[ok.l], lwd = lwd[ok.l], col = col[ok.l])
      # if (!merge)
      xt <- xt + (seg.len+x.off) * xchar
    }
    if (has.pch) {			#- draw points -------------------
      pch   <- rep_len(pch, n.leg)
      pt.bg <- rep_len(pt.bg, n.leg)
      pt.cex<- rep_len(pt.cex, n.leg)
      pt.lwd<- rep_len(pt.lwd, n.leg)
      ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
      x1 <- (if(merge && do.lines) xt-(seg.len/2)*xchar else xt)[ok]
      y1 <- yt[ok]
      if(trace)
        catn("  points2(", x1,",", y1,", pch=", pch[ok],", ...)")
      if(plot)
        points2(x1, y1, pch = pch[ok], col = col[ok],
                cex = pt.cex[ok], bg = pt.bg[ok], lwd = pt.lwd[ok])
      ##D	if (!merge) xt <- xt + dx.pch
    }

    xt <- xt + x.intersp * xchar
    if(plot) {
      if (!is.null(title))
        text2(left + w*title.adj, top - ymax, labels = title,
              adj = c(title.adj, 0), cex = cex, col = title.col)

      text2(xt, yt, labels = legend, adj = adj, cex = cex,
            col = text.col, font = text.font)
    }
    invisible(list(rect = list(w = w, h = h, left = left, top = top),
                   text = list(x = xt, y = yt)))
  }

if (!exists("rep_len")) {
  rep_len <- function(x, length.out) rep(x,length.out)
}


##################################################################################################
## Rescale vertex igraph   #######################################################################
## Sourced from 'netdiffuseR' package (https://github.com/USCCANA/netdiffuseR)

#' Rescale Vertex Igraph
#'
#' @description Rescale vertex size to be used in \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' This function rescales a vertex size before passing it to \code{\link[igraph:plot.igraph]{plot.igraph}}
#' so that the resulting vertices have the desired size relative to the x-axis.
#'
#' @param vertex.size Numeric vector of unscaled vertices' sizes. This is unit-free.
#' @param par.usr Integer vector of length 4 with the coordinates of plotting region.
#'  by default uses \code{par("usr")}.
#' @param minmax.relative.size A numeric vector of length 2. Represents the
#'  desired min and max vertex sizes relative to the x-axis in terms of percentage
#'  (see details).
#' @param adjust Numeric scalar. Adjustment made to the resulting adjusted size
#'  (see details).
#' @author George G. Vega Yon
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.rescale_vertex_igraph <- function(vertex.size, par.usr=par("usr"), minmax.relative.size=c(0.01, 0.04), adjust=200) {

  if (!length(vertex.size)) return(.rescale_vertex_igraph(1, par.usr, minmax.relative.size, adjust))

  # adjusting x
  xrange <- range(vertex.size)
  xscale <- (par.usr[2] - par.usr[1])*minmax.relative.size
  vertex.size <- (vertex.size-xrange[1]+1e-15)/(xrange[2]-xrange[1]+1e-15)*(xscale[2]-xscale[1])+xscale[1]
  return(vertex.size*adjust/2)
}



##################################################################################################
## Set layout matrix   ###########################################################################

#' Set Layout for Grouped Plots
#'
#' @description
#' This function creates a layout matrix for plotting several groups of individuals,
#' with dividers (empty space) between each group. It calculates the required number of rows
#' and columns to accommodate all plots, including optional space for a legend.
#'
#' @param n_cols An integer specifying the number of columns to use for the layout.
#' @param id.groups A list of character vectors, where each vector represents the identifiers
#' of individuals in each group to be plotted. Each group can contain a variable number of individuals.
#' @param plots.height A numeric value indicating the height of each plot. Default is 6.
#' @param dividers.height A numeric value indicating the height of the dividers (empty space)
#' between groups. Default is 1.
#' @param legend A logical value indicating whether to include space for a legend. Default is FALSE.
#' @param min.legend.plots An integer specifying the minimum number of contiguous plots required for the legend. Default is 2.
#' @param expand.legend A logical value determining whether the space allocated for the legend should occupy
#' all remaining plots until the end of the matrix (TRUE) or only the number of plots specified by min.legend.plots (FALSE).
#' Default is TRUE.
#'
#' @details
#' The function calculates the total number of plots required, along with the number of rows
#' needed to fit all plots and the specified dividers. It then constructs a layout matrix that
#' can be used with R's `layout()` or `par(mfrow=)` functions for grouped plotting. The function
#' also calculates the positions of each group within the layout for later reference.
#'
#' @return A list containing:
#' \item{matrix}{A matrix representing the layout of plots, with NA values indicating empty spaces.}
#' \item{heights}{A numeric vector containing the heights of each row in the layout.}
#' \item{group_positions}{A numeric vector containing the calculated vertical positions for each group.}
#'
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.setLayout <- function(n_cols,
                       id.groups,
                       plots.height=6,
                       dividers.height=1,
                       legend,
                       min.legend.plots=2,
                       expand.legend=TRUE) {

  #####################################################################
  # set main plots grid ###############################################

  # number of groups
  n_groups <- length(id.groups)
  # length of each group
  group_length <- sapply(id.groups, length)
  # total number of plots
  total_plots <- sum(group_length)
  # calculate the number of rows needed, including space for dividers
  n_rows <- ceiling(total_plots/n_cols) + (n_groups - 1)

  # initialize the layout matrix and variables to track the current plot number, row height and row index
  layout_mat <- matrix(NA, nrow=n_rows, ncol=n_cols)
  group_positions <- numeric(n_groups)
  row_heights <- numeric(n_rows)
  current_plot <- 1
  row_index <- 1

  # loop through each group
  for (group in seq_along(id.groups)) {
    # reset column index to start from the first column in the next group
    col_index <- 1
    # loop through each individual in the current group
    for (i in seq_along(id.groups[[group]])) {
      # ensure we have enough rows
      if (row_index > nrow(layout_mat)) {
        layout_mat <- rbind(layout_mat, rep(NA, n_cols))
        row_heights <- c(row_heights, 0)
      }
      # place the current plot number in the layout matrix
      layout_mat[row_index, col_index] <- current_plot
      # set the row height for the current row
      row_heights[row_index] <- plots.height
      # move to the next column
      col_index <- col_index + 1
      # if the end of a row is reached, move to the next row
      if (col_index > n_cols) {
        col_index <- 1
        row_index <- row_index + 1
      }
      # increment the current plot number
      current_plot <- current_plot + 1
    }
    # if the current row is not completely filled, move to the next row
    if (col_index != 1) row_index <- row_index + 1
    # add an empty row between groups, except after the last group
    if (group < n_groups) {
      if (row_index > nrow(layout_mat)) {
        layout_mat <- rbind(layout_mat, rep(NA, n_cols))
        row_heights <- c(row_heights, dividers.height)
      } else {
        row_heights[row_index] <- dividers.height
      }
      row_index <- row_index + 1
    }
  }

  # Calculate cumulative heights
  cumulative_heights <- cumsum(row_heights)
  current_plot_index <- 1

  for (group in seq_along(id.groups)) {
    group_length <- length(id.groups[[group]])
    start_row <- which(layout_mat == current_plot_index, arr.ind = TRUE)[1, 1]
    end_row <- which(layout_mat == (current_plot_index + group_length - 1), arr.ind = TRUE)[1, 1]

    # Calculate the midpoint for the current group
    if (start_row == end_row) {
      group_positions[group] <- cumulative_heights[start_row] - row_heights[end_row] / 2
    } else {
      group_positions[group] <- (cumulative_heights[start_row] + cumulative_heights[end_row] - row_heights[end_row]) / 2
    }

    # Update the current_plot_index for the next group
    current_plot_index <- current_plot_index + group_length
  }


  #####################################################################
  # should additional space be added to accommodate for a legend? ####
  if (legend == TRUE) {
    # find the last row that has plots in it
    filled_positions <- which(!is.na(layout_mat), arr.ind = TRUE)
    last_filled_row <- max(filled_positions[, 1])
    last_filled_col <- max(filled_positions[filled_positions[, 1] == last_filled_row, 2])

    # check if there is enough space for the legend
    if ((last_filled_col <= n_cols - min.legend.plots) || (nrow(layout_mat) - last_filled_row > 1)) {
      if (expand.legend) {
        # occupy all remaining plots in the last row for the legend
        layout_mat[last_filled_row, (last_filled_col + 1):n_cols] <- current_plot
      } else {
        # only occupy the specified number of plots for the legend
        if (last_filled_col + min.legend.plots <= n_cols) {
          layout_mat[last_filled_row, (last_filled_col + 1):(last_filled_col + min.legend.plots)] <- current_plot
        } else {
          layout_mat <- rbind(layout_mat, rep(current_plot, n_cols))
          row_heights <- c(row_heights, plots.height)
        }
      }
    } else {
      # not enough space, add a new row
      layout_mat <- rbind(layout_mat, rep(NA, n_cols))
      if(expand.legend)  layout_mat[nrow(layout_mat), seq_len(n_cols)] <- current_plot
      else layout_mat[nrow(layout_mat), 1:min.legend.plots] <- current_plot
      #layout_mat <- rbind(layout_mat, rep(current_plot, n_cols))
      row_heights <- c(row_heights, plots.height)
    }
  }

  # return the layout matrix and row heights
  list(matrix=layout_mat, heights=row_heights, group_positions=group_positions)
}
