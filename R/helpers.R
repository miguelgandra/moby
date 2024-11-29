#######################################################################################################
#######################################################################################################
# HELPER FUNCTIONS ####################################################################################
#######################################################################################################
#######################################################################################################

# Set of functions for internal use within the 'moby' package


##################################################################################################
## Import namespaces   ###########################################################################

#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
NULL


##################################################################################################
## Print to console   ############################################################################

#' Print to console
#'
#' @description Prints a string to the console with a specific formatting.
#' @param string A character string to be printed to the console.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal
#' @noRd

.printConsole <- function(string){
  wrapped_text <- strwrap(string, width=getOption("width")*1.2)
  cat(paste0("\033[0;", 1, "m", wrapped_text, "\033[0m", "\n"))
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
## Get position function  ########################################################################
## - get plot coordinates by keyword (adapted from graphics::legend) #############################

#' Get Position
#'
#' @description Calculates the position of a keyword on a plot with optional inset adjustments.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.getPosition <- function(keyword, inset, bar.width = 0) {
  usr <- par("usr")

  # Ensure `inset` has exactly two values
  if (length(inset) == 1) inset <- rep_len(inset, 2)

  # Calculate inset offsets for x and y
  insetx <- inset[1] * (usr[2] - usr[1])
  insety <- inset[2] * (usr[4] - usr[3])

  # Determine left and top positions based on keyword
  left <- switch(keyword,
                 bottomright = usr[2] - insetx - bar.width,
                 topright = usr[2] - insetx - bar.width,
                 right = usr[2] - insetx - bar.width,
                 bottomleft = usr[1] + insetx,
                 left = usr[1] + insetx,
                 topleft = usr[1] + insetx,
                 bottom = (usr[1] + usr[2]) / 2 - bar.width / 2,
                 top = (usr[1] + usr[2]) / 2 - bar.width / 2,
                 center = (usr[1] + usr[2]) / 2 - bar.width / 2,
                 stop("Invalid keyword for position: ", keyword)
  )

  top <- switch(keyword,
                bottomright = usr[3] + insety,
                bottom = usr[3] + insety,
                bottomleft = usr[3] + insety,
                topleft = usr[4] - insety,
                top = usr[4] - insety,
                topright = usr[4] - insety,
                left = (usr[3] + usr[4]) / 2,
                right = (usr[3] + usr[4]) / 2,
                center = (usr[3] + usr[4]) / 2,
                stop("Invalid keyword for position: ", keyword)
  )

  return(c(left, top))
}



##################################################################################################
## Scale Bar function ############################################################################
## Adapted from raster::scalebar #################################################################
##  - added new 'bar.lwd' argument to set the border line width of the scalebars and
##  - added new 'bar.height' argument to control the scalebar thickness
##  - added new 'label.color' argument to control the color of the text labels
##  - added new 'label.offset' argument to control the spacing of the text labels

#' Scale Bar
#'
#' @description Adds a scale bar to a plot, supporting both line and bar types.
#'
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.scalebar <- function (d, xy=NULL, type="line", divs=2, below="", bar.lwd=0.4, bar.height=1.5,
                      lonlat=NULL, label, adj=c(0.5, -0.5), lwd=2, label.color="black", label.offset=0.6, ...){

  # check if scale bar type is valid
  stopifnot(type %in% c("line", "bar"))

  # retrieve current plot parameters
  pr <- graphics::par()

  # determine if coordinates are in longitude/latitude format based on plot range
  if (is.null(lonlat)) {
    lonlat <- pr$usr[1] > -181 & pr$usr[2] < 181 & pr$yaxp[1] > -200 & pr$yaxp[2] < 200
  }

  # if in longitude/latitude, calculate bar distance based on latitude midpoint
  if (lonlat) {
    lat <- mean(pr$yaxp[1:2])
    if (missing(d)) {
      # Estimate distance d if missing, based on plot range
      dx <- (pr$usr[2] - pr$usr[1]) / 10
      d <- raster::pointDistance(cbind(0, lat), cbind(dx, lat), TRUE)
      d <- signif(d / 1000, 2)  # Convert meters to kilometers and round
      label <- NULL
    }
    # calculate end point of the scale bar
    p <- cbind(0, lat)
    dd <- .destPoint(p, d * 1000)[1, 1]  # Convert km to meters
  } else {
    # in projected units, estimate distance d if missing, based on plot width
    if (missing(d)) d <- round(10 * (pr$usr[2] - pr$usr[1]) / 10) / 10
    dd <- d
  }

  # set default scale bar position if `xy` is not provided
  if (is.null(xy)) {
    padding <- c(5, 5) / 100  # Padding as a percentage of plot range
    parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
    xy <- c(pr$usr[1] + padding[1] * parrange[1], pr$usr[3] + padding[2] * parrange[2])
  }

  # adjust `adj` for label offset (distance from scale bar)
  adj <- c(0.5, -label.offset)

  if (type == "line") {
    # draw a line for the scale bar
    lines(matrix(c(xy[1], xy[2], xy[1] + dd, xy[2]), byrow = TRUE, nrow = 2), lwd = bar.lwd, ...)
    # set label text if not provided, defaulting to distance `d`
    label <- if (missing(label)) paste(d) else label
    text(xy[1] + 0.5 * dd, xy[2], labels = label, adj = adj, col = label.color, ...)

  } else if (type == "bar") {
    # for segmented bar, calculate bar height based on `dd` and `bar.height`
    stopifnot(divs > 0)
    lwd <- dd / 25 * bar.height

    # if 2 divisions, draw a two-color bar and label start, middle, and end
    if (divs == 2) {
      half <- xy[1] + dd / 2
      graphics::polygon(c(xy[1], xy[1], half, half), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]), col = "white", lwd = bar.lwd)
      graphics::polygon(c(half, half, xy[1] + dd, xy[1] + dd), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]), col = "black", lwd = bar.lwd)
      label <- if (missing(label)) c("0", "", d) else label
      text(xy[1], xy[2] + lwd, labels = label[1], adj = adj, col = label.color, ...)
      text(xy[1] + 0.5 * dd, xy[2] + lwd, labels = label[2], adj = adj, col = label.color, ...)
      text(xy[1] + dd, xy[2] + lwd, labels = label[3], adj = adj, col = label.color, ...)

    }

    # if `below` text provided, place it below the scale bar with adjusted `adj`
    if (below != "") {
      adj[2] <- label.offset  # Adjust adj to place below text further away
      text(xy[1] + 0.5 * dd, xy[2]-lwd, labels = below, adj = adj, col = label.color, ...)
    }
  }
}


.destPoint <- function (p, d, b=90, r=6378137) {
  toRad <- pi/180
  lon1 <- p[, 1] * toRad
  lat1 <- p[, 2] * toRad
  b <- b * toRad
  lat2 <- asin(sin(lat1) * cos(d/r) + cos(lat1) * sin(d/r) * cos(b))
  lon2 <- lon1 + atan2(sin(b) * sin(d/r) * cos(lat1), cos(d/r) - sin(lat1) * sin(lat2))
  lon2 <- (lon2 + pi)%%(2 * pi) - pi
  cbind(lon2, lat2)/toRad
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
## Rescale function  #############################################################################
## Sourced from 'scales' package (https://CRAN.R-project.org/package=scales)

#' Rescale
#'
#' @description Rescales a numeric vector to a specified range.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.rescale <- function (x, to=c(0, 1), from=range(x, na.rm=TRUE, finite=TRUE), ...) {
  if (.zero_range(from) || .zero_range(to)) return(ifelse(is.na(x), NA, mean(to)))
  (x - from[1])/diff(from) * diff(to) + to[1]
}

.zero_range <- function(x, tol=1000*.Machine$double.eps) {
  if (length(x) == 1) return(TRUE)
  if (length(x) != 2) stop("'x' must be length 1 or 2", call.=FALSE)
  if (any(is.na(x))) return(NA)
  if (x[1] == x[2]) return(TRUE)
  if (all(is.infinite(x))) return(FALSE)
  m <- min(abs(x))
  if (m == 0) return(FALSE)
  abs((x[1] - x[2]) / m) < tol
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
      if(expand.legend)  layout_mat[nrow(layout_mat), 1:n_cols] <- current_plot
      else layout_mat[nrow(layout_mat), 1:min.legend.plots] <- current_plot
      #layout_mat <- rbind(layout_mat, rep(current_plot, n_cols))
      row_heights <- c(row_heights, plots.height)
    }
  }

  # return the layout matrix and row heights
  list(matrix=layout_mat, heights=row_heights, group_positions=group_positions)
}


##################################################################################################
## Projection check function #####################################################################

#' Check Projection
#'
#' @description Determines whether the supplied spatial object is projected or unprojected (geographic).
#' It checks the CRS of the spatial object and determines its projection status.
#' If the CRS is not defined, it verifies if the spatial object contains valid geographic coordinates.
#' @param spatial.object An 'sf' or 'Raster' object.
#' @return A string indicating whether the spatial object is "projected", "unprojected", or "no CRS".
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.checkProjection <- function(spatial.object) {

  # check if the input is an 'sf' object
  if (inherits(spatial.object, "sf")) {

    # extract the CRS
    coords_crs <- sf::st_crs(spatial.object)

    # check if CRS is NULL
    if (is.na(coords_crs)) {
      # get coordinates and check if they fall within valid geographic ranges
      coords <- sf::st_coordinates(spatial.object)
      lon_is_geographic <- all(coords[,1] >= -180 & coords[,1] <= 180, na.rm = TRUE)
      lat_is_geographic <- all(coords[,2] >= -90 & coords[,2] <= 90, na.rm = TRUE)
      # if both longitude and latitude fall within these ranges, it's unprojected
      if (lon_is_geographic && lat_is_geographic) return("geographic")
      else return("projected")
    } else if (sf::st_is_longlat(spatial.object)){
      return("geographic")
    } else {
      return("projected")
    }
  }

  # check if the input is a 'Raster' object
  else if (inherits(spatial.object, "Raster")) {

    # extract the CRS
    coords_crs <- raster::crs(spatial.object)

    # verify if the CRS is undefined (NULL)
    if (is.na(coords_crs)) {
      stop("The Raster layer does not contain CRS information.", call.=FALSE)
      # check if the CRS corresponds to geographic coordinates (WGS84)
    } else if (grepl("+proj=longlat +datum=WGS84", coords_crs@projargs, fixed=TRUE)) {
      return("geographic")
      # if the CRS does not match geographic, it's considered projected
    } else {
      return("projected")
    }
  }

  # If neither an 'sf' nor 'Raster' object
  stop("Input must be an 'sf' or 'Raster' object.", call.=FALSE)
}


##################################################################################################
## CRS management function #######################################################################

#' Manage Spatial Coordinate Reference Systems
#'
#' @description This function manages the coordinate reference systems (CRS) for spatial coordinates
#' and an optional spatial layer. It determines whether the provided coordinates and spatial
#' object are projected or geographic and processes them accordingly, based on a supplied EPSG code.
#'
#' @param coords A spatial object of class `sf`, containing longitude and latitude values.
#' @param spatial.layer An optional spatial object of class `sf` or `RasterLayer`.
#' @param epsg.code An optional numeric EPSG code for the desired coordinate reference system.
#'
#' @return A list containing:
#' \item{coords}{The transformed coordinates in the specified CRS.}
#' \item{spatial.layer}{The transformed spatial layer (if provided).}
#' \item{epsg}{The EPSG code used for transformations.}
#'
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd


.processSpatial <- function(coords, spatial.layer, epsg.code) {

  ##############################################################################
  # Determine the coordinate reference system of the provided spatial objects ##

  coords_crs <- .checkProjection(coords)
  if(!is.null(spatial.layer)) layer_crs <- .checkProjection(spatial.layer)
  epsg_supplied <- !is.null(epsg.code)
  if(inherits(spatial.layer, "sf")) layer_epsg <- sf::st_crs(spatial.layer)
  else if(inherits(spatial.layer, "Raster")) layer_epsg <- sf::st_crs(raster::crs(spatial.layer)@projargs)

  # check if epsg.code is a projected 'crs' object or an integer
  if(epsg_supplied){
    if(inherits(epsg.code, "numeric")) {
      epsg.code <- sf::st_crs(epsg.code)
    }else if(!inherits(epsg.code, "crs")) {
      stop("Invalid EPSG code format. Must be either a numeric value or a valid 'crs' object.", call.=FALSE)
    }
    if(!is.na(epsg.code$epsg) && epsg.code$epsg==4326){
      stop("Invalid EPSG code. The supplied code corresponds to WGS84 (EPSG:4326), a geographic coordinate system. Please provide a projected coordinate system instead.", call. = FALSE)
    }
  }

  # if the 'spatial.layer' variable is not NULL
  if(!is.null(spatial.layer)){
    # check whether 'spatial.layer' is either an 'sf' object or a Raster object
    if(!inherits(spatial.layer, c("sf","Raster"))) stop("Spatial.layer must be an 'sf' or 'Raster' object.", call.=FALSE)
    # if 'spatial.layer' is an 'sf' object, remove all attributes
    if(inherits(spatial.layer, "sf")) spatial.layer[] <- list(geometry=sf::st_geometry(spatial.layer))
  }


  # initialize string to return warning message
  warning_message <- c()


  ##############################################################################
  # Handle cases when spatial.layer is not provided ###############################

  if(is.null(spatial.layer)){

    # 1. Coordinates (geographic), EPSG (missing)
    if (coords_crs=="geographic" && !epsg_supplied) {
      stop("Longitudes/latitudes seem to be in a geographic (unprojected) format. Please either project them to a suitable CRS or provide an EPSG code for proper processing.", call.=FALSE)

      # 2. Coordinates (geographic) and EPSG (supplied)
    } else if (coords_crs=="geographic" && epsg_supplied) {
      sf::st_crs(coords) <- 4326
      coords <- sf::st_transform(coords, epsg.code)

      # 3. Coordinates (projected), EPSG (missing)
    } else if (coords_crs=="projected" && !epsg_supplied) {
      stop("Longitudes/latitude values seem to be projected but no 'epsg.code' has been supplied. Please provide the corresponding EPSG using the 'epsg.code' argument.", call.=FALSE)

      # 4. Coordinates (projected), EPSG (supplied)
    } else if (coords_crs=="projected" && epsg_supplied) {
      sf::st_crs(coords) <- epsg.code
    }

    ##############################################################################
    # Handle cases when spatial.layer is provided ###################################

  } else {

    # 1. Coordinates (geographic), spatial.layer (geographic), EPSG (missing)
    if (coords_crs=="geographic" && layer_crs=="geographic" && !epsg_supplied) {
      stop("Both spatial.layer and longitudes/latitudes seem to be in a geographic (unprojected) format. Please either project them to a suitable CRS or provide an EPSG code for proper processing.", call.=FALSE)

      # 2. Coordinates (geographic), spatial.layer (geographic), EPSG (supplied)
    } else if (coords_crs=="geographic" && layer_crs=="geographic" && epsg_supplied) {
      sf::st_crs(coords) <- 4326
      coords <- sf::st_transform(coords, epsg.code)
      if(inherits(spatial.layer, "sf")) spatial.layer <- sf::st_transform(spatial.layer, epsg.code)
      if(inherits(spatial.layer, "Raster")) spatial.layer <- raster::projectRaster(spatial.layer, crs=sf::st_crs(epsg.code$epsg)$proj4string, method="ngb")

      # 3. Coordinates (geographic), spatial.layer (projected), EPSG (missing)
    } else if (coords_crs=="geographic" && layer_crs=="projected" && !epsg_supplied) {
      sf::st_crs(coords) <- 4326
      epsg.code <- layer_epsg
      coords <- sf::st_transform(coords, epsg.code)
      if(!is.na(epsg.code$epsg)){
        warning_message <- paste0("No EPSG code supplied. Coordinates projected assuming CRS projection with EPSG:", epsg.code$epsg,
                                  " based on the provided spatial.layer.")
      }else{
        warning_message <- "No EPSG code supplied. Coordinates projected assuming the CRS projection from the provided spatial.layer."
      }

      # 4. Coordinates (geographic), spatial.layer (projected), EPSG (supplied)
    } else if (coords_crs=="geographic" && layer_crs=="projected" && epsg_supplied) {
      if (layer_epsg!=epsg.code) {
        if(inherits(spatial.layer, "sf")) spatial.layer <- sf::st_transform(spatial.layer, epsg.code)
        if(inherits(spatial.layer, "Raster")) spatial.layer <- raster::projectRaster(spatial.layer, crs=sf::st_crs(epsg.code$epsg)$proj4string, method="ngb")
        warning_message <- paste("The spatial layer has been projected to the supplied EPSG code:", epsg.code$epsg)
      }
      sf::st_crs(coords) <- 4326
      coords <- sf::st_transform(coords, epsg.code)

      # 5. Coordinates (projected), spatial.layer (geographic), EPSG (missing)
    } else if (coords_crs=="projected" && layer_crs=="geographic" && !epsg_supplied) {
      stop("Longitudes/latitude values seem to be projected but no EPSG code has been supplied. Please provide the corresponding epsg.code or supply longitude and latitude in a geographic CRS / unprojected format (WGS84).", call.=FALSE)

      # 6. Coordinates (projected), spatial.layer (geographic), EPSG (supplied)
    } else if (coords_crs=="projected" && layer_crs=="geographic" && epsg_supplied) {
      sf::st_crs(coords) <- epsg.code
      if(inherits(spatial.layer, "sf")) spatial.layer <- sf::st_transform(spatial.layer, epsg.code)
      if(inherits(spatial.layer, "Raster")) spatial.layer <- raster::projectRaster(spatial.layer, crs=epsg.code$wkt, method="ngb")

      # 7. Coordinates (projected), spatial.layer (projected), EPSG (missing)
    } else if (coords_crs=="projected" && layer_crs=="projected" && !epsg_supplied) {
      epsg.code <- layer_epsg
      sf::st_crs(coords) <- epsg.code
      if(!is.na(epsg.code$epsg)){
        warning_message <- paste0("No EPSG code supplied. Assuming CRS projection with EPSG:", epsg.code$epsg,
                                  " based on the provided spatial.layer.")
      }else{
        warning_message <- "No EPSG code supplied. Assuming CRS projection from the provided spatial.layer."
      }

      # 8. Coordinates (projected), spatial.layer (projected), EPSG (supplied)
    } else if (coords_crs=="projected" && layer_crs=="projected" && epsg_supplied) {
      if (layer_epsg!=epsg.code) {
        if(inherits(spatial.layer, "sf")) spatial.layer <- sf::st_transform(spatial.layer, epsg.code)
        if(inherits(spatial.layer, "Raster")) spatial.layer <- raster::projectRaster(spatial.layer, crs=epsg.code$wkt, method="ngb")
        warning_message <- paste("The spatial layer has been reprojected to the supplied EPSG code:", epsg.code$epsg)
      }
      sf::st_crs(coords) <- epsg.code
    }
  }

  ##############################################################################
  # Print warnings #############################################################

  if (length(warning_message)>0){
    warning_message <- sapply(warning_message, function(x) paste("-", x))
    warning_message <- sapply(warning_message, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    sapply(warning_message, function(x) warning(x, call.=FALSE))
  }



  ##############################################################################
  # Return transformed coordinates and spatial.layer (if applicable) ###########

  return(list(coords=coords, spatial.layer=spatial.layer, epsg.code=epsg.code))
}



##################################################################################################
## Decimal Places   ##############################################################################
## Sourced from https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r

#' Decimal Places
#'
#' @description Determines the number of decimal places in a numeric value.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.decimalPlaces <- function(x) {
  if(is.na(x)){return(NA)}
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', format(x, scientific=FALSE)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
.decimalPlaces <- Vectorize(.decimalPlaces)


##################################################################################################
## Economist color palette #######################################################################
## Based in ggthemes::economist_pal()

#' Economist color palette
#'
#' @description This function generates a color palette using the Economist color scheme,
#' inspired by the color palette provided in the ggthemes::economist_pal() function.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.economist_pal <- function(n){
  economist_colors <- c("#6794a7", "#014d64", "#01a2d9", "#7ad2f6","#00887d",
                        "#76c0c1", "#7c260b", "#ee8f71", "#a18376", "#adadad")
  if(n==1) economist_colors[2]
  else if(n==2) return(economist_colors[c(3,2)])
  else if(n==3) return(economist_colors[c(1,2,3)])
  else if(n==4) return(economist_colors[c(1,2,3,9)])
  else if(n==5) return(economist_colors[c(1,2,4,3,6)])
  else if(n==6) return(economist_colors[c(1,2,4,3,6,5)])
  else if(n==7) return(economist_colors[c(1:6,9)])
  else if(n==8) return(economist_colors[c(1:6,8:9)])
  else if(n==9) return(economist_colors[1:9])
  else if(n==10) return(economist_colors[1:10])
  else stop("Invalid number of colors requested for the economist palette. Please choose a number between 1 and 10.")
}


##################################################################################################
## Viridis color palette #########################################################################
## Sourced from viridis::viridis(100)

#' Viridis color palette
#'
#' @description This function generates a color palette using the Viridis color scheme,
#' known for its perceptually uniform properties. The palette is adapted from the viridis package.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.viridis_pal <- function(n) {
  viridis_colors <- c(
    "#440154FF", "#450558FF", "#46085CFF", "#470D60FF", "#471063FF",
    "#481467FF", "#481769FF", "#481B6DFF", "#481E70FF", "#482173FF",
    "#482576FF", "#482878FF", "#472C7AFF", "#472F7CFF", "#46327EFF",
    "#453581FF", "#453882FF", "#443B84FF", "#433E85FF", "#424186FF",
    "#404587FF", "#3F4788FF", "#3E4A89FF", "#3D4D8AFF", "#3C508BFF",
    "#3B528BFF", "#39558CFF", "#38598CFF", "#375B8DFF", "#355E8DFF",
    "#34608DFF", "#33638DFF", "#32658EFF", "#31688EFF", "#2F6B8EFF",
    "#2E6D8EFF", "#2D708EFF", "#2C718EFF", "#2B748EFF", "#2A768EFF",
    "#29798EFF", "#287C8EFF", "#277E8EFF", "#26818EFF", "#26828EFF",
    "#25858EFF", "#24878EFF", "#238A8DFF", "#228D8DFF", "#218F8DFF",
    "#20928CFF", "#20938CFF", "#1F968BFF", "#1F998AFF", "#1E9B8AFF",
    "#1F9E89FF", "#1FA088FF", "#1FA287FF", "#20A486FF", "#22A785FF",
    "#24AA83FF", "#25AC82FF", "#28AE80FF", "#2BB07FFF", "#2EB37CFF",
    "#31B67BFF", "#35B779FF", "#39BA76FF", "#3DBC74FF", "#41BE71FF",
    "#47C06FFF", "#4CC26CFF", "#51C56AFF", "#56C667FF", "#5BC863FF",
    "#61CA60FF", "#67CC5CFF", "#6DCD59FF", "#73D056FF", "#78D152FF",
    "#7FD34EFF", "#85D54AFF", "#8CD646FF", "#92D741FF", "#99D83DFF",
    "#A0DA39FF", "#A7DB35FF", "#ADDC30FF", "#B4DE2CFF", "#BBDE28FF",
    "#C2DF23FF", "#C9E020FF", "#D0E11CFF", "#D7E219FF", "#DDE318FF",
    "#E4E419FF", "#EBE51AFF", "#F1E51DFF", "#F7E620FF", "#FDE725FF"
  )
  return(colorRampPalette(viridis_colors)(n))
}



##################################################################################################
## Palr bathy deep color palette #########################################################################
## Sourced from palr::bathy_deep_pal

#' Bathy deep color palette
#'
#' @description This function generates a color palette using the Bathy Deep color scheme,
#' sourced from palr::bathy_deep_pal function.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.bathy_deep_pal <- function (x, palette = FALSE, alpha = 1) {
  breaks <- c(-5500, seq(-5000, -1000, by = 1000), -500, 0)
  breaks <- seq(-5500, 0, length = 255)
  cols <- colorRampPalette(rgb(c(0,18,60,103,141,194,255), c(0,18,60,111,163,216,255),
                               c(0,26,85,135,173,216,255), maxColorValue=255))(256)
  hexalpha <- as.hexmode(round(255 * alpha))
  if (nchar(hexalpha) == 1L) hexalpha <- paste(rep(hexalpha, 2L), collapse = "")
  cols <- paste0(cols, hexalpha)
  if (palette)  return(list(breaks = breaks, cols = cols))
  if (missing(x)) return(colorRampPalette(cols))
  if (length(x) == 1L) {return(paste0(colorRampPalette(cols)(x), hexalpha))}
  else {return(cols[findInterval(x, breaks)])}
}



##################################################################################################
## Jet color palette #############################################################################
## Sourced from pals::jet

#' Jet color palette
#'
#' @description This function generates a color palette using the Jet color scheme,
#' sourced from palr::jet function.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.jet_pal <- function(n=25){
  colors <- c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F", "#ffff00",
              "#FF7F00", "#ff0000", "#7F0000")
  colorRampPalette(colors)(n)
}





##################################################################################################
##################################################################################################
##################################################################################################
