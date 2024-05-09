
#######################################################################################################
# Helper functions ####################################################################################
#######################################################################################################

##################################################################################################
## Updated filled.contour function  ##############################################################
## adapted from https://gist.github.com/epijim/6514388 ###########################################

filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
            col = color.palette(length(levels) - 1), plot.title, plot.axes,
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
            axes = TRUE, frame.plot = axes, mar, invert.scale=FALSE, ...)
  {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page

    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2]) * par("csi") * 2.54
    # par(las = las)
    # mar <- mar.orig
    plot.new()

    if(invert.scale){
      zscale <- rev(range(levels))
    }else{
      zscale <- range(levels)
    }

    # par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1)
      stop("no proper 'z' matrix specified")
    if (!is.double(z))
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
                    col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot)
      box()
    if (missing(plot.title))
      title(...)
    else plot.title
    invisible()
  }



##################################################################################################
## Updated shape::colorlegend function  ##########################################################
## adapted from https://rdrr.io/cran/shape/src/R/colorlegend.R  ##################################
## added main.adj + support for scientific notation


colorlegend2 <- function(col=femmecol(100), zlim, zlevels=5,
                         dz=NULL, zval=NULL, log=FALSE, posx=c(0.9,0.93), posy=c(0.05,0.9),
                         main=NULL, main.cex=1.0, main.col="black", main.adj=0.5, lab.col="black",
                         digit=0, left=FALSE, lab.scientific=FALSE, ...) {

  ncol   <- length (col)
  par (new=TRUE)
  omar <- nmar <- par("mar")
  nmar[c(2,4)]<-0
  par (mar = nmar)

  shape::emptyplot()
  pars   <- par("usr")

  ## Rectangle positions on x and y-axis
  dx     <- pars[2]-pars[1]
  xmin   <- pars[1]+posx[1]*dx
  xmax   <- pars[1]+posx[2]*dx

  dy     <- pars[4]-pars[3]
  ymin   <- pars[3]+posy[1]*dy
  ymax   <- pars[3]+posy[2]*dy

  ## z-values
  if (!is.null(zval)) {
    zz<-zval
    dz<-NULL
  }

  if (is.null(dz)&is.null(zval))
    if (! is.null(zlevels)) {
      if (log) {
        zz <- 10^(pretty(log10(zlim),n=(zlevels+1)))
      } else
        zz <-     pretty(zlim,n=(zlevels+1))
    } else zz <- NULL
  if (!is.null(dz)) {
    if (log)
      zz <- 10^(seq(log10(zlim[1]),log10(zlim[2]),by=dz))
    if (!log)
      zz <- seq(zlim[1],zlim[2],by=dz)
  }

  if (log) {
    zlim <- log10(zlim)
    if (! is.null(zz))
      zz   <- log10(zz)
  }

  zmin   <- zlim[1]
  zmax   <- zlim[2]

  ## colors
  Y <- seq(ymin,ymax,length.out=ncol+1)
  rect(xmin,Y[-(ncol+1)],xmax,Y[-1],col=col,border=NA)
  rect(xmin,ymin,xmax,ymax,border=lab.col)

  if (! is.null(zz)) {
    ## labels
    dx     <- (xmax-xmin)
    dy     <- (ymax-ymin)

    if (left) {
      Dx  <-  -dx  # labels on left..
      pos <-   2
      xpos <- xmin+Dx*0.5
    } else {
      Dx  <- +dx  # labels on right..
      pos <- 4
      xpos <- xmax+Dx*0.5
    }

    Ypos <- ymin+(zz-zmin)/(zmax-zmin)*dy
    segments(xmin,Ypos,xmax,Ypos,col=lab.col)
    segments(xpos+Dx*0.25,Ypos,xmin,Ypos,col=lab.col)
    if(lab.scientific==T){
      labels <- format(zz, scientific=TRUE)
    }else{
      labels <- formatC(zz,digits=digit,format="f")
    }
    text (xpos,Ypos,labels,pos=pos,col=lab.col,...)
  }

  if  (!is.null(main)) {
    for (i in length(main):1)
      if(main.adj==0){main_pos<-xmin}else if(main.adj==1){main_pos<-xmax}else{main_pos<-mean(c(xmin,xmax))}
    text (x=main_pos,y=ymax+0.05*(length(main)-i+1),
          labels=main[i], adj=c(main.adj, 0.5), cex=main.cex, col=main.col)
  }
  par (new=FALSE)
  par (mar=omar)

}



##################################################################################################
# Auxiliary function I - get plot coordinates by keyword #########################################
# adapted from legend function ###################################################################

getPosition <- function(keyword, inset) {
  usr <- par("usr")
  inset <- rep_len(inset, 2)
  insetx <- inset[1] * (usr[2] - usr[1]) * 2
  left <- switch(keyword, bottomright=, topright=, right=usr[2]-insetx,
                 bottomleft=, left=, topleft=usr[1]+insetx, bottom=,
                 top=, center=(usr[1] + usr[2L])/2)
  insety <- inset[2] * (usr[4] - usr[3])
  top <- switch(keyword, bottomright=, bottom=, bottomleft=usr[3]+insety,
                topleft=, top=, topright=usr[4]-insety, left=, right=,
                center=(usr[3]+usr[4])/2)
  return(c(left, top))
}


##################################################################################################
## Updated raster::scalebar function  ############################################################
## added new 'bar.lwd' argument to set the border line width of the scalebars and
## added new 'bar.height' argument to control the scalebar thickness

scalebar2 <- function (d, xy=NULL, type="line", divs=2, below="", bar.lwd=0.4, bar.height=1.5,
                      lonlat=NULL, label, adj=c(0.5, -0.5), lwd=2, ...)
{
  stopifnot(type %in% c("line", "bar"))
  pr <- graphics::par()
  if (is.null(lonlat)) {
    if (pr$usr[1] > -181 & pr$usr[2] < 181 & pr$yaxp[1] >
        -200 & pr$yaxp[2] < 200) {
      lonlat <- TRUE
    }
    else {
      lonlat <- FALSE
    }
  }
  if (lonlat) {
    lat <- mean(pr$yaxp[1:2])
    if (missing(d)) {
      dx <- (pr$usr[2] - pr$usr[1])/10
      d <- pointDistance(cbind(0, lat), cbind(dx, lat),
                         TRUE)
      d <- signif(d/1000, 2)
      label <- NULL
    }
    p <- cbind(0, lat)
    dd <- .destPoint(p, d * 1000)
    dd <- dd[1, 1]
  }
  else {
    if (missing(d)) {
      d <- round(10 * (pr$usr[2] - pr$usr[1])/10)/10
      label <- NULL
    }
    dd <- d
  }
  if (is.null(xy)) {
    padding = c(5, 5)/100
    parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
    xy <- c(pr$usr[1] + (padding[1] * parrange[1]), pr$usr[3] +
              (padding[2] * parrange[2]))
  }
  if (type == "line") {
    lines(matrix(c(xy[1], xy[2], xy[1] + dd, xy[2]), byrow = T,
                 nrow = 2), lwd = lwd, ...)
    if (missing(label)) {
      label <- paste(d)
    }
    if (is.null(label)) {
      label <- paste(d)
    }
    if (missing(adj)) {
      adj <- c(0.5, -0.2 - lwd/20)
    }
    text(xy[1] + (0.5 * dd), xy[2], labels = label, adj = adj,
         ...)
  }
  else if (type == "bar") {
    stopifnot(divs > 0)
    if (missing(adj)) {
      adj <- c(0.5, -1)
    }

    lwd <- dd/25 * bar.height
    if (divs == 2) {
      half <- xy[1] + dd/2
      graphics::polygon(c(xy[1], xy[1], half, half), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]),
                        col = "white", lwd=bar.lwd)
      graphics::polygon(c(half, half, xy[1] + dd, xy[1] + dd), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]),
                        col = "black", lwd=bar.lwd)
      if (missing(label)) {
        label <- c("0", "", d)
      }
      if (is.null(label)) {
        label <- c("0", "", d)
      }
      text(xy[1], xy[2], labels = label[1], adj = adj,
           ...)
      text(xy[1] + 0.5 * dd, xy[2], labels = label[2],
           adj = adj, ...)
      text(xy[1] + dd, xy[2], labels = label[3], adj = adj,
           ...)
    }
    else {
      q1 <- xy[1] + dd/4
      half <- xy[1] + dd/2
      q3 <- xy[1] + 3 * dd/4
      end <- xy[1] + dd
      graphics::polygon(c(xy[1], xy[1], q1, q1), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]),
                        col = "white", lwd=bar.lwd)
      graphics::polygon(c(q1, q1, half, half), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]),
                        col = "black", lwd=bar.lwd)
      graphics::polygon(c(half, half, q3, q3), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]),
                        col = "white", lwd=bar.lwd)
      graphics::polygon(c(q3, q3, end, end), c(xy[2], xy[2] + lwd, xy[2] + lwd, xy[2]),
                        col = "black", lwd=bar.lwd)
      if (missing(label)) {
        label <- c("0", round(0.5 * d), d)
      }
      if (is.null(label)) {
        label <- c("0", round(0.5 * d), d)
      }
      text(xy[1], xy[2], labels = label[1], adj = adj,
           ...)
      text(half, xy[2], labels = label[2], adj = adj, ...)
      text(end, xy[2], labels = label[3], adj = adj, ...)
    }
    if (below != "") {
      adj[2] <- -adj[2]
      text(xy[1] + (0.5 * dd), xy[2], labels = below, adj = adj,
           ...)
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
## Updated legend function  ######################################################################
## Adapted by Ben Nolker (http://ms.mcmaster.ca/~bolker/R/misc/legendx.R)

legend2 <- function(x, y = NULL, legend, fill = NULL, col = par("col"), border="black",
           lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"),
           box.lwd = par("lwd"), box.lty = par("lty"), box.col = par("fg"),
           box.cex = c(0.8,0.5),
           pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd,
           xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5),
           text.width = NULL, text.col = par("col"), text.font = NULL,
           merge = do.lines && has.pch, trace = FALSE,
           plot = TRUE, ncol = 1, horiz = FALSE, title = NULL,
           inset = 0, xpd, title.col = text.col, title.adj = 0.5,
           seg.len = 2)
  {
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
# Helper function I - return number of decimal places ############################################
# Retrieved from https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r

decimalPlaces <- function(x) {
  if(is.na(x)){return(NA)}
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
decimalPlaces <- Vectorize(decimalPlaces)


##################################################################################################
##################################################################################################
##################################################################################################