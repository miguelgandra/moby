#######################################################################################################
## Internal helpers: spatial #########################################################################
## CRS detection/projection handling, scale-bar drawing and geodesic helpers.
#######################################################################################################


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
  if (!type %in% c("line", "bar")) .mobyAbort("'type' must be either 'line' or 'bar'.")

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
      d <- geosphere::distGeo(cbind(0, lat), cbind(dx, lat))
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
    if (divs <= 0) .mobyAbort("'divs' must be a positive integer.")
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
## Projection check function #####################################################################

#' Check Projection
#'
#' @description Determines whether the supplied spatial object is projected or unprojected (geographic).
#' It checks the CRS of the spatial object and determines its projection status.
#' If the CRS is not defined, it verifies if the spatial object contains valid geographic coordinates.
#' @param spatial.object An 'sf' or 'Raster' object.
#' @return A string indicating whether the spatial object is "geographic" or "projected".
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

  # check if the input is a terra 'SpatRaster'
  else if (inherits(spatial.object, "SpatRaster")) {

    # extract the CRS (terra::crs returns a WKT string; "" when undefined - sf::st_crs("") errors,
    # so guard the empty string before handing it to sf)
    crs_wkt <- terra::crs(spatial.object)

    # verify if the CRS is undefined
    if (is.na(crs_wkt) || !nzchar(crs_wkt)) {
      stop("The raster does not contain CRS information.", call.=FALSE)
      # check if the CRS corresponds to geographic coordinates (WGS84)
    } else if (isTRUE(sf::st_is_longlat(sf::st_crs(crs_wkt)))) {
      return("geographic")
      # if the CRS is not geographic (any datum, not just WGS84), it's projected
    } else {
      return("projected")
    }
  }

  # If neither an 'sf' nor 'SpatRaster' object
  stop("Input must be an 'sf' or 'SpatRaster' object.", call.=FALSE)
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
  else if(inherits(spatial.layer, "SpatRaster")) layer_epsg <- sf::st_crs(terra::crs(spatial.layer))

  # check if epsg.code is a projected 'crs' object or an integer
  if(epsg_supplied){
    if(inherits(epsg.code, "numeric")) {
      epsg.code <- sf::st_crs(epsg.code)
    }else if(!inherits(epsg.code, "crs")) {
      stop("Invalid EPSG code format. Must be either a numeric value or a valid 'crs' object.", call.=FALSE)
    }
    if(!is.na(epsg.code$epsg) && epsg.code$epsg==4326){
      stop("Invalid EPSG code. The supplied code corresponds to WGS84 (EPSG:4326), a geographic coordinate system. Please provide a projected coordinate system instead.", call. = FALSE)
    }else if(sf::st_is_longlat(epsg.code)){
      stop("Invalid EPSG code. The supplied code is a geographic coordinate system. Please provide a projected coordinate system with meter units.", call. = FALSE)
    }
  }

  # if the 'spatial.layer' variable is not NULL
  if(!is.null(spatial.layer)){
    # check whether 'spatial.layer' is either an 'sf' object or a Raster object
    if(!inherits(spatial.layer, c("sf","SpatRaster"))) stop("Spatial.layer must be an 'sf' or 'SpatRaster' object.", call.=FALSE)
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
      if(inherits(spatial.layer, "SpatRaster")) spatial.layer <- terra::project(spatial.layer, epsg.code$wkt, method="near")

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
        if(inherits(spatial.layer, "SpatRaster")) spatial.layer <- terra::project(spatial.layer, epsg.code$wkt, method="near")
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
      if(inherits(spatial.layer, "SpatRaster")) spatial.layer <- terra::project(spatial.layer, epsg.code$wkt, method="near")

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
        if(inherits(spatial.layer, "SpatRaster")) spatial.layer <- terra::project(spatial.layer, epsg.code$wkt, method="near")
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
## Shared basemap helpers ########################################################################
## Extracted from plotMaps and plotMovements so the two share ONE spatial pipeline (CRS handling,
## bounding box, canvas, land overlay and scale bar) instead of diverging.

#' Reconcile the CRS of coordinates, land shape AND background raster
#'
#' @description Extends `.processSpatial()` to also reproject the optional `background.layer` onto the
#' resolved map CRS. `.processSpatial` handles only coords + land.shape, which was the source of a
#' silent crash in `plotMaps` when a background raster was supplied in a different CRS. Reprojection
#' uses the WKT string (not the deprecated/lossy proj4string) so datum differences are respected.
#' @return list(coords, land.shape, background.layer, epsg.code) all on a common CRS.
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd
.prepareMapLayers <- function(coords, land.shape, background.layer, epsg.code) {
  sp <- .processSpatial(coords, land.shape, epsg.code)
  coords <- sp$coords; land.shape <- sp$spatial.layer; epsg.code <- sp$epsg.code
  target <- if (inherits(epsg.code, "crs")) epsg.code else sf::st_crs(epsg.code)

  if (!is.null(background.layer)) {
    # terra::crs returns "" (not NA) for a CRS-less raster, and sf::st_crs("") errors - guard it
    bg_wkt <- terra::crs(background.layer)
    if (is.na(bg_wkt) || !nzchar(bg_wkt)) {
      warning("- The 'background.layer' has no CRS; assuming it matches the map projection.", call. = FALSE)
    } else if (!isTRUE(sf::st_crs(bg_wkt) == target)) {
      # categorical layers must use nearest-neighbour; continuous layers bilinear
      method <- if (isTRUE(terra::is.factor(background.layer)[1])) "near" else "bilinear"
      background.layer <- terra::project(background.layer, target$wkt, method = method)
      warning("- The 'background.layer' was reprojected to match the map's coordinate system.", call. = FALSE)
    }
  }
  list(coords = coords, land.shape = land.shape, background.layer = background.layer, epsg.code = epsg.code)
}

#' Expand an sf bounding box by a factor in all directions (preserves names)
#' @keywords internal
#' @noRd
.expandBbox <- function(bbox, extent.factor) {
  dx <- bbox["xmax"] - bbox["xmin"]; dy <- bbox["ymax"] - bbox["ymin"]
  bbox["xmin"] <- bbox["xmin"] - dx * (extent.factor - 1)
  bbox["xmax"] <- bbox["xmax"] + dx * (extent.factor - 1)
  bbox["ymin"] <- bbox["ymin"] - dy * (extent.factor - 1)
  bbox["ymax"] <- bbox["ymax"] + dy * (extent.factor - 1)
  bbox
}

#' Convert an sf-order bbox (xmin,ymin,xmax,ymax) to a terra SpatExtent in the CORRECT axis order
#'
#' @description `terra::ext()`/`terra::crop()` read a plain numeric POSITIONALLY as
#' (xmin, xmax, ymin, ymax) and ignore names, so passing an sf-order bbox scrambles the axes
#' (mis-crop, or an "xmin >= xmax" crash for near-equatorial/local grids). This centralises the
#' explicit reorder so every raster call is axis-correct.
#' @keywords internal
#' @noRd
.bboxToExtent <- function(bbox) {
  terra::ext(as.numeric(bbox[c("xmin", "xmax", "ymin", "ymax")]))
}

#' Open an empty, square (asp = 1) map canvas spanning a bounding box
#' @keywords internal
#' @noRd
.newMapCanvas <- function(bbox) {
  graphics::plot(x = bbox[c("xmin", "xmax")], y = bbox[c("ymin", "ymax")], type = "n",
                 xlab = "", ylab = "", main = "", xlim = bbox[c("xmin", "xmax")], ylim = bbox[c("ymin", "ymax")],
                 axes = FALSE, asp = 1, xaxs = "i", yaxs = "i")
}

#' Overlay a coastline/land polygon on the current map
#' @keywords internal
#' @noRd
.drawLandOverlay <- function(land.shape, land.color) {
  if (!is.null(land.shape)) graphics::plot(sf::st_geometry(land.shape), col = land.color, border = NA, add = TRUE)
}

#' Fetch a default coastline for the study area when the user supplies no `land.shape`
#'
#' @description Returns land polygons (an `sfc` in the map CRS) cropped to the projected `bbox`, so the
#' plotting functions can draw a coastline without the user having to supply one. Sources, in order of
#' preference: the (Suggests) `rnaturalearth` package (better resolution) then `maps`; if neither is
#' installed it warns once per session and returns `NULL` (maps are then drawn without land, as before).
#' `coastline` may be `TRUE` (auto-pick a scale from the extent) or an explicit rnaturalearth scale
#' (`"small"`/`"medium"`/`"large"`). Regional windows entirely over water return `NULL`.
#' @keywords internal
#' @noRd
.defaultCoastline <- function(bbox, epsg.code, coastline = TRUE, verbose = TRUE) {
  # geographic (EPSG:4326) window matching the projected study area, to crop a global coastline to it
  win <- tryCatch({
    poly <- sf::st_as_sfc(sf::st_bbox(c(xmin = bbox[["xmin"]], ymin = bbox[["ymin"]],
                                        xmax = bbox[["xmax"]], ymax = bbox[["ymax"]]),
                                      crs = sf::st_crs(epsg.code)))
    sf::st_bbox(sf::st_transform(poly, 4326))
  }, error = function(e) NULL)
  if (is.null(win)) return(NULL)
  span  <- max(win[["xmax"]] - win[["xmin"]], win[["ymax"]] - win[["ymin"]])
  scale <- if (is.character(coastline)) coastline else if (span > 20) "small" else "medium"

  land <- NULL
  if (requireNamespace("rnaturalearth", quietly = TRUE) && requireNamespace("rnaturalearthdata", quietly = TRUE))
    land <- tryCatch(sf::st_geometry(rnaturalearth::ne_countries(scale = scale, returnclass = "sf")),
                     error = function(e) NULL)
  if (is.null(land) && requireNamespace("maps", quietly = TRUE))
    land <- tryCatch(sf::st_geometry(sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))),
                     error = function(e) NULL)
  if (is.null(land)) {
    if (is.null(mobyEnv$coastline_warned)) {
      .mobyInform("No 'land.shape' supplied and no coastline source is installed; drawing without land. ",
                  "Install 'rnaturalearth' (+ 'rnaturalearthdata') or 'maps' for an automatic coastline, ",
                  "or pass your own 'land.shape'.", verbose = verbose)
      mobyEnv$coastline_warned <- TRUE
    }
    return(NULL)
  }
  # normalise the source to WGS84 before cropping: 'maps' land is unset or Clarke-1866 longlat, which
  # would otherwise mismatch the 4326 window; rnaturalearth is already 4326 (transform is then a no-op).
  if (is.na(sf::st_crs(land))) land <- sf::st_set_crs(land, 4326)
  land <- tryCatch(suppressWarnings(sf::st_transform(land, 4326)), error = function(e) land)
  land <- tryCatch(suppressWarnings(sf::st_crop(sf::st_make_valid(land), win)), error = function(e) NULL)
  if (is.null(land) || length(land) == 0) return(NULL)   # study window entirely over water
  tryCatch(suppressWarnings(sf::st_transform(land, sf::st_crs(epsg.code))), error = function(e) NULL)
}

#' Draw a projected (metric) scale bar with a single, consistent placement convention
#'
#' @description Both spatial plots draw the same two-division km scale bar. `scale.km = NULL`
#' auto-sizes to ~20 percent of the map width. Coordinates are projected by construction (they pass
#' through `.prepareMapLayers()`), so `lonlat = FALSE` is forced - avoiding the geographic/projected
#' auto-detection misfire in `.scalebar()`. Returns the km length used (for the diagnostic summary).
#' @keywords internal
#' @noRd
.drawScaleBar <- function(bbox, scale.km = NULL, scale.pos = "bottomright", scale.inset = 0.05,
                          height = 1.5, cex = 0.7, color = "black") {
  if (is.null(scale.km)) scale.km <- pretty((bbox["xmax"] - bbox["xmin"]) * 0.2 / 1000)[2]
  scale.m <- scale.km * 1000
  xy <- .getPosition(scale.pos, inset = scale.inset, bar.width = scale.m)
  .scalebar(d = scale.m, xy = xy, type = "bar", divs = 2, below = "km", bar.height = height,
            label = c(0, scale.km / 2, scale.km), label.color = color, lonlat = FALSE,
            cex = cex, bar.lwd = 0.2, lwd = 0.2)
  invisible(as.numeric(scale.km))
}

##################################################################################################
##################################################################################################
##################################################################################################
