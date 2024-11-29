##################################################################################################
## Plot maps with movement trajectories + home ranges ############################################
##################################################################################################

#' Plot maps with movement trajectories and home range areas.
#'
#' @description This function generates maps that display inferred animal movement trajectories
#' (shortest in-water paths) alongside kernel density home range estimates. It can optionally
#' overlay additional spatial information such as bathymetry or sea surface temperature using a
#' background raster layer.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal positions (centers of activity, COAs) with longitude
#' and latitude coordinates.
#' @param kernel.densities Kernel utilization areas ('estUD' class), as returned by \code{\link{calculateKUDs}}.
#' @param animal.tracks The output of \code{\link{calculateTracks}}.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param land.shape Optional. A shapefile containing the coastlines or landmasses to be plotted.
#' @param epsg.code The EPSG code (integer) representing the coordinate reference system (CRS) to be used
#' for projecting the positions. If not specified, the function will attempt to use the CRS from the
#' provided land.shape (if available).
#' @param land.color Color used to represent land areas on the map. Default to "gray50".
#' @param background.layer Optional. A projected raster object representing an environmental
#' variable (e.g., bathymetry, temperature) to be displayed in the background.
#' @param background.pal Optional. Color palette for the background raster layer.
#' @param background.title Optional. The title or label for the background raster layer legend,
#' indicating the type of environmental variable being displayed (e.g., "Bathymetry" or "Temperature").
#' If not provided, the default title is "Layer values".
#' @param discard.missing Logical, indicating whether individuals with fewer than 5 detections
#' should be excluded from the plots. Defaults to TRUE. If set to FALSE, empty maps will be drawn
#' for these individuals, containing only the land.shape and background layer (if provided).
#' @param color.pal Optional. A color palette used to represent the kernel utilization distribution (KUD)
#' density values on the map. If NULL (default), If NULL (default), the \code{terrain.colors} palette is used.
#' @param kud.transparency.threshold Numeric threshold (between 0 and 1) for kernel utilization density (KUD) values.
#' Areas with KUD values below this threshold will be set to NA, effectively removing lower-density regions
#' (only relevant when a background layer is displayed). This threshold corresponds to a contour level;
#' for example, a value of 0.95 will exclude areas outside the 95% kernel utilization contour. Defaults to 0.95.
#' @param tracks.color Optional. The color used for the animal movement trajectories
#' on the map. Defaults to "black".
#' @param tracks.lty Optional. The line type for the animal movement trajectories, as per \code{lty}
#' in R graphics. Defaults to 2 (dashed line).
#' @param tracks.lwd Optional. The line width for the animal movement trajectories. Defaults to 0.2.
#' @param plot.detections Logical. If TRUE (default), individual detection points will be plotted on the map.
#' If FALSE, only trajectories and kernel density estimates (if provided) will be shown without displaying individual locations.
#' @param pts.color Optional. Specifies the color of the individual detection points on the map.
#' It can be a single color or a vector of three colors: the first color corresponds to the earliest detection,
#' the second color applies to all detections in between, and the third color represents the latest detection.
#' Defaults to c("blue", "white", "red").
#' @param pts.cex Optional. Specifies the size of the individual detection points. It can be a single
#' size value or a vector of three sizes: the first size applies to the earliest detection,
#' the second size is for detections in between, and the third size is for the latest detection.
#' Defaults to c(0.8, 0.7, 0.8).
#' @param pts.color Optional. Color for the individual detection points on the plot. Defaults to "white".
#' @param pts.cex Optional. Size of the individual detection points. Defaults to 0.7.
#' @param same.scale Logical, whether to enforce the same scale for kernel density (KUD) areas across
#' all individuals for comparison. When set to FALSE, each plot will display its own y-axis.
#' If TRUE, only the y-axis of the left-most plots will be shown to conserve space. Defaults to FALSE.
#' @param kud.legend Logical, whether to display a legend for the KUD density values. Defaults to TRUE.
#' @param legend.cex Size of the legend text for KUD density and background layer values (if available).
#' Defaults to 0.8.
#' @param title.color Optional. The color of the title (animal ID) displayed on each map. Defaults to "black".
#' @param title.pos Position of the title (animal ID) on the plot. This can be one of the following
#' predefined keywords: `"bottomright"`, `"bottomleft"`, `"topright"`, `"topleft"`, `"bottom"`,
#' `"top"`,  `"left"`, `"right"`, or `"center"`. Defaults to `"bottomright"`. Defaults to `"topleft"`.
#' @param title.inset Controls how far the title is placed from the plot edges.
#' Can be specified as a single value (applied to both x and y directions) or as a two-element vector
#' (with x and y inset values specified separately).
#' - **Positive values** move the scale bar *inward*, towards the center of the plot.
#' - **Negative values** move the scale bar *outward*, potentially placing it outside the visible plot area.
#'
#' Defaults to c(-0.08, 0).
#' @param title.cex Determines the size of the plot title (animal ID). Defaults to 1.5.
#' @param scale.meters The length of the scale bar in meters. If set to NULL (default), it is
#' automatically set to 20% of the plot width.
#' @param scale.color Optional. The color of the scale bar legend displayed on each map. Defaults to "black".
#' @param scale.pos The position of the scale bar on the map. This can be one of the following
#' predefined keywords: `"bottomright"`, `"bottomleft"`, `"topright"`, `"topleft"`, `"bottom"`,
#' `"top"`,  `"left"`, `"right"`, or `"center"`. Defaults to `"bottomright"`.
#' @param scale.inset Controls how far the scale bar is placed from the plot edges.
#' Can be specified as a single value (applied to both x and y directions) or as a two-element vector
#' (with x and y inset values specified separately).
#' - **Positive values** move the scale bar *inward*, towards the center of the plot.
#' - **Negative values** move the scale bar *outward*, potentially placing it outside the visible plot area.
#'
#' Defaults to `0.2`.
#' @param extent.factor Numeric. Factor by which to adjust the extent of the plotting region,
#' defined based on the bounding box around animal positions/detections. A value of 1 keeps
#' the original bounding box, values greater than 1 increase the extent, and values less
#' than 1 decrease it. Defaults to 1.1 (10% increase).
#' @param cols Number of columns in the plot panel layout (used in the 'mfrow' argument). Defaults to 3.
#' @export
#'
#' @details This function is designed to create multi-panel plots of animal movement and home ranges.
#' It is flexible enough to visualize multiple individuals or animal groups (e.g. species),
#' while also supporting the overlay of environmental layers such as bathymetry or sea surface temperature.
#' The function processes input data (animal positions, kernel utilization distributions, and optional
#' movement tracks) and outputs a grid of maps. Each map represents an individual's home range,
#' along with optional trajectories and background layers. When multiple individuals are included,
#' it allows the option of enforcing a common scale for kernel density areas across all maps for
#' comparability.
#'
#' @examples
#' \dontrun{
#' # Example usage (adjust PDF dimensions as needed)
#' pdf("./home-range-maps.pdf", width=12, height=20)
#' plotMaps(data = animal_coas,
#'          kernel.densities = kud_output,
#'          animal.tracks = trajectories,
#'          land.shape = coastline_shapefile,
#'          id.col = "ID",
#'          lat.col = "lat",
#'          lon.col = "lon")
#'  dev.off()
#' }

plotMaps <- function(data,
                     kernel.densities = NULL,
                     animal.tracks = NULL,
                     id.groups = NULL,
                     id.col = getDefaults("ID"),
                     lon.col = getDefaults("lon"),
                     lat.col = getDefaults("lat"),
                     land.shape = NULL,
                     epsg.code = getDefaults("epsg"),
                     land.color = "gray50",
                     background.layer = NULL,
                     background.pal = NULL,
                     background.title = NULL,
                     discard.missing = TRUE,
                     color.pal = NULL,
                     kud.transparency.threshold = 0.95,
                     tracks.color = "black",
                     tracks.lty = 2,
                     tracks.lwd = 0.2,
                     plot.detections = TRUE,
                     pts.color = "white",
                     pts.cex = 0.7,
                     same.scale = FALSE,
                     kud.legend = TRUE,
                     legend.cex = 0.8,
                     title.color = "black",
                     title.pos = "topleft",
                     title.inset = c(-0.08, 0),
                     title.cex = 1.5,
                     scale.meters = NULL,
                     scale.color = "black",
                     scale.pos = "bottomright",
                     scale.inset = 0.05,
                     extent.factor = 1.1,
                     cols = 3) {



  ##############################################################################
  ## Initial setup and argument validation #####################################
  ##############################################################################

  # print message indicating the type of map being generated
  if (!is.null(kernel.densities) && !is.null(animal.tracks)) {
    .printConsole("Generating maps with kernel densities and animal movement trajectories")
  } else if (!is.null(kernel.densities)) {
    .printConsole("Generating maps with kernel density estimates")
  } else if (!is.null(animal.tracks)) {
    .printConsole("Generating maps with animal movement trajectories")
  } else {
    .printConsole("Generating maps with animal positions only")
  }

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  land.shape <- reviewed_params$land.shape


  # check the format of the kernel.densities input (if provided)
  if(!is.null(kernel.densities)){
    # validate the kernel density input format
    if(any(grepl("kernel_density", names(kernel.densities), fixed=TRUE))){
      kernel.densities <- kernel.densities$kernel_density
    }
    # check if the kernel density object has the correct class (second attempt)
    if(!inherits(kernel.densities, "estUDm")){
      if(inherits(kernel.densities[[1]], "estUDm")){
        ids <- unlist(lapply(kernel.densities, names))
        kernel.densities <- unlist(kernel.densities)
        names(kernel.densities) <- ids
        class(kernel.densities) <- "estUDm"
      }else{
        stop("Supplied kernel.densities appear to be in a wrong format. Please provide the output of the 'calculateKUDs()' function.", call.=FALSE)
      }
    }
    # ensure the data contains matching IDs for kernel densities
    if(all(!unique(data[,id.col]) %in% names(kernel.densities))){
      stop("Names of kernel.densities do not match the supplied data IDs", call.=FALSE)
    }
  }

  # check the format of the animal.tracks input (if provided)
  if(!is.null(animal.tracks)){
    if(!is.null(animal.tracks$trajectories)){
      animal.tracks <- animal.tracks$trajectories
    }else{
      animal.tracks <- animal.tracks
      warning("- A 'trajectories' list was not found in the supplied animal.tracks. Function will assume that animal.tracks already contains all linestring objects", call.=FALSE)
    }
  }


  ##############################################################################
  ## Cleanup data ##############################################################
  ##############################################################################

  # define single id.group if needed
  if(is.null(id.groups)){
    id.groups <- list(levels(data[,id.col]))
  }

  # optionally discard individuals with fewer than 5 detections (as required for KUD estimation)
  if(discard.missing) {
    missing_individuals <- names(table(data[,id.col])[table(data[,id.col])<=5])
    data <- data[!data[,id.col] %in% missing_individuals,]
    data[,id.col] <- droplevels(data[,id.col])
    id.groups <- lapply(id.groups, function(x) x[!x %in% missing_individuals])
  }

  # manage spatial objects
  coords <- sf::st_as_sf(data, coords=c(lon.col, lat.col))
  spatial_data <- .processSpatial(coords, land.shape, epsg.code)
  coords <- spatial_data$coords
  land.shape <- spatial_data$spatial.layer
  epsg.code <- spatial_data$epsg.code


  ##############################################################################
  ## Prepare plot parameters ###################################################
  ##############################################################################

  # set the background palette and title based on whether a background layer is provided
  if(!is.null(background.layer)){
    if(is.null(background.pal)) background.pal <- rev(.bathy_deep_pal(100)[c(30:100)])
    if(is.null(background.title)) background.title <- "Layer values"
  } else{
    if(is.null(background.pal)) background.pal <- "#F3F7F7"
  }

  # set the color palette for plotting
  if(is.null(color.pal)){
    #color.pal <- colorRampPalette(c("#00032F", "blue2", "blue", "#8754A2", "#FD4CB3", "#F9A35D", "yellow", "#faf873"))(100)
    color.pal <- rev(terrain.colors(100))
  }

  # ensure pts.color and pts.cex can be used as vectors of three elements
  if(length(pts.color)==1) pts.color <- rep(pts.color, 3)
  if(length(pts.cex)==1) pts.cex <- rep(pts.cex, 3)

  # define the plot boundaries (bounding box) based on the animal positions
  bbox <- sf::st_bbox(coords)

  # expand bounding box by a given % in all directions
  dx <- bbox["xmax"]-bbox["xmin"]
  dy <- bbox["ymax"]-bbox["ymin"]
  bbox <- bbox + c(-dx, -dy, dx, dy) * (extent.factor - 1)

  # crop the land.shape and background.layer (if provided)
  if(!is.null(land.shape)) suppressWarnings(land.shape <- sf::st_crop(land.shape, bbox))
  if(!is.null(background.layer)) background.layer <- raster::crop(background.layer, bbox)

  # automatically define scale bar length if not specified
  if(is.null(scale.meters)) {
    scale_km <- pretty((bbox[3]-bbox[1])*0.2/1000)[2]
    scale.meters <- scale_km*1000
  }else{
    scale_km <- scale.meters/1000
  }

  # calculate the range of KUD densities for scaling, used if same.scale == TRUE
  if(!is.null(kernel.densities)){
    kud_range <- range(unlist(lapply(kernel.densities, function(x) range(x$ud))))
  }

  # handle space for a common legend, if applicable
  if((!is.null(kernel.densities) && same.scale && kud.legend) || !is.null(background.layer)) add_legend_space <- TRUE
  else add_legend_space <- FALSE

  # set layout grid
  layout_params <- .setLayout(cols, id.groups, plots.height=6, dividers.height=0.6,
                              legend=add_legend_space, min.legend.plots=1,
                              expand.legend=FALSE)
  nplots <- max(layout_params$matrix, na.rm=TRUE)
  layout_params$matrix[is.na(layout_params$matrix)] <- nplots + 1
  layout(mat=layout_params$matrix, heights=layout_params$heights)

  # set outer margin and inner margin parameters (bottom, left, top, right) )
  oma <- c(0, 3, 0, 1)
  mar <- c(0.4, 0.4, 0.4, 0.4)
  # add additional space to the right inner margin to accommodate the legend (if needed)
  if(kud.legend && !same.scale) mar[4] <- 4.5
  # if only one group (id.groups) is being plotted, remove all outer margins
  else if(length(id.groups)==1) oma <- c(0,0,0,0)
  # apply the updated margin settings
  par(mar=mar, oma=oma)

  ##############################################################################
  ## Generate maps #############################################################
  ##############################################################################

  # split the data by individual animals
  data_individual <- split(data, f=data[,id.col], drop=FALSE)

  # set progress bar
  pb <- txtProgressBar(min=1, max=length(data_individual), initial=0, style=3)


  # iterate through each plot index
  for (i in 1:nplots) {

    #################################################################
    # generate individual maps ######################################
    if(i %in% c(1:nlevels(data[,id.col]))){

      # retrieve individual ID
      id <- as.character(levels(data[,id.col])[i])

      ###########################################################
      # create an empty plot with specified bounding box ########
      plot(x=c(bbox[1], bbox[3]), c(bbox[2], bbox[4]), type='n', xlab="", ylab="", main="",
           xlim=c(bbox[1], bbox[3]), ylim=c(bbox[2], bbox[4]), axes=FALSE, asp=1, xaxs="i", yaxs="i")

      ###########################################################
      # add background layer to the plot if provided ############
      if(!is.null(background.layer)){
        raster::image(background.layer, col=background.pal, legend=FALSE, axes=FALSE, add=TRUE)
      }else{
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.pal[1], border=NA)
      }

      ###########################################################
      # check if kernel utilization densities (KUD) are available
      plot_kernel <- FALSE
      if(!is.null(kernel.densities)){
        if(nrow(data_individual[[i]])>0 && id %in% names(kernel.densities)){
          plot_kernel <- TRUE
          # convert the KUD densities for the current individual to a raster format
          r <- raster::raster(methods::as(kernel.densities[[id]], "SpatialPixelsDataFrame"))
          # save original density values
          kud_values <- raster::values(r)
          # set density values above threshold to NA for transparency
          if(!is.null(kud.transparency.threshold)) {
            # calculate the cumulative distribution
            sorted_kud_values <- sort(kud_values[!is.na(kud_values)])
            cumulative_distribution <- cumsum(sorted_kud_values)/sum(sorted_kud_values)
            # find the threshold value that corresponds to the desired percentage
            threshold_index <- which(cumulative_distribution >= (1-kud.transparency.threshold))[1]
            density_threshold <- sort(kud_values)[threshold_index]
            # apply the transparency threshold
            raster::values(r)[raster::values(r)<density_threshold] <- NA
          }
          # determine the KUD scale based on the same.scale option
          if(same.scale) kud_scale <- kud_range
          else kud_scale <- range(kud_values, na.rm=TRUE)
          # crop
          r <- raster::crop(r, raster::extent(bbox))
          # overlay KUD density on the plot
          suppressWarnings(raster::image(r, zlim=kud_scale, axes=FALSE, col=color.pal, asp=1, add=TRUE))
        }
      }

      ###########################################################
      # overlay land shape to the plot ##########################
      if(!is.null(land.shape)){
        plot(sf::st_geometry(land.shape), col=land.color, border=NA, main="", add=TRUE)
      }

      ###########################################################
      # if animal tracks are provided, plot them along with the individual detection points
      if(!is.null(animal.tracks) && !is.null(animal.tracks[[id]])){
        lines(sf::st_coordinates(animal.tracks[[id]]), col=tracks.color, lwd=tracks.lwd, lty=tracks.lty)
      }

      ###########################################################
      # overlay animal positions  ###############################
      if(plot.detections && nrow(data_individual[[i]])>0) {
        id_coords <- coords[coords[[id.col]]==id,]
        points(sf::st_coordinates(id_coords), pch=21, bg=pts.color[2], cex=pts.cex[2], col="gray15", lwd=0.25)
        points(sf::st_coordinates(id_coords[1,]), pch=21, bg=pts.color[1], cex=pts.cex[1], lwd=0.25)
        points(sf::st_coordinates(id_coords[nrow(id_coords),]), pch=21, bg=pts.color[3], cex=pts.cex[3], lwd=0.25)
      }

      ###########################################################
      # add individual ID as a label on the top left corner of the plot
      legend(title.pos, legend=id, inset=title.inset , cex=title.cex, text.col=title.color, bty="n", text.font=2)

      ###########################################################
      # position the scale bar based on the specified scale position and inset
      scale_xy <- .getPosition(scale.pos, inset=scale.inset, bar.width = scale.meters)
      .scalebar(d=scale.meters, xy=scale_xy, type="bar", divs=2, below="km", lwd=0.2,
                label=c(0, scale_km/2, scale_km), label.color=scale.color, cex=0.7)

      ###########################################################
      # add a box around the plot with specified line type and width
      box(lty="solid", lwd=0.75, cex=0.5)

      ###########################################################
      # if same.scale is FALSE and kud.legend is TRUE, add a KUD density color legend
      if(!is.null(kernel.densities) && !same.scale && kud.legend && plot_kernel){
        kud_ticks <- quantile(kud_scale, probs=seq(0, 1, by=0.2), na.rm=TRUE)
        kud_labels <- paste0(seq(0, 100, by=20), "%")
        # get relative plot area dimensions
        plt <- par("plt")
        plot_width <- plt[2] - plt[1]
        # define fixed offset and bar width as a fraction of plot width
        fixed_offset <- 0.025
        bar_width <- 0.025
        # calculate posx relative to the plot area, ensuring a consistent distance from the right
        posx <- c(plt[2] + fixed_offset, plt[2] + fixed_offset + bar_width)
        # plot color scale
        .colorlegend(col=color.pal, zlim=kud_scale, zval=kud_ticks, zlab=kud_labels,
                     posx=posx, posy=c(0.2, 0.7), main="KUD\ndensity\n\n",
                     main.cex=legend.cex+0.1, main.adj=0, cex=legend.cex)
      }
    }

    #################################################################
    # add common legends ############################################
    if(i==nplots && add_legend_space){

      # add empty plot
      plot.new()

      # initialize a variable to check if the first legend was plotted
      first_legend <- FALSE

      # add color legend for uniformized KUD densities if applicable
      if(!is.null(kernel.densities) && same.scale && kud.legend){
        kud_ticks <- quantile(kud_range, probs=seq(0, 1, by=0.2), na.rm=TRUE)
        kud_labels <- paste0(seq(0, 100, by=20), "%")
        first_legend <- TRUE
        .colorlegend(col=color.pal, zlim=kud_range, zval=kud_ticks, zlab=kud_labels,
                     posx=c(0.15, 0.85), posy=c(0.75, 0.8), main="KUD density",
                     main.cex=legend.cex+0.2, lab.scientific=FALSE,
                     cex=legend.cex, horizontal=TRUE)
      }

      # add color legend for the background layer if available
      if(!is.null(background.layer)) {
        layer_range <- range(raster::values(background.layer), na.rm=TRUE)
        layer_labs <- pretty(raster::values(background.layer), min.n=4)
        layer_labs <- layer_labs[layer_labs>=layer_range[1] & layer_labs<=layer_range[2]]
        display_digits <- max(.decimalPlaces(layer_labs))
        posy_vals <- if (first_legend) c(0.4, 0.45) else c(0.75, 0.8)
        .colorlegend(col=background.pal, zlim=layer_range, zval=layer_labs,
                     main=background.title, posx=c(0.15, 0.85), posy=posy_vals,
                     main.cex=legend.cex+0.2, digit=display_digits,
                     cex=legend.cex, horizontal=TRUE)
      }
    }

    #################################################################
    # update the progress bar #######################################
    setTxtProgressBar(pb,i)
  }


  #################################################################
  # add id.group labels ###########################################
  if (length(id.groups)>1) {
    label_pos <- layout_params$group_positions
    layout_height <- sum(layout_params$heights)
    label_pos <- grconvertY(1 - (label_pos / layout_height), "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels = names(id.groups),
         srt = 90, cex = title.cex + 0.2, font = 2, xpd = NA, adj = c(0.5, 0.5))
  }


  ##############################################################################
  ## Wrap it up ################################################################
  ##############################################################################

  # close the progress bar
  close(pb)

  # print a summary message to the console if individuals were discarded
  if(discard.missing) warning(paste0("- ", length(missing_individuals), " individual(s) with < 5 detections discarded."), call.=FALSE)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
