#######################################################################################################
# Function to plot spatial network(s) #################################################################
#######################################################################################################

#' Plot migrations network
#'
#' @description Plots a spatial network representing animal transitions/movements between receivers/sites.
#' If tagging information is supplied, the number of tagged individuals is indicated for each location/site, and
#' the release site is included in the calculated transitions.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections.
#' @param lon.col Name of the column containing longitude values (either projected or unprojected). Defaults to 'lon'.
#' @param lat.col Name of the column containing latitude values (either projected or unprojected). Defaults to 'lat'.
#' @param spatial.col Name of the column in both `data` and `id.metadata` (if supplied)
#' that contains spatial information for calculating transitions. This can include receiver IDs
#' in case all movements/transitions are of interest, or contain (for example) location/habitat
#' classes for broader-scale analyses.
#' @param id.groups Optional. A list containing ID groups (e.g. different species or size-classes).
#' If supplied, a separate map will be plotted for each group.
#' @param id.metadata Optional. A data frame containing information on the tagging location of each animal.
#' Used to indicate the number of tagged individuals at each site.
#' @param land.shape Optional. A shape file containing coastlines.
#' @param land.color Optional. Color of land areas. Defaults to "gray50".
#' @param background.layer Optional. A projected raster containing a variable to be displayed in the background
#' (e.g. bathymetry, temperature, etc.).
#' @param background.pal Color palette for the background layer. If a background layer is not supplied,
#' the first color is used as background. Defaults to "#F3F7F7".
#' @param title.cex  Size of the plot title. Defaults to 1.
#' @param color.nodes.by Variable used to color nodes (locations). It can be set to any variable in the dataset.
#' Alternatively, it can be set to 'detections' (default) to distinguish between nodes with and without detections,
#' or to 'group' to assign different colors to different ID groups (one color per map). Defaults to "detection".
#' @param nodes.color Color(s) for the nodes. Defaults to c("darkblue", "black").
#' @param nodes.alpha Numeric. Transparency level for the nodes (0 to 1): 0=fully transparent, 1=fully opaque. Defaults to 0.8.
#' @param nodes.size A numeric vector of length 2. Represents the desired min and max vertex sizes relative
#' to the x-axis in terms of percentage (see details). See \code{\link[netdiffuseR]{rescale_vertex_igraph}}. Defaults to c(0.04, 0.08).
#' @param nodes.label.wrap Logical. If TRUE, splits node labels into multiple lines (to better fit within the nodes). Defaults to FALSE.
#' @param nodes.label.cex The font size for vertex labels. Defaults to 0.5.
#' @param nodes.label.color The font color for vertex labels. Defaults to white.
#' @param repel.nodes Logical. If TRUE, nodes are plotted using a repulsion algorithm to avoid overlap. Defaults to TRUE
#' @param repel.buffer Controls the amount of space between nodes if repel.nodes is set to TRUE. Defaults to 1.1.
#' @param edge.type A character string indicating the metric to be used to calculate the network edges.
#' It can be either "movements" to represent the number of movements between sites or "individuals"
#' to represent the number of individual animals moving between sites. Defaults to "movements".
#' @param edge.color Color(s) for the edges. Defaults to "darkblue".
#' @param edge.curved Specifies whether to draw curved edges, or not. This can be a logical or a numeric vector or scalar.
#' A numeric value specifies the curvature of the edge; zero curvature means straight edges, negative values mean the edge bends
#' clockwise, positive values the opposite. TRUE means curvature 0.5, FALSE means curvature zero. Defaults to 0.5.
#' @param edge.width A numeric vector of length 2. Represents the desired min and max
#' width/thickness of the edges. Defaults to c(0.4, 3.5).
#' @param edge.arrow.size The size of the arrows. Defaults to 0.5.
#' @param edge.arrow.width The width of the arrows. Defaults to 1.5.
#' @param edge.label.cex The font size for edge labels. Defaults to 0.6.
#' @param edge.label.color The color of the edge labels. Defaults to black.
#' @param edge.label.font The font for the edge labels. It is interpreted the same way
#' as the font graphical parameter: 1 is plain text, 2 is bold face, 3 is italic,
#' 4 is bold and italic and 5 specifies the symbol font. Defaults to 1.
#' @param scale.km Distance covered by the scale bar, in kilometers. If NULL, it is
#' automatically defined as 20% of the plot region.
#' @param scale.pos Position of the map scale, specified by keyword.
#' See \code{\link[grDevices]{xy.coords}}. Defaults to "bottom".
#' @param scale.inset Controls how far the scale bar is placed from the plot edges.
#' Can be specified as a single value (applied to both x and y directions) or as a two-element vector
#' (with x and y inset values specified separately).
#' - **Positive values** move the scale bar *inward*, towards the center of the plot.
#' - **Negative values** move the scale bar *outward*, potentially placing it outside the visible plot area.
#'
#' Defaults to c(0, 0.05).
#' @param scale.height Controls the thickness of the scale bar. Defaults to 1.5.
#' @param scale.cex Size of the scale bar values. Defaults to 0.6.
#' @param extent.factor Numeric. Factor by which to adjust the extent of the plotting region,
#' defined based on the bounding box around animal positions/detections. A value of 1 keeps
#' the original bounding box, values greater than 1 increase the extent, and values less
#' than 1 decrease it. Defaults to 1.1 (10% increase).
#' @param cols Number of columns in the plot panel layout (used in the 'mfrow' argument). Defaults to 1.
#' @param ... Arguments passed to the \code{\link[igraph]{plot.igraph}} function.
#' See \code{\link[igraph]{igraph.plotting}} for the complete parameters list.
#'
#' @return A plot displaying the migrations network. Nodes correspond to different sites,
#' with coordinates reflecting their approximate location in geographic space.
#' Arrows illustrate movements between sites. The size of each node indicates the number of
#' individuals detected at that site, and labels on the arrows show either the number of
#' transition movements or the the number of individuals transiting between sites.
#' Optionally, the number of tagged individuals can be indicated in parentheses
#' within the nodes of the corresponding sites.
#'
#' @examples
#' \dontrun{
#'pdf("./spatial-network.pdf", width=6, height=8)
#' plotMigrations(data=data_filtered,
#'                id.groups=species_ids,
#'                spatial.col="site",
#'                land.shape=coastline,
#'                edge.type="individuals",
#'                edge.color=c("darkblue", "darkred"),
#'                edge.width=c(0.2, 4),
#'                color.nodes.by="group",
#'                nodes.color=c("darkblue", "darkred"),
#'                nodes.size=c(0.05,0.1),
#'                nodes.label.wrap=T,
#'                repel.nodes=T)
#' dev.off()
#' }

#' @export

plotMigrations <- function(data,
                           id.col = getDefaults("ID"),
                           lon.col = getDefaults("lon"),
                           lat.col = getDefaults("lat"),
                           datetime.col = getDefaults("datetime"),
                           spatial.col,
                           id.groups = NULL,
                           id.metadata = NULL,
                           land.shape = NULL,
                           land.color = "gray50",
                           epsg.code = getDefaults("epsg"),
                           background.layer = NULL,
                           background.pal = "#F3F7F7",
                           title.cex = 1,
                           color.nodes.by = "detection",
                           nodes.color = c("darkblue", "black"),
                           nodes.alpha = 0.8,
                           nodes.size = c(0.05, 0.10),
                           nodes.label.wrap = FALSE,
                           nodes.label.cex = 0.5,
                           nodes.label.color = "white",
                           repel.nodes = TRUE,
                           repel.buffer = 1.1,
                           edge.type = "movements",
                           edge.color = "darkblue",
                           edge.curved = 0.5,
                           edge.width = c(0.4, 3.5),
                           edge.arrow.size = 0.5,
                           edge.arrow.width = 1.5,
                           edge.label.cex = 0.6,
                           edge.label.color = "black",
                           edge.label.font = 1,
                           scale.km = NULL,
                           scale.pos = "bottom",
                           scale.inset = c(0, 0.05),
                           scale.height = 1.5,
                           scale.cex = 0.6,
                           extent.factor = 1.1,
                           cols = 1,
                           ...) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  # print to console
  .printConsole("Generating migrations network")

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  land.shape <- reviewed_params$land.shape

  # validate additional parameters
  errors <- c()
  # check if id.metadata contains the required columns
  if(!is.null(id.metadata)){
    if(!id.col %in% colnames(id.metadata)) errors <- errors("'id.col' variable not found in the supplied id.metadata")
    if(!spatial.col %in% colnames(id.metadata)) errors <- errors("'spatial.col' variable not found in the supplied id.metadata")
  }
  if(!requireNamespace("igraph", quietly=TRUE)) errors <- errors("The 'igraph' package is required but is not installed. Please install 'igraph' using install.packages('igraph') and try again.")
  if(repel.nodes && !requireNamespace("packcircles", quietly=TRUE)) errors <- errors("The 'packcircles' package is required for the repel.nodes feature but is not installed.
                                                                                     Please install it by running install.packages('packcircles'), and then try again.")
  if(!edge.type %in% c("movements", "individuals")) errors <- errors("'edge.type' must be either 'movements' or 'individuals'")
  if(!color.nodes.by %in% c("detection", "group", colnames(data))) errors <- errors("Please set a valid 'color.nodes.by' argument")
  if(length(nodes.size)!=2) errors <- errors("Please supply a min and max value for the node sizes")
  if(length(edge.width)!=2) errors <- errors("Please supply a min and max value for the edge width")
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # save the current par settings and ensure they are restored upon function exit
  original_par <- par(no.readonly=TRUE)
  on.exit(par(original_par))

  # manage spatial objects
  site_coords <- stats::aggregate(data[,c(lon.col, lat.col)], by=list(data[,spatial.col]), mean)
  colnames(site_coords)[1] <- "site"
  coords <- sf::st_as_sf(site_coords, coords=c(lon.col, lat.col))
  spatial_data <- .processSpatial(coords, land.shape, epsg.code)
  coords <- spatial_data$coords
  land.shape <- spatial_data$spatial.layer
  epsg.code <- spatial_data$epsg.code

  # check layer projections
  if(!is.null(background.layer) && !is.null(land.shape)){
    if(!raster::compareCRS(raster::projection(background.layer), epsg.code$proj4string)){
      background.layer <- raster::projectRaster(background.layer, crs=epsg.code$proj4string)
      warning_message <- "- The 'background.layer' has been reprojected to match 'land.shape'."
      warning(paste(strwrap(warning_message, width=getOption("width")), collapse="\n"), call.=FALSE)
    }
  }

  # check if there color palette is
  if(!is.null(background.layer) && length(background.pal)==1){
    warning_message <- "- Using a new color palette for the background layer. To override, please supply > 1 color."
    warning(paste(strwrap(warning_message, width=getOption("width")), collapse="\n"), call.=FALSE)
    if(is.null(levels(background.layer))) background.pal <- adjustcolor(terrain.colors(100), 0.75)
    if(!is.null(levels(background.layer))) background.pal <- terrain.colors(length(levels(background.layer)[[1]][,2]))
  }

  # check color nodes argument
  if(color.nodes.by=="group" && is.null(id.groups)){
    warning_message <- "Parameter 'color.nodes.by' is set to 'group' but id.groups were not supplied. Only the first node color will be used."
    warning(paste(strwrap(warning_message, getOption("width")), collapse="\n"), call.=FALSE)
  }



  #####################################################################################
  # Prepare data ######################################################################
  #####################################################################################

  # set id groups
  if(is.null(id.groups)) id.groups <- list(levels(data[,id.col]))

  # convert to factor
  if(!is.factor(data[,spatial.col])) data[,spatial.col] <- as.factor(data[,spatial.col])

  # split data
  n_groups <- length(id.groups)
  data_groupped <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])
  groupped_tags <- lapply(id.groups, function(x) id.metadata[id.metadata[,id.col] %in% x,])
  tagged_stats <- list()


  #####################################################################################
  # Calculate animal movements (transitions) ##########################################
  #####################################################################################

  transition_movements <- list()
  transition_animals <- list()
  node_detections <- list()
  node_individuals <- list()
  node_tagged <- list()

  # calculate metrics independently for each group
  for(g in 1:n_groups){

    #  data by individual
    data_individual <- split(data_groupped[[g]], f=data_groupped[[g]][,id.col], drop=TRUE)

    # calculate transitions
    transitions_list <- list()
    for(i in 1:length(data_individual)){
      transition_data <- data_individual[[i]][,c(datetime.col, spatial.col)]
      transition_data <- as.character(transition_data[order(transition_data[,datetime.col]),-1])
      transition_data <- factor(transition_data, levels=levels(data_groupped[[g]][,spatial.col]))
      # include transitions from the release site if tagging info is available
      if(!is.null(id.metadata)){
        id <- unique(data_individual[[i]][,id.col])
        tagging_site <- as.character(id.metadata[id.metadata[,id.col]==id, spatial.col])
        tagging_site <- factor(tagging_site, levels=levels(transition_data))
        transition_data <- c(tagging_site, transition_data)
      }
      transitions_list[[i]] <- table(head(transition_data,-1), tail(transition_data,-1))
    }
    transition_movements[[g]] <- unclass(Reduce("+", transitions_list))
    transitions_list <- lapply(transitions_list, function(x) {x[x>0]<-1; return(x)})
    transition_animals[[g]] <- unclass(Reduce("+", transitions_list))

    # calculate nº detections on each site
    node_detections[[g]] <- table(data_groupped[[g]][,spatial.col])

    # calculate nº individuals on each site
    nindividuals <- stats::aggregate(data_groupped[[g]][,id.col], by=list(data_groupped[[g]][,spatial.col]),
                             function(x) length(unique(x)), drop=FALSE)
    colnames(nindividuals) <- c(spatial.col, "nids")
    nindividuals$nids[is.na(nindividuals$nids)] <- 0
    node_individuals[[g]]  <- nindividuals

    # calculate proportion of tagged/visiting individuals on each site
    if(!is.null(id.metadata)){
      detected_individuals <- split(data_groupped[[g]][,id.col], f=data_groupped[[g]][,spatial.col])
      detected_individuals <- lapply(detected_individuals, function(x) as.character(unique(x)))
      tagged_individuals <- split(groupped_tags[[g]][,id.col], f=groupped_tags[[g]][,spatial.col])
      tagged_individuals <- lapply(tagged_individuals, function(x) as.character(unique(x)))
      tagged <- mapply(function(group, ids){length(which(ids %in% tagged_individuals[[group]]))},
                       group=names(detected_individuals), ids=detected_individuals)
      detected <- unlist(lapply(detected_individuals, length))
      visitors <- mapply(function(group, ids){length(which(!ids %in% tagged_individuals[[group]]))},
                         group=names(detected_individuals), ids=detected_individuals)
      tag_stats <- data.frame("tagged"=tagged, "detected"=detected, "visiting"=visitors)
      tag_stats$site <- rownames(tag_stats)
      rownames(tag_stats) <- NULL
      tag_stats <- tag_stats[,c(4,1:3)]
      node_tagged[[g]] <- tag_stats
    }
  }


  ##############################################################################
  ## Prepare plot parameters ###################################################
  ##############################################################################

  # define the plot boundaries (bounding box) based on the animal positions
  bbox <- sf::st_bbox(coords)

  # expand bounding box by a given % in all directions
  dx <- bbox["xmax"]-bbox["xmin"]
  dy <- bbox["ymax"]-bbox["ymin"]
  bbox <- bbox + c(-dx, -dy, dx, dy) * (extent.factor - 1)

  # crop 'land.shape' and 'background.layer' to the bounding box
  if(!is.null(land.shape)) land.shape <- sf::st_crop(sf::st_geometry(land.shape), bbox)
  if(!is.null(background.layer)) background.layer <- raster::crop(background.layer, raster::extent(bbox))

  # configure layout dimensions for multi-panel plots
  rows <- ceiling(n_groups/cols)
  par(mfrow=c(rows, cols), oma=c(0,0,0,0))

  # set layout variables
  if(!is.null(background.layer)) par(mar=c(1,2,2,6))
  else par(mar=c(1,2,2,1))

  # set diagonal values to 0 in movement transition matrices (self-transitions not relevant)
  transition_movements <- lapply(transition_movements, function(x) {diag(x)<-0; return(x)})
  transition_animals <- lapply(transition_animals, function(x) {diag(x)<-0; return(x)})

  # define edge values based on the specified 'edge.type' parameter
  if(edge.type=="movements") all_movements <- reshape2::melt(transition_movements)$value
  else if(edge.type=="individuals") all_movements <- reshape2::melt(transition_animals)$value

  # if a single edge color is provided, replicate it to match the number of groups
  if(length(edge.color)==1){
    edge.color <- rep(edge.color, n_groups)
  }


  #####################################################################################
  # Plot network(s) ###################################################################
  #####################################################################################

  # loop through each one of the id.groups
  for(g in 1:n_groups){

    # create transition matrix (edges)
    if(edge.type=="movements") edges <- reshape2::melt(transition_movements[[g]])
    if(edge.type=="individuals") edges <- reshape2::melt(transition_animals[[g]])
    colnames(edges) <- c("site1", "site2", "value")
    edges <- edges[edges$value>0,]
    edges <- edges[edges$site1!=edges$site2,]

    # create an empty plot with specified bounding box
    plot(x= bbox[c("xmin", "xmax")], y=bbox[c("ymin", "ymax")], type='n',
         xlab="", ylab="", main="", axes=FALSE, asp=1, xaxs="i", yaxs="i")

    # add background layer or background color
    if(is.null(background.layer)){
      rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4], col=background.pal[1], border="black")
    }else{
      if(is.null(levels(background.layer))) {
        raster::plot(background.layer, col=background.pal, legend=TRUE, axes=FALSE, add=TRUE,
             axis.args=list(cex.axis=0.7), legend.width=1.2, legend.shrink=0.6)
      }else{
        raster::plot(background.layer, col=background.pal, legend=FALSE, axes=FALSE, add=TRUE)
        legend("right", legend=rev(levels(background.layer)[[1]][,2]), xpd=TRUE,
               fill=rev(background.pal), inset=-0.18, bty="n", cex=0.7)
      }
      box()
    }

    # add land.shape
    if(!is.null(land.shape)) plot(sf::st_geometry(land.shape), col=land.color, border=NA, add=TRUE)

    # add network title
    if(n_groups>1) title(main=names(id.groups)[g], cex.main=title.cex, line=1.2, xpd=TRUE)
    else title(main="Spatial Network", cex.main=title.cex, line=1.2, xpd=TRUE)
    metric <- ifelse(edge.type=="movements", "n\u00ba of movements", "n\u00ba of transiting individuals")
    title(main=tools::toTitleCase(metric), cex.main=title.cex-0.3, line=0.4, font.main=1, xpd=TRUE)

    # add scale bar
    if(is.null(scale.km)) scale.km <- round(min(pretty((bbox["xmax"]-bbox["xmin"])*0.15))/1000)
    scale_xy <- .getPosition(scale.pos, inset=scale.inset)
    scale.meters <- scale.km*1000
    # adjusting the x position based on alignment
    if (grepl("right", scale.pos)) scale_xy[1] <- scale_xy[1] - scale.meters
    else if (grepl("center", scale.pos) || scale.pos %in% c("top", "bottom")) scale_xy[1] <- scale_xy[1] - (scale.meters/2)
    .scalebar(d=scale.meters, xy=scale_xy, type="bar", divs=2, below="km", bar.height=scale.height,
              label=c(0, scale.km/2, scale.km), lwd=0.2, cex=scale.cex, bar.lwd=0.2)


    #set network variables
    nanimals <- node_individuals[[g]]
    colnames(nanimals)[1] <- "site"
    nanimals <- nanimals[match(coords$site, nanimals$site), "nids"]
    map_bbox <- raster::extent(bbox)
    node_size <- .rescale_vertex_igraph(nanimals, minmax.relative.size=c(nodes.size[1], nodes.size[2]))

    # set node names
    if(!is.null(id.metadata)){
      node_titles <- paste0(node_individuals[[g]]$site, " (", node_tagged[[g]]$tagged, ")",
      "\n", "n=", node_individuals[[g]]$nids)
      node_titles <- gsub(" (0)", "", node_titles, fixed=TRUE)
    }else{
      node_titles <- paste0(node_individuals[[g]]$site, "\n", "n=", node_individuals[[g]]$nids)
    }
    if(nodes.label.wrap==TRUE) node_titles <- gsub(" ", "\n", node_titles, fixed=TRUE)

    # set node colors
    if(color.nodes.by=="detection"){
      vertex_color <- rep(nodes.color[1], length.out=length(nanimals))
      vertex_color[nanimals==0] <- nodes.color[2]
    } else if(color.nodes.by=="group"){
      vertex_color <- rep(nodes.color[g], length.out=length(nanimals))
    }else if(is.numeric(data[,color.nodes.by])){
      node_classes <- stats::aggregate(data[,color.nodes.by], by=list(data[,spatial.col]), mean, na.rm=TRUE)
      colnames(node_classes) <- c("node", "class")
      node_classes$class <- round(.rescale(node_classes$class, to=c(1,100)))
      nodes.color <- colorRampPalette(nodes.color)(100)
      vertex_color <- nodes.color[node_classes$class]
    } else {
      node_classes <- stats::aggregate(data[,color.nodes.by], by=list(data[,spatial.col]), function(x) names(table(x))[which.max(table(x))])
      colnames(node_classes) <- c("node", "class")
      node_classes$class <- factor(node_classes, levels=levels(data[,color.nodes.by]))
      vertex_color <- nodes.color[node_classes$class]
    }
    vertex_color <- adjustcolor(vertex_color, alpha.f=nodes.alpha)

    # set nodes position
    if(repel.nodes){
      radii <- (par("usr")[2]-par("usr")[1])*c(0.035, 0.07)
      radii <- plotrix::rescale(node_size, newrange=radii)
      repel_data <- cbind(sf::st_coordinates(coords), "radius"=radii*repel.buffer)
      new_coords <- packcircles::circleRepelLayout(repel_data, xlim=map_bbox[1:2], ylim=map_bbox[3:4],
                                                   sizetype=c("radius"), wrap=FALSE, weights=1)$layout
      new_coords <- cbind(site_coords[,1], new_coords[,c("x","y")])
      colnames(new_coords) <- c(spatial.col, "x", "y")
      coords <- sf::st_as_sf(new_coords, coords=c("x", "y"), crs=epsg.code)
    }

    # plot network
    network <- igraph::graph_from_data_frame(edges, directed=TRUE, vertices=coords)
    igraph::V(network)$name <- node_titles
    igraph::E(network)$width <- .rescale(edges$value, from=range(all_movements), to=edge.width)
    igraph::plot.igraph(network,
                        layout = sf::st_coordinates(coords),
                        add = TRUE,
                        rescale = FALSE,
                        edge.color = edge.color[g],
                        edge.curved = edge.curved,
                        vertex.size = node_size,
                        vertex.color = vertex_color,
                        vertex.label.family = "Helvetica",
                        vertex.label.color = nodes.label.color,
                        vertex.label.cex = nodes.label.cex,
                        edge.arrow.size = edge.arrow.size,
                        edge.arrow.width = edge.arrow.width,
                        edge.label = paste0("(",edges$value,")\n"),
                        edge.label.family = "Helvetica",
                        edge.label.cex = edge.label.cex,
                        edge.label.color = edge.label.color,
                        edge.label.font = edge.label.font,
                        ...)
  }
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
