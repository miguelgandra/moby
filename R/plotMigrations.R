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
#' @param spatial.col Name of the variable/column containing the spatial information used to calculate transitions.
#' It can include receiver IDs in case all movements/transitions are of interest,
#' or contain (for example) location/habitat classes for broader-scale analyses.
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
#' @param repel.nodes Logical. If TRUE, nodes are plotted using a repulsion algorithm to avoid overlap. Defaults to FALSE.
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
#' @param scale.inset Inset distance(s) of the map scale from the margins as a fraction
#' of the plot region. If a single value is given, it is used for both margins; if two values
#' are given, the first is used for x- distance, the second for y-distance. Defaults to c(0, 0.05).
#' @param scale.height Controls the thickness of the scale bar. Defaults to 1.5.
#' @param scale.cex Size of the scale bar values. Defaults to 0.6.
#' @param bbox.ext Amount of extension of the bounding box around the coordinates.
#' It can be numeric (same extension into all four directions), vector of two (first x, then y directional extension)
#' or vector of four (xmin, xmax, ymin, ymax extension). Default is 0.2 (extends the bounding box by 20%).
#' Only considered if a background layer is not supplied. Defaults to 0.2.
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
                           id.col = getDefaults("id"),
                           lon.col = getDefaults("lon"),
                           lat.col = getDefaults("lat"),
                           datetime.col = getDefaults("datetime"),
                           spatial.col,
                           id.groups = NULL,
                           id.metadata = NULL,
                           land.shape = NULL,
                           land.color = "gray50",
                           background.layer = NULL,
                           background.pal = "#F3F7F7",
                           title.cex = 1,
                           color.nodes.by = "detection",
                           nodes.color = c("darkblue", "black"),
                           nodes.alpha = 0.8,
                           nodes.size = c(0.04, 0.08),
                           nodes.label.wrap = FALSE,
                           nodes.label.cex = 0.5,
                           nodes.label.color = "white",
                           repel.nodes = FALSE,
                           repel.buffer = 1.1,
                           edge.type = c("movements", "individuals"),
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
                           bbox.ext = 0.2,
                           ...) {

  #####################################################################################
  # Initial checks ####################################################################
  #####################################################################################

  # print to console
  .printConsole("Generating migrations network")

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # check if id.metadata contains the required columns
  if(!is.null(id.metadata)){
    if(!id.col %in% colnames(id.metadata)) stop("'id.col' variable not found in the supplied id.metadata", call.=FALSE)
    if(!spatial.col %in% colnames(id.metadata)) stop("'spatial.col' variable not found in the supplied id.metadata", call.=FALSE)
  }

  # check if edge.type is valid
  if(!edge.type %in% c("movements", "individuals")) stop("'edge.type' must be either 'movements' or 'individuals'", call. = FALSE)

  # check color.nodes.by argument
  if(!color.nodes.by %in% c("detection", "group", colnames(data))) stop("Please set a valid 'color.nodes.by' argument", call.=FALSE)
  if(color.nodes.by=="group" & is.null(id.groups)) warning("'color.nodes.by' set to group but id.groups were not supplied. Only the first node color will be used", call.=FALSE)

  # check if node size is of length 2
  if(length(nodes.size)!=2) stop("Please supply a min and max value for the node sizes", call.=FALSE)

  # check if edge edge.width is of length 2
  if(length(edge.width)!=2) stop("Please supply a min and max value for the edge width", call.=FALSE)

  # check if dataset coordinates are projected (meters)
  geographic_coords <- all(data[, lon.col]>=-180 & data[,lon.col]<=180 & data[,lat.col]>=-90 & data[,lat.col]<=90)
  if(geographic_coords) {
    if(!is.null(land.shape) || !is.null(background.layer)) {
      if(!is.null(land.shape)) epsg_code <- raster::proj4string(land.shape)
      else if(!is.null(background.layer)) epsg_code <- raster::proj4string(background.layer)
      if(!raster::isLonLat(epsg_code)){
        warning("Projecting coordinates based on the CRS of the supplied land shape or background layer", call.=FALSE)
        coords <- sp::SpatialPoints(data[,c(lon.col, lat.col)])
        raster::projection(coords) <- sp::CRS("+proj=longlat +datum=WGS84")
        coords <- sp::spTransform(coords, epsg_code)
        data[,lon.col] <- coords@coords[,1]
        data[,lat.col] <- coords@coords[,2]
        geographic_coords <- FALSE
      }
    }
  }

  # check if scale bar can be plotted with geographic coordinates
  if(!is.null(scale.km) && geographic_coords) warning("Projected coordinates are required to plot a scale bar. Please provide coordinates in a projected coordinate system.", call.=FALSE)

  # check layer projections
  if(!is.null(background.layer) & !is.null(land.shape)){
    if(!raster::compareCRS(land.shape@proj4string, background.layer@crs))
      warning("land.shape and background.layer seem to have different CRS projections", call.=FALSE)
  }

  # check if there color palette is
  if(!is.null(background.layer) & length(background.pal)==1){
    warning("Using a new color palette for the background layer. To override, please supply > 1 color\n")
    if(is.null(levels(background.layer))) background.pal <- adjustcolor(terrain.colors(100), 0.75)
    if(!is.null(levels(background.layer))) background.pal <- terrain.colors(length(levels(background.layer)[[1]][,2]))
  }

  # check if the igraph package is installed
  if(!requireNamespace("igraph", quietly=TRUE)){
    stop("The 'igraph' package is required but is not installed. Please install 'igraph' using install.packages('igraph') and try again.", call.=FALSE)
  }


  #####################################################################################
  # Prepare data ######################################################################
  #####################################################################################

  # set id groups
  if(is.null(id.groups)) {
    id.groups<-list(levels(data[,id.col]))
  }

  # convert to factor
  if(!is.factor(data[,spatial.col])){
    data[,spatial.col] <- as.factor(data[,spatial.col])
  }

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
    data_individual <- split(data_groupped[[g]], f=data_groupped[[g]][,id.col], drop=T)

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
        transition_data <- c(tagging_site, transition_data)
        transition_data <- factor(transition_data, levels=levels(data[,spatial.col]))
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
                             function(x) length(unique(x)), drop=F)
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


  #####################################################################################
  # Set geographic coordinates ########################################################
  #####################################################################################

  # convert node coords to spatial points
  site_coords <- stats::aggregate(data[,c(lon.col, lat.col)], by=list(data[,spatial.col]), mean)
  colnames(site_coords)[1] <- "site"
  site_coords <- sp::SpatialPointsDataFrame(site_coords[,2:3], data=site_coords)

  # project coords if needed
  if(geographic_coords==T){
    sp::proj4string(site_coords) <- sp::CRS("+proj=longlat +datum=WGS84")
  }else{
    sp::proj4string(site_coords) <- sp::CRS(epsg_code)
  }

  # format bbox.ext multiplier
  if(length(bbox.ext)==1) bbox.ext <- rep(bbox.ext, 4)
  if(length(bbox.ext)==2) bbox.ext <- rep(bbox.ext/2, each=2)

  # set bounding box polygon
  if(!is.null(background.layer)){
    map_extent <- extent(background.layer)
  }else{
    map_extent <- extent(site_coords)
    map_extent[1] <- map_extent[1] - (map_extent[2]-map_extent[1])*bbox.ext[1]
    map_extent[2] <- map_extent[2] + (map_extent[2]-map_extent[1])*bbox.ext[2]
    map_extent[3] <- map_extent[3] - (map_extent[4]-map_extent[3])*bbox.ext[3]
    map_extent[4] <- map_extent[4] + (map_extent[4]-map_extent[3])*bbox.ext[4]
  }

  if(!is.null(land.shape)) {land.shape <- crop(land.shape, map_extent)}
  map_extent <- SpatialPolygons(list(Polygons(list(Polygon(map_extent)),"1")))


  #####################################################################################
  # Plot network(s) ###################################################################
  #####################################################################################

  # set layout variables
  if(!is.null(background.layer)) {
    par(mfrow=c(n_groups,1), mar=c(1,2,2,6), oma=c(0,0,0,0))
  }else{
    par(mfrow=c(n_groups,1), mar=c(1,2,2,1), oma=c(0,0,0,0))
  }
  transition_movements <- lapply(transition_movements, function(x) {diag(x)<-0; return(x)})
  transition_animals <- lapply(transition_animals, function(x) {diag(x)<-0; return(x)})

  if(edge.type=="movements") all_movements <- reshape2::melt(transition_movements)$value
  else if(edge.type=="individuals") all_movements <- reshape2::melt(transition_animals)$value

  # loop through each one of the id.groups
  for(g in 1:n_groups){

    # create transition matrix (edges)
    if(edge.type=="movements") edges <- reshape2::melt(transition_movements[[g]])
    if(edge.type=="individuals") edges <- reshape2::melt(transition_animals[[g]])
    colnames(edges) <- c("site1", "site2", "value")
    edges <- edges[edges$value>0,]
    edges <- edges[edges$site1!=edges$site2,]

    # plot bounding box
    plot(map_extent, col=NA, border=F, axes=F, xaxs="i", yaxs="i")
    rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4], col=background.pal[1], border="black")

    # add background layer or background color
    if(!is.null(background.layer)){
      if(is.null(levels(background.layer))) {
        plot(background.layer, col=background.pal, legend=T, axes=F, add=T,
             axis.args=list(cex.axis=0.7))
      }else{
        plot(background.layer, col=background.pal, legend=F, axes=F, add=T)
        legend("right", legend=rev(levels(background.layer)[[1]][,2]), xpd=T,
               fill=rev(background.pal), inset=-0.18, bty="n", cex=0.7)
      }
    }

    if(!is.null(land.shape)){
      plot(land.shape, col=land.color, border=NA, add=T)
    }

    # add network title
    if(n_groups>1) title(main=names(id.groups)[g], cex.main=title.cex, line=1.2, xpd=T)
    else title(main="Spatial Network", cex.main=title.cex, line=1.2, xpd=T)
    metric <- ifelse(edge.type=="movements", "n\u00ba of movements", "n\u00ba of transiting individuals")
    title(main=tools::toTitleCase(metric), cex.main=title.cex-0.3, line=0.4, font.main=1, xpd=T)

    # # add network title
    # metric <- ifelse(edge.type=="movements", "nº of movements", "nº of transiting individuals")
    # if(n_groups>1) title(main=paste(names(id.groups)[g], "-", metric), cex.main=title.cex, line=1, xpd=T)
    # else title(main=tools::toTitleCase(metric), cex.main=title.cex, line=1, xpd=T)

    # add scale bar
    if(geographic_coords==FALSE){
      if(is.null(scale.km)) scale.km <- min(pretty((bbox[2]-bbox[1])*0.15))/1000
      scale_xy <- .getPosition(scale.pos, inset=scale.inset)
      scale.meters <- scale.km*1000
      .scalebar(d=scale.meters, xy=scale_xy, type="bar", divs=2, below="km", bar.height=scale.height,
                       label=c(0, scale.km/2, scale.km), lwd=0.2, cex=scale.cex, bar.lwd=0.2)
    }

    #set network variables
    nanimals <- node_individuals[[g]]
    colnames(nanimals)[1] <- "site"
    nanimals <- nanimals[match(site_coords@data$site, nanimals$site), "nids"]
    map_bbox <- raster::extent(map_extent)
    node_size <- .rescale_vertex_igraph(nanimals, minmax.relative.size=c(nodes.size[1], nodes.size[2]))

    # set node names
    if(!is.null(id.metadata)){
      node_titles <- paste0(node_individuals[[g]]$site, " (", node_tagged[[g]]$tagged, ")",
      "\n", "n=", node_individuals[[g]]$nids)
      node_titles <- gsub(" (0)", "", node_titles, fixed=T)
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
    } else {
      node_classes <- stats::aggregate(data[,color.nodes.by], by=list(data[,spatial.col]), function(x) names(table(x))[which.max(table(x))])
      colnames(node_classes) <- c("node", "class")
      node_classes$class <- factor(node_classes, levels=levels(data[,color.nodes.by]))
      vertex_color <- nodes.color[node_classes$class]
    }
    vertex_color <- adjustcolor(vertex_color, alpha.f=nodes.alpha)

    # set nodes position
    if(repel.nodes==TRUE){
      radii <- (par("usr")[2]-par("usr")[1])*c(0.035, 0.07)
      radii <- plotrix::rescale(node_size, newrange=radii)
      repel_data <- cbind(site_coords@coords, "radius"=radii*repel.buffer)
      new_coords <- packcircles::circleRepelLayout(repel_data, xlim=map_bbox[1:2], ylim=map_bbox[3:4],
                                                   sizetype=c("radius"), wrap=F, weights=1)$layout
      site_coords@coords <- as.matrix(new_coords[,1:2])
      site_coords@data[,2:3] <- new_coords[,1:2]
    }


    # plot network
    network <- igraph::graph_from_data_frame(edges, directed=T, vertices=site_coords@data)
    igraph::V(network)$name <- node_titles
    igraph::E(network)$width <- .rescale(edges$value, from=range(all_movements), to=edge.width)
    igraph::plot.igraph(network, layout=as.matrix(site_coords@data[,2:3]), add=T, rescale=F,
                        edge.color=edge.color[g], edge.curved=edge.curved,
                        vertex.size=node_size, vertex.color=vertex_color,
                        vertex.label.family="Helvetica", vertex.label.color=nodes.label.color, vertex.label.cex=nodes.label.cex,
                        edge.arrow.size=edge.arrow.size, edge.arrow.width=edge.arrow.width,
                        edge.label=paste0("(",edges$value,")\n"), edge.label.family="Helvetica",
                        edge.label.cex=edge.label.cex, edge.label.color=edge.label.color, edge.label.font=edge.label.font, ...)
  }

  #reset par
  par(mfrow=c(1,1), mar=c(5,4,4,2) + 0.1)

}

#######################################################################################################
#######################################################################################################
#######################################################################################################
