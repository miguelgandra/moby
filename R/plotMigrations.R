#######################################################################################################
# Function to plot migrations network ##############################################################
#######################################################################################################

#' Plot map representing animal transitions/movements between receivers/sites.
#' If tagging information is supplied, the number of tagged individuals is indicated for
#' each location/site and the release site is included in the calculated transitions.
#'
#' @param data A data frame containing binned animal detections.
#' @param id.groups Optional. A list containing ID groups (e.g. different species or size-classes).
#' If supplied, a separate map will be plotted for each group.
#' @param aggregate.by Name of the variable/column containing the spatial information used to calculate transitions.
#' It can include receiver IDs in case all movements/transitions are of interest,
#' or contain (for example) location/habitat classes for broader-scale analyses.
#' @param land.shape Optional. A shape file containing coastlines.
#' @param land.color Optional. Color of land areas.
#' @param tags.info Optional. A data frame containing information on the tagging location of each animal.
#' Used to xxx.
#' @param background.layer A projected raster containing a variable to be displayed in the background
#' (e.g. bathymetry, temperature, etc.).
#' @param background.pal Color palette for the background layer. If a background layer is not supplied,
#' the first color is used as background.
#' @param color.nodes.by Variable used to color nodes (locations). Alternatively, it can be
#' set to 'detectios' (default) to distinguish between nodes with detections and nodes without detections or
#' to 'group' to assign different colors to different ID groups (one color per map).
#' @param nodes.color Color of the nodes where animals have been detected.
#' @param nodes.alpha Opacity of the nodes: 0=fully transparent, 1=fully opaque.
#' @param nodes.size A numeric vector of length 2. Represents the desired min and max vertex sizes relative
#' to the x-axis in terms of percentage (see details). See \code{\link[netdiffuseR]{rescale_vertex_igraph}}.
#' @param repel.nodes Used to automatically separate nodes if they overlap or are too clustered in space.
#' Defaults to False.
#' @param repel.buffer Controls the amount of space between nodes if repel.nodes is set to True.
#' @param scale.km Distance covered by the scale bar, in kilometers. If null, it is
#' automatically defined as 20% of the plot region.
#' @param scale.pos Position of the map scale, specified by keyword.
#' See \code{\link[base]{xy.coords}}.
#' @param scale.inset Inset distance(s) of the map scale from the margins
#' as a fraction of the plot region.
#' @param scale.height Controls the thickness of the scale bar. Defaults to 1.5.
#' @param id.col Name of the column containing animal IDs Defaults to 'ID'.
#' @param lon.col Name of the column containing longitude values. Defaults to 'lon'.
#' @param lat.col Name of the column containing latitude values. Defaults to 'lat'.
#' @param bbox.ext Amount of extension of the bounding box around the coordinates.
#' It can be numeric (same extension into all four directions), vector of two (first x, then y directional extension)
#' or vector of four (xmin, xmax, ymin, ymax extension). Default is .2 (extends the bounding box by 20%).
#' Only considered if a background layer is not supplied.
#' @param ... Arguments passed to the \code{\link[igraph]{plot.igraph}} function.
#' @export

plotMigrations <- function(data, id.groups=NULL, tag.info=NULL, aggregate.by,
                           land.shape=NULL, land.color="gray50", background.layer=NULL, background.pal="#F3F7F7",
                           color.nodes.by="detection", nodes.color=c("darkblue", "black"), nodes.alpha=0.8, nodes.size=c(0.035, 0.08),
                           repel.nodes=F, repel.buffer=1.1, scale.km=NULL, scale.pos="bottomright", scale.inset=0.05,
                           scale.height=1.5, id.col="ID", lon.col="lon", lat.col="lat", bbox.ext=0.2, ...) {


  #####################################################################################
  # First checks ######################################################################

  # print to console
  cat("Generating migrations network\n")

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # check if tagging data contains the aggregate.by column
  if(!is.null(tag.info)){
    if(!id.col %in% colnames(tag.info)){
      stop("'id.col' variable not found in the supplied tag.info")
    }
    if(!aggregate.by %in% colnames(tag.info)){
      stop("'aggregate.by' variable not found in the supplied tag.info")
    }
  }

  # check if data contains lon.col and lat.col
  if(any(!c(lon.col, lat.col) %in% colnames(data))){
    stop("Longitude/latitude columns not found. Please assign the correct column names with 'lon.col' and 'lat.col'")
  }

  # check if dataset coordinates are projected (meters)
  geographic_coords <-  all(data[,lon.col]>=(-180) & data[,lon.col]<=180 & data[,lat.col]>=(-90) & data[,lat.col]<=90)
  if(geographic_coords==T){
    if(is.null(land.shape) & is.null(background.layer)){
      stop("Coordinates seem to be in degrees/geographic format (unprojected).")
    }else{
      cat("Warning: Projecting coordinates  based on the CRS of the supplied land shape or background layer\n")
      if(!is.null(land.shape)) epsg_code <- land.shape@proj4string
      else if (!is.null(background.layer)) epsg_code <- background.layer@crs
    }
  }

  # check color.nodes.by argument
  if(!color.nodes.by %in% c("detection", "group", colnames(data))){
    stop("Please set a valid 'color.nodes.by' argument")
  }
  if(color.nodes.by=="group" & is.null(id.groups)){
    warning("'color.nodes.by' set to group but id.groups were not supplied. Only the first node color will be used")
  }

  # check if node size is of length 2
  if(length(nodes.size)!=2){
    stop("Please supply a min and max value for the node sizes.")
  }

  # check layer projections
  if(!is.null(background.layer) & !is.null(land.shape)){
    if(!raster::compareCRS(land.shape@proj4string, background.layer@crs))
      warning("Land.shape and background.layer seem to have different CRS projections\n")
  }

  # check if there color pallete is
  if(!is.null(background.layer) & length(background.pal)==1){
    cat("Warning: using a new color palette for the background layer. To override, please supply > 1 color\n")
    if(is.null(levels(background.layer))) background.pal <- adjustcolor(pals::ocean.deep(100), 0.75)
    if(!is.null(levels(background.layer))) background.pal <- pals::parula(length(levels(background.layer)[[1]][,2]))
  }


  #####################################################################################
  # Prepare data ######################################################################

  # set id groups
  if(is.null(id.groups)) {
    id.groups<-list(levels(data[,id.col]))
  } else {
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {"Some of the ID(s) in id.groups don't match the IDs in the data"}
  }

  # convert to factor
  if(!is.factor(data[,aggregate.by])){
    data[,aggregate.by] <- as.factor(data[,aggregate.by])
  }

  # split data
  n_groups <- length(id.groups)
  data_groupped <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])
  groupped_tags <- lapply(id.groups, function(x) tag.info[tag.info[,id.col] %in% x,])
  tagged_stats <- list()


  #####################################################################################
  # Calculate animal movements (transitions) ##########################################

  group_transitions <- list()
  group_detections <- list()
  group_individuals <- list()
  group_tagged <- list()

  # calculate metrics independently for each group
  for(g in 1:n_groups){

    #  data by individual
    data_individual <- split(data_groupped[[g]], f=data_groupped[[g]][,id.col], drop=T)

    # retrieve transitions
    transitions_list <- list()
    for(i in 1:length(data_individual)){
      transition_data <- data_individual[[i]][,c("timebin", aggregate.by)]
      transition_data <- as.character(transition_data[order(transition_data$timebin),-1])
      # include transitions from the release site if tagging info is available
      if(!is.null(tag.info)){
        id <- unique(data_individual[[i]][,id.col])
        tagging_site <- as.character(tag.info[tag.info[,id.col]==id, aggregate.by])
        transition_data <- c(tagging_site, transition_data)
        transition_data <- factor(transition_data, levels=levels(data[,aggregate.by]))
      }
      transition_matrix <- table(head(transition_data,-1), tail(transition_data,-1))
      transitions_list[[i]] <- transition_matrix
    }
    group_transitions[[g]] <- unclass(Reduce("+", transitions_list))

    # calculate nº detections on each site
    group_detections[[g]] <- table(data_groupped[[g]][,aggregate.by])

    # calculate nº individuals on each site
    nindividuals <- aggregate(data_groupped[[g]][,id.col], by=list(data_groupped[[g]][,aggregate.by]),
                             function(x) length(unique(x)), drop=F)
    colnames(nindividuals) <- c(aggregate.by, "nids")
    nindividuals$nids[is.na(nindividuals$nids)] <- 0
    group_individuals[[g]]  <- nindividuals

    # calculate proportion of tagged/visiting individuals on each site
    if(!is.null(tag.info)){
      detected_individuals <- split(data_groupped[[g]][,id.col], f=data_groupped[[g]][,aggregate.by])
      detected_individuals <- lapply(detected_individuals, function(x) as.character(unique(x)))
      tagged_individuals <- split(groupped_tags[[g]][,id.col], f=groupped_tags[[g]][,aggregate.by])
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
      group_tagged[[g]] <- tag_stats
    }
  }


  #####################################################################################
  # Setting geographic coordinates ####################################################

  # convert node coords to spatial points
  site_coords <- aggregate(data[,c(lon.col, lat.col)], by=list(data[,aggregate.by]), mean)
  colnames(site_coords)[1] <- "site"
  site_coords <- sp::SpatialPointsDataFrame(site_coords[,2:3], data=site_coords)

  # project coords if needed
  if(geographic_coords==T){
    projection(site_coords) <- CRS("+proj=longlat +datum=WGS84")
    site_coords <- sp::spTransform(site_coords, epsg_code)
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


  ####################################################################################
  # Plot network(s) ################################################################

  # set layout variables
  if(!is.null(background.layer)) {
    par(mfrow=c(n_groups,1), mar=c(1,2,1,6), oma=c(0,0,0,0))
  }else{
    par(mfrow=c(n_groups,1), mar=c(1,2,1,1), oma=c(0,0,0,0))
  }
  group_transitions <- lapply(group_transitions, function(x) {diag(x)<-0; return(x)})
  all_movements <- reshape2::melt(group_transitions)$value

  # loop through each one of the id.groups
  for(g in 1:n_groups){

    # create transition matrix (edges)
    edges <- reshape2::melt(group_transitions[[g]])
    colnames(edges) <- c("site1", "site2", "movements")
    edges <- edges[edges$movements>0,]
    edges <- edges[edges$site1!=edges$site2,]

    # plot bounding box
    plot(map_extent, col=background.pal[1], border=F, axes=F, xaxs="i", yaxs="i")
    box()

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

    # add map title and scale
    if(n_groups>1) {title(main=names(id.groups)[g], cex.main=0.9)}
    if(is.null(scale.km)){scale.km <- min(pretty((par("usr")[2]-par("usr")[1])*0.15))*1000}
    scale_xy <- moby:::getPosition(scale.pos, inset=scale.inset)
    scale.meters <- scale.km*1000
    moby:::scalebar2(d=scale.meters, xy=scale_xy, type="bar", divs=2, below="km", bar.height=scale.height,
                         label=c(0, scale.km/2, scale.km), lwd=0.2, cex=0.5, bar.lwd=0.2)

    #set network variables
    nanimals <- group_individuals[[g]]
    colnames(nanimals)[1] <- "site"
    nanimals <- nanimals[match(site_coords@data$site, nanimals$site), "nids"]
    map_bbox <- extent(template_raster)
    node_size <- netdiffuseR::rescale_vertex_igraph(nanimals, minmax.relative.size=c(nodes.size[1], nodes.size[2]))

    # set node names
    if(!is.null(tag.info)){
      node_titles <- paste0(group_tagged[[g]]$site, " (", group_tagged[[g]]$tagged, ")")
      node_titles <- gsub(" (0)", "", node_titles, fixed=T)
    }else{
      node_titles <- paste0(group_tagged[[g]]$site)
    }
    node_titles <- paste0(node_titles, "\n", "n=", group_tagged[[g]]$detected)


    # set node colors
    if(color.nodes.by=="detection"){
      vertex_color <- rep(nodes.color[1], length.out=length(nanimals))
      vertex_color[nanimals==0] <- nodes.color[2]
    } else if(color.nodes.by=="group"){
      vertex_color <- rep(nodes.color[g], length.out=length(nanimals))
    } else {
      node_classes <- aggregate(data[,color.nodes.by], by=list(data[,aggregate.by]), function(x) names(table(x))[which.max(table(x))])
      colnames(node_classes) <- c("node", "class")
      node_classes$class <- factor(node_classes, levels=levels(data[,color.nodes.by]))
      vertex_color <- nodes.color[node_classes$class]
    }
    vertex_color <- adjustcolor(vertex_color, alpha.f=nodes.alpha)

    # set nodes position
    if(repel.nodes==T){
      radii <- (par("usr")[2]-par("usr")[1])*c(0.035, 0.07)
      radii <- plotrix::rescale(node_size, newrange=radii)
      repel_data <- cbind(site_coords@coords, "radius"=radii*repel.buffer)
      new_coords <- packcircles::circleRepelLayout(repel_data, xlim=map_bbox[1:2], ylim=map_bbox[3:4],
                                                   sizetype=c("radius"), wrap=F, weights=1)$layout
      site_coords@coords <- as.matrix(new_coords[,1:2])
      site_coords@data[,2:3] <- new_coords[,1:2]
    }

    # plot network
    network <- igraph::graph.data.frame(edges, directed=T, vertices=site_coords@data)
    igraph::V(network)$name <- node_titles
    igraph::E(network)$width <- rescale2(edges$movements, range(all_movements), newrange=c(0.4, 3))
    igraph::plot.igraph(network, layout=as.matrix(site_coords@data[,2:3]), add=T, rescale=F,
                        edge.color="darkblue", edge.curved=0.5,
                        vertex.size=node_size, vertex.color=vertex_color,
                        vertex.label.family="Helvetica", vertex.label.color="white", vertex.label.cex=0.5,
                        edge.arrow.size=0.3, edge.arrow.width=1.5,
                        edge.label=paste0("(",edges$movements,")\n"), edge.label.family="Helvetica",
                        edge.label.cex=0.45,  edge.label.color="black", edge.label.font=2, ...)
  }

  #reset par
  par(mfrow=c(1,1), mar=c(5,4,4,2) + 0.1)


}



rescale2 <- function(x, oldrange, newrange){
  newrange[1] + (x-oldrange[1])/(oldrange[2]-oldrange[1])*(newrange[2]-newrange[1])
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
