##################################################################################################
## Plot maps with movement trajectories + home ranges ############################################
##################################################################################################

#' Plot maps with movement trajectories + home ranges
#'
#' @description Produces maps with inferred movement trajectories
#' (shortest paths in water) and estimated kernel home-ranges.
#' It can also overlay information such as bathymetry or
#' sea surface temperature, by plotting a raster in the background.
#'
#' @param data A data frame with animal positions (COAs), containing 'longitude' and 'latitude' columns.
#' @param kud.densities Kernel utilization areas ('estUD' class), as returned by \code{\link{calculateKUDs}}.
#' @param animal.tracks The output of \code{\link{calculateDistances}}.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param id.col Name of the column containing animal IDs. Defaults to "ID".
#' @param lon.col Name of the column containing projected longitudes. Defaults to "lon".
#' @param lat.col Name of the column containing projected latitudes. Defaults to "lat".
#' @param land.shape A projected shape file containing coastlines.
#' @param background.layer A projected raster containing a variable to be displayed in the background
#' (e.g. bathymetry, temperature, etc.).
#' @param background.pal Color palette for the background layer.
#' @param land.color Color of land areas.
#' @param discard.missing If true, only individuals with detections are included.
#' If false, plots of missing individuals are drawn only with land.shape and provided background.
#' @param kud.transparency Transparency of the kernel utilization areas (0-1).
#' @param scale.meters Distance covered by the scale bar, in meters. If null, it is
#' automatically defined as 20% of the plot region.
#' @param scale.pos Position of the scale bar, specified by keyword.
#' @param scale.inset Inset distance(s) from the margins as a fraction of the
#' plot region, relatively to the scale bar position.
#' @param same.scale Forces same KUD scale across all plots, allowing for
#' density comparison between individuals. If set to false, all y-axis are displayed.
#' Otherwise they are only displayed in the left-most plots to save space.
#' @param cols Number of columns in the panel (passed to the mfrow argument). Defaults to 3.
#' @export


################################################################################
# Main function - plot maps ####################################################

plotMaps <- function(data, kud.densities, animal.tracks=NULL, id.groups=NULL, id.col="ID",
                     lon.col="lon", lat.col="lat", land.shape, background.layer=NULL,
                     background.pal=NULL, background.color="#F3F7F7", land.color="gray50",
                     discard.missing=T, kud.transparency=0, scale.meters=NULL,  scale.pos="bottomright",
                     scale.inset=0.1, same.scale=F, kud.legend=T, cols=3) {


  ############################################################################
  ## Initial checks ##########################################################
  ############################################################################

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # check ID groups
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {"Some of the ID(s) in id.groups don't match the IDs in the data"}
    data <- data[data[,id.col] %in% unlist(id.groups),]
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
  }

  # check if data contains lon.col and lat.col
  if(any(!c(lon.col, lat.col) %in% colnames(data))){
    stop("Longitude/latitude columns not found. Please assign the correct column names with 'lon.col' and 'lat.col'")
  }

  # check if data contains lon.col and lat.col
  if(any(!c(lon.col, lat.col) %in% colnames(data))){
    stop("Longitude/latitude columns not found. Please assign the correct column names with 'lon.col' and 'lat.col'")
  }

  # check animal.tracks format
  if(!is.null(animal.tracks)){
    if(!is.null(animal.tracks$trajectories)){
      animal.tracks <- animal.tracks$trajectories
    }else{
      animal.tracks <- animal.tracks
      cat("Warning: a 'trajectories' list was not found in the supplied animal.tracks. Function will assume that animal.tracks already contains all linestring objects\n\n")
    }
  }

  # check kernel densities format (first attempt)
  if(any(grepl("kernel_density", names(kud.densities), fixed=T))){
    kud.densities <- kud.densities$kernel_density
  }

  # check kernel densities format (second attempt)
  if(class(kud.densities)!="estUDm"){
     if(class(kud.densities[[1]])=="estUDm"){
      ids <- unlist(lapply(kud.densities, names))
      kud.densities <- unlist(kud.densities)
      names(kud.densities) <- ids
      class(kud.densities) <- "estUDm"
    }else{
      stop("Supplied kud.densities appear to be in a wrong format. Please check calculateKUDs function")
    }
  }

  # check if data contains lon.col and lat.col
  if(all(!unique(data[,id.col]) %in% names(kud.densities))){
    stop("names of kud.densities do not match the supplied data ids")
  }

   # print to console
  cat("Plotting utilization area maps...\n")


  ############################################################################
  ## Set plot variables ######################################################
  ############################################################################

  if(discard.missing==T) {
    missing_individuals <- names(table(data[,id.col])[table(data[,id.col])<=5])
    data[,id.col] <- droplevels(data[,id.col])
  }

  # set background and scale variables
  if(is.null(background.pal)) {background.pal <- rev(palr::bathy_deep_pal(100)[c(30:100)])}
  bbox <- sf::st_bbox(sf::st_multipoint(cbind(data[,lon.col], data[,lat.col])))
  if(is.null(scale.meters)) {scale.meters <- pretty((bbox[3]-bbox[1])*0.2)[2]}
  scale_km <- scale.meters/1000

  # get range of KUD densities across all individuals (only used if same.scale=T)
  kud_range <- range(unlist(lapply(kud.densities, function(x) range(x$ud))))

  # set layout variables
  if(is.null(id.groups)){
    nindividuals <- nlevels(data[,id.col])
    rows <- ceiling(nindividuals/cols)
    nplots <- rows*cols
    if(nindividuals<nplots){
      animal_indexes <- c(1:nindividuals, rep(NA, nplots-nindividuals))
    }else{
      animal_indexes <- 1:nindividuals
    }
  }else{
    group_ids_selected <- lapply(id.groups, function(x) x[x %in% levels(data$ID)])
    group_numbers <- lapply(group_ids_selected,  length)
    group_rows <- lapply(group_numbers, function(x) ceiling(x/cols))
    rows <- do.call("sum", group_rows)
    group_plots <- lapply(group_rows, function(x) x*cols)
    plot_ids <- mapply(function(nids, nplots) {if(nids<nplots){c(1:nids, rep(NA, nplots-nids))}else{1:nids}},
                       nids=group_numbers, nplots=group_plots, SIMPLIFY=F)
    for(i in 2:length(plot_ids)){plot_ids[[i]]<-plot_ids[[i]]+max(plot_ids[[i-1]], na.rm=T)}
    plot_ids <- unlist(plot_ids, use.names=F)
  }
  plot_layout <- matrix(plot_ids, nrow=rows, ncol=cols, byrow=T)


  # arrange space for the common legends if required
  if((same.scale==T & kud.legend==T) | !is.null(background.layer)){
    empty_plot <- max(which(is.na(plot_ids)))
    if(is.na(empty_plot)){
      plot_layout <- rbind(plot_layout, NA)
      plot_ids <- unlist(plot_layout)
      empty_plot <- min(which(is.na(plot_ids)))
      rows <- rows+1
    }
    plot_ids[empty_plot] <- 9999
  }



  # set par
  if(kud.legend==T & same.scale==F){
    par(mfrow=c(rows, cols),  mar=c(0.4,0.4,0.4,4.5), oma=c(0,3,0,0))
  }else{
    par(mfrow=c(rows, cols),  mar=c(0.4,0.4,0.4,0.4), oma=c(0,3,0,0))
  }


  ############################################################################
  ## Generate maps ###########################################################
  ############################################################################

  # split data by individual
  data_individual <- split(data, f=data[,id.col])

  # set progress bar
  pb <- txtProgressBar(min=1, max=length(data_individual), initial=0, style=3)

  # iterate through each plot index
  for (i in plot_ids) {

    # if NA, add empty plot and go to next
    if(is.na(i)) {plot.new(); next}


    # plot legends
    if(i == 9999){

      plot.new()

      # add color legend of uniformized kud densities (if applicable)
      if(same.scale==T & kud.legend==T){
        kud_labs <- pretty(kud_range, min.n=4)
        kud_labs <- kud_labs[kud_labs>=min(kud_range) & kud_labs<=max(kud_range)]
        .colorlegend(col=color_pal, zlim=kud_range, zval=kud_labs,
                     posx=c(0.130, 0.155), posy=c(0.1, 0.8), main="KUD\ndensity\n\n",
                     main.cex=0.8, main.adj=0, lab.scientific=T, cex=0.8)
      }

      # add color legend of background layer (if available)
      if(!is.null(background.layer)) {
        layer_labs <- pretty(values(background.layer), min.n=4)
        layer_labs <- layer_labs[layer_labs>=min(kud_range) & layer_labs<=max(kud_range)]
        display_digits <- max(.decimalPlaces(layer_labs))
        .colorlegend(col=color_pal, zlim=kud_range, zval=layer_labs,
                     posx=c(0.200, 0.225), posy=c(0.1, 0.8), main="Layer values",
                     main.cex=0.8, digit=display_digits, main.adj=0, cex=0.8)
      }

      next
    }


    # update progress bar
    setTxtProgressBar(pb,i)

    # convert KUD to raster
    id <- as.character(levels(data[,id.col])[i])
    if(nrow(data_individual[[i]])==0 | !id %in% names(kud.densities)){
      r <- raster::raster(as(kud.densities[[1]],"SpatialPixelsDataFrame"))
      r[] <- NA
    } else {
      coords <- SpatialPoints(cbind(data_individual[[i]][,lon.col], data_individual[[i]][,lat.col]))
      r <- raster::raster(as(kud.densities[[id]],"SpatialPixelsDataFrame"))
    }

    # set color palette and zlims
    color_pal <- c(adjustcolor("gray96", alpha.f=0), rev(adjustcolor(terrain.colors(254), alpha.f=(1-kud.transparency))))
    if(same.scale==T){kud_scale <- kud_range}else{kud_scale <- range(values(r))}

    # plot maps
    suppressWarnings(image(r,  zlim=kud_scale, axes=F, col=adjustcolor("gray96", alpha.f=0), asp=1))
    if(!is.null(background.layer)){plot(background.layer, col=background.pal, legend=F, axes=F, add=T)}
    if(is.null(background.layer)){rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.color, border=NA)}
    suppressWarnings(image(r,  axes=F, col=color_pal, asp=1, add=T))
    legend("topleft", legend=id, inset=c(-0.08, 0) , cex=1.5, bty="n", text.font=2)
    plot(land.shape, col=land.color, border=NA, add=T)
    if(!is.null(animal.tracks) & !is.null(animal.tracks[[id]])){
      lines(animal.tracks[[id]], lwd=0.2, lty=2)
      points(coords, pch=21, bg="white", cex=0.7, col="gray15", lwd=0.25)
      points(coords[1,], pch=21, cex=0.8, bg=adjustcolor("blue", alpha.f=0.7) , lwd=0.25)
      points(coords[nrow(coords@coords),], pch=21, cex=0.8, bg=adjustcolor("red", alpha.f=0.7), lwd=0.25)
    }
    scale_xy <- .getPosition(scale.pos, inset=scale.inset)
    .scalebar(d=scale.meters, xy=scale_xy, type="bar", divs=2, below="km", label=c(0, scale_km/2, scale_km), lwd=0.2, cex=0.7)
    box(lty="solid", lwd=0.75, cex=0.5)

    if(same.scale==F & kud.legend==T){
      kud_labs <- pretty(kud_scale, min.n=4)
      kud_labs <- kud_labs[kud_labs>=min(kud_scale) & kud_labs<=max(kud_scale)]
      display_digits <- max(.decimalPlaces(kud_labs))
      .colorlegend(col=color_pal, zlim=kud_scale, zval=kud_labs,
                         posx=c(0.830, 0.855), posy = c(0.2, 0.7), main="KUD\ndensity\n\n",
                         main.cex=0.8, digit=display_digits, main.adj=0, lab.scientific=T, cex=0.6)
    }


  }

  # add id.group labels
  if(!is.null(id.groups)){
    label_pos <- rev(unlist(lapply(group_rows, function(x) x/2)))
    for(i in 2:length(group_rows)){label_pos[i]<-label_pos[i]+rev(group_rows)[[i-1]]}
    label_pos <- grconvertY(label_pos/rows, "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels=rev(names(id.groups)),
         srt=90, cex=2, font=2, xpd=NA)
  }

  #reset par and close progress bar
  par(mar=c(5, 4, 4, 2) + 0.1)
  close(pb)

  # print to console
  if(discard.missing==T){cat(paste0(length(missing_individuals), " individual(s) with < 5 detections discarded\n"))}

}


