#######################################################################################################
# Function to plot retention/attrition curve(s) #################################################################
#######################################################################################################

#' Plot retention curve(s)
#'
#' @description Function used to plot retention rate(s) at each one of the tagging sites.
#' It calculates the percentage of individuals that occur within their original tagging location over
#' a given time frame, relatively to the total number of released individuals. The algorithm considers the last detection
#' of each animal in that particular site (independently of the presences/absence periods or
#' migrations in-between).
#'
#' @param data A data frame containing binned animal detections.
#' @param site.col Should match the name of a variable in the supplied data containing
#' similar categories as those used in the tagging sites (e.g. pointing to specific receivers or broader areas/locations).
#' Defaults to "habitat".
#' @param id.metadata A data frame containing the tagging location of each animal ('tagging_location' column),
#'  as well as an optional (continuous) variable to generate a color scale.
#' @param id.col Name of the column containing animal identifications. Defaults to "ID".
#' @param tagdate.col Name of the column in id.metadata containing animal tagging/release dates in
#' POSIXct format. Defaults to 'tagging_date'.
#' @param tagsite.col Name of the column in id.metadata containing animal tagging/release locations.
#' Defaults to 'tagging_location'.
#' @param color.by Optional. Variable contained in the id.metadata table, used to generate
#' a color gradient (e.g. length).
#' @param color.pal Optional. Color palette used for the color gradient.
#' @param aggregate.fun Function used to aggregate numeric values if a color.by argument is supplied.
#' Defaults to "mean".
#' @param same.scale Logical. If TRUE, uses the same color scale for all panels. Defaults to FALSE.
#' @param background.col Background color of the plot. Defaults to "gray96".
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @param von.bertalanffy Optional. In case 'color.by' represents individual lengths measured at tagging,
#' it predicts new lengths at each time frame based on a Von Bertalanffy Growth model.
#' @param VBGF.params A list containing the required parameters to fit the  Von Bertalanffy Growth model.
#' e.g.: 'Linf', 'K' and 't0'. Check \code{\link[TropFishR]{VBGF}} function for further details.
#' @export


plotRetention <- function(data, id.metadata, site.col="habitat", id.col="ID", tagdate.col="tagging_date",
                          tagsite.col="tagging_location", color.by=NULL, color.pal=NULL, aggregate.fun="mean", same.scale=F,
                          background.col="grey96", cols=1, von.bertalanffy=F, VBGF.params=NULL) {


  ############################################################################
  ## Initial checks ##########################################################
  ############################################################################

  cat("Generating retention curve(s)\n")

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # check if data contains site.col
  if(!site.col %in% colnames(data)){
      stop("site.col column not found")
  }

  # check if id.metadata contains id.col
  if(!id.col %in% colnames(id.metadata)){
    stop("ID column not found in id.metadata. Please assign the correct column name with 'id.col'")
  }

  # check id.col format
  if(class(id.metadata[,id.col])!="factor"){
    cat("Converting metadata ids to factor\n")
    id.metadata[,id.col] <- as.factor(id.metadata[,id.col])
  }else{
    id.metadata[,id.col] <- droplevels(id.metadata[,id.col])
  }

  # get unique IDs
  ids <- levels(id.metadata[,id.col])
  data <- data[data[,id.col] %in% ids,]
  data[,id.col] <- factor(data[,id.col], levels=ids)

  # check if metadata contains tagdate.col
  if(!tagdate.col %in% colnames(id.metadata)){
    stop("tagdate.col column not found in id.metadata")
  }

  # check if tagdate.col is in the right format
  if(class(id.metadata[,tagdate.col])[1]!="POSIXct"){
    stop("Please convert the tagging dates to POSIXct format")
  }

  # check if metadata contains tagsite.col
  if(!tagsite.col %in% colnames(id.metadata)){
    stop("tagsite.col column not found in id.metadata")
  }

  # check tagsite.col format
  if(class(id.metadata[,tagsite.col])!="factor"){
    id.metadata[,tagsite.col] <- as.factor(id.metadata[,tagsite.col])
  }

  # check von bertalanffy arguments
  if(von.bertalanffy==T){
    if(is.null(VBGF.params)){
      stop("Please supply VBGF.params when von.bertalanffy is set to true")
    }
    cat(paste0("Applying Von Bertalanffy Growth curve to predict lengths at each time frame, based on the n\u00ba days passed since tagging. Assuming that the 'color.by' variable represents fish lengths\n"))
    if(VBGF.params$Linf<max(id.metadata$length, na.rm=T)){
      cat("Warning: Some individuals are larger than the supplied Linf (infinite length for investigated species in cm)\n")
    }
  }

  # set a color palette if not defined
  if(is.null(color.pal)){
    color.pal <- adjustcolor(rev(.viridis_pal(100)), alpha.f=0.9)
  }


  ##############################################################################
  # Prepare data ###############################################################

  data_individual <- split(data, f=data[,id.col], drop=F)
  data_individual <- mapply(getDaysAtLiberty, data=data_individual, tagdate=id.metadata[,tagdate.col], SIMPLIFY=F)
  data <- do.call("rbind", data_individual)

  ids_by_site <- split(id.metadata[,id.col], f=id.metadata[,tagsite.col])
  sample_by_site <- unlist(lapply(ids_by_site, length))
  data_attrition <- lapply(ids_by_site, function(x) data[data[,id.col] %in% x,])
  data_attrition <- mapply(function(data,site){data[data[,site.col]==site,]}, data=data_attrition, site=names(data_attrition), SIMPLIFY=F)
  data_attrition <- lapply(data_attrition, function(x) stats::aggregate(x$days_post_tag, by=list(x[,id.col]), max))
  data_attrition <- lapply(data_attrition, function(x) {colnames(x)<-c("ID", "days_post_tag"); return(x)})
  if(!is.null(color.by)){
    data_attrition <- lapply(data_attrition, function(x) plyr::join(x, id.metadata[,c("ID", color.by)], by="ID", type="left"))
  }

  data_attrition <- lapply(data_attrition, function(x) x[order(x$days_post_tag),])
  max_days <- max(unlist(lapply(data_attrition, function(x) max(x$days_post_tag))))
  day_seq <- seq(0, max_days, by=1)
  days_seq <- lapply(data_attrition, function(x) seq(0, max(x$days_post_tag), by=1))


  ##############################################################################
  # Calculate dispersal rate(s) ##################################################

  attritions <- list()
  for(i in 1:length(data_attrition)){
    attrition <- sapply(days_seq[[i]], function(d) length(which(data_attrition[[i]]$days_post_tag>=d)))
    attrition <- data.frame("days_post_tag"=days_seq[[i]], "individuals"=attrition)
    attrition$percentage <- attrition$individuals / sample_by_site[i] * 100
    attrition$site <- names(data_attrition)[i]
    if(!is.null(color.by)){
      agg.fun <- match.fun(aggregate.fun)
      if(von.bertalanffy==T){
        lengths_at_tagging <-  sapply(days_seq[[i]], function(d) data_attrition[[i]][,color.by][data_attrition[[i]]$days_post_tag>=d])
        ages_at_tagging <- sapply(lengths_at_tagging, function(l) TropFishR::VBGF(param=VBGF.params, L=l, na.rm=F))
        ages_at_period <- mapply(function(ages, days_after){ages+(days_after/365)}, ages=ages_at_tagging, days_after=days_seq[[i]], SIMPLIFY=F)
        ages_at_period <- lapply(ages_at_period, function(x) {x[is.nan(x)]<-100; return(x)})
        lengths_at_period <- sapply(ages_at_period, function(a) TropFishR::VBGF(param=VBGF.params, t=a, na.rm=F))
        attrition$aggr_stat <-  sapply(lengths_at_period, function(l) agg.fun(l, na.rm=T))
      }else{
        attrition$aggr_stat <-  sapply(days_seq[[i]], function(d) agg.fun(data_attrition[[i]][,color.by][data_attrition[[i]]$days_post_tag>=d], na.rm=T))
      }
    }
    attritions[[i]] <- attrition
  }

  # retrieve min and max values for the color.by var to uniformize color scale
  if(!is.null(color.by) & same.scale==T){
    var_range <- do.call("rbind", attritions)
    var_range <- range(var_range$aggr_stat, na.rm=T)
  }

  # assign colors
  for(i in 1:length(attritions)){
    if(!is.null(color.by)){
      if(same.scale==T){breaks <- seq(var_range[1], var_range[2], length.out=100)
      }else{breaks <- 100}
      attritions[[i]]$color <- color.pal[as.numeric(cut(attritions[[i]]$aggr_stat, breaks=breaks, include.lowest=T))]
    }else{
      attritions[[i]]$color <- color.pal
    }
  }


  ##############################################################################
  # Plot curves ###############################################################

  # set layout variables
  rows <- ceiling(length(data_attrition)/cols)
  par(mfrow=c(rows,cols), mar=c(4,4,2,6), mgp=c(3,0.6,0))
  max_percent <- max(unlist(lapply(attritions, function(x) max(x$percentage))))

  # generate a plot per site
  for(i in 1:length(attritions)){
    attrition <- attritions[[i]]
    plot(x=attrition$days_post_tag, y=attrition$percentage, type="n", axes=F, xlim=c(0,max_days), xlab="", ylab="", ylim=c(0, max_percent+2))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=background.col, border=NA)
    abline(h=seq(0,100, by=10), lwd=0.05, col="grey70")
    points(attrition$days_post_tag, y=attrition$percentage, type="l", lwd=0.6)
    points(attrition$days_post_tag, y=attrition$percentage, pch=16, cex=0.8, col=attrition$color)
    title(main=paste0(names(data_attrition)[i], " (n=", sample_by_site[i], ")"), cex.main=1)
    title(xlab="Days after tagging", cex.lab=1, line=2)
    title(ylab="N\u00ba individuals (%)", cex.lab=1, line=2.5)
    axis(1, at=pretty(c(0, max_days)), labels=pretty(c(0, max_days)), cex.axis=0.9)
    axis(2, at=pretty(c(0, max_percent)), labels=pretty(c(0, max_percent)), cex.axis=0.9, las=1)
    box()
    legend.title <- ifelse(von.bertalanffy==F, paste(aggregate.fun, color.by), paste("predicted", aggregate.fun, color.by))
    if(same.scale==T){
      length_labs <- pretty(var_range, min.n=4)
      length_labs <- length_labs[length_labs>=min(var_range) & length_labs<=max(var_range)]
      .colorlegend(col=color.pal, zlim=var_range, zval=length_labs,
                   posx=c(0.85, 0.88), posy = c(0.1, 0.9), main=legend.title,
                   main.cex=0.8, digit=1, cex=0.8)
    }else {
      length_labs <- pretty(c(min(attrition$aggr_stat), max(attrition$aggr_stat)), min.n=4)
      length_labs <- length_labs[length_labs>=min(var_range) & length_labs<=max(var_range)]
      .colorlegend(col=color.pal, zlim=range(attrition$aggr_stat), zval=length_labs,
                  posx=c(0.85, 0.88), posy = c(0.1, 0.9), main=legend.title,
                  main.cex=0.8, digit=1, cex=0.8)
    }

  }
}


################################################################################
# Auxiliary function ###########################################################

getDaysAtLiberty <- function(data, tagdate) {
  data$days_post_tag <- as.numeric(difftime(data$timebin, tagdate, units="day"))
  return(data)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
