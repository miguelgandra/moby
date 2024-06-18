#######################################################################################################
## Plot receiver statistics ###########################################################################
#######################################################################################################

#' Calculate and plot receivers' statistics

#' @description Calculates several receivers' statistics and generates corresponding barplots.
#' Implemented statistics include total detections, average detection frequency (per individual),
#' number of animals detected and number of co-occurrences.
#'
#' @param data A data frame containing animal detections with corresponding time-bins.
#' @param overlaps List containing counts of stations where animals co-occurred, as returned by
#' \code{\link{calculateOverlap}}. Only required if type contains "co-occurrences".
#' @param type Type of statistic to calculate/plot. If more than one type is supplied a grouped barplot
#' is generated. Possible values are "detections", "average detections", "individuals" and "co-occurrences".
#' @param id.col Name of the column containing animal identifications. Defaults to "ID".
#' @param site.col Should match the name of a variable in the supplied data containing the
#' spatial categories to plot (e.g. receivers or broader areas/locations). Defaults to "station".
#' @param number.stations Replace receiver names with numbers? Useful if receiver names are
#' too long to be displayed correctly.
#' @param rotate.labels Another alternative to address long receiver names.
#' @param show.percentage If true, percentages are displayed on top of each bar.
#' @param lab.threshold Optional. If defined, text labels are only only displayed to
#' frequencies above this threshold (e.g., if lab.threshold=0.05 count values are only displayed
#' in bars with frequencies >= 5%).
#' @param legend.pos Legend position.
#' @param legend.style Legend style (horizontal vs vertical).
#' @param id.groups Optional. A list containing ID groups, used to
#' generate different plots to each class (e.g. different species).
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' One of "within", "between" or "all". Defaults to "all".
#' @param color.pal Color palette for the bars. Should match the number of types defined (1 or 2).
#' @param cols Number of columns in the final panel (passed to the mfrow argument).

#' @export


#######################################################################################
## Main function ######################################################################
#######################################################################################

plotStationStats <- function(data, overlap=NULL, type=c("detections", "co-occurrences"), id.col="ID",
                             site.col="station", number.stations=F, rotate.labels=F, show.percentage=T,
                             lab.threshold=NULL, legend.pos="top", legend.style="horizontal",
                             id.groups=NULL, color.pal=NULL, group.comparisons="all", cols=1) {


  ######################################################################################
  # Initial checks #####################################################################
  ######################################################################################

  # check type argument
  if(length(type)>4) {
    stop("Currently only a maximum of four 'types' can be defined")}

  # check type argument
  if(any(!type %in% c("detections", "average detections", "individuals", "co-occurrences"))) {
    stop("Wrong 'type' argument. Please choose up to two of: 'detections', 'average detections', 'individuals' or 'co-occurrences'")}

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please specify the correct column using 'id.col'")}

  # check if id.col is of class factor (if not, convert)
  if(class(data[,id.col])!="factor"){
    cat("Converting ids to factor\n")
    data[,id.col] <- as.factor(data[,id.col])}

  # check if data contains site.col
  if(!site.col %in% colnames(data)){
    stop("Site column not found. Please specify the correct column using 'site.col'")
    }

  # check if site.col is of class factor (if not, convert)
  if(class(data[,site.col])!="factor"){
    cat("Converting site.col to factor\n")
    data[,site.col] <- as.factor(data[,site.col])}

  # check if site.col levels agree with the supplied overlap results
  if(!is.null(overlap)){
    unique(unlist(overlap$station_counts))
    overlap_values <- unique(unlist(lapply(overlap$station_counts, function(x) names(x))))
    if(!all(overlap_values %in% levels(data[,site.col]))){
      stop("Overlap values do not match with the levels in the supplied 'site.col' variable")
    }
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {"Some of the ID(s) in id.groups don't match the IDs in the data"}
    data <- data[data[,id.col] %in% unlist(id.groups),]
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))

    # check group.comparisons argument
    if(!is.null(group.comparisons)){
      if(any(!group.comparisons %in% c("within", "between", "all"))) {
        stop("Wrong 'group.comparisons' argument. Please choose one of: 'within', 'between' or 'all'")}
    }
  }



  ######################################################################################
  # Convert to list ####################################################################
  ######################################################################################

  if(!is.null(id.groups)){
    ids_table <- reshape2::melt(id.groups)
    colnames(ids_table) <- c("ID", "group")
    if(group.comparisons=="within"){
      data_list <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])
      data_list <- lapply(data_list, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
    }else if(group.comparisons=="between"){
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse="<->"))
      comb_ids <- apply(group_comb, 1, function(x) ids_table$ID[ids_table$group %in% x], simplify=F)
      data_list <- lapply(comb_ids, function(x) data[data[,id.col] %in% x,])
      data_list <- lapply(data_list, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
      names(data_list) <- group_comb_names
    }else{
      data_list <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])
      data_list <- lapply(data_list, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse="<->"))
      comb_ids <- apply(group_comb, 1, function(x) ids_table$ID[ids_table$group %in% x], simplify=F)
      data_group <- lapply(comb_ids, function(x) data[data[,id.col] %in% x,])
      data_group <- lapply(data_group, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
      names(data_group) <- group_comb_names
      data_list <- c(data_list, data_group)
    }
  }else{
    data_list <- list(data)
    ids_table <- data.frame("ID"=levels(data[,id.col]), "group"=1)
  }


  ######################################################################################
  # Calculate stats ####################################################################
  ######################################################################################

  group_results <- list()
  for(i in 1:length(data_list)){

    # initialize list
    counts_list <- list()
    freqs_list <- list()

    # grab subsetted data
    group_data <- data_list[[i]]

    # calculate nº detections per station
    if(any(type=="detections")) {
      station_detections <- data.frame(rbind(table(group_data[,site.col])), check.names=F)
      counts_list <- append(counts_list, list("detections"=station_detections))
      freqs_list <- append(freqs_list, list("detections"=station_detections/sum(station_detections)))
    }

    # calculate average detection frequency per station
    if(any(type=="average detections")) {
      station_mean_detections <- by(group_data, group_data[,id.col], FUN=function(x) table(x[,site.col]))
      station_mean_detections <- do.call("rbind", station_mean_detections)
      station_mean_detections <- station_mean_detections/rowSums(station_mean_detections)
      station_mean_detections <- colMeans(station_mean_detections)
      station_mean_detections <- data.frame(t(station_mean_detections), check.names=F)
      counts_list <- append(counts_list, list("average detections"=station_mean_detections))
      freqs_list <- append(freqs_list, list("average detections"=station_mean_detections))
    }

    # calculate number of individuals per station
    if(any(type=="individuals")) {
      nindividuals <- nlevels(group_data[,id.col])
      station_individuals <- stats::aggregate(group_data[,id.col], by=list(group_data[,site.col]), function(x) length(unique(x)), drop=F)
      colnames(station_individuals) <- c("station", "individuals")
      station_individuals <- data.frame(t(station_individuals), check.names=F)
      colnames(station_individuals) <- station_individuals[1,]
      station_individuals <- station_individuals[-1,]
      station_individuals[1,] <- as.numeric(station_individuals[1,])
      station_individuals[] <- lapply(station_individuals, function(x) as.numeric(x))
      counts_list <- append(counts_list, list("individuals"=station_individuals))
      freqs_list <- append(freqs_list, list("individuals"=station_individuals/nindividuals))
    }

    # calculate nº co-occurrences per location
    if(any(type=="co-occurrences")) {
      if(is.null(overlap)) {stop("Overlap results are required to calculate co-occurrences")}
      station_co_occurrences <- overlap$station_counts
      station_co_occurrences <- station_co_occurrences[!sapply(station_co_occurrences, is.null)]
      pair_ids <- sapply(names(station_co_occurrences), function(x) strsplit(x, split="<->", fixed=T))
      pair_ids <- lapply(pair_ids, function(x) plyr::mapvalues(x, ids_table$ID, ids_table$group, warn_missing=F))

      if(group.comparisons=="within" | (group.comparisons=="all" & !grepl("<->", names(data_list)[i], fixed=T))){
        pair_ids <- lapply(pair_ids, function(x) all(x==names(data_list)[i]))
      }else if(group.comparisons=="between" | (group.comparisons=="all" & grepl("<->", names(data_list)[i], fixed=T))){
        comp_groups <- unique(unlist(strsplit(names(data_list)[i], split="<->", fixed=T)))
        pair_ids <- lapply(pair_ids, function(x) all(comp_groups %in% x))
      }
      station_co_occurrences <- station_co_occurrences[which(unlist(pair_ids))]
      station_co_occurrences <- plyr::rbind.fill(lapply(station_co_occurrences, function(x) as.data.frame.matrix(t(x))))
      station_co_occurrences <- colSums(station_co_occurrences, na.rm=T)
      station_co_occurrences <- data.frame(t(station_co_occurrences), check.names=F)
      missing_stations <- as.character(levels(group_data[,site.col])[!levels(group_data[,site.col]) %in% colnames(station_co_occurrences)])
      if(length(missing_stations)>0) {
        missing_cols <- as.data.frame(matrix(0, ncol=length(missing_stations)))
        colnames(missing_cols) <- missing_stations
        station_co_occurrences <- cbind(station_co_occurrences, missing_cols)
        stations_order <- match(levels(group_data[,site.col]), colnames(station_co_occurrences))
        station_co_occurrences <- station_co_occurrences[,stations_order]
      }
      counts_list <- append(counts_list, list("co-occurrences"=station_co_occurrences))
      freqs_list <- append(freqs_list, list("co-occurrences"=station_co_occurrences/sum(station_co_occurrences)))
    }
      counts_list <- do.call("rbind", counts_list)
      freqs_list <- do.call("rbind", freqs_list)
      group_results[[i]] <- list("Counts"=counts_list, "Freqs"=freqs_list)
  }

  names(group_results) <- names(data_list)



  ######################################################################################
  # Prepare plot matrices ################################################################
  ######################################################################################

  # reorder rows if needed to match the supplied type order
  group_results <- lapply(group_results, function(x) lapply(x, function(y) y[match(type, rownames(y)),]))

  # format decimal places
  formatVals <- function(x) {if(all(x%%1==0)){as.character(x)}else{sprintf("%.2f", x)}}
  group_results <- lapply(group_results, function(x) lapply(x, function(y){y[is.na(y)]<-0; return(y)}))
  group_results <- lapply(group_results, function(x) lapply(x, function(y) t(apply(y, 1, formatVals))))
  group_results <- lapply(group_results, function(x) lapply(x, function(y) t(apply(y, 1, as.numeric))))

  # assign variable display names
  var_names <- data.frame(type=c("detections", "average detections", "individuals", "co-occurrences"),
                          display=c("Detections", "Average Detection Frequency", "Nº individuals", "Co-occurrences"))

  group_results <- lapply(group_results, function(x) lapply(x, function(y){
    rownames(y) <- plyr::mapvalues(rownames(y), var_names$type, var_names$display, warn_missing=F); return(y)}))


  # convert to matrix
  group_results <- lapply(group_results, function(x) lapply(x, as.matrix))



  ######################################################################################
  # Prepare plot and layout variables ##################################################
  ######################################################################################

  # set color palette
  if(is.null(color.pal)){
    if(length(type)==1){
      color.pal <- c("grey45")
    }else{
      color.pal <- rev(grey.colors(length(type)))
    }
  }

  # set legend orientation
  legend.style <- switch(legend.style, "horizontal"=T, "vertical"=F)


  # set layout
  nplots <- length(group_results)
  rows <- ceiling(nplots/cols)
  par(mfrow=c(rows, cols))

  x_las <- 1
  if(number.stations==T){
    station_labs <- 1:nlevels(data[,site.col])
    if(length(station_labs)<14){
      par(mar=c(4,4,2,8))
      legend_cols <- 1}
    if(length(station_labs)>=14){
      par(mar=c(4,4,2,14))
      legend_cols <- 2}}
  if(number.stations==F){
    par(mar=c(4,4,2,2))
    legend_cols <- 1
    station_labs <- levels(data[,site.col])}
  if(rotate.labels==T){
    par(mar=par("mar")+c(3,0,0,0))
    x_las <- 2}


  ######################################################################################
  # Plot barplot(s) ####################################################################

  par(mgp=c(3,0.75,0), lwd=0.6)
  for(i in 1:length(group_results)){

    station_freqs <- group_results[[i]]$Freqs
    station_counts <- group_results[[i]]$Counts

    # generate barplots
    station_bars <- barplot(station_freqs, beside=T, col=NA, border=NA, ylim=c(0, 1.1), space=c(0,0.4), axes=F,
                            width=0.4, legend.text=F, names.arg=rep("",ncol(station_freqs)))
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray98", border=NULL)
    station_bars <- barplot(station_freqs, beside=T, ylim=c(0, 1.1), space=c(0,0.4), axes=F, col=color.pal,
                            width=0.4, cex.names=0.8, legend.text=F,  axisnames=F, lwd=0.5,
                            names.arg=rep("",ncol(station_freqs)), add=T)

    # add legend if more than two types
     if(length(type)>=2){
      .legend(legend.pos, legend=rownames(station_freqs), horiz=legend.style, fill=color.pal, bty="n",
                         y.intersp=1.4, cex=0.7, box.cex=c(1.2, 1))
     }

    # add plot titles
    if(!is.null(id.groups)){
      if(length(type)==1){title(main=paste0(names(data_list)[i], "\n", rownames(station_freqs)), line=1, font=2, cex.main=0.9)
      }else{title(main=names(data_list)[i], line=1, font=2, cex.main=0.9)}
    }else{
      if(length(type)==1){title(main=rownames(station_freqs), line=1, font=2, cex.main=0.8)}
    }

    # add axis and labels
    axis(1, at=colMeans(station_bars), labels=station_labs, cex.axis=0.75, las=x_las)
    axis(2, at=seq(0, 1, by=0.2), labels=paste0(seq(0, 100, by=20),"%"), las=1, cex.axis=0.8)
    mtext(side=2, text="Frequency", line=3, cex=0.8)
    if(rotate.labels==F) {mtext(side=1, text=tools::toTitleCase(site.col), line=2.5, cex=0.8)}
    station_counts <- paste0("n=", station_counts)
    display_freqs <- paste0("(", round(station_freqs*100), "%)")
    if(!is.null(lab.threshold)){
      station_counts[station_freqs<lab.threshold] <- ""
      display_freqs[station_freqs<lab.threshold] <- ""
    }
    if(show.percentage==F){
      text(x=station_bars, y=station_freqs+0.04, labels=station_counts, cex=0.6)}
    if(show.percentage==T){
      text(x=station_bars, y=station_freqs+0.085, labels=station_counts, cex=0.6)
      text(x=station_bars, y=station_freqs+0.035, labels=display_freqs, cex=0.6)}
    if(number.stations==T) {
      usr <- par("usr")
      legend(x=usr[2]*1.00, y=usr[4], legend=paste0(station_labs, ". ", colnames(station_freqs)),
             bty="n", y.intersp=1.3, cex=0.6, xpd=T, ncol=legend_cols, text.width=usr[2]*0.16 )}
    par(lwd=1)
    box()
  }

}


##################################################################################################
##################################################################################################
##################################################################################################
