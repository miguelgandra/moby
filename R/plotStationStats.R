#######################################################################################################
## Plot receiver statistics ###########################################################################
#######################################################################################################

#' Calculate and plot receivers' statistics

#' @description Calculates and plots receiver-based statistics, such as total and average detections,
#' the number of unique individuals detected, and co-occurrence events. This function generates grouped bar plots,
#' enabling side-by-side comparisons of specified statistics across stations or other broader spatial categories.
#' Additionally, users can compute and plot statistics separately for animal groups (e.g., different species or sexes)
#' and adjust a range of plot aesthetic options.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections with corresponding time-bins.
#' @param type Type of statistic to calculate/plot. If more than one type is supplied a grouped barplot
#' is generated. Possible values are "detections", "average detections", "individuals" and "co-occurrences".
#' @param id.groups Optional. A list containing ID groups, used to generate different
#' plots to each class (e.g. different species).
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' Only relevant to estimate co-occurrences. One of "within" (intra-group), "between" (inter-group) or "all" (both). Defaults to "all".
#' @param aggregate.by Column name for aggregating statistics (defaults to `station.col`).
#' @param color.pal A color palette for the bars; the length should match the number of `type` values provided.
#' @param number.locations Logical. If `TRUE`, replaces location names with numbers (for clearer display when names are long).
#' @param rotate.names Logical. If `TRUE`, rotates location names on the x-axis.
#' @param show.percentage Logical. If `TRUE`, percentages are displayed above each bar.
#' @param show.counts Logical. If `TRUE`, shows the count above each bar.
#' @param annot.threshold Numeric (optional). Minimum count required to display
#' labels on bars (e.g., `1` only labels bars with counts greater than 1).
#' Useful for reducing visual clutter in plots with many bars.
#' @param annot.cex Numeric. Font size for bar annotation text. Defaults to 0.7.
#' @param title.cex Numeric. Font size for the plot title. Defaults to 1.1.
#' @param title.pos Position of the title on the plot (by keyword). Defaults to "top".
#' @param title.inset Inset distance of the title from the specified position, as a fraction of the
#' plot region. Defaults to c(0, 0.02)
#' @param label.cex Numeric. Font size for axis labels. Defaults to 1.
#' @param axis.cex Numeric. Font size for axis text. Defaults to 0.9.
#' @param legend.pos Legend position.
#' @param legend.style Orientation of the legend, either "horizontal" or "vertical".
#' @param legend.cex Numeric. Font size for legend text. Defaults to 0.7.
#' @param cols Number of columns in the final panel (passed to the mfrow argument).

#' @export


#######################################################################################
## Main function ######################################################################
#######################################################################################

plotStationStats <- function(data,
                             type = c("detections", "co-occurrences"),
                             id.groups = NULL,
                             group.comparisons = "all",
                             id.col = getDefaults("ID"),
                             timebin.col = getDefaults("timebin"),
                             station.col = getDefaults("station"),
                             aggregate.by = NULL,
                             color.pal = NULL,
                             number.locations = FALSE,
                             rotate.names = FALSE,
                             show.percentage = TRUE,
                             show.counts = TRUE,
                             annot.threshold = NULL,
                             annot.cex = 0.7,
                             title.pos = "top",
                             title.inset = c(0, 0.02),
                             title.cex = 1.1,
                             label.cex = 1,
                             axis.cex = 0.9,
                             legend.pos = "topright",
                             legend.style = "vertical",
                             legend.cex = 0.7,
                             cols = 1) {


  ######################################################################################
  # Initial checks #####################################################################
  ######################################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # validate additional parameters
  errors <- c()
  # check type argument
  if(length(type)>4) errors <- c(errors, "Currently only a maximum of four 'types' can be defined.")
  if(any(!type %in% c("detections", "average detections", "individuals", "co-occurrences"))) errors <- c(errors, "Wrong 'type' argument. Please choose up to two of: 'detections', 'average detections', 'individuals' or 'co-occurrences'.")
  if(length(group.comparisons) > 1) errors <- c(errors, "The 'group.comparisons' argument must be a single value.")
  if(!is.null(id.groups) && !group.comparisons %in% c("within", "between", "all")) errors <- c(errors, "Wrong 'group.comparisons' argument. Please choose one of: 'within', 'between' or 'all'.")
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # save the current par settings and ensure they are restored upon function exit
  original_par <- par(no.readonly=TRUE)
  on.exit(par(original_par))

  # convert station column to a factor if not already
  if(!is.factor(data[,station.col])){
    warning("- Converting 'station.col' to factor.", call.=FALSE)
    data[,station.col] <- as.factor(data[,station.col])
  }

  # check for and set grouping variable
  if(is.null(aggregate.by)){
    aggregate.by <- station.col
  }else{
    # check for stations associated with more than one group
    station_group_counts <- table(data[[station.col]], data[[aggregate.by]])
    # find stations with more than one associated group
    multi_group_stations <- rownames(station_group_counts)[rowSums(station_group_counts>0)>1]
    # if there are any such stations, warn the user or stop the process
    if(length(multi_group_stations) > 0) {
      stop("Each station must be associated with a unique 'aggregate.by' group. The following stations are linked to multiple groups, which may cause errors:\n",
           paste(multi_group_stations, collapse = ", "), call. = FALSE)
    }
    # create station_habitat_df by subsetting data based on the specified columns
    station_groups <- unique(data[,c(station.col, aggregate.by)])
  }

  # convert aggregation column to factor if not already
  if(!is.factor(data[,aggregate.by])){
    warning("Converting 'aggregate.by' to factor.", call.=FALSE)
    data[,aggregate.by] <- as.factor(data[,aggregate.by])
  }


  ######################################################################################
  # Process data into a list based on groups ###########################################
  ######################################################################################

  if(!is.null(id.groups)){
    ids_table <- reshape2::melt(id.groups)
    colnames(ids_table) <- c("ID", "group")
    if(group.comparisons=="within"){
      data_list <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])
      data_list <- lapply(data_list, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
    }else if(group.comparisons=="between"){
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse=" <-> "))
      comb_ids <- apply(group_comb, 1, function(x) ids_table$ID[ids_table$group %in% x], simplify=F)
      data_list <- lapply(comb_ids, function(x) data[data[,id.col] %in% x,])
      data_list <- lapply(data_list, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
      names(data_list) <- group_comb_names
    }else{
      data_list <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])
      data_list <- lapply(data_list, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse=" <-> "))
      comb_ids <- apply(group_comb, 1, function(x) ids_table$ID[ids_table$group %in% x], simplify=F)
      data_group <- lapply(comb_ids, function(x) data[data[,id.col] %in% x,])
      data_group <- lapply(data_group, function(x) {x[,id.col] <- droplevels(x[,id.col]); return(x)})
      names(data_group) <- group_comb_names
      data_list <- c(data_list, data_group)
    }
  }else{
    data_list <- list(data)
    names(data_list) <- "All"
    ids_table <- data.frame("ID"=levels(data[,id.col]), "group"=1)
  }


  ######################################################################################
  # Calculate stats ####################################################################
  ######################################################################################

  # initialize list to hold group results
  group_results <- list()

  # iterate over each group
  for(i in 1:length(data_list)){

    # initialize lists to hold counts and frequencies for each type
    counts_list <- list()
    freqs_list <- list()

    # extract current subset of data
    group_data <- data_list[[i]]

    # count detections per location
    if(any(type=="detections")) {
      station_detections <- data.frame(rbind(table(group_data[,aggregate.by])), check.names=F)
      counts_list <- append(counts_list, list("detections"=station_detections))
      freqs_list <- append(freqs_list, list("detections"=station_detections/sum(station_detections)))
    }

    # compute average detection frequency per location
    if(any(type=="average detections")) {
      station_mean_detections <- by(group_data, group_data[,id.col], FUN=function(x) table(x[,aggregate.by]))
      station_mean_detections <- do.call("rbind", station_mean_detections)
      station_mean_detections <- station_mean_detections/rowSums(station_mean_detections)
      station_mean_detections <- colMeans(station_mean_detections)
      station_mean_detections <- data.frame(t(station_mean_detections), check.names=F)
      counts_list <- append(counts_list, list("average detections"=station_mean_detections))
      freqs_list <- append(freqs_list, list("average detections"=station_mean_detections))
    }

    # calculate number of individuals per location
    if(any(type=="individuals")) {
      nindividuals <- nlevels(group_data[,id.col])
      station_individuals <- stats::aggregate(group_data[,id.col], by=list(group_data[,aggregate.by]), function(x) length(unique(x)), drop=FALSE)
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
      # convert subsetted data to wide format
      group_table <- suppressWarnings(createWideTable(group_data, id.col=id.col, timebin.col=timebin.col, value.col=station.col, verbose=FALSE))
      group_table <- group_table[,-1]
      # replace IDs by group ID
      colnames(group_table) <- plyr::mapvalues(colnames(group_table), ids_table$ID, ids_table$group, warn_missing=FALSE)
      # calculate co-occurring group sizes for each timebin
      if(group.comparisons=="within" || (group.comparisons=="all" && !grepl("<->", names(data_list)[i], fixed=T))){
        co_occurrence_events <- unlist(apply(group_table, 1, .countJointDetections, comparison="within", groups=colnames(group_table)))
      }else if(group.comparisons=="between" || (group.comparisons=="all" && grepl("<->", names(data_list)[i], fixed=T))){
        co_occurrence_events <- unlist(apply(group_table, 1, .countJointDetections, comparison="between", groups=colnames(group_table)))
      }
      co_occurrence_events <- data.frame("station"=co_occurrence_events)
      # join co-occurrence data with grouping information
      if(aggregate.by!=station.col){
        co_occurrence_events <- merge(co_occurrence_events, station_groups, by.x="station", by.y=station.col)
        co_occurrence_events[,2] <- as.character(co_occurrence_events[,2])
      }
      co_occurrences <- table(co_occurrence_events[,aggregate.by])
      co_occurrences <- as.data.frame(co_occurrences)
      colnames(co_occurrences) <- c(aggregate.by, "count")
      # add stations/groups with zero co-occurrences to the results
      missing_levels <- as.character(levels(group_data[,aggregate.by])[!levels(group_data[,aggregate.by]) %in% co_occurrences[,aggregate.by]])
      if(length(missing_levels)>0) {
        missing_rows <- data.frame("levels"=missing_levels, "count"=rep(0, length(missing_levels)))
        colnames(missing_rows)[1] <- aggregate.by
        co_occurrences <- rbind(co_occurrences, missing_rows)
      }
      # reorder co_occurrences to match original levels in group_data for consistency
      co_occurrences <- co_occurrences[match(levels(group_data[,aggregate.by]), co_occurrences[,aggregate.by]),]
      # append final counts and frequencies to lists
      counts_list <- append(counts_list, list("co-occurrences"=co_occurrences$count))
      freqs_list <- append(freqs_list, list("co-occurrences"=co_occurrences$count/sum(co_occurrences$count)))
    }

    # consolidate results
    counts_list <- do.call("rbind", counts_list)
    freqs_list <- do.call("rbind", freqs_list)
    group_results[[i]] <- list("Counts"=counts_list, "Freqs"=freqs_list)
  }

  # assign group names
  names(group_results) <- names(data_list)



  ######################################################################################
  # Prepare plot matrices ##############################################################
  ######################################################################################

  # reorder rows if needed to match the supplied type order
  group_results <- lapply(group_results, function(x) lapply(x, function(y) {
    y <- as.data.frame(y)
    y[match(type, rownames(y)), , drop = FALSE]
    colnames(y) <- NULL
    return(y)
  }))

  # format decimal places
  formatVals <- function(x) {if(all(x%%1==0)){as.numeric(x)}else{as.numeric(sprintf("%.2f", x))}}
  group_results <- lapply(group_results, function(x) lapply(x, function(y){y[is.na(y)]<-0; return(y)}))
  group_results <- lapply(group_results, function(x) lapply(x, function(y) t(apply(y, 1, formatVals))))


  # assign variable display names
  var_names <- data.frame(type=c("detections", "average detections", "individuals", "co-occurrences"),
                          display=c("Detections", "Average Detection Frequency", "N\u00ba Individuals", "Co-Occurrence Events"))
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

  # configure layout dimensions for multi-panel plots
  nplots <- length(group_results)
  rows <- ceiling(nplots/cols)
  par(mfrow=c(rows, cols))
  mat <- matrix(1:nplots, nrow=rows)
  topright <- mat[1, ncol(mat)]
  bottom <- mat[nrow(mat),]

  # set default label rotation for x-axis labels
  x_las <- 1

  # adjust margins and determine the number of columns in the legend
  oma <- c(3,0,0,0)
  if(number.locations){
    station_labs <- 1:nlevels(data[,aggregate.by])
    if(length(station_labs)<14) {
      oma[4] <- 8
      legend_cols <- 1
    }else{
      oma[4] <- 14
      legend_cols <- 2
    }
  }else{
    # use factor levels as station labels if `number.locations` is FALSE
    legend_cols <- 1
    station_labs <- levels(data[,aggregate.by])
  }
  # rotate x-axis labels
  if(rotate.names==T){
    oma[1] <- 9
    x_las <- 2
  }


  ######################################################################################
  # Generate barplot(s) ################################################################
  ######################################################################################

  # set graphical parameters for barplots
  par(mar=c(1,5,1,2), oma=oma, mgp=c(3,0.75,0), lwd=0.6)

  # iterate over each group
  for(i in 1:length(group_results)){

    # access frequency and count data for current group
    station_freqs <- group_results[[i]]$Freqs
    colnames(station_freqs) <- levels(data[,aggregate.by])
    station_counts <- group_results[[i]]$Counts
    colnames(station_counts) <- levels(data[,aggregate.by])

    # generate blank barplot layout
    station_bars <- barplot(station_freqs, beside=TRUE, col=NA, border=NA, ylim=c(0, 1.2), space=c(0,0.4), axes=FALSE,
                            width=0.4, legend.text=FALSE, names.arg=rep("",ncol(station_freqs)))
    # add background color
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray98", border=NULL)

    # draw barplot
    station_bars <- barplot(station_freqs, beside=T, ylim=c(0, 1.2), space=c(0,0.4), axes=FALSE, col=color.pal,
                            width=0.4, cex.names=0.8, legend.text=FALSE,  axisnames=FALSE, lwd=0.5,
                            names.arg=rep("",ncol(station_freqs)), add=TRUE)

    # add legend if there is more than one type:
    if(length(type)>=2){
      .legend(legend.pos, legend=rownames(station_freqs), horiz=legend.style, fill=color.pal, bty="n",
                         y.intersp=1.4, x.intersp=0.8, cex=legend.cex, box.cex=c(1.2, 1))
     }

    # add plot titles
    if(!is.null(id.groups)){
      if(length(type)==1){legend(title.pos, inset=title.inset, legend=paste0(names(data_list)[i], "\n", rownames(station_freqs)), text.font=2, cex=title.cex, bty="n")
      }else{legend(title.pos, inset=title.inset, legend=names(data_list)[i], text.font=2, cex=title.cex, bty="n")}
    }else{
      if(length(type)==1){
        legend(title.pos, inset=title.inset, legend=rownames(station_freqs), text.font=2, cex=title.cex, bty="n")
      }
    }

    # add axes
    axis(2, at=seq(0, 1, by=0.2), labels=paste0(seq(0, 100, by=20),"%"), las=1, cex.axis=axis.cex)
    mtext(side=2, text="Frequency", line=3.5, cex=label.cex)
    # x-axis only included in the last plot
    if(i %in% bottom){
      axis(1, at=colMeans(station_bars), labels=station_labs, cex.axis=axis.cex, las=x_las)
      if(rotate.names==F) {mtext(side=1, text=tools::toTitleCase(aggregate.by), line=2.5, cex=label.cex)}
    }

    # prepare and format location frequencies and counts for display:
    display_counts <- paste0("n=", station_counts)
    display_freqs <- paste0("(", round(station_freqs*100), "%)")
    if(!is.null(annot.threshold)){
      display_counts[station_counts<annot.threshold] <- ""
      display_freqs[station_counts<annot.threshold] <- ""
    }

    # add frequency and count labels
    if(show.percentage){
      if(show.counts){
        text(x=station_bars, y=station_freqs+0.085, labels=display_counts, cex=annot.cex)
        text(x=station_bars, y=station_freqs+0.035, labels=display_freqs, cex=annot.cex)
      }else{text(x=station_bars, y=station_freqs+0.04, labels=display_freqs, cex=annot.cex)}
    }else{
      if(show.counts){text(x=station_bars, y=station_freqs+0.04, labels=display_counts, cex=annot.cex)}
    }

    # display legend with station labels if `number.locations` is TRUE:
    if(number.locations && i %in% topright) {
      usr <- par("usr")
      legend(x=usr[2], y=usr[4], legend=paste0(station_labs, ". ", colnames(station_freqs)),
             bty="n", y.intersp=1.2, cex=legend.cex, xpd=NA, ncol=legend_cols, text.width=usr[2]*0.16)
    }

    # draw box
    par(lwd=1)
    box()
  }

}


################################################################################
# Auxiliary function ###########################################################
################################################################################

#' Count nº of co-occurring groups in each timebin
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.countJointDetections <- function(row, comparison, groups) {
  row <- as.character(row)
  names(row) <- groups
  counts <- table(row)
  if(any(counts>1)){
    stations <- names(counts)[counts>1]
    groups <- sapply(stations, function(x) names(row)[which(row==x)], simplify=FALSE)
    if(comparison == "between"){
      valid_indices <- lapply(groups, function(x) any(length(unique(x))>1))
    }else if(comparison == "within"){
      valid_indices <- lapply(groups, function(x) any(table(x)>1))
    }
    return(stations[unlist(valid_indices)])
  }else{
    return(NULL)
  }
}


##################################################################################################
##################################################################################################
##################################################################################################
