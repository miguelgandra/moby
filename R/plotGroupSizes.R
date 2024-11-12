#######################################################################################################
## Plot frequency distributions for co-occurring group sizes ##########################################
#######################################################################################################

#' Plot frequency distributions for co-occurring group sizes

#' @description Function to calculate the frequency distribution of the
#' nº of animals co-occurring at each time-bin, either for the
#' entire study duration or split by a given variable (e.g. diel phase).
#'
#' @param table A data frame containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding to the receiver/station with the
#' highest number of detections), as returned by \code{\link{createWideTable}}.
#' @param id.groups Optional. A list containing ID groups, used to generate different
#' plots to each class (e.g. different species).
#' @param split.by Optional. If defined, frequencies are calculated separately for each level
#' of this variable. If NULL, frequencies are calculated for the entire study duration.
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' One of "within", "between" or "all". Defaults to "all".
#' @param levels.order Numeric vector with the preferred order of the levels defined by the 'split.by' variable
#' (in the stacked barplot).
#' @param bar.colors A vector of colors for the bars. Nº of colors should match
#' the nº of levels in the 'split.by' variable.
#' @param title.cex Numeric. Font size for the plot title. Defaults to 1.1.
#' @param label.cex Numeric. Font size for axis labels. Defaults to 1.
#' @param axis.cex Numeric. Font size for axis text. Defaults to 0.9.
#' @param annot.cex Numeric. Font size for bar annotation text. Defaults to 0.7.
#' @param legend.pos Legend position.
#' @param legend.style Orientation of the legend, either "horizontal" or "vertical".
#' @param legend.cex Numeric. Font size for legend text. Defaults to 0.7.
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @export


######################################################################################
# Main function ######################################################################
######################################################################################

plotGroupSizes <- function(table,
                           id.groups = NULL,
                           split.by = NULL,
                           group.comparisons = "all",
                           bar.colors = NULL,
                           levels.order = NULL,
                           title.cex = 1.1,
                           label.cex = 1,
                           axis.cex = 0.9,
                           annot.cex = 0.7,
                           legend.pos = "topright",
                           legend.style = "vertical",
                           legend.cex = 0.7,
                           cols = 1) {


  ######################################################################################
  # Initial checks #####################################################################
  ######################################################################################

  # validate parameters
  errors <- c()
  if (!is.data.frame(table)) errors <- c(errors, "The 'table' argument must be a data frame.")
  if(!c("ids") %in% names(attributes(table))) errors <- c(errors, "The supplied table does not seem to be in the wide format. Please use the output of the 'createWideTable' function.")
  if(!is.null(split.by) && !split.by %in% colnames(table)) errors <- c(errors, paste0("The specified split.by variable '", split.by, "' was not found in the provided table."))
  if(!is.null(id.groups)){
    if(!group.comparisons %in% c("within", "between", "all")) errors <- c(errors, "Wrong 'group.comparisons' argument. Please choose one of: 'within', 'between' or 'all'.")
    if(any(duplicated(unlist(id.groups)))) errors <- c(errors, "Repeated ID(s) in id.groups.")
  }
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # print message
  .printConsole("Generating plot of frequency distributions for co-occurring group sizes")

  # save the current par settings and ensure they are restored upon function exit
  original_par <- par(no.readonly=TRUE)
  on.exit(par(original_par))

  # extracts the ID attributes from the table to specify detection columns
  ids <- attributes(table)$ids
  id_cols <- which(colnames(table) %in% ids)

  # select only the columns in 'table' that correspond to detections
  table_detections <- table[,id_cols]

  # checks if the detection columns are character or factor (corresponding to the receiver ID)
  if(all(!sapply(id_cols, function(x) inherits(table[,x], c("character","factor"))))){
    stop("Please supply a table with values corresponding to the receiver with the highest number of detections.", call.=FALSE)
  }

  # reorder ID levels based on user-defined 'id.groups', if provided
  if(!is.null(id.groups)){
    if(any(!unlist(id.groups) %in% ids)) {warning("Some of the ID(s) in id.groups don't match the IDs in the data.", call.=FALSE)}
    if(any(!ids %in% unlist(id.groups))) {warning("Some of the ID(s) in the supplied table were not found within id.groups and will be discarded", call.=FALSE)}
    table_detections <- table_detections[,colnames(table_detections) %in% unlist(id.groups)]
  }

  # check if 'split.by' column is specified for data splitting
  if(!is.null(split.by)){
    # print to console
    cat(paste0("Splitting data by '", split.by, "'\n"))
    # convert column to factor if it isn't one already
    if(!is.factor(table[,split.by])){
      warning("Converting split.by to factor.", call.=FALSE)
      table[,split.by] <- as.factor(table[,split.by])
    }
  } else{
    # if 'levels.order' is provided without 'split.by', issue a warning and reset 'levels.order' to default
    if(!is.null(levels.order)) {
      warning("'levels.order' is defined but will be ignored, as no 'split.by' variable is specified.", call.=FALSE)
      levels.order <- 1
    }
  }

  # validate the 'bar.colors' palette based on 'split.by' settings
  if(!is.null(bar.colors)){
    # if 'split.by' is specified, check if color count matches the number of levels
    if(!is.null(split.by)){
      if(length(bar.colors)!=nlevels(table[,split.by])) warning("The number of colors in 'bar.colors' does not match the number of levels in 'split.by'.", call.=FALSE)
    }else{
    # if 'split.by' is not specified, issue a warning if more than one color is supplied
      if(length(bar.colors)>1) warning("Multiple colors provided in 'bar.colors', but 'split.by' is not specified. Only the first color will be used.", call.=FALSE)
    }
  }

  ######################################################################################
  # Process data into a list based on groups ###########################################
  ######################################################################################

  if(!is.null(id.groups)){
    ids_table <- reshape2::melt(id.groups)
    colnames(ids_table) <- c("ID", "group")
    if(group.comparisons=="within"){
      data_list <- lapply(id.groups, function(x) table_detections[,colnames(table_detections) %in% x])
      names(data_list) <- names(id.groups)
    }else if(group.comparisons=="between"){
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse=" <-> "))
      comb_ids <- apply(group_comb, 1, function(x) ids_table$ID[ids_table$group %in% x], simplify=F)
      data_list <- lapply(comb_ids, function(x) table_detections[,colnames(table_detections) %in% x])
      names(data_list) <- group_comb_names
    }else{
      data_list <- lapply(id.groups, function(x) table_detections[,colnames(table_detections) %in% x])
      names(data_list) <- names(id.groups)
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse=" <-> "))
      comb_ids <- apply(group_comb, 1, function(x) ids_table$ID[ids_table$group %in% x], simplify=F)
      data_group <- lapply(comb_ids, function(x) table_detections[,colnames(table_detections) %in% x])
      names(data_group) <- group_comb_names
      data_list <- c(data_list, data_group)
    }
  }else{
    data_list <- list(table_detections)
    names(data_list) <- "complete"
    ids_table <- data.frame("ID"=ids, "group"=1)
  }


  ######################################################################################
  # Calculate and plot frequencies #####################################################
  ######################################################################################

  # initialize variable
  list_results <- list()

  # iterate through each group
   for (i in 1:length(data_list)){

     # extract current subset of data
     group_table <- data_list[[i]]


     ####################################################################
     # calculate overall frequencies of co-occurring group sizes ########

     if(is.null(split.by)) {

       # replace IDs by group ID
       colnames(group_table) <- plyr::mapvalues(colnames(group_table), from=ids_table$ID, to=ids_table$group, warn_missing=FALSE)

       # calculate co-occurring group sizes for each timebin
       if(group.comparisons=="within" | (group.comparisons=="all" & !grepl("<->", names(data_list)[i], fixed=TRUE))){
         overlapping_animals <- unlist(apply(group_table, 1, .countOverlappingIDs, type="within", groups=colnames(group_table)))
       }else if(group.comparisons=="between" | (group.comparisons=="all" & grepl("<->", names(data_list)[i], fixed=TRUE))){
         overlapping_animals <- unlist(apply(group_table, 1, .countOverlappingIDs, type="between", groups=colnames(group_table)))
       }
       overlapping_animals <- overlapping_animals[overlapping_animals>0]
       overlapping_animals <- table(overlapping_animals)
       overlapping_animals <-  as.data.frame(overlapping_animals)
       overlapping_animals <- data.frame("total"=overlapping_animals[,-1], row.names=overlapping_animals[,1])
       overlapping_freqs <- overlapping_animals/sum(overlapping_animals)*100
       list_results[[i]] <- list("freqs"=overlapping_freqs, "counts"=overlapping_animals)
     }

     ####################################################################
     # calculate frequency of co-occurring animals results per split.by variable
     if(!is.null(split.by)) {

       # recover remaining columns
       group_table <- cbind(group_table, table[,split.by])
       colnames(group_table)[ncol(group_table)] <- split.by

       # subset data by 'split.by' variable
       data_subset <- split(group_table, f=group_table[,split.by])
       data_subset <- lapply(data_subset, function(x) x[,-ncol(x)])
       data_subset <- data_subset[order(match(names(data_subset), levels(table[,split.by])))]

        # replace IDs by group ID
       data_subset <- lapply(data_subset, function(x){colnames(x) <- plyr::mapvalues(colnames(x), from=ids_table$ID, to=ids_table$group, warn_missing=FALSE); return(x)})

       # calculate co-occurring group sizes for each timebin
       if(group.comparisons=="within" | (group.comparisons=="all" & !grepl("<->", names(data_list)[i], fixed=TRUE))){
         overlapping_animals <- reshape2::melt(lapply(data_subset, function(x) unlist(apply(x, 1, .countOverlappingIDs, type="within", groups=colnames(x)))))
       }else if(group.comparisons=="between" | (group.comparisons=="all" & grepl("<->", names(data_list)[i], fixed=TRUE))){
         overlapping_animals <- reshape2::melt(lapply(data_subset, function(x) unlist(apply(x, 1, .countOverlappingIDs, type="between", groups=colnames(x)))))
       }
       colnames(overlapping_animals) <- c("overlapping_animals", "group")
       overlapping_animals <- overlapping_animals[overlapping_animals$overlapping_animals>0,]
       overlapping_animals <- as.data.frame.matrix(table(overlapping_animals))
       missing_levels <- levels(table[,split.by])[!levels(table[,split.by]) %in% colnames(overlapping_animals)]
       if(length(missing_levels)>0){
         overlapping_animals[, missing_levels] <- 0
       }
       overlapping_animals <- overlapping_animals[,match(colnames(overlapping_animals), levels(table[,split.by]))]
       overlapping_freqs <- (overlapping_animals/sum(overlapping_animals))*100
       list_results[[i]] <- list("freqs"=overlapping_freqs, "counts"=overlapping_animals)
     }
   }


  ######################################################################################
  # Set layout variables ###############################################################
  ######################################################################################

  # set layout
  nplots <- length(data_list)
  rows <- ceiling(nplots/cols)
  par(mfrow=c(rows, cols), mgp=c(3,0.75,0))

  # set legend orientation
  legend.style <- switch(legend.style, "horizontal"=TRUE, "vertical"=FALSE)

  # level names
  if(!is.null(split.by)){
    group_labels <-  unique(unlist(lapply(list_results, function(x) colnames(x$counts))))
    group_labels <- group_labels[match(levels(table[,split.by]), group_labels)]
    group_labels <- group_labels[!is.na(group_labels)]
  }

  # set color palette
  if(is.null(bar.colors)) {
    if(!is.null(split.by)){bar.colors <- colorRampPalette(c("gray30", "gray90"))(length(group_labels))
    }else{bar.colors <- "grey55"}
  }

  # group sizes
  size_labels <- as.integer(unique(unlist(lapply(list_results, function(x) rownames(x$counts)))))
  size_labels <- size_labels[order(size_labels)]


  ######################################################################################
  # Generate barplots ##################################################################
  ######################################################################################

  for (i in 1:length(data_list)){

    # retrieve plot data
    overlapping_freqs <- list_results[[i]]$freqs
    overlapping_animals <- list_results[[i]]$counts

    # define plot levels
    group_levels <- 1
    if(!is.null(split.by)){
      plot_levels <- colnames(overlapping_animals)
      group_levels <- which(group_labels %in% plot_levels)
      if(!is.null(levels.order)){
        overlapping_freqs <- overlapping_freqs[,levels.order]
        overlapping_animals <- overlapping_animals[,levels.order]
      }
    }

    # standardize group sizes (needed if multiple plots are being generated - per id.group)
    group_sizes <- as.numeric(rownames(overlapping_animals))
    missing_sizes <- size_labels[!size_labels %in% group_sizes]
    if(length(missing_sizes)>0){
      overlapping_animals[nrow(overlapping_animals)+length(missing_sizes),] <- NA
      overlapping_animals[is.na(overlapping_animals)] <- 0
      overlapping_freqs[nrow(overlapping_freqs)+length(missing_sizes),] <- NA
      overlapping_freqs[is.na(overlapping_freqs)] <- 0
      rownames(overlapping_animals)[tail(1:nrow(overlapping_animals), 2)] <- missing_sizes
      rownames(overlapping_freqs)[tail(1:nrow(overlapping_animals), 2)] <- missing_sizes
    }

    # calculate bar heights
    bar_heights <- rowSums(overlapping_freqs)
    total_counts <- rowSums(overlapping_animals)

    # generate blank plot
    barplot(t(overlapping_freqs), col=NA, border=NA, xpd=FALSE, axes=FALSE, axisnames=FALSE, ylim=c(0,110),
            width=0.75, xlab="Co-occurring group sizes", ylab="Frequency", cex.lab=label.cex)

    # add background color
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray98", border=NULL)

     # set color palette
    if(!is.null(levels.order)) color_pal <- bar.colors[levels.order][group_levels]
    else color_pal <- bar.colors[group_levels]

    # add stacked barplot
    p<-barplot(t(overlapping_freqs), col=color_pal, lwd=0.6, xpd=FALSE, axes=FALSE, axisnames=FALSE, ylim=c(0,110),
               width=0.75, xlab="", ylab="", add=TRUE)

    # add plot titles
    if(!is.null(id.groups)){title(main=names(data_list)[i], line=1, font=2, cex.main=title.cex)}

    # draw axes
    axis(1, at=p, labels=size_labels, cex.axis=axis.cex)
    axis(2, at=seq(0,100,by=20), paste0(seq(0,100,by=20),"%"), cex.axis=axis.cex, las=2)

    # add total count annotation above each column
    text(x=p, y=bar_heights+5, labels=paste0("n=", total_counts), cex=annot.cex)

    # add legend
    if(!is.null(split.by)) {
      .legend(legend.pos, legend=group_labels, fill=bar.colors, horiz=legend.style,
              box.cex=c(1.2, 1),  y.intersp=1.2, bty="n", cex=legend.cex)
    }

    # draw box
    box()
   }

}


################################################################################
# Auxiliary function ###########################################################
################################################################################

#' Count maximum nº of co-occurring animals in each timebin
#' @note This function is intended for internal use within the 'moby' package.
#' @keywords internal
#' @noRd

.countOverlappingIDs <- function(data, type, groups) {
  if(type=="between"){
    if(any(table(data)>1) & length(unique(groups[!is.na(data)]))>1){return(as.numeric(table(data)[table(data)>1]))
    }else {return(0)}
  }else if(type=="within"){
    if(any(table(data)>1)) {return(as.numeric(table(data)[table(data)>1]))
    }else {return(0)}
  }
}


##################################################################################################
##################################################################################################
##################################################################################################
