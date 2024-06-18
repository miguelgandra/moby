#######################################################################################################
## Plot co-occurrences by station and diel phase ######################################################
#######################################################################################################

#' Plot frequencies of co-occurring group sizes

#' @description Function to calculate the frequency distribution of the
#' nº of animals co-occurring at each time-bin, either for the
#' entire study duration or split by a given variable (e.g. diel phase).
#'
#' @param table A data frame object containing binned detections in the wide format
#' (time bin x individual matrix, with values corresponding to the receiver with the
#' highest number of detections), as returned by \code{\link{createWideTable}}.
#' @param n.individuals Number of tagged individuals.
#' @param split.by Optional. If defined, frequencies are calculated separately for each level
#' of this variable. If NULL, frequencies are calculated for the entire study duration.
#' @param id.groups Optional. A list containing ID groups, used to generate different
#' plots to each class (e.g. different species).
#' @param group.comparisons Controls the type of comparisons to be run, when id.groups are defined.
#' One of "within", "between" or "all". Defaults to "all".
#' @param levels.order Numeric vector with the preferred order of the levels defined by the 'split.by' variable
#' (in the stacked barplot).
#' @param bar.colors A vector of colors for the bars. Nº of colors should match
#' the nº of levels in the 'split.by' variable.
#' @param legend.pos Legend position.
#' @param legend.style Legend style (horizontal vs vertical).
#' @param cols Number of columns in the final panel (passed to the mfrow argument).
#' @seealso \code{\link{calculateOverlap}}
#' @export


######################################################################################
# Main function - plot frequencies of co-occurring group sizes #######################
######################################################################################

plotOverlapFrequences <- function(table, n.individuals, split.by=NULL, id.groups=NULL,
                                  group.comparisons="all", legend.pos="topright", legend.style="horizontal",
                                  bar.colors=NULL, levels.order=NULL, cols=1) {


  ######################################################################################
  # Initial checks #####################################################################
  ######################################################################################

  # set right data format
  out_cols <- setdiff(1:ncol(table), 1:n.individuals+1)
  id_cols <- (1:ncol(table))[-out_cols]
  ids <- colnames(table[,id_cols])
  table_detections <- table[,id_cols]

  if(all(!sapply(id_cols, function(x) class(table[,x])) %in% c("character", "factor"))){
    stop("Please supply a table with values corresponding to the receiver with the highest number of detections")
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% ids)) {cat("Warning: Some of the ID(s) in id.groups don't match the IDs in the data")}
    if(any(!ids %in% unlist(id.groups))) {cat("Warning: Some of the ID(s) in the supplied table were not found within id.groups and will be discarded")}
    table_detections <- table_detections[,colnames(table_detections) %in% unlist(id.groups)]
  }

  # check split.by variable
  if(!is.null(split.by)){
    if(!split.by %in% colnames(table)) {
      stop("split.by variable not found in the supplied table")
    }
    if(class(table[,split.by])!="factor"){
      cat("Warning: converting split.by variable to factor\n")
      table[,split.by] <- as.factor(table[,split.by])
    }
    # print to console
    cat(paste0("Calculating frequency of co-occurring group sizes by ", split.by, "\n"))
  }else{
    if(!is.null(levels.order)) {
      cat("Warning: levels.order supplied but will have no effect - no split.by variable defined")
      levels.order <- 1
    }
    # print to console
    cat("Calculating overall frequencies of co-occurring group sizes\n")
  }

  # check color palette
  if(!is.null(bar.colors)){
    if(!is.null(split.by)){
      if(length(bar.colors)!=nlevels(table[,split.by])){
        cat("Warning: the number of supplied colors does not match the number of levels in split.by variable\n")}
    }else{
      if(length(bar.colors)>1){
        cat("Warning: multiple colors supplied but only the first one will be used\n")}
    }
  }



  # check group.comparisons argument
  if(any(!group.comparisons %in% c("within", "between", "all"))) {
    stop("Wrong 'group.comparisons' argument. Please choose one of: 'within', 'between' or 'all'")}


  ######################################################################################
  # Split data (if required) ###########################################################
  ######################################################################################

  if(!is.null(id.groups)){
    ids_table <- reshape2::melt(id.groups)
    colnames(ids_table) <- c("ID", "group")
    if(group.comparisons=="within"){
      data_list <- lapply(id.groups, function(x) table_detections[,colnames(table_detections) %in% x])
      names(data_list) <- names(id.groups)
    }else if(group.comparisons=="between"){
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse="<->"))
      comb_ids <- apply(group_comb, 1, function(x) ids_table$ID[ids_table$group %in% x], simplify=F)
      data_list <- lapply(comb_ids, function(x) table_detections[,colnames(table_detections) %in% x])
      names(data_list) <- group_comb_names
    }else{
      data_list <- lapply(id.groups, function(x) table_detections[,colnames(table_detections) %in% x])
      names(data_list) <- names(id.groups)
      group_comb <- t(combn(names(id.groups), 2))
      group_comb_names <- apply(group_comb, 1, function(x) paste0(x, collapse="<->"))
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

     table_detections <- data_list[[i]]

     ####################################################################
     # calculate overall frequencies of co-occurring group sizes ########

     if(is.null(split.by)) {

       # replace IDs by group ID
       colnames(table_detections) <- plyr::mapvalues(colnames(table_detections), ids_table$ID, ids_table$group, warn_missing=F)

       # calculate co-occurring group sizes for each timebin
       if(group.comparisons=="within" | (group.comparisons=="all" & !grepl("<->", names(data_list)[i], fixed=T))){
         overlapping_animals <- unlist(apply(table_detections, 1, countOverlappingIDs, type="within", groups=colnames(table_detections)))
       }else if(group.comparisons=="between" | (group.comparisons=="all" & grepl("<->", names(data_list)[i], fixed=T))){
         overlapping_animals <- unlist(apply(table_detections, 1, countOverlappingIDs, type="between", groups=colnames(table_detections)))
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
       table_detections <- cbind(table_detections, table[,split.by])
       colnames(table_detections)[ncol(table_detections)] <- split.by

       # subset data by 'split.by' variable
       data_subset <- split(table_detections, f=table_detections[,split.by])
       data_subset <- lapply(data_subset, function(x) x[,-ncol(x)])
       data_subset <- data_subset[order(match(names(data_subset), levels(table[,split.by])))]

        # replace IDs by group ID
       data_subset <- lapply(data_subset, function(x){colnames(x) <- plyr::mapvalues(colnames(x), ids_table$ID, ids_table$group, warn_missing=F); return(x)})

       # calculate co-occurring group sizes for each timebin
       if(group.comparisons=="within" | (group.comparisons=="all" & !grepl("<->", names(data_list)[i], fixed=T))){
         overlapping_animals <- reshape2::melt(lapply(data_subset, function(x) unlist(apply(x, 1, countOverlappingIDs, type="within", groups=colnames(x)))))
       }else if(group.comparisons=="between" | (group.comparisons=="all" & grepl("<->", names(data_list)[i], fixed=T))){
         overlapping_animals <- reshape2::melt(lapply(data_subset, function(x) unlist(apply(x, 1, countOverlappingIDs, type="between", groups=colnames(x)))))
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
  par(mfrow=c(rows, cols))

  # set legend orientation
  legend.style <- switch(legend.style, "horizontal"=T, "vertical"=F)

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


    barplot(t(overlapping_freqs), col=NA, border=NA, xpd=F, axes=F, axisnames=F, ylim=c(0,110),
            width=0.75, xlab="Co-occurring group sizes", ylab="Frequency", cex.lab=0.95)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray98", border=NULL)
    if(!is.null(levels.order)){
      color_pal <- bar.colors[levels.order][group_levels]
    }else{
      color_pal <- bar.colors[group_levels]
    }
    p<-barplot(t(overlapping_freqs), col=color_pal, lwd=0.6, xpd=F, axes=F, axisnames=F, ylim=c(0,110),
               width=0.75, xlab="", ylab="", add=T)
    # add plot titles
    if(!is.null(id.groups)){title(main=names(data_list)[i], line=1, font=2, cex.main=0.95)}
    axis(1, at=p, labels=size_labels, cex.axis=0.85)
    axis(2, at=seq(0,100,by=20), paste0(seq(0,100,by=20),"%"), cex.axis=0.85, las=2)
    text(x=p, y=bar_heights+5, labels=paste0("n=", total_counts), cex=0.75)
    if(!is.null(split.by)) {
      .legend(legend.pos, legend=group_labels, fill=bar.colors,
                         horiz=legend.style, box.cex=c(1.2, 1),  y.intersp=1.2, bty="n", cex=0.8)
    }
    box()


   }

}


################################################################################
# Auxiliary function - count maximum nº of co-occurring animals in each timebin
################################################################################

countOverlappingIDs <- function(data, type, groups) {
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
