#######################################################################################################
## Create summary table ###############################################################################
#######################################################################################################

#' Create summary table
#'
#' @description Creates a table containing information about the
#' tagged animals (tagging dates, last detections, etc.) as well
#' as residency metrics.
#'
#' @param data A data frame containing binned animal detections.
#' @param tagging.dates A POSIXct vector containing the tag/release date of each animal.
#' @param tags.info A data frame containing additional informations about the tagged animals,
#' such as length, sex or transmitter type. If more than one row exists per animal,
#' variables are collapsed before being merged with the remaining stats.
#' @param residency.by Variable used to calculate partial residencies (e.g. array or habitat).
#' Defaults to NULL.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species or life stages).
#' @param error.stat Statistic to present with the mean. Either standard deviation (sd) or standard error (se).
#' @importFrom plotrix std.error
#' @export


summaryTable <- function(data, tagging.dates, tags.info=NULL, residency.by=NULL,
                         id.groups=NULL, error.stat="sd") {


  ####################################################################
  ## Set error function ##############################################

  if(!error.stat %in% c("sd", "se")){
    stop("Wrong error.stat argument, please choose either 'sd' or 'se'.")
  }

  getError <- function(x) {
    if(error.stat=="sd"){return(sd(x, na.rm=T))}
    if(error.stat=="se"){return(plotrix::std.error(x))}
  }

  ####################################################################
  ## Prepare data ####################################################

  # set id groups
  if(is.null(id.groups)) {
    id.groups<-list(levels(data$ID))
  } else {
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% colnames(table))) {"Some of the ID(s) in id.groups don't match the IDs in the data"}
  }

  # format tagging dates
  tag_dates <- strftime(tagging.dates, format="%d/%m/%Y", tz="UTC")

  # retrieve last detections dates
  last_detections <- tapply(X=data$timebin, INDEX=data$ID, FUN=max)
  last_detections <- as.POSIXct(last_detections, origin='1970-01-01', tz="UTC")
  last_dates <- strftime(last_detections, format="%d/%m/%Y", tz="UTC")

  # number of detections
  if("detections" %in% colnames(data)){
    detections <- as.integer(aggregate(detections~ID, data=data, FUN=sum, drop=F)$detections)
  }else{
    print("Warning: No 'detections' column found, assuming one detection per row\n")
    detections <- as.integer(table(data$ID))
  }
  detections[detections==0] <- NA

  #number of receivers
  receivers <- aggregate(data$station, by=list(data$ID), function(x) length(unique(x)), drop=F)$x

  # days between 1st and last detection (Di)
  getTimeSeqs <- function(start, end) {if(is.na(end)) {return(NA)}else{seq.POSIXt(start, end, by="day")}}
  timeseqs <- mapply(getTimeSeqs, start=tagging.dates, end=last_detections, SIMPLIFY=F)
  Di <- as.integer(lapply(timeseqs, function(x) length(x[!is.na(x)])))
  Di[Di==0] <- NA

  # days with detections (Dd)
  data$date <- strftime(data$timebin, format="%d-%m-%Y", tz="UTC")
  Dd <- aggregate(data$date, by=list(data$ID), function(x) length(unique(x)), drop=F)$x

  # complete IR
  Ir <- round(Dd/Di, 2)

  # calculate partial residencies if a residency.by variable is supplied
  if(!is.null(residency.by)){
    data_groupped <- split(data, f=data[,residency.by])
    Ir_partial <- list()
    for(i in 1:length(data_groupped)){
      data_subset <- data_groupped[[i]]
      Dd_partial <-  aggregate(data_subset$date, by=list(data_subset$ID), function(x) length(unique(x)), drop=F)$x
      Ir_partial[[i]] <- round(Dd_partial/Di, 2)
    }
  }

  # aggregate stats
  stats <- data.frame("ID"=levels(data$ID), "Tagging date"=tag_dates, "Last detection"=last_dates,
                      "N Detect"=detections, "N Receiv"=receivers, "Detection span (d)"=Di,
                      "N days detected"=Dd, "IR"=Ir, row.names=NULL, check.names=F)

  # add partial residencies (if available)
  if(!is.null(residency.by)){
    for(i in 1:length(data_groupped)){
      stats[,names(data_groupped)[i]] <- Ir_partial[i]
    }}

  # format additional tag info (if available) and merge
  if(!is.null(tags.info)) {
    column_types <- sapply(1:ncol(tags.info),function(c) class(tags.info[,c]))
    numeric_cols <- which(column_types %in% c("numeric", "integer"))
    animal_info <- aggregate(tags.info[,-1], by=list(tags.info$ID), function(x) paste(unique(x), collapse="/"), drop=F)
    animal_info[animal_info=="NA"] <- NA
    animal_info[,numeric_cols] <- as.numeric(animal_info[,numeric_cols] )
    colnames(animal_info)[1] <- "ID"
    stats <- plyr::join(animal_info, stats, by="ID", type="left")
  }

  # calculate means ± se and format missing values
  stats$ID <- as.character(stats$ID)
  column_types <- sapply(1:ncol(stats),function(c) class(stats[,c]))
  numeric_cols <- which(column_types %in% c("numeric", "integer"))
  decimal_digits <- apply(stats[,numeric_cols], 2, moby:::decimalPlaces)
  decimal_digits <- apply(decimal_digits, 2, max, na.rm=T)

  n_groups <- length(id.groups)
  group_stats <- lapply(id.groups, function(x) stats[stats$ID %in% x,])
  for(i in 1:n_groups){
    group <- group_stats[[i]]
    group[nrow(group)+1,] <- NA
    group$ID[nrow(group)] <- "mean"
    group[nrow(group), numeric_cols] <- sprintf(paste0("%.", decimal_digits, "f"), colMeans(group[,numeric_cols], na.rm=T))
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(group[,numeric_cols], 2, getError)))
    group[nrow(group), numeric_cols] <- paste(group[nrow(group), numeric_cols], "±", errors)
    for(c in 1:length(numeric_cols)){
      col_number <- numeric_cols[c]
      group[-nrow(group), col_number] <- sprintf(paste0("%.", decimal_digits[c], "f"), as.numeric(group[-nrow(group), col_number]))
    }
    #Ir_cols <- which(grepl("IR", colnames(group), fixed=T))
    #if(!is.null(residency.by)){Ir_cols <- c(Ir_cols, which(colnames(group) %in% names(data_groupped)))}
    #for(c in Ir_cols){group[-nrow(group), c] <- sprintf("%.2f", as.numeric(group[-nrow(group), c]))}
    group[group=="NA"] <- "-"
    group[group=="NaN ± NA"] <- "-"
    group[is.na(group)] <- "-"
    group_stats[[i]]<-group
  }

  # add group names
  if(n_groups>1){
    group_labels <- stats[0,]
    group_labels[1:n_groups,] <- ""
    group_labels$ID <- names(id.groups)
    group_labels <- split(group_labels, f=group_labels$ID)
    group_labels <- group_labels[match(names(group_labels), names(id.groups))]
    group_stats <- mapply(function(label, stats) {rbind(label, stats)}, label=group_labels, stats=group_stats, SIMPLIFY=F)
  }

  # return table
  stats <- do.call("rbind", group_stats)
  rownames(stats) <- NULL
  return(stats)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
