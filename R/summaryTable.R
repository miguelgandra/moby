#######################################################################################################
## Create summary table ###############################################################################
#######################################################################################################

#' Generate summary table for tagged animals
#'
#' @description This function generates a summary table with information about tagged animals,
#' including tagging dates, last detection dates, and various monitoring and residency metrics. It allows for the
#' inclusion of additional metadata and the calculation of overall means and error metrics.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections. Each row should represent an individual detection event,
#' unless a 'detections' column is included to indicate the number of detections for each row.
#' @param id.metadata A data frame containing metadata about the tagged animals, such as their length,
#' sex, or transmitter type. If there are multiple rows per animal, variables will be collapsed before merging with other statistics.
#' @param residency.by Optional. Variable used to calculate partial residencies (e.g. array or habitat).
#' Defaults to NULL.
#' @param id.groups Optional. A list where each element is a group of IDs, used for visually aggregating
#' animals belonging to the same class (e.g., different species or life stages). If supplied, averages
#' will be calculated independently for each group.
#' @param error.stat The statistic to use for variability/error calculation, either 'sd' (standard deviation)
#' or 'se' (standard error). Defaults to 'sd'.
#' @importFrom plotrix std.error
#' @export


summaryTable <- function(data, tagging.dates=getDefaults("tagdates"), id.metadata=NULL, id.col=getDefaults("id"),
                         datetime.col=getDefaults("datetime"), residency.by=NULL, id.groups=NULL, error.stat="sd") {


  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  # check if data contains id.col
  if(!id.col %in% colnames(data)) stop("ID column not found. Please specify the correct column using 'id.col'")
  # check if data contains id.col
  if(!id.col %in% colnames(id.metadata)) stop("ID column not found in id.metadata. Please specify the correct column using 'id.col'")
  # check if data contains datetime.col
  if(!datetime.col %in% colnames(data)) stop("Datetime column not found. Please specify the correct column using 'datetime.col'")
  # check if datetimes are in the right format
  if(!grepl("POSIXct", paste(class(data[, datetime.col]), collapse = " "))) stop("Datetimes must be provided in POSIXct format")
  # check if tagging.dates were supplied
  if(is.null(tagging.dates)) stop("Tagging dates not specified. Please provide the required dates via the 'tagging.dates' argument")
  # check if tagging.dates are in the right format
  if(!grepl("POSIXct", paste(class(tagging.dates), collapse = " "))) stop("'tagging.dates' must be provided in POSIXct format")
  # check if data contains datetime.col
  if(!residency.by %in% colnames(data)) stop("Variable for partial residencies not found in the supplied data. Please specify the correct column using 'residency.by'")
   # check error function
  if(!error.stat %in% c("sd", "se")) stop("Wrong error.stat argument, please choose between 'sd' and 'se'.")

  # convert IDs to factor
  if(class(data[,id.col])!="factor") {
    data[,id.col] <- as.factor(data[,id.col])
    cat("Warning: 'id.col' converted to factor\n")
  }

  # check number of tagging dates
  if(length(tagging.dates)>1 & length(tagging.dates)!=nlevels(data[,id.col])){
    stop("Incorrect number of tagging.dates. Must be either a single value or
           a vector containing a tagging date for each individual")
  }else if(length(tagging.dates)==1){
    tagging.dates <- rep(tagging.dates, nlevels(data[,id.col]))
  }

  getErrorFun <- function(x) {
    if(error.stat=="sd"){return(sd(x, na.rm=T))}
    if(error.stat=="se"){return(plotrix::std.error(x))}
  }

  # reorder ID levels if ID groups are defined
  if(!is.null(id.groups)){
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {cat("Warning: Some of the ID(s) in id.groups don't match the IDs in the data")}
    data <- data[data[,id.col] %in% unlist(id.groups),]
    tagging.dates <- tagging.dates[match(unlist(id.groups), levels(data[,id.col]))]
    data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
  }else{
    id.groups <- list(levels(data[,id.col]))
  }


  ##############################################################################
  ## Generate table ############################################################
  ##############################################################################

  # format tagging dates
  tag_dates <- strftime(tagging.dates, format="%d/%m/%Y", tz="UTC")

  # retrieve last detections dates
  last_detections <- tapply(X=data[,datetime.col], INDEX=data[,id.col], FUN=max)
  last_detections <- as.POSIXct(last_detections, origin='1970-01-01', tz="UTC")
  last_dates <- strftime(last_detections, format="%d/%m/%Y", tz="UTC")

  # number of detections
  if("detections" %in% colnames(data)){
    detections <- as.integer(graphics::aggregate(as.formula(paste0("detections~", id.col)), data=data, FUN=sum, drop=F)$detections)
  }else{
    warnin
    print("Warning: No 'detections' column found, assuming one detection per row\n")
    detections <- as.integer(table(data[,id.col]))
  }
  detections[detections==0] <- NA

  #number of receivers
  receivers <- graphics::aggregate(data$station, by=list(data[,id.col]), function(x) length(unique(x)), drop=F)$x

  # days between 1st and last detection (Di)
  getTimeSeqs <- function(start, end) {if(is.na(end)) {return(NA)}else{seq.POSIXt(start, end, by="day")}}
  timeseqs <- mapply(getTimeSeqs, start=tagging.dates, end=last_detections, SIMPLIFY=F)
  Di <- as.integer(lapply(timeseqs, function(x) length(x[!is.na(x)])))
  Di[Di==0] <- NA

  # days with detections (Dd)
  data$date <- strftime(data[,datetime.col], format="%d-%m-%Y", tz="UTC")
  Dd <- graphics::aggregate(data$date, by=list(data[,id.col]), function(x) length(unique(x)), drop=F)$x

  # complete IR
  Ir <- round(Dd/Di, 2)

  # calculate partial residencies if a residency.by variable is supplied
  if(!is.null(residency.by)){
    data_groupped <- split(data, f=data[,residency.by])
    Ir_partial <- list()
    for(i in 1:length(data_groupped)){
      data_subset <- data_groupped[[i]]
      Dd_partial <-  graphics::aggregate(data_subset$date, by=list(data_subset[,id.col]), function(x) length(unique(x)), drop=F)$x
      Ir_partial[[i]] <- round(Dd_partial/Di, 2)
    }
  }

  # aggregate stats
  stats <- data.frame("ID"=levels(data[,id.col]), "Tagging date"=tag_dates, "Last detection"=last_dates,
                      "N Detect"=detections, "N Receiv"=receivers, "Detection span (d)"=Di,
                      "N days detected"=Dd, "IR"=Ir, row.names=NULL, check.names=F)

  # add partial residencies (if available)
  if(!is.null(residency.by)){
    for(i in 1:length(data_groupped)){
      stats[,names(data_groupped)[i]] <- Ir_partial[i]
    }}

  # format additional tag info (if available) and merge
  if(!is.null(id.metadata)) {
    column_types <- sapply(1:ncol(id.metadata),function(c) class(id.metadata[,c]))
    numeric_cols <- which(column_types %in% c("numeric", "integer"))
    animal_info <- graphics::aggregate(id.metadata[,-1], by=list(id.metadata[,id.col]), function(x) paste(unique(x), collapse="/"), drop=F)
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
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(group[,numeric_cols], 2, getErrorFun)))
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
