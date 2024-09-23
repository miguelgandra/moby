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
#' @param tag.durations Optional. A numeric vector containing the estimated battery
#' duration of the deployed tags (in days). If a single value is provided, it will be applied to all IDs.
#' @param residency.by Optional. Variable used to calculate partial residencies (e.g. array or habitat).
#' Defaults to NULL.
#' @param residency.index A character string specifying the type of residency index to calculate.
#' Options include:
#'  - "IR1": Residency Index 1, calculated as the number of days the animal was detected (Dd) divided by the detection interval (Di), i.e., the number of days between tagging and last detection (days at liberty). This represents a maximum residency value, considering only the period for which the animal was known to be alive and the tag operational.
#'  - "IR2": Residency Index 2, calculated as the number of days the animal was detected (Dd) divided by the study interval (Dt), i.e., the total number of days between release and last data download or tag expiration date. This approach provides a minimum residency value, assuming the animal was alive and detectable throughout the study period.
#'  - "IWR": Weighted Residency Index, which combines two ratios: Dd/Dt (days detected over study interval) weighted by Di/Dt (detection interval over study interval). It provides a more nuanced measure of residency that balances the time an animal is detected with the total study period.
#'
#' The choice of index can affect the interpretation of residency patterns, so it's important to select the one that best fits the study objectives.
#' Defaults to "IR1". If `tag.durations` are provided, the default changes to "IR2".
#' Further information on residency estimation methods can be found in Kraft et al. (2023) - see below in the references section.
#' @param id.groups Optional. A list where each element is a group of IDs, used for visually aggregating
#' animals belonging to the same class (e.g., different species or life stages). If supplied, averages
#' will be calculated independently for each group.
#' @param error.stat The statistic to use for variability/error calculation, either 'sd' (standard deviation)
#' or 'se' (standard error). Defaults to 'sd'.
#' @references
#' Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., & Abecasis, D. (2023).
#' Residency and space use estimation methods based on passive acoustic telemetry data.
#' Movement Ecology, 11(1), 12.
#' https://doi.org/10.1186/s40462-023-00349-y
#' @importFrom plotrix std.error
#' @export


summaryTable2 <- function(data, tagging.dates=getDefaults("tagdates"), tag.durations=NULL,
                         id.metadata=NULL, id.col=getDefaults("id"), datetime.col=getDefaults("datetime"),
                         station.col=getDefaults("station"), residency.index="IR1",
                         residency.by=NULL, id.groups=NULL, error.stat="sd") {


  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates

  # validate additional parameters
  errors <- c()
  # check if data contains residency.by
  if(!is.null(residency.by) && !residency.by %in% colnames(data)) errors <- c(errors, "Variable used to calculate partial residencies not found in the supplied data. Please specify the correct column using 'residency.by'.")
  # check error function
  error.stat <- tolower(error.stat)
  if(!error.stat %in% c("sd", "se")) errors <- c(errors, "Wrong error.stat argument, please choose between 'sd' and 'se'.")
  # check if id.metadata contains id.col
  if(!is.null(id.metadata) && !id.col %in% colnames(id.metadata))   errors <- c(errors, "The specified ID column ('id.col') does not exist in 'id.metadata'. Please ensure that the column name in 'id.metadata' matches the 'id.col' specified.")
  # check residency index
  if(!residency.index %in% c("IR1", "IR2", "IWR"))   errors <- c(errors, "Invalid 'residency.index' argument. Please select one of the following options: 'IR1' (detection interval), 'IR2' (study interval), or 'IWR' (weighted residency index).")
  # print all errors
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # define error function
  getErrorFun <- function(x) {
    if(error.stat=="sd"){return(sd(x, na.rm=T))}
    if(error.stat=="se"){return(plotrix::std.error(x))}
  }




  ##############################################################################
  ## Generate table ############################################################
  ##############################################################################

  # define single id.group if needed
  if(is.null(id.groups)){
    id.groups <- list(levels(data[,id.col]))
  }

  # format tagging dates
  tag_dates <- strftime(tagging.dates, format="%d/%m/%Y", tz="UTC")

  # retrieve last detections dates
  last_detections <- tapply(X=data[,datetime.col], INDEX=data[,id.col], FUN=max)
  last_detections <- as.POSIXct(last_detections, origin='1970-01-01', tz="UTC")
  last_dates <- strftime(last_detections, format="%d/%m/%Y", tz="UTC")

  # number of detections
  if("detections" %in% colnames(data)){
    detections <- as.integer(stats::aggregate(as.formula(paste0("detections~", id.col)), data=data, FUN=sum, drop=F)$detections)
  }else{
    warning("No 'detections' column found, assuming one detection per row.", call.=FALSE)
    detections <- as.integer(table(data[,id.col]))
  }
  detections[detections==0] <- NA

  #number of receivers
  receivers <- stats::aggregate(data[,station.col], by=list(data[,id.col]), function(x) length(unique(x)), drop=F)$x

  # days between 1st and last detection (Di)
  getTimeSeqs <- function(start, end) {if(is.na(end)) {return(NA)}else{seq.POSIXt(start, end, by="day")}}
  timeseqs <- mapply(getTimeSeqs, start=tagging.dates, end=last_detections, SIMPLIFY=F)
  Di <- as.integer(lapply(timeseqs, function(x) length(x[!is.na(x)])))
  Di[Di==0] <- NA

  # days with detections (Dd)
  data$date <- strftime(data[,datetime.col], format="%d-%m-%Y", tz="UTC")
  Dd <- stats::aggregate(data$date, by=list(data[,id.col]), function(x) length(unique(x)), drop=F)$x

  # complete IR
  Ir <- round(Dd/Di, 2)

  # calculate partial residencies if a residency.by variable is supplied
  if(!is.null(residency.by)){
    data_groupped <- split(data, f=data[,residency.by])
    Ir_partial <- list()
    for(i in 1:length(data_groupped)){
      data_subset <- data_groupped[[i]]
      Dd_partial <-  stats::aggregate(data_subset$date, by=list(data_subset[,id.col]), function(x) length(unique(x)), drop=F)$x
      Ir_partial[[i]] <- round(Dd_partial/Di, 2)
    }
  }

  # aggregate stats
  stats <- data.frame("ID"=levels(data[,id.col]), "Tagging date"=tag_dates, "Last detection"=last_dates,
                      "N Detect"=detections, "N Receiv"=receivers, "Detection span (d)"=Di,
                      "N days detected"=Dd, "IR"=Ir, row.names=NULL, check.names=F)

  # add residency index
  colnames(stats)[colnames(stats)=="IR"] <- residency.index

  # add partial residencies (if available)
  if(!is.null(residency.by)){
    for(i in 1:length(data_groupped)){
      stats[,names(data_groupped)[i]] <- Ir_partial[i]
    }}

  # format additional tag info (if available) and merge
  if(!is.null(id.metadata)) {
    column_types <- sapply(1:ncol(id.metadata),function(c) class(id.metadata[,c]))
    numeric_cols <- which(column_types %in% c("numeric", "integer"))
    animal_info <- stats::aggregate(id.metadata[,-1], by=list(id.metadata[,id.col]), function(x) paste(unique(x), collapse="/"), drop=F)
    animal_info[animal_info=="NA"] <- NA
    animal_info[,numeric_cols] <- as.numeric(animal_info[,numeric_cols] )
    colnames(animal_info)[1] <- "ID"
    stats <- plyr::join(animal_info, stats, by="ID", type="left")
  }

  # calculate means ± se and format missing values
  stats$ID <- as.character(stats$ID)
  column_types <- sapply(1:ncol(stats),function(c) class(stats[,c]))
  numeric_cols <- which(column_types %in% c("numeric", "integer"))
  decimal_digits <- apply(stats[,numeric_cols], 2, .decimalPlaces)
  decimal_digits <- apply(decimal_digits, 2, max, na.rm=T)

  n_groups <- length(id.groups)
  group_stats <- lapply(id.groups, function(x) stats[stats$ID %in% x,])
  for(i in 1:n_groups){
    group <- group_stats[[i]]
    group[nrow(group)+1,] <- NA
    group$ID[nrow(group)] <- "mean"
    group[nrow(group), numeric_cols] <- sprintf(paste0("%.", decimal_digits, "f"), colMeans(group[,numeric_cols], na.rm=T))
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(group[,numeric_cols], 2, getErrorFun)))
    group[nrow(group), numeric_cols] <- paste(group[nrow(group), numeric_cols], "\u00b1", errors)
    for(c in 1:length(numeric_cols)){
      col_number <- numeric_cols[c]
      group[-nrow(group), col_number] <- sprintf(paste0("%.", decimal_digits[c], "f"), as.numeric(group[-nrow(group), col_number]))
    }
    #Ir_cols <- which(grepl("IR", colnames(group), fixed=T))
    #if(!is.null(residency.by)){Ir_cols <- c(Ir_cols, which(colnames(group) %in% names(data_groupped)))}
    #for(c in Ir_cols){group[-nrow(group), c] <- sprintf("%.2f", as.numeric(group[-nrow(group), c]))}
    group[group=="NA"] <- "-"
    group[group=="NaN \u00b1 NA"] <- "-"
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
