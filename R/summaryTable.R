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
#' duration of the deployed tags (in days). The length of this vector should match the number of
#' unique animal IDs, and the values must be in the same order as the ID levels.
#' Alternatively, if a single value is provided, it will be applied to all IDs.
#' @param residency.index A character string specifying the type of residency index to calculate.
#' Options include:
#'  - "IR1": Residency Index 1, calculated as the number of days the animal was detected (Dd) divided by the detection interval (Di), i.e., the number of days between release/first detection and last detection (days at liberty).
#'  This represents a maximum residency value, considering only the period for which the animal was known to be alive and the tag operational.
#'  - "IR2": Residency Index 2, calculated as the number of days the animal was detected (Dd) divided by the study interval (Dt), i.e., the total number of days between release/first detection and last data download or tag expiration date.
#'  This approach provides a minimum residency value, assuming the animal was alive and detectable throughout the study period.
#'  - "IWR": Weighted Residency Index, which corresponds to the IR2 index weighted by the ratio between the detection interval (Di, the number of days between the first and last detection) and the study interval (Dt, the total monitoring period).
#'  This accounts for the number of days detected and the spread of detections within the monitoring period, providing a measure of residency that balances the frequency of detections with their temporal distribution.
#'
#' The choice of index can affect the interpretation of residency patterns, so it's important to select the one(s) that best fits the study objectives.
#' Defaults to "IR1". If `tag.durations` are provided, the default changes to "IR2".
#' Further information on residency estimation methods can be found in Kraft et al. (2023) and Appert et al. (2023) - see below in the references section.
#' @param start.point A character string specifying the starting point for calculating days at liberty.
#' Options include:
#' - `"release"`: The release day is used as the starting point for calculating days at liberty.
#' - `"first.detection"`: The first detection after tagging is used as the starting point. This option can help account for potential post-release behavior effects.
#'
#' Defaults to `"release"`.
#' @param residency.by Optional. Variable used to calculate partial residencies (e.g. array or habitat).
#' Defaults to NULL.
#' @param id.groups Optional. A list where each element is a group of IDs, used for visually aggregating
#' animals belonging to the same class (e.g., different species or life stages). If supplied, averages
#' will be calculated independently for each group.
#' @param error.stat The statistic to use for variability/error calculation, either 'sd' (standard deviation)
#' or 'se' (standard error). Defaults to 'sd'.
#'
#' @return A data frame summarizing information on tagged animals, with the following columns:
#' - `ID`: Unique identifier for each tagged animal.
#' - Any additional metadata columns from `id.metadata` if provided.
#' - `Tagging date`: The date when the animal was tagged.
#' - `Tag duration (d)`: Duration of the tag in days (if provided).
#' - `Last detection`: The date of the last detection.
#' - `N Detect`: Total number of detections for the animal.
#' - `N Receiv`: Number of unique receivers that detected the animal.
#' - `Monitoring duration (d)`: Total duration of monitoring in days.
#' - `Detection span (d)`: Number of days between release/first detection and last detection (days at liberty)
#' - `N days detected`: Total number of days the animal was detected.
#' - Additional columns for each residency index specified in the `residency.index` parameter.
#' - If `residency.by` is specified, additional columns for partial residency metrics will be included.
#'
#' @references
#' Kraft, S., Gandra, M., Lennox, R. J., Mourier, J., Winkler, A. C., & Abecasis, D. (2023).
#' Residency and space use estimation methods based on passive acoustic telemetry data.
#' Movement Ecology, 11(1), 12.
#' https://doi.org/10.1186/s40462-023-00349-y
#'
#' Appert, C., Udyawer, V., Simpfendorfer, C. A., Heupel, M. R., Scott, M., Currey-Randall, L. M., ... & Chin, A. (2023).
#' Use, misuse, and ambiguity of indices of residence in acoustic telemetry studies.
#' Marine Ecology Progress Series, 714, 27-44.
#' https://doi.org/10.3354/meps14300
#'
#' @importFrom plotrix std.error
#' @export


summaryTable <- function(data,
                         tagging.dates = getDefaults("tagdates"),
                         tag.durations = NULL,
                         id.metadata = NULL,
                         id.col = getDefaults("id"),
                         datetime.col = getDefaults("datetime"),
                         station.col = getDefaults("station"),
                         residency.index = c("IR1", "IR2", "IWR"),
                         start.point = "release",
                         residency.by = NULL,
                         id.groups = NULL,
                         error.stat = "sd") {


  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates
  tag.durations <- reviewed_params$tag.durations


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
  if(!all(residency.index %in% c("IR1", "IR2", "IWR"))) errors <- c(errors, "Invalid 'residency.index' argument. Please select one of the following options: 'IR1' (detection interval), 'IR2' (study interval), or 'IWR' (weighted residency index).")
  if (!is.character(start.point) || length(start.point) !=1) errors <- c(errors, "Invalid 'start.point' parameter: it must be a single character string.")
  if (!start.point %in% c("first.detection", "release")) errors <- c(errors, "Invalid 'start.point' parameter: must be one of 'release' or 'first.detection'.")
  # print all errors
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }

  # print warning in case residency.index included 'IR2' or 'IWR' and no tag.durations were provided
  if (any(residency.index %in% c("IR2", "IWR")) && is.null(tag.durations)) {
    warning("The residency index includes 'IR2' or 'IWR', but no tag durations have been supplied. The function will default to assuming the last detection date in the dataset as the end of the monitoring period. This may lead to spurious index values, particularly if the estimated expiration date of a tag is earlier than the last detection in the dataset.", call.=FALSE)
  }

  # prompt the user for multiple residency.by levels and residency.index combinations
  if(!is.null(residency.by)) {
    nlevels <- length(unique(data[[residency.by]]))
    if((nlevels*length(residency.index))>=9 && length(residency.index)>1) {
      # print prompt
      proceed <- menu(c("Yes", "No"), title="Residency.by has multiple levels and more than one residency index will be calculated. This can result in a table with many columns. Do you want to proceed?")
      # if user selects 'No', stop the function
      if(proceed==2) {
        message("Operation cancelled. No summary table was generated.")
        return(NULL)
      }
    }
  }

  # estimate tag lifetime dates
  if(!is.null(tag.durations)){
    end_dates <- as.POSIXct(rep(NA, length(tagging.dates)), tz="UTC")
    for(e in 1:length(tagging.dates)){
      end_dates[e] <- tagging.dates[e] + tag.durations[e]*60*60*24
    }
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

  # calculate number of detections
  if("detections" %in% colnames(data)){
    detections <- as.integer(stats::aggregate(as.formula(paste0("detections~", id.col)), data=data, FUN=sum, drop=F)$detections)
  }else{
    warning("No 'detections' column found, assuming one detection per row.", call.=FALSE)
    detections <- as.integer(table(data[,id.col]))
  }
  detections[detections==0] <- NA

  # calculate number of receivers
  receivers <- stats::aggregate(data[,station.col], by=list(data[,id.col]), function(x) length(unique(x)), drop=F)$x


  # determine start dates based on the 'start.point' parameter
  if (start.point=="first.detection") {
    first_detections <- tapply(X=data[, datetime.col], INDEX=data[, id.col], FUN=min, na.rm=TRUE)
    start_dates <- as.POSIXct(first_detections, origin='1970-01-01', tz="UTC")
  } else {
    start_dates <- tagging.dates
  }

  # determine end dates (if tag.duration are not available, use the last dataset detection)
  if (is.null(tag.durations)) {
    end_dates <- rep(max(data[,datetime.col], na.rm=T), length(tagging.dates))
  }

  # days with detections (Dd)
  data$date <- strftime(data[,datetime.col], format="%d-%m-%Y", tz="UTC")
  Dd <- stats::aggregate(data$date, by=list(data[,id.col]), function(x) length(unique(x)), drop=F)$x

  # calculate days between 1st and last detection (Di) - days at liberty
  # + 1 accounts for the inclusive nature of days (i.e., if an animal is tagged on a specific day and last detected on the same day = 1 day at liberty)
  Di <- as.integer(difftime(last_detections, start_dates, units="days")) + 1
  Di[Di<0] <- NA

  # study interval (Dt) - time between first and last monitoring day
  Dt <- as.integer(difftime(end_dates, start_dates, units="days"))
  Dt[Dt==0] <- NA

  # helper function to estimate residency metrics
  calculateResidencyIndex <- function(Dd, Di, Dt, index) {
    if (index=="IR1") return(round(Dd/Di, 2))
    if (index=="IR2") return(round(Dd/Dt, 2))
    if (index=="IWR") return(round(Dd/Dt*Di/Dt, 2))
  }

  # aggregate stats
  stats <- data.frame("ID"=levels(data[,id.col]), "Tagging date"=tag_dates, "Last detection"=last_dates,
                      "N Detect"=detections, "N Receiv"=receivers, "Monitoring duration (d)"=Dt, "Detection span (d)"=Di,
                      "N days detected"=Dd, row.names=NULL, check.names=F)

  # add tag durations if available
  if(!is.null(tag.durations)){
    stats$`Tag duration (d)` <- round(tag.durations, 2)
    stats <- stats[, c("ID", "Tagging date", "Tag duration (d)", "Last detection", "N Detect",
                       "N Receiv", "Monitoring duration (d)", "Detection span (d)", "N days detected")]
  }

  #  calculate the requested residency indexes
  for (index in residency.index) {
    stats[[index]] <- calculateResidencyIndex(Dd, Di, Dt, index)
  }

  # calculate partial residencies if a residency.by variable is supplied
  if (!is.null(residency.by)) {
    data_groupped <- split(data, f=data[,residency.by])
    for (index in residency.index) {
      for (i in seq_along(data_groupped)) {
        data_subset <- data_groupped[[i]]
        Dd_partial <- stats::aggregate(data_subset$date, by=list(data_subset[,id.col]),  function(x) length(unique(x)), drop=F)$x
        stats[paste(index, names(data_groupped)[i])] <- calculateResidencyIndex(Dd_partial, Di, Dt, index)
      }
    }
  }

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
  IR_cols <- which(grepl("IR|IWR", colnames(stats), fixed=FALSE))
  numeric_cols <- which(column_types %in% c("numeric", "integer"))
  numeric_cols <- numeric_cols[!numeric_cols %in% IR_cols]
  decimal_digits <- apply(stats[,numeric_cols], 2, .decimalPlaces)
  decimal_digits <- apply(decimal_digits, 2, max, na.rm=T)
  n_groups <- length(id.groups)
  group_stats <- lapply(id.groups, function(x) stats[stats$ID %in% x,])
  for(i in 1:n_groups){
    group <- group_stats[[i]]
    group[nrow(group)+1,] <- NA
    group$ID[nrow(group)] <- "mean"
    # format numeric columns
    group[nrow(group), numeric_cols] <- sprintf(paste0("%.", decimal_digits, "f"), colMeans(group[,numeric_cols], na.rm=T))
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(group[,numeric_cols], 2, getErrorFun)))
    group[nrow(group), numeric_cols] <- paste(group[nrow(group), numeric_cols], "\u00b1", errors)
    for(c in numeric_cols) group[-nrow(group), c] <- sprintf(paste0("%.", decimal_digits[which(numeric_cols==c)], "f"),  as.numeric(group[-nrow(group), c]))
    # format residency columns
    group[nrow(group), IR_cols] <- sprintf("%.2f", colMeans(group[,IR_cols], na.rm=T))
    errors <- sprintf("%.2f", unlist(apply(group[,IR_cols], 2, getErrorFun)))
    group[nrow(group), IR_cols] <- paste(group[nrow(group), IR_cols], "\u00b1", errors)
    for(c in IR_cols) group[-nrow(group), c] <- sprintf("%.2f", as.numeric(group[-nrow(group), c]))
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

  # aggregate table
  stats <- do.call("rbind", group_stats)
  rownames(stats) <- NULL

  # create new attributes to save relevant params
  attr(stats, 'residency.index') <- residency.index
  attr(stats, 'start.point') <- start.point
  attr(stats, 'residency.by') <- residency.by
  attr(stats, 'id.groups') <- id.groups
  attr(stats, 'processing.date') <- Sys.time()

  # return table
  return(stats)
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
