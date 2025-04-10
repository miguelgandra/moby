#######################################################################################################
## Create summary table ###############################################################################
#######################################################################################################

#' Generate summary table for tagged animals
#'
#' @description This function generates a comprehensive summary table for tagged animals,
#' including tagging and last detection dates, along with various monitoring and residency metrics.
#' It also allows the inclusion of additional sensor data (e.g., depth, temperature) and provides
#' an option to summarize these readings by calculating the mean, minimum, and maximum values for each specified sensor column.
#' Additionally, the function can incorporate additional metadata and automatically
#' computes overall means and error metrics.
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections. Each row should represent an individual detection event,
#' unless a 'detections' column is included to indicate the number of detections for each row.
#' @param id.metadata A data frame containing metadata about the tagged animals, such as their length,
#' sex, or transmitter type. All columns in this data frame will be summarized and included
#' in the final table. If there are multiple rows per animal, variables will be collapsed
#' before merging with other statistics.
#' @param tag.durations Optional. A numeric vector containing the estimated battery
#' duration of the deployed tags (in days). This parameter must be either:
#' - A single numeric value, which will be applied to all unique animal IDs; or
#' - A named numeric vector, where the names correspond to the animal IDs in the `id.col` column.
#' If multiple tag durations are provided, the vector must include all IDs and will be reordered to align with the levels of `id.col`.
#' @param sensor.cols Optional. A character vector specifying column names in `data` that contain
#' additional sensor readings (e.g., depth, temperature). For each column specified, the mean,
#' minimum, and maximum values will be calculated and incorporated into the summary.
#' @param sensor.titles Optional. Titles for the sensor readings in `sensor.cols`.
#' If not provided, defaults to the column names.
#' @param residency.index A character string specifying the type of residency index to calculate.
#' Options include:
#'  - "IR1": Residency Index 1, calculated as the number of days the animal was detected (Dd) divided by the detection interval (Di), i.e., the number of days between release/first detection and last detection (days at liberty).
#'  This represents a maximum residency value, considering only the period for which the animal was known to be alive and the tag operational.
#'  - "IR2": Residency Index 2, calculated as the number of days the animal was detected (Dd) divided by the study interval (Dt), i.e., the total number of days between release/first detection and last data download or tag expiration date.
#'  This approach provides a minimum residency value, assuming the animal was alive and detectable throughout the study period.
#'  - "IWR": Weighted Residency Index, which corresponds to the IR2 index weighted by the ratio between the detection interval (Di, the number of days between the first and last detection) and the study interval (Dt, the total monitoring period).
#'  This accounts for the number of days detected and the spread of detections within the monitoring period, providing a measure of residency that balances the frequency of detections with their temporal distribution.
#'  - "IR2/IR1": The ratio of IR2 to IR1, providing a measure of the gap between the last detection and the end of the monitoring period.
#'
#' The choice of index can affect the interpretation of residency patterns, so it's important to select the one(s) that best fits the study objectives.
#' Further information on residency estimation methods can be found in Kraft et al. (2023) and Appert et al. (2023) - see below in the references section.
#' @param start.point A character string specifying the starting point for calculating days at liberty.
#' Options include:
#' - `"release"`: The release day is used as the starting point for calculating days at liberty.
#' - `"first.detection"`: The first detection after tagging is used as the starting point.
#'
#' Defaults to `"release"`.
#' @param last.monitoring.date Optional. A POSIXct object or a named vector of POSIXct objects specifying
#' the last timestamp when data could be retrieved, typically corresponding to the last data download date
#' or the final day receivers were operational. If a single value is provided, it will be applied to all
#' individuals. If a named vector is provided, the names should correspond to individual IDs, allowing
#' for unique timestamps per individual. When `tag.durations` are also supplied, the total monitoring
#' duration for each individual will be estimated based on the shortest of the two values: the tag expiration
#' date or the last monitoring day.
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
#' - `Last detection`: The date of the last detection.
#' - `N Detect`: Total number of detections for the animal.
#' - `N Receiv`: Number of unique receivers that detected the animal.
#' - `Monitoring duration (d)`: Total duration of monitoring in days. This is determined
#' by the tag duration, if provided, or alternatively calculated as the time
#' between release and the last detection in the dataset (assumed to represent the
#' final data download).
#' - `Detection span (d)`: Number of days between release/first detection and last detection (days at liberty)
#' - `N days detected`: Total number of days the animal was detected.
#' - Additional columns for each residency index specified in the `residency.index` parameter.
#' - If `residency.by` is specified, additional columns for partial residency metrics will be included.
#' - Sensor data metrics: For each column in `sensor.cols`, the mean, minimum, and maximum values, using titles specified in `sensor.titles` (if provided).
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
                         id.col = getDefaults("ID"),
                         datetime.col = getDefaults("datetime"),
                         station.col = getDefaults("station"),
                         sensor.cols = NULL,
                         sensor.titles = NULL,
                         residency.index = c("IR1", "IR2", "IR2/IR1"),
                         start.point = "release",
                         last.monitoring.date = NULL,
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
  last.monitoring.date <- reviewed_params$last.monitoring.date

  # validate additional parameters
  errors <- c()
  # check if data contains residency.by
  if(!is.null(residency.by) && !residency.by %in% colnames(data)) {
    errors <- c(errors, "Variable used to calculate partial residencies not found in the supplied data. Please specify the correct column using 'residency.by'.")
  }
  # check error function
  error.stat <- tolower(error.stat)
  if(!error.stat %in% c("sd", "se")) {
    errors <- c(errors, "Wrong error.stat argument, please choose between 'sd' and 'se'.")
  }
  # check if id.metadata contains id.col
  if(!is.null(id.metadata) && !id.col %in% colnames(id.metadata)) {
    errors <- c(errors, "The specified ID column ('id.col') does not exist in 'id.metadata'. Please ensure that the column name in 'id.metadata' matches the 'id.col' specified.")
  }
  # check if data contains sensor.cols
  if(!is.null(sensor.cols) && !all(sensor.cols %in% colnames(data))) {
    errors <- c(errors, "One or more specified sensor columns ('sensor.cols') were not found in the supplied data. Please check the column names and ensure they exist in the data.")
  }
  # check residency index
  if(!all(residency.index %in% c("IR1", "IR2", "IWR", "IR2/IR1"))) {
    errors <- c(errors, "Invalid 'residency.index' argument. Please select one of the following options: 'IR1' (detection interval), 'IR2' (study interval), 'IWR' (weighted residency index) or 'IR2/IR1' (quotient).")
  }
  # check start.point
  if(!is.character(start.point) || length(start.point) !=1) {
    errors <- c(errors, "Invalid 'start.point' parameter: it must be a single character string.")
  }
  if(!start.point %in% c("first.detection", "release")) {
    errors <- c(errors, "Invalid 'start.point' parameter: must be one of 'release' or 'first.detection'.")
  }
  # ensure required durations are provided for IR2 and IWR indices
  if (any(residency.index %in% c("IR2", "IWR")) && is.null(tag.durations) && is.null(last.monitoring.date)) {
    errors <- c(errors, "The residency index includes 'IR2' or 'IWR', but no tag durations or last monitoring dates have been provided")
  }
  # print all errors
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
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

  # define error function
  getErrorFun <- function(x) {
    if(error.stat=="sd"){return(sd(x, na.rm=TRUE))}
    if(error.stat=="se"){return(plotrix::std.error(x))}
  }

  # check for NA values in the datetime column
  if (any(is.na(data[, datetime.col]))) {
    warning(paste("- NA values detected in the", datetime.col, "column."), call. = FALSE)
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
  last_detections <- tapply(X=data[,datetime.col], INDEX=data[,id.col], FUN=max, na.rm=TRUE)
  last_detections <- as.POSIXct(last_detections, origin='1970-01-01', tz="UTC")
  last_dates <- strftime(last_detections, format="%d/%m/%Y", tz="UTC")

  # calculate number of detections
  if("detections" %in% colnames(data)){
    detections <- as.integer(stats::aggregate(as.formula(paste0("detections~", id.col)), data=data, FUN=sum, drop=FALSE)$detections)
  }else{
    warning("- No 'detections' column found, assuming one detection per row.", call.=FALSE)
    detections <- as.integer(table(data[,id.col]))
  }
  detections[detections==0] <- NA

  # calculate number of receivers
  receivers <- stats::aggregate(data[,station.col], by=list(data[,id.col]), function(x) length(unique(x)), drop=FALSE)$x

  # determine start dates based on the 'start.point' parameter
  if (start.point=="first.detection") {
    first_detections <- tapply(X=data[, datetime.col], INDEX=data[, id.col], FUN=min, na.rm=TRUE)
    start_dates <- as.POSIXct(first_detections, origin='1970-01-01', tz="UTC")
  } else {
    start_dates <- tagging.dates
  }

  # determine end dates (using the shortest duration between tag expiration and last monitoring date)
  if (!is.null(tag.durations)) {
    # calculate tag expiration dates based on tagging dates and durations
    end_dates <- as.POSIXct(rep(NA, length(tagging.dates)), tz = "UTC")
    tag_expiration_dates <- tagging.dates + tag.durations * 60 * 60 * 24
    # use the shortest value: tag expiration date or last monitoring date (if available)
    for (e in 1:length(tagging.dates)) {
      if (!is.null(last.monitoring.date)) end_dates[e] <- min(tag_expiration_dates[e], last.monitoring.date[e], na.rm = TRUE)
      else end_dates[e] <- tag_expiration_dates[e]
    }
  } else if (!is.null(last.monitoring.date)) {
    # use last monitoring dates if no tag durations are provided
    end_dates <- last.monitoring.date
  } else {
    # fallback to the maximum detection time in the dataset
    end_dates <- rep(max(data[, datetime.col], na.rm = TRUE), length(tagging.dates))
    warning_msg <- paste("- Neither `tag.durations` nor `last.monitoring.date` were provided.",
                  "The function will default to assuming the last detection date",
                  "in the dataset as the end of the monitoring period for all individuals.")
    warning(paste(strwrap(warning_msg, width=getOption("width")), collapse="\n"), call.=FALSE)
  }

  # days with detections (Dd)
  data$date <- strftime(data[,datetime.col], format="%d-%m-%Y", tz="UTC")
  Dd <- stats::aggregate(data$date, by=list(data[,id.col]), function(x) length(unique(x)), drop=FALSE)$x

  # calculate days between 1st and last detection (Di) - days at liberty
  # + 1 accounts for the inclusive nature of days (i.e., if an animal is tagged on a specific day and last detected on the same day = 1 day at liberty)
  Di <- as.integer(difftime(last_detections, start_dates, units="days")) + 1
  Di[Di<0] <- NA

  # study interval (Dt) - time between release and last monitoring day
  Dt <- as.integer(difftime(end_dates, tagging.dates, units="days"))
  Dt[Dt==0] <- NA

  # helper function to estimate residency metrics
  calculateResidencyIndex <- function(Dd, Di, Dt, metric) {
    formulas <- list(IR1 = Dd / Di,
                     IR2 = Dd / Dt,
                     IWR = (Dd / Dt) * (Di / Dt),
                     `IR2/IR1` = (Dd / Dt) / (Dd / Di))
    result <- formulas[[metric]]
    return(pmin(round(result, 2), 1, na.rm = FALSE))
  }

  # aggregate stats
  stats <- data.frame("ID" = levels(data[,id.col]),
                      "Tagging date" = tag_dates,
                      "Last detection" = last_dates,
                      "N Detect" = detections,
                      "N Receiv" = receivers,
                      "Monitoring duration (d)" = Dt,
                      "Detection span (d)" = Di,
                      "N days detected" = Dd,
                      row.names = NULL,
                      check.names = FALSE)

  # add sensor measurement summaries, if available
  if(!is.null(sensor.cols)){
    if(is.null(sensor.titles)) sensor.titles <- tools::toTitleCase(sensor.cols)
    sensor_mean <- stats::aggregate(data[,sensor.cols], by=list(data[,id.col]), function(x){
      if(all(is.na(x))) return(NA) else round(mean(x, na.rm=TRUE), 1)})
    colnames(sensor_mean) <- c("ID", paste(sensor.titles, "- mean"))
    sensor_min <- stats::aggregate(data[,sensor.cols], by=list(data[,id.col]), function(x){
      if(all(is.na(x))) return(NA) else round(min(x, na.rm=TRUE), 1)})
    colnames(sensor_min) <- c("ID", paste(sensor.titles, "- min"))
    sensor_max <- stats::aggregate(data[,sensor.cols], by=list(data[,id.col]), function(x){
      if(all(is.na(x))) return(NA) else round(max(x, na.rm=TRUE), 1)})
    colnames(sensor_max) <- c("ID", paste(sensor.titles, "- max"))
    sensor_stats <- Reduce(function(x, y) plyr::join(x, y, by="ID", type="left"), list(sensor_mean, sensor_min, sensor_max))
    # merge and reorder columns
    ordered_cols <- c("ID", unlist(lapply(sensor.titles, function(x) c(paste(x, "- mean"), paste(x, "- min"), paste(x, "- max")))))
    sensor_stats <- sensor_stats[,ordered_cols]
    col_index <- which(colnames(stats)=="Last detection")
    first_cols <- names(stats)[1:col_index]
    last_cols <- names(stats)[(col_index+1):ncol(stats)]
    stats <- plyr::join(stats, sensor_stats, by="ID", type="left")
    stats <- stats[,c(first_cols, ordered_cols[-1], last_cols)]
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
        Dd_partial <- stats::aggregate(data_subset$date, by=list(data_subset[,id.col]),  function(x) length(unique(x)), drop=FALSE)$x
        stats[paste(index, names(data_groupped)[i])] <- calculateResidencyIndex(Dd_partial, Di, Dt, index)
      }
    }
  }

  # format additional tag info (if available) and merge
  if(!is.null(id.metadata)) {
    column_types <- sapply(1:ncol(id.metadata),function(c) class(id.metadata[,c]))
    numeric_cols <- which(column_types %in% c("numeric", "integer"))
    animal_info <- stats::aggregate(id.metadata[,-1], by=list(id.metadata[,id.col]), function(x) paste(unique(x), collapse="/"), drop=FALSE)
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
  decimal_digits <- apply(decimal_digits, 2, max, na.rm=TRUE)
  n_groups <- length(id.groups)
  group_stats <- lapply(id.groups, function(x) stats[stats$ID %in% x,])
  for(i in 1:n_groups){
    group <- group_stats[[i]]
    group[nrow(group)+1,] <- NA
    group$ID[nrow(group)] <- "mean"
    # format numeric columns
    group[nrow(group), numeric_cols] <- sprintf(paste0("%.", decimal_digits, "f"), colMeans(group[,numeric_cols, drop=FALSE], na.rm=TRUE))
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(group[,numeric_cols, drop=FALSE], 2, getErrorFun)))
    group[nrow(group), numeric_cols] <- paste(group[nrow(group), numeric_cols], "\u00b1", errors)
    for(c in numeric_cols) group[-nrow(group), c] <- sprintf(paste0("%.", decimal_digits[which(numeric_cols==c)], "f"),  as.numeric(group[-nrow(group), c]))
    # format residency columns
    group[nrow(group), IR_cols] <- sprintf("%.2f", colMeans(group[,IR_cols, drop=FALSE], na.rm=TRUE))
    errors <- sprintf("%.2f", unlist(apply(group[,IR_cols, drop=FALSE], 2, getErrorFun)))
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
    group_stats <- mapply(function(label, stats) {rbind(label, stats)}, label=group_labels, stats=group_stats, SIMPLIFY=FALSE)
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
