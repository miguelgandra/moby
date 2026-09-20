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
#' @inheritParams as_moby
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
#' @param residency.index A character vector specifying the residency indices to calculate.
#' Options include:
#'  - "IR1": Residency Index 1, calculated as the number of days the animal was detected (Dd) divided by the detection interval (Di), i.e., the number of days between release/first detection and last detection (days at liberty).
#'  This represents a maximum residency value, considering only the period for which the animal was known to be alive and the tag operational.
#'  - "IR2": Residency Index 2, calculated as the number of days the animal was detected (Dd) divided by the study interval (Dt), i.e., the calendar dates with monitoring opportunity from tagging to the end of receiver monitoring or tag life.
#'  This approach provides a minimum residency value, assuming the animal was alive and detectable throughout the study period.
#'  - "IWR": Weighted Residency Index, which corresponds to the IR2 index weighted by the ratio between the detection interval (Di, from the chosen `start.point` to last detection) and the study interval (Dt, the total monitoring period).
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
#' the end of receiver monitoring as a cutoff timestamp. A cutoff exactly at midnight does not
#' count the new date as a monitoring day (to include all of December 31, use January 1 at
#' midnight). If a single value is provided, it will be applied to all
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
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#'
#' @return A \code{mobyTable}: a TYPED data frame (counts stay integer, indices numeric, dates
#' POSIXct), one row per tagged animal, so the result can be computed on directly. Presentation - fixed
#' precision, the display-only `mean +/- error` row, group headings - is applied by
#' \code{\link[=format.mobyTable]{format}} and \code{\link[=print.mobyTable]{print}}. Export the
#' rendered version with \code{write.csv(format(x), file, row.names = FALSE)}.
#'
#' Columns use stable snake_case names, which is what `format(decimals=)` and `format(group.by=)`
#' are keyed on; the publication headers live in `format(style = "report")`.
#' \itemize{
#'   \item `ID` (or your `id.col`), and a `group` factor when `id.groups` names more than one group
#'   \item `tagging_date`, `last_detection` (POSIXct)
#'   \item `n_detections`, `n_receivers`, `monitoring_duration_d`, `detection_span_d`,
#'     `n_days_detected`
#'   \item one column per requested `residency.index` (`IR1`, `IR2`, `IR2/IR1`, ...), plus the
#'     `residency.by` variants when supplied
#'   \item `<sensor>_mean`, `<sensor>_min`, `<sensor>_max` for each of `sensor.cols`
#'   \item any columns carried over from `id.metadata`
#' }
#' `monitoring_duration_d` is `NA` when only IR1 is requested without a monitoring cutoff.
#'
#' @seealso \code{\link{calculateResidency}} for the underlying numeric residency indices
#' (useful for plotting or downstream statistical analysis).
#' @examples
#' data(rays)
#' # per-animal summary with residency indices, grouped by species
#' summaryTable(rays,
#'              last.monitoring.date = as.POSIXct("2023-12-31", tz = "UTC"),
#'              id.groups = mobyMeta(rays)$id.groups)
#' @export


summaryTable <- function(data,
                         id.metadata = NULL,
                         id.col = NULL,
                         datetime.col = NULL,
                         station.col = NULL,
                         id.groups = NULL,
                         tagging.dates = NULL,
                         tag.durations = NULL,
                         sensor.cols = NULL,
                         sensor.titles = NULL,
                         residency.index = c("IR1", "IR2", "IR2/IR1"),
                         start.point = "release",
                         last.monitoring.date = NULL,
                         residency.by = NULL,
                         error.stat = "sd",
                         verbose = getOption("moby.verbose", TRUE)) {


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
  # all indices involving Dt need a defined monitoring endpoint
  if (any(residency.index %in% c("IR2", "IWR", "IR2/IR1")) && is.null(tag.durations) && is.null(last.monitoring.date)) {
    errors <- c(errors, "The requested index uses the monitoring duration, but no tag durations or last monitoring dates have been provided")
  }
  # print all errors
  if(length(errors)>0){
    stop_message <- c("\n", paste0("- ", errors, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }


  # prompt the user for multiple residency.by levels and residency.index combinations
  if(!is.null(residency.by) && interactive()) {
    nlevels <- length(unique(data[[residency.by]]))
    if((nlevels*length(residency.index))>=9 && length(residency.index)>1) {
      # print prompt
      proceed <- utils::menu(c("Yes", "No"), title="Residency.by has multiple levels and more than one residency index will be calculated. This can result in a table with many columns. Do you want to proceed?")
      # if user selects 'No', stop the function
      # NOTE: deliberately UNGATED (not .mobyNote). This is the direct reply to an action the user
      # just took at an interactive prompt, not verbose narration - declining must never be met with
      # silence, even under verbose = FALSE.
      if(proceed==2) {
        message("Operation cancelled. No summary table was generated.")
        return(NULL)
      }
    }
  }

  # ---- header ---------------------------------------------------------------------------------
  # Placed after the argument checks and after the cancel prompt, so a call that is going to error -
  # or that the user declines - never prints a banner first.
  # This function prints a header and nothing else: the returned table IS the summary, so any
  # closing line would just restate what the user is about to look at. Criteria carry only the
  # choices that change the NUMBERS: which residency index/indices are computed, where the days at
  # liberty start, and which statistic every reported error is (sd vs se).
  crit <- c("residency index" = paste(residency.index, collapse=paste0(" ", .mobyGlyph("mid"), " ")),
            "start point"     = if(start.point=="release") "release date" else "first detection",
            "error"           = if(error.stat=="sd") "standard deviation (sd)" else "standard error (se)")

  .mobyHeader("summaryTable()", "Summarising monitoring and residency metrics per individual",
              input = paste0(.fmtCount(.nDetections(data), "detection"), " ", .mobyGlyph("mid"), " ",
                             .fmtCount(.nObserved(data[,id.col]), "individual")),
              criteria = crit, verbose = verbose)
  .noteUndetected(data[,id.col], verbose, "row retained, values shown as '-'")


  # define error function
  getErrorFun <- function(x) {
    if(error.stat=="sd"){return(sd(x, na.rm=TRUE))}
    if(error.stat=="se"){return(.stdError(x))}
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

  # determine the timezone of the supplied timestamps, so that derived dates are
  # reported in the same timezone as the input data (rather than forced to UTC)
  tz <- .dataTZ(data[,datetime.col])

  # format tagging dates
  # kept as POSIXct: the table is typed, and format() renders the dates on the way out

  # retrieve last detections dates
  last_detections <- tapply(X=data[,datetime.col], INDEX=data[,id.col], FUN=max, na.rm=TRUE)
  last_detections <- as.POSIXct(last_detections, origin='1970-01-01', tz=tz)

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

  # compute residency indices and their temporal building blocks via the dedicated
  # numeric function (single source of truth; this function formats them for display).
  # Name the per-individual vectors so they are matched by ID rather than by position.
  id_levels <- levels(data[,id.col])
  tagging.dates_named <- tagging.dates
  if (length(tagging.dates_named) == length(id_levels)) names(tagging.dates_named) <- id_levels
  tag.durations_named <- tag.durations
  if (!is.null(tag.durations_named) && length(tag.durations_named) == length(id_levels)) names(tag.durations_named) <- id_levels
  last.monitoring.date_named <- last.monitoring.date
  if (!is.null(last.monitoring.date_named) && length(last.monitoring.date_named) == length(id_levels)) names(last.monitoring.date_named) <- id_levels

  residency <- calculateResidency(data, tagging.dates=tagging.dates_named, tag.durations=tag.durations_named,
                                  id.col=id.col, datetime.col=datetime.col,
                                  last.monitoring.date=last.monitoring.date_named,
                                  residency.index=residency.index, start.point=start.point,
                                  residency.by=residency.by, cap=TRUE, verbose=FALSE)
  # align residency rows with the table's ID order
  residency <- residency[match(id_levels, as.character(residency[[id.col]])), , drop=FALSE]
  Dd <- residency$days_detected
  Di <- residency$detection_span
  Dt <- residency$monitoring_duration

  # aggregate stats
  stats <- data.frame("ID" = levels(data[,id.col]),
                      "tagging_date" = tagging.dates,
                      "last_detection" = last_detections,
                      "n_detections" = detections,
                      "n_receivers" = receivers,
                      "monitoring_duration_d" = Dt,
                      "detection_span_d" = Di,
                      "n_days_detected" = Dd,
                      row.names = NULL,
                      check.names = FALSE)

  # add sensor measurement summaries, if available
  if(!is.null(sensor.cols)){
    if(is.null(sensor.titles)) sensor.titles <- tools::toTitleCase(sensor.cols)
    sensor_mean <- stats::aggregate(data[,sensor.cols], by=list(data[,id.col]), function(x){
      if(all(is.na(x))) return(NA) else round(mean(x, na.rm=TRUE), 1)})
    colnames(sensor_mean) <- c("ID", paste0(sensor.titles, "_mean"))
    sensor_min <- stats::aggregate(data[,sensor.cols], by=list(data[,id.col]), function(x){
      if(all(is.na(x))) return(NA) else round(min(x, na.rm=TRUE), 1)})
    colnames(sensor_min) <- c("ID", paste0(sensor.titles, "_min"))
    sensor_max <- stats::aggregate(data[,sensor.cols], by=list(data[,id.col]), function(x){
      if(all(is.na(x))) return(NA) else round(max(x, na.rm=TRUE), 1)})
    colnames(sensor_max) <- c("ID", paste0(sensor.titles, "_max"))
    sensor_stats <- Reduce(function(x, y) .joinKeep(x, y, by="ID", type="left"), list(sensor_mean, sensor_min, sensor_max))
    # merge and reorder columns
    ordered_cols <- c("ID", unlist(lapply(sensor.titles, function(x) paste0(x, c("_mean", "_min", "_max")))))
    sensor_stats <- sensor_stats[,ordered_cols]
    col_index <- which(colnames(stats)=="last_detection")
    first_cols <- names(stats)[seq_len(col_index)]
    last_cols <- names(stats)[(col_index+1):ncol(stats)]
    stats <- .joinKeep(stats, sensor_stats, by="ID", type="left")
    stats <- stats[,c(first_cols, ordered_cols[-1], last_cols)]
  }

  #  add the requested residency indexes (computed above by calculateResidency)
  for (index in residency.index) {
    stats[[index]] <- residency[[index]]
  }

  # add partial (spatially structured) residencies, if a residency.by variable was supplied
  if (!is.null(residency.by)) {
    building_blocks <- c(id.col, "tagging_date", "first_detection", "last_detection",
                         "monitoring_end", "days_detected", "detection_span", "monitoring_duration")
    partial_cols <- setdiff(colnames(residency), c(building_blocks, residency.index))
    for (pc in partial_cols) stats[[pc]] <- residency[[pc]]
  }

  # format additional tag info (if available) and merge
  if(!is.null(id.metadata)) {
    column_types <- sapply(seq_len(ncol(id.metadata)),function(c) class(id.metadata[,c]))
    numeric_cols <- which(column_types %in% c("numeric", "integer"))
    # drop=FALSE keeps id.metadata[,-1] a data.frame when only one non-ID column remains; otherwise it
    # collapses to a vector and aggregate() names the result "x", silently losing the real column name.
    animal_info <- stats::aggregate(id.metadata[,-1, drop=FALSE], by=list(id.metadata[,id.col]), function(x) paste(unique(x), collapse="/"), drop=FALSE)
    animal_info[animal_info=="NA"] <- NA
    # coerce column-wise: as.numeric() on a multi-column data.frame subset would error ("'list'
    # object cannot be coerced to type 'double'") when id.metadata has 2+ numeric non-ID columns.
    animal_info[numeric_cols] <- lapply(animal_info[numeric_cols], as.numeric)
    colnames(animal_info)[1] <- "ID"
    stats <- .joinKeep(animal_info, stats, by="ID", type="left")
  }

  # ---- typed output ---------------------------------------------------------------------------
  # The table returns its NUMBERS. Precision, the mean row, the group headings and the "-" for
  # missing all belong to format()/print(), so the object stays something you can compute on.
  stats$ID <- as.character(stats$ID)
  n_groups <- length(id.groups)

  # A `group` column when the caller declared named ID groups. Previously the split was rendered as
  # blank heading rows spliced into the table, which exported to CSV as empty records; as a column it
  # is both analysable and what format(group.by=) groups on.
  if (n_groups > 1) {
    grp <- rep(NA_character_, nrow(stats))
    for (i in seq_len(n_groups)) {
      grp[stats$ID %in% as.character(id.groups[[i]])] <- names(id.groups)[i]
    }
    stats$group <- factor(grp, levels = names(id.groups))
    # rows are emitted in group order, as they were when each group was its own block
    stats <- stats[order(stats$group, seq_len(nrow(stats))), , drop = FALSE]
  }

  # the ID column carries the caller's own name
  if (!identical(id.col, "ID")) colnames(stats)[colnames(stats) == "ID"] <- id.col

  .newMobyTable(stats, kind = "summary", error.stat = error.stat, label.col = id.col,
                extra = list(residency.index = residency.index,
                             start.point = start.point,
                             residency.by = residency.by,
                             id.groups = id.groups,
                             processing.date = Sys.time()))
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
