#######################################################################################################
## Validate Arguments #################################################################################
#######################################################################################################

#' Validate Arguments
#'
#' This internal function performs a comprehensive set of validation checks on the arguments
#' provided to functions within the `moby` package. It ensures that the data and various parameters
#' are in the correct format and meet the expected criteria, issuing errors and warnings as necessary.
#' The function tries to issue as many errors and warnings at the same time as possible, to speed up the
#' completion of required changes and reduce the time it takes to get the function call right.
#'
#' @return A list containing potentially modified versions of `data`, `tagging.dates`, and `tag.durations`.
#'         If validation fails, the function will stop and return an error message.
#'
#' @details
#' The `.validateArguments` function performs rigorous checks to ensure that all input arguments
#' are in the correct format, of the correct type, and consistent with the requirements of the `moby` package.
#' This includes checks for data structure, column specifications, and additional parameters. The validation
#' aims to provide detailed error and warning messages to guide users in correcting their input.
#'
#'#' Specifically, the following checks are performed:
#' \describe{
#'   \item{\code{data}}{Ensures the first argument is a data frame.}
#'   \item{\code{id.col}}{Validates that the ID column exists in the data frame, is unique for each individual, and is converted to a factor if necessary.}
#'   \item{\code{id.groups}}{Ensures ID groups are unique and correctly referenced in the data.}
#'   \item{\code{datetime.col}}{Checks that the datetime column exists and is in POSIXct format.}
#'   \item{\code{timebin.col}}{Validates the existence of the time-bin column and checks that it is in POSIXct format.}
#'   \item{\code{lon.col} and \code{lat.col}}{Ensures the longitude and latitude columns exist and contain numeric values.}
#'   \item{\code{station.col}}{Checks for the existence of the station column.}
#'   \item{\code{spatial.col}}{Validates spatial column(s), ensuring they exist and match the required format for spatial analyses.}
#'   \item{\code{dist.col}}{Checks the distance column, ensuring it contains numeric values.}
#'   \item{\code{split.by}}{Validates the grouping variable(s) used to split the data.}
#'   \item{\code{subset}}{Checks that subset criteria are valid and appropriately referenced in the data.}
#'   \item{\code{variable}}{Ensures a single variable is specified and exists in the data.}
#'   \item{\code{variables}}{Validates that all specified variables exist in the data.}
#'   \item{\code{color.by}}{Checks that the color-by column exists, contains no missing values, and is converted to a factor if needed.}
#'   \item{\code{tagging.dates}}{Ensures tagging dates are in POSIXct format, match the unique IDs, and are appropriately named if provided as a vector.}
#'   \item{\code{tag.durations}}{Validates tag durations, ensuring they are numeric and consistent with the unique IDs.}
#'   \item{\code{start.dates} and \code{end.dates}}{Checks start and end dates, ensuring they are in POSIXct format and consistent with the data.}
#'   \item{\code{cutoff.dates}}{Validates cutoff dates, ensuring they are in POSIXct format and logically set.}
#'   \item{\code{last.monitoring.date}}{Ensures the last monitoring date is in POSIXct format and appropriately set relative to the data.}
#'   \item{\code{land.shape}}{Checks that the provided spatial object is either an \code{sf} object or convertible to one.}
#'   \item{\code{color.pal}}{Validates the color palette, ensuring it is a character vector or a function.}
#'   \item{\code{diel.lines}}{Checks that diel line values are valid (e.g., 0, 2, or 4).}
#'   \item{\code{polygons}}{Ensures polygons are validly set to \code{'diel'}, \code{'season'}, or \code{FALSE}.}
#'   \item{\code{background.color}}{Checks that the background color is a valid named color or hexadecimal code.}
#'   \item{\code{style}}{Validates that the style parameter is either \code{'raster'} or \code{'points'}.}
#'   \item{\code{cores}}{Ensures that the number of cores for parallel processing is specified as a positive integer. If this parameter is missing or invalid, a default value of 1 is used.}
#'   \item{\code{cols}}{Checks that the \code{cols} parameter is correctly formatted as a single positive integer.}
#' }
#'
#' @note
#' This function is intended for internal use within the `moby` package. It helps streamline
#' argument validation across multiple functions, ensuring consistency and reducing redundancy.
#' @keywords internal

.validateArguments <- function() {


  ##############################################################################
  # Retrieve arguments #########################################################
  ##############################################################################

  calling_fun <- deparse(sys.call(-1)[[1]])
  args <- mget(names(formals(sys.function(sys.parent()))), sys.frame(sys.nframe()-1L))
  data <- args[["data"]]

  valid_ids <- FALSE
  color.by <- NULL
  tagging.dates <- NULL
  tag.durations <- NULL
  start.dates <- NULL
  end.dates <- NULL
  cutoff.dates <- NULL
  last.monitoring.date <- NULL
  land.shape <- NULL

  # initialize strings to return all error and warning messages at once
  errors <- c()
  warnings <- c()


  ##############################################################################
  # Perform checks ############################################################
  ##############################################################################

  ##############################################################################
  # validate data ##############################################################
  if (!inherits(data, "data.frame")) {
    stop(paste0(names(args)[1], " must be a data frame."), call.=FALSE)
  }

  ##############################################################################
  # validate id.col ############################################################
  if ("id.col" %in% names(args)) {
    id.col <- args$id.col
    errors <- c(errors, .checkColumn(id.col, "ID column", data))
    if (is.null(errors)) valid_ids <- TRUE
    # convert to factor
    if (valid_ids && !inherits(data[, id.col], "factor")){
      data[, id.col] <- as.factor(data[, id.col])
      warnings <- c(warnings, "'id.col' converted to factor.")
    }
  }

  ##############################################################################
  # validate id.groups #########################################################
  if ("id.groups" %in% names(args) && valid_ids) {
    id.groups <- args$id.groups
    if(!is.null(id.groups)){
      if(!inherits(id.groups, "list") || is.null(names(id.groups)) ) errors <- c(errors, "id.groups should be supplied as a named list.")
      else{
        empty_groups <- sapply(id.groups, function(group) length(group) == 0 || is.null(group))
        if(any(empty_groups)) {
          empty_ids <- names(id.groups)[empty_groups]
          errors <- c(errors, paste("The following ID groups are empty or NULL:", paste(empty_ids, collapse = ", ")))
        }else if(any(duplicated(unlist(id.groups)))){
          errors <- c(errors, "Repeated ID(s) in id.groups.")
        }else{
          if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {warnings <- c(warnings, "Some of the ID(s) in id.groups don't match the IDs in the data.")}
          data <- data[data[,id.col] %in% unlist(id.groups),]
          data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
        }
      }
    }
  }

  ##############################################################################
  # validate datetime.col ######################################################
  if ("datetime.col" %in% names(args)) {
    datetime.col <- args$datetime.col
    datetime_msg <- .checkColumn(datetime.col, "datetime column", data)
    if(!is.null(datetime_msg))  errors <- c(errors, datetime_msg)
    else if (!inherits(data[, datetime.col], "POSIXct")) errors <- c(errors, "Datetimes must be provided in POSIXct format.")
  }

  ##############################################################################
  # validate timebin.col   #####################################################
  if ("timebin.col" %in% names(args)) {
    timebin.col <- args$timebin.col
    timebin_msg <- .checkColumn(timebin.col, "time-bin column", data)
    if(!is.null(timebin_msg))  errors <- c(errors, timebin_msg)
    else if (!inherits(data[, timebin.col], "POSIXct")) errors <- c(errors, "Time-bins must be provided in POSIXct format.")
  }

  ##############################################################################
  # validate lon.col   #####################################################
  if ("lon.col" %in% names(args)) {
    lon.col <- args$lon.col
    lon_msg <- .checkColumn(lon.col, "longitude column", data)
    if(!is.null(lon_msg))  errors <- c(errors, lon_msg)
    else if (!is.numeric(data[[lon.col]])) errors <- c(errors, paste("The column", lon.col, "must be numeric."))
  }

  ##############################################################################
  # validate lat.col   #####################################################
  if ("lat.col" %in% names(args)) {
    lat.col <- args$lat.col
    lat_msg <- .checkColumn(lat.col, "latitude column", data)
    if(!is.null(lat_msg))  errors <- c(errors, lat_msg)
    else if (!is.numeric(data[[lat.col]])) errors <- c(errors, paste("The column", lat.col, "must be numeric."))
  }


  ##############################################################################
  # validate station.col #######################################################
  if ("station.col" %in% names(args)) {
    station.col <- args$station.col
    if(!is.null(station.col)){
      errors <- c(errors, .checkColumn(station.col, "station column", data))
    }
  }

  ##############################################################################
  # validate spatial.col #######################################################
  if ("spatial.col" %in% names(args)) {
    spatial.col <- args$spatial.col
    errors <- c(errors, .checkColumn(spatial.col, "spatial column", data))
  }

  ##############################################################################
  # validate dist.col ##########################################################
  if ("dist.col" %in% names(args)) {
    dist.col <- args$dist.col
    if(!is.null(dist.col)) {
      errors <- c(errors, .checkColumn(dist.col, "distance column", data))
    }
  }


  ##############################################################################
  # validate split.by ##########################################################
  if ("split.by" %in% names(args)) {
    split.by <- args$split.by
    if(!is.null(split.by)) {
      for(s in 1:length(split.by)){
        errors <- c(errors, .checkColumn(split.by[s], "grouping variable", data))
      }
    }
  }

  ##############################################################################
  # validate subset ############################################################
  if ("subset" %in% names(args)) {
    subset <- args$subset
    if(!is.null(subset)) errors <- c(errors, .checkColumn(subset, "subset column", data))
  }

  ##############################################################################
  # validate variable ##########################################################
  if ("variable" %in% names(args)) {
    variable <- args$variable
    if(!is.null(variable)) errors <- c(errors, .checkColumn(variable, "variable", data))
  }

  ##############################################################################
  # validate variables ##########################################################
  if ("variables" %in% names(args)) {
    all_variables <- args$variables
    if(!is.null(all_variables) && calling_fun!="plotChronogram"){
      for(v in 1:length(all_variables)){
        variables <- all_variables[v]
        errors <- c(errors, .checkColumn(variables, "variable", data))
      }
    }
  }

  ##############################################################################
  # validate color.by ##########################################################
  if ("color.by" %in% names(args)) {
    color.by <- args$color.by
    if(!is.null(color.by)){
      colorby_msg <- .checkColumn(color.by, "color variable", data)
      if(!is.null(colorby_msg)) errors <- c(errors, colorby_msg)
      else if(any(is.na(data[, color.by]))) errors <- c(errors, "Missing values in color.by variable.")
      else if(inherits(data[,color.by], "character")) {
        data[,color.by] <- as.factor(data[,color.by])
        warnings <- c(warnings, "'color.by' variable converted to factor.")}
    }
  }

  ##############################################################################
  # validate tagging.dates   ###################################################
  if ("tagging.dates" %in% names(args)) {
    tagging.dates <- args$tagging.dates
    if (is.null(tagging.dates)){
      errors <- c(errors, "Tagging dates not specified. Please provide the required POSIXct values via the 'tagging.dates' argument.")
    } else {
      check_result <- .checkAnimalParams(tagging.dates, "Tagging dates", expected_class="POSIXct", data, id.col)
      tagging.dates <- check_result$vector
      errors <- c(errors, check_result$errors)
    }
  }

  ##############################################################################
  # validate tag durations #####################################################
  if ("tag.durations" %in% names(args)) {
    tag.durations <- args$tag.durations
    if(!is.null(tag.durations)){
      check_result <- .checkAnimalParams(tag.durations, "Tag durations", expected_class="numeric", data, id.col)
      tag.durations <- check_result$vector
      errors <- c(errors, check_result$errors)
    }
  }

  ##############################################################################
  # validate start dates #######################################################
  if ("start.dates" %in% names(args)) {
    start.dates <- args$start.dates
    if(!is.null(start.dates)){
      check_result <- .checkAnimalParams(start.dates, "Start dates", expected_class="POSIXct", data, id.col)
      start.dates <- check_result$vector
      errors <- c(errors, check_result$errors)
    }
  }

  ##############################################################################
  # validate end dates #########################################################
  if ("end.dates" %in% names(args)) {
    end.dates <- args$end.dates
    if(!is.null(end.dates)){
      check_result <- .checkAnimalParams(end.dates, "Start dates", expected_class="POSIXct", data, id.col)
      end.dates <- check_result$vector
      errors <- c(errors, check_result$errors)
    }
  }

  ##############################################################################
  # validate cutoff dates #########################################################
  if ("cutoff.dates" %in% names(args)) {
    cutoff.dates <- args$cutoff.dates
    if(!is.null(cutoff.dates)){
      check_result <- .checkAnimalParams(cutoff.dates, "Cut-off dates", expected_class="POSIXct", data, id.col)
      cutoff.dates <- check_result$vector
      errors <- c(errors, check_result$errors)
    }
  }


  ##############################################################################
  # validate last monitoring date ##############################################
  if ("last.monitoring.date" %in% names(args)) {
    last.monitoring.date <- args$last.monitoring.date
    if(!is.null(last.monitoring.date)){
      check_result <- .checkAnimalParams(last.monitoring.date, "Last monitoring date", expected_class="POSIXct", data, id.col)
      last.monitoring.date <- check_result$vector
      errors <- c(errors, check_result$errors)
    }
  }


  ##############################################################################
  # validate land.shape ########################################################
  if ("land.shape" %in% names(args)) {
    land.shape <- args$land.shape

    # ensure 'land.shape' is not NULL
    if(!is.null(land.shape)){
      # check if the object is already of class 'sf'
      if(!inherits(land.shape, "sf")) {
        # if 'land.shape' is not of class 'sf', check if it can be converted
        if(!inherits(land.shape, c("SpatialPolygonsDataFrame", "SpatialPolygons"))) {
          errors <- c(errors, "The provided 'land.shape' could not be converted to class 'sf'. Please provide a valid spatial object.")
        }else{
          # attempt to convert 'land.shape' to an 'sf' object
          tryCatch({
            land.shape <- sf::st_as_sf(land.shape)
            warnings <- c(warnings, "The provided 'land.shape' was internally converted to 'sf' for spatial processing.")
          }, error = function(e) {
            errors <- c(errors, "The provided 'land.shape' could not be converted to class 'sf'. Please provide a valid spatial object.")
          })
        }
      }
    }
  }

  ##############################################################################
  # validate color.pal   #######################################################
  if ("color.pal" %in% names(args)) {
    color.pal <- args$color.pal
    if(!is.null(color.pal) && !inherits(color.pal, c("character", "function")))
      errors <- c(errors, "Invalid color palette. Please supply either a character vector w/ the colors or a function that generates colors.")
  }


  ##############################################################################
  # validate diel.lines ########################################################
  if ("diel.lines" %in% names(args)) {
    diel.lines <- args$diel.lines
    if(!diel.lines %in% c(0, 2, 4)) errors <- c(errors, "Invalid 'diel.lines' specified. Accepted values are '0' (no lines),
                                                  2' (sunrise and sunset) or '4' (dawn, sunrise, sunset, dusk).")
  }

  ##############################################################################
  # validate polygons ##########################################################
  if ("polygons" %in% names(args)) {
    polygons <- args$polygons
    if(!polygons %in% c(FALSE, "diel","season")) errors <- c(errors,"Invalid 'polygons' specified. Accepted values are: 'diel', 'season' or 'FALSE' (disabled).")
  }

  ##############################################################################
  # Validate background.color ##################################################
  if ("background.color" %in% names(args)) {
    background.color <- args$background.color
    if (length(background.color) != 1) {
      errors <- c(errors, "Invalid background color specified. The background color must be a single value.")
    } else if (!is.character(background.color)) {
      errors <- c(errors, "Invalid background color specified. The background color must be a string.")
    } else if (!grepl("^#[0-9A-Fa-f]{6}$", background.color) && !background.color %in% colors()) {
      errors <- c(errors, "Invalid background color specified. The background color must be a valid hex code or named color.")
    }
  }

  ##############################################################################
  # Validate 'style' argument###################################################
  if ("style" %in% names(args)) {
    style <- args$style
    if(!style %in% c("raster", "points")) errors <- c(errors, "Invalid style specified. Accepted values are: 'raster' or 'points'.")
  }


  ##############################################################################
  # Validate 'cores' argument for parallel computing ###########################
  if ("cores" %in% names(args)) {
    cores <- args$cores
    if (cores>1){
      if (!is.numeric(cores) || cores < 1 || cores %% 1 != 0) errors <- c(errors, "The 'cores' parameter must be a positive integer.")
      if (!requireNamespace("foreach", quietly=TRUE)) errors <- c(errors, "The 'foreach' package is required for parallel computing but is not installed. Please install 'foreach' using install.packages('foreach') and try again.")
      if (!requireNamespace("doSNOW", quietly=TRUE)) errors <- c(errors, "The 'doSNOW' package is required for parallel computing but is not installed. Please install 'doSNOW' using install.packages('doSNOW') and try again.")
      if (!requireNamespace("parallel", quietly=TRUE)){
        errors <- c(errors, "The 'parallel' package is required for parallel computing but is not installed. Please install 'parallel' using install.packages('parallel') and try again.")
      }else{
          if(is.numeric(cores) && parallel::detectCores()<cores) errors <- c(errors, paste("Please choose a different number of cores for parallel computing (only", parallel::detectCores(), "available)."))
      }
    }
  }


  ##############################################################################
  # Validate 'cols' argument ###################################################
  if ("cols" %in% names(args)) {
    cols <- args$cols
    if (!is.null(cols)){
      if (!is.numeric(cols) || length(cols) > 1 || cols %% 1 != 0) errors <- c(errors, "The 'cols' parameter (number of columns) must be a single positive integer.")
    }
  }

  ##############################################################################
  # Return errors and/or warnings ##############################################
  if (length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }
  if (length(warnings)>0){
    warning_message <- sapply(warnings, function(x) paste("-", x))
    warning_message <- sapply(warning_message, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    sapply(warning_message, function(x) warning(x, call.=FALSE))
  }


  ##############################################################################
  # return potentially modified data   #########################################
  return(list("data" = data,
              "tagging.dates" = tagging.dates,
              "tag.durations" = tag.durations,
              "start.dates" = start.dates,
              "end.dates" = end.dates,
              "cutoff.dates" = cutoff.dates,
              "last.monitoring.date" = last.monitoring.date,
              "land.shape" = land.shape))
}


################################################################################
# Helper function for column argument checks ###################################
################################################################################

#' Validate Column Presence and Format
#'
#' This function checks whether a specified column name is valid, exists in the provided dataset,
#' and is correctly formatted as a single character string. It returns informative error messages if
#' the validation fails.
#'
#' @return A character string containing an error message if validation fails, or `NULL` if the column name is valid.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal


.checkColumn <- function(argument, col_label, data){
  col_name <- deparse(substitute(argument))
  if(is.null(argument)) {
    return(paste0("No ", col_label, " specified. Please provide the corresponding column name using the '", col_name, "' parameter."))
  }else if(length(argument)!=1) {
    return(paste0("The ", col_label, " name should be a single value. Please ensure you provide only one column name in the '", col_name, "' parameter."))
  }else if(!argument %in% colnames(data)){
    return(paste0("Specified ",  col_label, " ('", argument, "') not found. Please provide the correct column name using the '", col_name, "' parameter."))
  }
}


################################################################################
# Helper function to display a warning only once per session ###################
################################################################################

#' Display a warning only once per session
#'
#' This function checks if a warning has already been shown in the current session.
#' If not, it displays the warning message and ensures that the warning will not
#' be shown again during the same session.
#'
#' @param message A character string containing the warning message to be displayed.
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal

# Function to display the warning once per session
.showWarningOnce <- function(arg_name, message) {

  # create a unique identifier for the warning message
  warning_id <- paste0("warn_", arg_name)

  # check if the warning flag for this specific message is set in mobyEnv
  if (!exists(warning_id, envir = mobyEnv) || !mobyEnv[[warning_id]]) {
    # set the flag to TRUE to prevent future warnings for this message
    assign(warning_id, TRUE, envir=mobyEnv)
    # display the warning
    message <- paste(strwrap(message, width=getOption("width")), collapse="\n")
    warning(message, call.=FALSE)
  }
}


################################################################################
# Helper function for animal-specific parameters checks ########################
################################################################################

#' Validate and Process Animal-Specific Parameters
#'
#' This function validates and processes an argument intended for animal-specific parameters,
#' ensuring it conforms to the required structure and class. It handles single values,
#' unnamed vectors, and named vectors while issuing appropriate warnings or errors as needed.
#'
#' @return A list with three elements:
#' \itemize{
#'   \item `vector`: The processed argument (replicated, reordered, or unchanged as appropriate).
#'   \item `errors`: A character vector of error messages, or `NULL` if no errors are found.
#'   \item `warnings`: A character vector of warning messages, or `NULL` if no warnings are found.
#' }
#'
#' @note This function is intended for internal use within the `moby` package.
#' @keywords internal


.checkAnimalParams <- function(argument, arg_label, expected_class=c("POSIXct", "numeric"),
                               data, id_col) {

  # capture the name of the argument as a character string
  arg_name <- deparse(substitute(argument))

  # match the expected class to ensure it is either "POSIXct" or "numeric"
  expected_class <- match.arg(expected_class)

  # initialize empty lists for errors
  errors <- NULL
  processed_argument <- NULL

  # validate the class of the argument
  if (expected_class == "POSIXct" && !inherits(argument, "POSIXct")) {
    # argument must be of class POSIXct
    errors <- c(errors, paste(arg_label, "must be provided in POSIXct format."))
  } else if (expected_class == "numeric" && !(is.numeric(argument) || is.integer(argument))) {
    # argument must be numeric or integer
    errors <- c(errors, paste(arg_label, "must be provided in numeric format."))
  }

  # case 1: single value provided, replicate for all unique IDs
  if (length(argument) == 1) {
    processed_argument <- rep(argument, nlevels(data[[id_col]]))

    # case 2: multiple values provided
  } else if (length(argument) > 1) {

    # ensure the vector is named
    if (is.null(names(argument))) {
      # issue a warning for unnamed vector and keep the original order
      processed_argument <- argument
      .showWarningOnce(arg_name, paste("- The vector for",  arg_name, "is unnamed and couldn't be auto-matched to specific individuals.",
                                       "Please double-check that the supplied values match the order of the animal ID levels.",
                                       "This warning will be displayed only once during this session."))
    } else {
      # identify IDs in the data that are missing from the argument
      missing_ids <- setdiff(levels(data[[id_col]]), names(argument))
      if (length(missing_ids) > 0) {
        errors <- c(errors, paste0("The following IDs are missing in ", arg_name, ": ",
                                   paste(missing_ids, collapse = ", ")))
      } else {
        # reorder the argument to match the ID levels
        processed_argument <- argument[levels(data[[id_col]])]
      }
    }
  }

  # return the result: a list with the processed argument and errors
  return(list(vector=processed_argument, errors=errors))
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
