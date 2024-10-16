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
#' @param ... Arguments to validate. This function dynamically retrieves and checks the arguments
#'        passed to the calling function. The first argument should be a data frame containing the telemetry data.
#'
#' @return A list containing potentially modified versions of `data`, `tagging.dates`, and `tag.durations`.
#'         If validation fails, the function will stop and return an error message.
#'
#' @details
#' The `validateArguments` function checks for the following:
#' \describe{
#'   \item{data}{Ensures the first argument is a data frame.}
#'   \item{id.col}{Verifies the ID column exists and converts it to a factor if necessary.}
#'   \item{datetime.col}{Ensures the datetime column exists and is in POSIXct format.}
#'   \item{timebin.col}{Checks the time-bin column for existence and correct format.}
#'   \item{lon.col}{Ensures the longitude column exists and is numeric.}
#'   \item{lat.col}{Ensures the latitude column exists and is numeric.}
#'   \item{station.col}{Checks for the existence of the station column.}
#'   \item{split.by}{Validates the column(s) used to split the data.}
#'   \item{color.by}{Verifies the color-by column exists and contains no missing values; converts it to a factor if needed.}
#'   \item{tagging.dates}{Checks tagging dates are in POSIXct format and match the number of unique IDs if required.}
#'   \item{tag.durations}{Ensures tag durations are specified and correctly matched with unique IDs.}
#'   \item{id.groups}{Validates that ID groups are unique and present in the data.}
#'   \item{land.shape}{Checks if the provided spatial object can be converted to class `sf` if necessary.}
#'   \item{color.pal}{Ensures the color palette is either a character vector or a function.}
#'   \item{diel.lines}{Validates diel line values (0, 2, or 4).}
#'   \item{polygons}{Checks if the polygons parameter is set to a valid value ('diel', 'season', or FALSE).}
#'   \item{background.color}{Ensures the background color is a valid hex code or named color.}
#'   \item{style}{Verifies the style is either 'raster' or 'points'.}
#' }
#'
#' @note
#' This function is intended for internal use within the `moby` package. It helps streamline
#' argument validation across multiple functions, ensuring consistency and reducing redundancy.
#' @keywords internal

.validateArguments <- function() {


  ##############################################################################
  # Helper functions for common checks #########################################
  ##############################################################################

  checkColumn <- function(col_name, col_label){
    arg <- deparse(substitute(col_name))
    if(is.null(col_name)) return(paste0("No ", col_label, " specified. Please provide the corresponding column name using the '", arg, "' parameter."))
    if(length(col_name)!=1) return(paste0("The ", col_label, " name should be a single value. Please ensure you provide only one column name in the '", arg, "' parameter."))
    if(!col_name %in% colnames(data)) return(paste0("Specified ",  col_label, " ('", col_name, "') not found. Please provide the correct column name using the '", arg, "' parameter."))
  }


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
    stop(paste0(names(args)[1], " must be a data frame"), call.=FALSE)
  }

  ##############################################################################
  # validate id.col ############################################################
  if ("id.col" %in% names(args)) {
    id.col <- args$id.col
    errors <- c(errors, checkColumn(id.col, "ID column"))
    if (is.null(errors)) valid_ids <- TRUE
    # convert to factor
    if (valid_ids && !inherits(data[, id.col], "factor")){
      data[, id.col] <- as.factor(data[, id.col])
      warnings <- c(warnings, "'id.col' converted to factor")
    }
  }

  ##############################################################################
  # validate datetime.col ######################################################
  if ("datetime.col" %in% names(args)) {
    datetime.col <- args$datetime.col
    datetime_msg <- checkColumn(datetime.col, "datetime column")
    if(!is.null(datetime_msg))  errors <- c(errors, datetime_msg)
    else if (!inherits(data[, datetime.col], "POSIXct")) errors <- c(errors, "Datetimes must be provided in POSIXct format.")
  }

  ##############################################################################
  # validate timebin.col   #####################################################
  if ("timebin.col" %in% names(args)) {
    timebin.col <- args$timebin.col
    timebin_msg <- checkColumn(timebin.col, "time-bin column")
    if(!is.null(timebin_msg))  errors <- c(errors, timebin_msg)
    else if (!inherits(data[, timebin.col], "POSIXct")) errors <- c(errors, "Time-bins must be provided in POSIXct format.")
  }

  ##############################################################################
  # validate lon.col   #####################################################
  if ("lon.col" %in% names(args)) {
    lon.col <- args$lon.col
    lon_msg <- checkColumn(lon.col, "longitude column")
    if(!is.null(lon_msg))  errors <- c(errors, lon_msg)
    else if (!is.numeric(data[[lon.col]])) errors <- c(errors, paste("The column", lon.col, "must be numeric."))
  }

  ##############################################################################
  # validate lat.col   #####################################################
  if ("lat.col" %in% names(args)) {
    lat.col <- args$lat.col
    lat_msg <- checkColumn(lat.col, "latitude column")
    if(!is.null(lat_msg))  errors <- c(errors, lat_msg)
    else if (!is.numeric(data[[lat.col]])) errors <- c(errors, paste("The column", lat.col, "must be numeric."))
  }


  ##############################################################################
  # validate station.col #######################################################
  if ("station.col" %in% names(args)) {
    station.col <- args$station.col
    errors <- c(errors, checkColumn(station.col, "station column"))
  }

  ##############################################################################
  # validate spatial.col #######################################################
  if ("spatial.col" %in% names(args)) {
    spatial.col <- args$spatial.col
    errors <- c(errors, checkColumn(spatial.col, "spatial column"))
  }

  ##############################################################################
  # validate dist.col ##########################################################
  if ("dist.col" %in% names(args)) {
    dist.col <- args$dist.col
    if(!is.null(dist.col)) {
      errors <- c(errors, checkColumn(dist.col, "distance column"))
    }
  }


  ##############################################################################
  # validate split.by ##########################################################
  if ("split.by" %in% names(args)) {
    split.by <- args$split.by
    if(!is.null(split.by)) {
      for(s in 1:length(split.by)){
        errors <- c(errors, checkColumn(split.by[s], "grouping variable"))
      }
    }
  }

  ##############################################################################
  # validate subset ############################################################
  if ("subset" %in% names(args)) {
    subset <- args$subset
    if(!is.null(subset)) errors <- c(errors, checkColumn(subset, "subset column"))
  }

  ##############################################################################
  # validate variable ##########################################################
  if ("variable" %in% names(args)) {
    variable <- args$variable
    if(!is.null(variable)) errors <- c(errors, checkColumn(variable, "variable"))
  }

  ##############################################################################
  # validate variables ##########################################################
  if ("variables" %in% names(args)) {
    variables <- args$variables
    if(!is.null(variables) && calling_fun!="plotChronogram"){
      for(v in 1:length(variables)){
        errors <- c(errors, checkColumn(variables[v], "variable"))
      }
    }
  }

  ##############################################################################
  # validate color.by ##########################################################
  if ("color.by" %in% names(args)) {
    color.by <- args$color.by
    if(!is.null(color.by)){
      colorby_msg <- checkColumn(color.by, "color variable")
      if(!is.null(colorby_msg)) errors <- c(errors, colorby_msg)
      else if(any(is.na(data[, color.by]))) errors <- c(errors, "Missing values in color.by variable.")
      else if(inherits(data[,color.by], "character")) {
        data[,color.by] <- as.factor(data[,color.by])
        warnings <- c(warnings, "'color.by' variable converted to factor")}
    }
  }

  ##############################################################################
  # validate tagging.dates   ###################################################
  if ("tagging.dates" %in% names(args)) {
    tagging.dates <- args$tagging.dates
    if (is.null(tagging.dates)) errors <- c(errors, "Tagging dates not specified. Please provide the required dates via the 'tagging.dates' argument.")
    else {
      if (!inherits(tagging.dates, "POSIXct")) errors <- c(errors, "Tagging dates must be provided in POSIXct format.")
      if (length(tagging.dates) > 1 && valid_ids  && length(tagging.dates)!= nlevels(data[, id.col]))
        errors <- c(errors, "Incorrect number of tagging.dates. Must be either a single value or a vector containing a tagging date for each individual.")
      if(length(tagging.dates)==1) tagging.dates <- rep(tagging.dates, nlevels(data[,id.col]))
    }
  }

  ##############################################################################
  # validate tag durations #####################################################
  if ("tag.durations" %in% names(args)) {
    tag.durations <- args$tag.durations
    if(!is.null(tag.durations)){
      if(length(tag.durations)>1 && valid_ids && length(tag.durations)!=nlevels(data[,id.col])){
        errors <- c(errors, "Incorrect number of tag.durations. Must be either a single value or
                             a vector containing the estimated tag duration for each individual.")
      }
      if(length(tag.durations)==1) tag.durations <- rep(tag.durations, nlevels(data[,id.col]))
    }
  }

  ##############################################################################
  # validate start.dates  ######################################################
  if ("start.dates" %in% names(args)) {
    start.dates <- args$start.dates
    if(!is.null(start.dates)) {
      if(!inherits(start.dates, "POSIXct")) errors <- c(errors, "Start dates must be provided in POSIXct format.")
      if(length(start.dates)>1 && valid_ids  && length(start.dates)!=nlevels(data[, id.col]))
        errors <- c(errors, "Incorrect number of start.dates. Must be either a single value or a vector containing a start date for each individual.")
    }
  }

  ##############################################################################
  # validate end.dates  ########################################################
  if ("end.dates" %in% names(args)) {
    end.dates <- args$end.dates
    if(!is.null(end.dates)) {
      if(!inherits(end.dates, "POSIXct")) errors <- c(errors, "End dates must be provided in POSIXct format.")
      if(length(end.dates)>1 && valid_ids  && length(end.dates)!=nlevels(data[, id.col]))
        errors <- c(errors, "Incorrect number of end.dates. Must be either a single value or a vector containing an end date for each individual.")
    }
  }


  ##############################################################################
  # validate cutoff.dates  #####################################################
  if ("cutoff.dates" %in% names(args)) {
    cutoff.dates <- args$cutoff.dates
    if(!is.null(cutoff.dates)) {
      if(!inherits(cutoff.dates, "POSIXct")) errors <- c(errors, "Cutoff dates must be provided in POSIXct format.")
      if(length(cutoff.dates)>1 && valid_ids  && length(cutoff.dates)!=nlevels(data[, id.col]))
        errors <- c(errors, "Incorrect number of cutoff.dates. Must be either a single value or a vector containing a cutoff date for each individual.")
    }
  }


  ##############################################################################
  # validate id.groups #########################################################
  if ("id.groups" %in% names(args) && valid_ids) {
    id.groups <- args$id.groups
    if(!is.null(id.groups)){
      if(!inherits(id.groups, "list") || is.null(names(id.groups)) ) errors <- c(errors, "id.groups should be supplied as a named list.")
      else{
        if(any(duplicated(unlist(id.groups)))) errors <- c(errors, "Repeated ID(s) in id.groups.")
        if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {warnings <- c(warnings, "Some of the ID(s) in id.groups don't match the IDs in the data")}
        data <- data[data[,id.col] %in% unlist(id.groups),]
        if(!is.null(tagging.dates)){
          tagging.dates <- tagging.dates[match(unlist(id.groups), levels(data[,id.col]))]
        }
        if(!is.null(tag.durations)){
          tag.durations <- tag.durations[match(unlist(id.groups), levels(data[,id.col]))]
        }
      }
      data[,id.col] <- factor(data[,id.col], levels=unlist(id.groups))
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
    if(!polygons %in% c(FALSE, "diel","season")) errors <- c(errors,"Invalid 'polygons' specified. Accepted values are: 'diel', 'season' or 'FALSE' (disabled)")
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
  # Return errors and/or warnings ##############################################
  if (length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }
  if (length(warnings)>0){
    warning_message <- sapply(warnings, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    warning_message <- c("\n", paste0("- ", warning_message, collapse="\n"))
    warning(warning_message, call.=FALSE)
  }


  ##############################################################################
  # return potentially modified data   #########################################
  return(list("data"=data, "tagging.dates"=tagging.dates, "tag.durations"=tag.durations, "land.shape"=land.shape))
}

#######################################################################################################
#######################################################################################################
#######################################################################################################
