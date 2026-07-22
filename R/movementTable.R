#######################################################################################################
## Create movement table ###############################################################################
#######################################################################################################

#' Create movement stats table
#'
#' @description Creates a publication-ready table of per-animal movement metrics (total distance
#' travelled, rate of movement, linearity index and home-range areas), with a summary mean +/- SE
#' row. This is a formatter: the underlying numeric values are computed by
#' \code{\link{calculateROM}} (total distance and rate of movement) and
#' \code{\link{calculateLinearityIndex}} (movement directness); use those functions directly when
#' you need the raw values rather than a formatted table.
#'
#' @inheritParams as_moby
#' @param data A data frame containing binned animal detections and distances traveled,
#' as returned by \code{\link{calculateStepDistances}}.
#' @param uds Output of \code{\link{calculateUDs}}.
#' @param land.shape Optional. A projected shape file containing coastlines, used (when supplied) to
#' compute net displacements along the shortest in-water path for the linearity index.
#' @param epsg.code Coordinate reference system used to project positions (class 'CRS').
#' If not supplied, CRS is assumed to be the same as in land.shape.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param dist.col Name of the column containing distance values (in meters). Defaults to 'dist_m'.
#' @param discard.missing If true, only individuals with detections are included.
#' @param ... Additional arguments passed to \code{\link{calculateLinearityIndex}} (and onwards to
#' \code{\link{calculateStepDistances}}), used to calculate distances between the first and last recorded
#' detections for each individual (e.g., `grid.resolution`, `mov.directions` and `cores`).
#' @seealso \code{\link{calculateROM}}, \code{\link{calculateLinearityIndex}}
#' @examples
#' \donttest{
#' data(rays)
#'
#' # build per-time-bin tracks with stepwise distances
#' coas <- calculateCOAs(rays)
#' tracks <- calculateStepDistances(coas, verbose = FALSE)
#'
#' if (requireNamespace("adehabitatHR", quietly = TRUE)) {
#'   # home-range areas (a coarse estimation grid keeps this example fast)
#'   grid <- terra::rast(terra::ext(-9.05, -8.95, 38.43, 38.48),
#'                       ncol = 60, nrow = 60, crs = "EPSG:4326")
#'   terra::values(grid) <- 0
#'   grid <- terra::project(grid, "EPSG:32629")
#'   kud <- calculateUDs(coas, method = "kde", bandwidth = 500,
#'                        spatial.grid = grid)
#'
#'   # publication-ready movement metrics table (one row per animal + mean +/- SE)
#'   movementTable(tracks, uds = kud)
#' }
#' }
#' @export


movementTable <- function(data,
                          uds,
                          id.col = NULL,
                          timebin.col = NULL,
                          lon.col = NULL,
                          lat.col = NULL,
                          dist.col = "dist_m",
                          id.groups = NULL,
                          land.shape = NULL,
                          epsg.code = NULL,
                          discard.missing = TRUE,
                          ...) {

  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  land.shape <- reviewed_params$land.shape

  # validate uds
  if(!c("bandwidth") %in% names(attributes(uds))) stop("The supplied uds do not seem to be in the right format. Please use the output of the 'calculateUDs' function.", call. = FALSE)


  ##############################################################################
  ## Compute numeric cores (delegated) #########################################
  ##############################################################################

  # total distance and rate of movement (handles interval detection / interpolation internally)
  metrics <- calculateROM(data, id.col=id.col, timebin.col=timebin.col, dist.col=dist.col)

  # movement linearity (net displacement / total distance)
  linearity <- calculateLinearityIndex(data, land.shape=land.shape, epsg.code=epsg.code, id.col=id.col,
                                       timebin.col=timebin.col, lon.col=lon.col, lat.col=lat.col,
                                       dist.col=dist.col, ...)

  # assemble a single per-individual numeric core; use the (possibly interpolated) total distance
  # from 'metrics' as the linearity denominator so the displayed distance and LI stay consistent
  core <- metrics
  core$net_distance_m <- linearity$net_distance_m[match(core[[id.col]], linearity[[id.col]])]
  core$linearity_index <- core$net_distance_m / core$total_distance_m
  core$linearity_index[!is.finite(core$linearity_index)] <- NA_real_

  # individuals with at least one detection
  detected <- as.character(unique(data[, id.col]))

  # define single id.group if needed
  if(is.null(id.groups)){
    id.groups <- list(levels(data[,id.col]))
  }

  # subset UD results per group
  summary_table <- uds$summary_table
  uds <- lapply(id.groups, function(x) summary_table[summary_table[[id.col]] %in% x, ])



  #####################################################################
  ## Format stats #####################################################

  movement_table <- list()

  # decide the rate-of-movement display unit ONCE for the whole table. Deciding it per id.group (as
  # before) gave group-specific column names ("ROM (m/h)" vs "ROM (km/h)") that broke the final
  # rbind() whenever groups differed in speed. Switch to km/h only when EVERY group is fast, so a slow
  # group is never shown in km/h at a precision that would collapse it to "0.0". For a single group
  # this reduces to the original per-dataset rule (m/h is the lossless base unit).
  group_fast <- vapply(id.groups, function(g){
    gi <- as.character(core[[id.col]]) %in% as.character(g)
    isTRUE(mean(core$mean_rom[gi], na.rm=TRUE) > 1000) && isTRUE(mean(core$max_rom[gi], na.rm=TRUE) > 1000)
  }, logical(1))
  if(length(group_fast) > 0 && all(group_fast)){
    rom_units <- "km/h"; rom_scale <- 1000
  }else{
    rom_units <- "m/h"; rom_scale <- 1
  }

  for(i in seq_along(id.groups)){

    # select this group's individuals (optionally dropping those without detections)
    group_ids <- as.character(id.groups[[i]])
    if(discard.missing) group_ids <- group_ids[group_ids %in% detected]
    sub <- core[match(group_ids, core[[id.col]]), , drop=FALSE]

    # total distance traveled (km)
    total_distance <- sprintf("%.1f", sub$total_distance_m/1000)

    # rate of movement (hourly distance); unit/scale decided once above so every group shares columns
    mean_rom <- sprintf("%.1f", sub$mean_rom/rom_scale)
    max_rom <- sprintf("%.1f", sub$max_rom/rom_scale)

    # linearity index
    li_index <- sprintf("%.2f", sub$linearity_index)
    li_index[li_index=="NaN" | li_index=="NA"] <- "NA"

    # overall movement stats
    movement_stats <- data.frame("ID"=group_ids, "Distance (km)"=total_distance,
                                 "ROM"=mean_rom, "Max ROM"=max_rom,
                                 "LI"=li_index, check.names=FALSE, row.names=NULL)
    colnames(movement_stats)[1] <- id.col
    colnames(movement_stats)[3] <- paste0(colnames(movement_stats)[3], " (", rom_units,")")
    colnames(movement_stats)[4] <- paste0(colnames(movement_stats)[4], " (", rom_units,")")
    movement_stats <- .joinKeep(movement_stats, uds[[i]], by=id.col, type="left")
    if("group" %in% colnames(movement_stats)) movement_stats <- .dropCols(movement_stats, "group")

    # calculate means ± se and format missing values
    values_chr <- as.matrix(movement_stats[,-1, drop=FALSE])
    values <- suppressWarnings(matrix(as.numeric(values_chr), nrow=nrow(values_chr),
                                      ncol=ncol(values_chr), dimnames=dimnames(values_chr)))
    # per-column maximum number of decimal places (robust to single columns / rows)
    decimal_digits <- apply(values, 2, function(col){
      dp <- .decimalPlaces(col)
      if(all(is.na(dp))) 0 else max(dp, na.rm=TRUE)
    })
    movement_stats[[id.col]] <- as.character(movement_stats[[id.col]])
    movement_stats[nrow(movement_stats)+1,] <- NA
    movement_stats[[id.col]][nrow(movement_stats)] <- "mean"
    movement_stats[nrow(movement_stats), -1] <- sprintf(paste0("%.", decimal_digits, "f"), colMeans(values, na.rm=TRUE))
    errors <- sprintf(paste0("%.", decimal_digits, "f"), unlist(apply(values, 2, .stdError)))
    movement_stats[nrow(movement_stats), -1] <- paste(movement_stats[nrow(movement_stats), -1], "\u00b1", errors)
    movement_stats[is.na(movement_stats)] <- "-"
    movement_stats[movement_stats=="NA"] <- "-"

    # add title
    if(length(id.groups)>1){
      table_title <- movement_stats[0,]
      table_title[1,] <- ""
      table_title[[id.col]] <- names(id.groups)[i]
      movement_stats <- rbind(table_title, movement_stats)
    }

    # save table
    movement_table[[i]] <- movement_stats
  }


  # aggregate group tables
  movement_table <- do.call("rbind", movement_table)

  # return table
  return(movement_table)
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
