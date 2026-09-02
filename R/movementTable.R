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
#' @param verbose Logical; print a summary of the operation. Defaults to
#' \code{getOption("moby.verbose", TRUE)}.
#' @param ... Additional arguments passed to \code{\link{calculateLinearityIndex}} (and onwards to
#' \code{\link{calculateStepDistances}}), used to calculate distances between the first and last recorded
#' detections for each individual (e.g., `grid.resolution`, `mov.directions` and `cores`).
#' @seealso \code{\link{calculateROM}}, \code{\link{calculateLinearityIndex}}
#' @return A \code{mobyTable}: a TYPED data frame (counts stay integer, indices numeric, dates
#' POSIXct), one row per individual, so the result can be computed on directly. Presentation - fixed
#' precision, the display-only `mean +/- error` row, group headings - is applied by
#' \code{\link[=format.mobyTable]{format}} and \code{\link[=print.mobyTable]{print}}. Export the
#' rendered version with \code{write.csv(format(x), file, row.names = FALSE)}.
#'
#' Columns use stable snake_case names, which is what `format(decimals=)` and `format(group.by=)`
#' are keyed on; the publication headers live in `format(style = "report")`.
#' \itemize{
#'   \item your `id.col`, and a `group` factor when `id.groups` names more than one group
#'   \item `distance_km`, `rom`, `rom_max`, `linearity_index`
#'   \item the home-range columns from \code{\link{calculateUDs}}
#' }
#'
#' `rom`/`rom_max` are in m/h unless every group is fast enough to warrant km/h, in which case they
#' are scaled and the unit is recorded on the table - so the column NAME never changes with the data,
#' and `format()` names the unit in the header.
#'
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
                          verbose = getOption("moby.verbose", TRUE),
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

  # ---- header ---------------------------------------------------------------------------------
  # The only methodological choice made here is how the net displacement behind the linearity index
  # is measured: supplying a land layer swaps straight lines for shortest in-water paths, which
  # changes what the index means. Everything else (the metrics themselves, the mean +/- SE row) is
  # visible in the returned table.
  crit <- c("net displacement" = if(is.null(land.shape)) "straight-line (great-circle)" else "least-cost around land")

  # calculateROM() rejects ambiguous steps, but it runs after this header, so the same guard is
  # applied here first: a call that is going to error must not print a banner ahead of the error.
  if(any(duplicated(data[, c(id.col, timebin.col)]))){
    stop("Duplicated detection timestamps found. Please check the data.", call.=FALSE)
  }

  .mobyHeader("movementTable()", "Summarising distance, rate of movement and space use per individual",
              input = paste0(.fmtCount(nrow(data), "position"), " ", .mobyGlyph("mid"), " ",
                             .fmtCount(.nObserved(data[, id.col]), "individual")),
              criteria = crit, verbose = verbose)


  ##############################################################################
  ## Compute numeric cores (delegated) #########################################
  ##############################################################################

  # total distance and rate of movement (handles interval detection / interpolation internally)
  # (children stay silent: this function's own header already describes the run)
  metrics <- calculateROM(data, id.col=id.col, timebin.col=timebin.col, dist.col=dist.col,
                          verbose=FALSE)

  # Silencing the child would otherwise swallow its one disclosure, and this is the documented entry
  # point for these metrics: when the series was irregular, every distance and rate below was computed
  # from interpolated rather than raw steps - a methodological fact the returned table cannot reveal.
  if(isTRUE(attr(metrics, "interpolated"))){
    .mobyBlank(verbose)
    .mobyNote("Irregular time-bin widths detected ", .mobyGlyph("mid"),
              " distances interpolated to a common interval", verbose = verbose)
  }

  # movement linearity (net displacement / total distance)
  linearity <- calculateLinearityIndex(data, land.shape=land.shape, epsg.code=epsg.code, id.col=id.col,
                                       timebin.col=timebin.col, lon.col=lon.col, lat.col=lat.col,
                                       dist.col=dist.col, ..., verbose=FALSE)

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

  # individuals that will be missing from the table altogether (discard.missing drops those without a
  # single detection). Counted here because the returned table cannot reveal it: their rows are
  # simply absent. as.character() guards against factor levels collapsing to integer codes.
  n_missing <- if(discard.missing) sum(!unlist(lapply(id.groups, as.character)) %in% detected) else 0

  # subset UD results per group
  summary_table <- uds$summary_table
  uds <- lapply(id.groups, function(x) summary_table[summary_table[[id.col]] %in% x, ])



  #####################################################################
  ## Format stats #####################################################

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

  # ---- typed output -----------------------------------------------------------------------------
  # One table of NUMBERS, in group order. Precision, the mean row and the group headings are applied
  # by format()/print(); the unit chosen above travels as metadata so the header can name it without
  # the column name changing with the data.
  n_groups <- length(id.groups)
  group_ids <- lapply(id.groups, function(g) {
    g <- as.character(g)
    if (discard.missing) g[g %in% detected] else g
  })
  ids <- unlist(group_ids, use.names = FALSE)
  sub <- core[match(ids, as.character(core[[id.col]])), , drop = FALSE]

  movement_table <- data.frame(ids,
                               distance_km      = sub$total_distance_m / 1000,
                               rom              = sub$mean_rom / rom_scale,
                               rom_max          = sub$max_rom / rom_scale,
                               linearity_index  = sub$linearity_index,
                               row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(movement_table)[1] <- id.col

  # home-range columns from calculateUDs(), stacked across groups and joined once
  ud_all <- .rbindFill(uds)
  if (!is.null(ud_all) && nrow(ud_all) > 0) {
    ud_all <- .dropCols(ud_all, "group")
    movement_table <- .joinKeep(movement_table, ud_all, by = id.col, type = "left")
  }

  if (n_groups > 1 && !is.null(names(id.groups))) {
    movement_table$group <- factor(rep(names(id.groups), lengths(group_ids)),
                                   levels = names(id.groups))
  }

  # ---- outcome --------------------------------------------------------------------------------
  # No completion line: the returned table IS the summary, and restating its size would only repeat
  # what the user is about to read. The one thing the table cannot show is who is not in it.
  if(n_missing > 0){
    .mobyBlank(verbose)
    .mobyNote(.fmtCount(n_missing, "individual"), " with no detections excluded from the table",
              verbose = verbose)
  }

  # return table
  .newMobyTable(movement_table, kind = "movement", error.stat = "se", label.col = id.col,
                units = list(rom = rom_units),
                extra = list(id.groups = id.groups, processing.date = Sys.time()))
}


#######################################################################################################
#######################################################################################################
#######################################################################################################
