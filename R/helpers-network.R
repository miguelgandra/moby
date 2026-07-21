##################################################################################################
## Internal helpers: association-network plots ###################################################
## Shared by plotAssociations and plotAssociationMatrix so the null-model validation, the id x id  ##
## matrix reshaping and the colour mixing live in ONE place.                                       ##
##################################################################################################


#' Validate a randomizeAssociations() result object
#'
#' @description Returns a character vector of error messages (empty if valid). Both association plots
#' call this so the "not the right object" check - and the correct function name in the message - is
#' defined once (the old inline checks referenced the non-existent `randomizeOverlaps`).
#' @keywords internal
#' @noRd
.validateRandomResults <- function(x) {
  errors <- character(0)
  if (!inherits(x, "list"))
    errors <- c(errors, "'random.results' must be a list, as returned by randomizeAssociations().")
  else if (!"metric" %in% names(attributes(x)))
    errors <- c(errors, "'random.results' not recognised. Please supply the output of randomizeAssociations().")
  errors
}

#' Reshape long pairwise results into an id x id matrix
#'
#' @description Pivots a long pairwise table (columns `id1`, `id2` and `value.var`) into a wide
#' id x id data frame with `id1` as row names. `drop = FALSE` in the dcast keeps every factor level so
#' the matrix stays square/padded (needed by the network plot); `discard.missing` optionally drops
#' rows/columns that are entirely NA, always with `drop = FALSE` so a single surviving row/column does
#' not silently collapse to a vector (which previously erased the ID labels).
#' @param pairwise Long data frame with `id1`, `id2` and the value column.
#' @param value.var Name of the value column to spread (e.g. "association" or "significance").
#' @param discard.missing Logical; drop all-NA rows and columns. Defaults to TRUE.
#' @return A data frame (id x id) with row names = id1 and column names = id2.
#' @keywords internal
#' @noRd
.associationMatrix <- function(pairwise, value.var, discard.missing = TRUE) {
  m <- .castWide(pairwise, "id1", "id2", value.var, fun.aggregate = function(z) z[1L], fill = NA)
  rownames(m) <- as.character(m$id1)
  m <- m[, -1, drop = FALSE]
  if (discard.missing) {
    m <- m[rowSums(!is.na(m)) > 0, colSums(!is.na(m)) > 0, drop = FALSE]
  }
  m
}

#' Average a set of colours in RGB space
#' @keywords internal
#' @noRd
.mixColors <- function(colors) {
  avg <- rowMeans(vapply(colors, grDevices::col2rgb, numeric(3)))
  grDevices::rgb(avg[1], avg[2], avg[3], maxColorValue = 255)
}

##################################################################################################
##################################################################################################
##################################################################################################
