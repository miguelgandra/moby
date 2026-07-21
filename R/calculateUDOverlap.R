#######################################################################################################
## Pairwise overlap between utilization distributions #################################################
#######################################################################################################

#' Pairwise overlap between utilization distributions
#'
#' @description Computes pairwise overlap between the individual utilization distributions (UDs)
#' produced by \code{\link{calculateUDs}}, returning a tidy table of dyadic overlap values (one row
#' per pair of animals). It is the home-range counterpart to \code{\link{calculateAssociations}}
#' (which measures spatiotemporal co-occurrence): here the question is how much two animals'
#' \emph{space use} overlaps, irrespective of timing.
#'
#' The overlap computation is delegated to the estimator that produced the UDs, so the result is
#' method-consistent: for autocorrelated KDE (AKDE, the \code{\link{calculateUDs}} default) it uses
#' \code{\link[ctmm]{overlap}}, which returns the Bhattacharyya coefficient \strong{with confidence
#' intervals}; for classic kernel density (\code{method = "kde"}) it uses
#' \code{\link[adehabitatHR]{kerneloverlaphr}}, which offers a wider set of indices but no CIs.
#'
#' @details Available `index` values depend on how the UDs were estimated:
#' \itemize{
#'   \item \strong{AKDE} (`ctmm`): only \code{"BA"} (Bhattacharyya coefficient), reported with a
#'   lower/upper confidence interval at the requested `level`. This autocorrelation-aware estimate is
#'   the recommended choice for tracking data.
#'   \item \strong{KDE} (`adehabitatHR`): \code{"BA"} (Bhattacharyya), \code{"UDOI"} (utilization
#'   distribution overlap index), \code{"HR"} (home-range area overlap), \code{"PHR"} (probability of
#'   finding animal j in animal i's range), \code{"VI"} (volume of intersection) or \code{"HD"}
#'   (Hellinger distance). See Fieberg & Kochanny (2005) for definitions.
#' }
#' \code{"BA"} is the default because it is the one index available under \emph{both} methods, so the
#' default works whether the UDs are AKDE or KDE. All indices are bounded in \code{[0, 1]} except
#' \code{UDOI}, which can exceed 1 when two ranges overlap and are both non-uniform. Overlap is
#' computed within each unit (e.g. species / \code{id.groups} block) that `ud` was estimated over,
#' never across units (whose UDs live on separate grids); pairs are formed among individuals of the
#' same unit only.
#'
#' \strong{Symmetric vs directional indices.} Most indices (\code{"BA"}, \code{"UDOI"}, \code{"VI"},
#' \code{"HD"}) are symmetric — the overlap of A with B equals that of B with A — and are reported once
#' per unordered pair. \code{"HR"} and \code{"PHR"} are \emph{directional}: \code{HR[i, j]} is the
#' proportion of animal \code{i}'s range covered by animal \code{j}, and \code{PHR[i, j]} the
#' probability of finding animal \code{j} inside animal \code{i}'s range, so \code{[i, j] != [j, i]}.
#' For these two indices the result therefore has one row per \emph{ordered} pair (both \code{id1->id2}
#' and \code{id2->id1}), the value in each row being the index evaluated for \code{(id1, id2)}.
#'
#' @param ud The output of \code{\link{calculateUDs}} (a list with a `$ud` element). The estimation
#' method (AKDE vs KDE) is detected automatically.
#' @param index Overlap index to compute (case-insensitive). One of \code{"BA"} (default), and for
#' KDE-estimated UDs additionally \code{"UDOI"}, \code{"HR"}, \code{"PHR"}, \code{"VI"}, \code{"HD"}.
#' @param contour KDE only: the home-range isopleth percentage the overlap is restricted to (a value
#' in (0, 100]). Defaults to 95. Ignored for AKDE (which integrates the full UD).
#' @param level AKDE only: the confidence level for the reported overlap CIs. Defaults to 0.95.
#' Ignored for KDE (which provides no CIs).
#' @param id.groups Optional named list of ID groups. When supplied, each pair is annotated with its
#' `group1`/`group2` membership and a `pair_type` of `"within"` or `"between"`.
#' @param verbose Logical; print progress messages. Defaults to TRUE.
#'
#' @return A tidy data frame with one row per pair of individuals (self-pairs excluded) — unordered
#' pairs for symmetric indices, ordered pairs for the directional indices \code{HR}/\code{PHR}:
#' \item{id1, id2}{the two animals in the pair.}
#' \item{<index>}{the overlap value (column named after `index`, e.g. `BA` or `UDOI`). For `HR`/`PHR`
#' this is the directional value for `(id1, id2)`.}
#' \item{<index>_lower, <index>_upper}{(AKDE only) the confidence-interval bounds.}
#' \item{group1, group2, pair_type}{(when `id.groups` is supplied or `ud` was grouped) the group
#' membership of each animal and whether the pair is within or between groups.}
#' When a single unit produces pairs, the estimator's own overlap matrix is attached as attribute
#' `"matrix"` (symmetric except for the directional `HR`/`PHR` indices, which stay asymmetric);
#' `"method"`, `"index"`, and `"contour"`/`"level"` record how it was computed.
#'
#' @references
#' Fieberg, J. & Kochanny, C. O. (2005). Quantifying home-range overlap: the importance of the
#' utilization distribution. Journal of Wildlife Management, 69(4), 1346-1359.
#'
#' Winner, K., Noonan, M. J., Fleming, C. H., et al. (2018). Statistical inference for home range
#' overlap. Methods in Ecology and Evolution, 9(7), 1679-1691.
#'
#' @seealso \code{\link{calculateUDs}}, \code{\link{calculateAssociations}}, \code{\link{plotMaps}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("adehabitatHR", quietly = TRUE)) {
#'   data(rays)
#'   # estimate UDs for a few animals (kde keeps the example fast), then measure pairwise overlap
#'   sub <- rays[rays$ID %in% head(levels(factor(rays$ID)), 3), ]
#'   ud  <- calculateUDs(sub, method = "kde", bandwidth = 500, verbose = FALSE)
#'   calculateUDOverlap(ud, index = "UDOI")
#' }
#' }
#' @export

calculateUDOverlap <- function(ud,
                               index = "BA",
                               contour = 95,
                               level = 0.95,
                               id.groups = NULL,
                               verbose = TRUE) {

  ##############################################################################
  ## Validate + detect method ##################################################
  ##############################################################################

  if (is.null(ud) || is.null(ud$ud))
    stop("'ud' must be the output of calculateUDs() (a list with a '$ud' element).", call. = FALSE)
  if (!is.numeric(contour) || length(contour) != 1 || contour <= 0 || contour > 100)
    stop("'contour' must be a single isopleth percentage in (0, 100].", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1)
    stop("'level' must be a single confidence level in (0, 1).", call. = FALSE)

  uds <- ud$ud
  method <- attr(ud, "method")
  if (is.null(method) || !method %in% c("akde", "kde")) method <- .detectUDMethod(uds)
  if (is.na(method)) stop("Could not determine the UD estimation method from 'kud'.", call. = FALSE)

  kde_idx <- c("BA", "UDOI", "HR", "PHR", "VI", "HD")
  index <- toupper(index)
  if (method == "akde" && index != "BA")
    stop("For AKDE utilization distributions, only index = 'BA' (Bhattacharyya, with CIs) is ",
         "available. Use index = 'BA', or estimate UDs with method = 'kde' for ",
         paste(setdiff(kde_idx, "BA"), collapse = "/"), ".", call. = FALSE)
  if (method == "kde" && !index %in% kde_idx)
    stop("Invalid 'index' for a KDE overlap. Choose from: ", paste(kde_idx, collapse = ", "), ".", call. = FALSE)

  pkg <- if (method == "akde") "ctmm" else "adehabitatHR"
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Computing overlap for method = '", method, "' requires the '", pkg, "' package.", call. = FALSE)

  ##############################################################################
  ## Compute overlap (per group), assemble tidy pairs ##########################
  ##############################################################################

  groups <- .udGroups(uds, method)   # named list of per-group UD collections ("all" if ungrouped)
  # HR and PHR are DIRECTIONAL (asymmetric) indices: kerneloverlaphr()[i, j] differs from [j, i], so
  # each ordered pair is reported on its own row. All other indices (BA/UDOI/VI/HD, and AKDE's BA) are
  # symmetric and reported once per unordered pair.
  directed <- method == "kde" && index %in% c("HR", "PHR")
  out <- list(); mats <- list()
  for (g in names(groups)) {
    coll <- groups[[g]]
    ids <- names(coll)
    if (length(ids) < 2) {
      if (verbose) message("- Group '", g, "' has < 2 individuals; no pairs computed.")
      next
    }
    if (method == "kde") {
      m  <- adehabitatHR::kerneloverlaphr(coll, method = index, percent = contour, conditional = TRUE)
      df <- .overlapPairs(m, index, ci = NULL, directed = directed)
    } else {
      ov  <- ctmm::overlap(coll, level = level)
      ci  <- ov$CI                                  # [n, n, 3]; 3rd dim = low / est / high
      n   <- length(ids)
      pull <- function(k) matrix(ci[, , k], n, n, dimnames = list(ids, ids))  # shape-agnostic slice
      m   <- pull(2L)
      df  <- .overlapPairs(m, index, ci = list(lower = pull(1L), upper = pull(3L)))
    }
    if (length(groups) > 1 && !identical(names(groups), "all")) df$group <- g
    out[[g]] <- df
    mats[[g]] <- m   # the estimator's own overlap matrix (asymmetric for HR/PHR), for the attribute
  }

  result <- if (length(out)) do.call(rbind, out)
            else .overlapPairs(matrix(numeric(0), 0, 0), index,
                               ci = if (method == "akde") list() else NULL, directed = directed)
  rownames(result) <- NULL

  ##############################################################################
  ## Annotate + attributes #####################################################
  ##############################################################################

  if (!is.null(id.groups) && nrow(result) > 0) {
    lut <- stats::setNames(rep(names(id.groups), lengths(id.groups)), unlist(id.groups))
    result$group1 <- unname(lut[result$id1])
    result$group2 <- unname(lut[result$id2])
    result$pair_type <- ifelse(is.na(result$group1) | is.na(result$group2), NA_character_,
                               ifelse(result$group1 == result$group2, "within", "between"))
  }

  attr(result, "method") <- method
  attr(result, "index") <- index
  if (method == "kde") attr(result, "contour") <- contour else attr(result, "level") <- level
  # attach the estimator's own overlap matrix when a single unit produced pairs (symmetric except HR/PHR)
  if (length(mats) == 1 && nrow(result) > 0) attr(result, "matrix") <- mats[[1]]
  attr(result, "processing.date") <- Sys.time()

  if (verbose) message("- Computed ", index, " overlap for ", nrow(result), " pair(s) (", method, ").")
  result
}


################################################################################
# Internal helpers #############################################################
################################################################################

#' Detect the UD estimation method from a calculateUDs() $ud element
#' @keywords internal
#' @noRd
.detectUDMethod <- function(ud) {
  if (inherits(ud, "estUDm")) return("kde")
  if (is.list(ud) && length(ud) > 0) {
    f <- ud[[1]]
    if (inherits(f, "estUDm")) return("kde")
    if (inherits(f, "UD")) return("akde")
    if (is.list(f) && length(f) > 0) {
      if (inherits(f[[1]], "estUDm")) return("kde")
      if (inherits(f[[1]], "UD")) return("akde")
    }
  }
  NA_character_
}

#' Normalise a $ud element into a named list of per-group UD collections
#' @keywords internal
#' @noRd
.udGroups <- function(ud, method) {
  if (method == "kde") {
    if (inherits(ud, "estUDm")) return(list(all = ud))
    if (is.list(ud) && all(vapply(ud, inherits, logical(1), "estUDm"))) return(ud)
    stop("Unrecognised KDE UD structure in 'kud$ud'.", call. = FALSE)
  }
  # akde: calculateUDs() returns a FLAT named list of ctmm UD objects keyed "unit::id" (each unit
  # is estimated on its own grid). Split by the unit prefix so overlap is only ever computed among
  # grid-aligned individuals of the same unit; a plain (no-"::") list is treated as a single group.
  if (all(vapply(ud, inherits, logical(1), "UD"))) {
    nms <- names(ud)
    if (is.null(nms) || !any(grepl("::", nms, fixed = TRUE)))
      return(list(all = stats::setNames(ud, nms)))
    # Keys are paste0(unit, "::", id); the id is the trailing segment, so split on the LAST "::"
    # (keeps unit names that themselves contain "::" intact, e.g. hierarchical id.groups labels).
    unit <- sub("::[^:]*$", "", nms)
    bare <- sub("^.*::", "", nms)
    ix   <- split(seq_along(ud), factor(unit, levels = unique(unit)))
    return(lapply(ix, function(i) {
      if (anyDuplicated(bare[i]))
        stop("Individual IDs collide within a unit after splitting UD keys; individual IDs must not ",
             "contain the reserved '::' delimiter.", call. = FALSE)
      stats::setNames(ud[i], bare[i])
    }))
  }
  # already nested as a list of per-unit UD lists
  if (is.list(ud) && all(vapply(ud, function(x) is.list(x) && length(x) > 0 &&
                                all(vapply(x, inherits, logical(1), "UD")), logical(1)))) return(ud)
  stop("Unrecognised AKDE UD structure in 'ud$ud'.", call. = FALSE)
}

#' Turn an n x n overlap matrix into tidy pairs.
#'
#' Symmetric indices give one row per unordered pair (i < j). Directional indices (`directed = TRUE`,
#' i.e. KDE HR/PHR) give one row per ordered pair (both i->j and j->i), since `m[i, j] != m[j, i]`.
#' @keywords internal
#' @noRd
.overlapPairs <- function(m, index, ci = NULL, directed = FALSE) {
  ids <- rownames(m)
  if (is.null(ids) || length(ids) < 2) {
    df <- data.frame(id1 = character(0), id2 = character(0), stringsAsFactors = FALSE)
    df[[index]] <- numeric(0)
    if (!is.null(ci)) { df[[paste0(index, "_lower")]] <- numeric(0); df[[paste0(index, "_upper")]] <- numeric(0) }
    return(df)
  }
  if (directed) {
    grid <- expand.grid(i = ids, j = ids, stringsAsFactors = FALSE)
    grid <- grid[grid$i != grid$j, , drop = FALSE]
    cmb  <- cbind(grid$i, grid$j)
  } else {
    cmb <- t(utils::combn(ids, 2))
  }
  df <- data.frame(id1 = cmb[, 1], id2 = cmb[, 2], stringsAsFactors = FALSE)
  df[[index]] <- as.numeric(m[cmb])
  if (!is.null(ci)) {
    df[[paste0(index, "_lower")]] <- as.numeric(ci$lower[cmb])
    df[[paste0(index, "_upper")]] <- as.numeric(ci$upper[cmb])
  }
  df
}
