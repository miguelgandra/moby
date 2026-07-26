#######################################################################################################
## Pipeline overview (documentation-only topic) ########################################################
#######################################################################################################

#' The moby pipeline: from a raw export to an analysis-ready dataset
#'
#' @description `moby`'s preparation functions look numerous, but they fill only **five roles**. Learn the
#' roles and the running order stops being something to memorise. Two things are true throughout, and
#' both simplify the picture:
#'
#' \itemize{
#'   \item Every function also accepts an **ordinary data frame**. A \code{\link{mobyData}} is a
#'     convenience that saves you repeating arguments - never a gate you must pass through.
#'   \item Only the \emph{enrich} and \emph{clean} roles change your data. The audit role cannot.
#' }
#'
#' @section The five roles:
#' \describe{
#'   \item{\strong{1. Read} - \code{\link{importDetections}}, \code{\link{importTags}},
#'     \code{\link{importDeployments}}}{Turn a vendor export into moby's canonical columns: renames
#'     columns, parses date-times, coerces coordinates. Each returns a plain data frame. Optional - skip
#'     them if your table already uses `ID`/`datetime`/`station`/`lon`/`lat` with a POSIXct datetime.
#'     See \code{\link{moby_import_schema}} for every field.}
#'   \item{\strong{2. Audit} - \code{\link{checkDeployments}}}{Validates the receiver log (coverage gaps,
#'     overlapping deployments, duplicate stations, invalid dates, implausible or on-land coordinates,
#'     short monitoring effort) and cross-checks detections against it. Returns a report and
#'     \strong{modifies nothing}. Recommended.}
#'   \item{\strong{3. Enrich} - \code{\link{assignAnimalIDs}}, \code{\link{matchDeployments}}}{Fold a side
#'     table into the detections. `assignAnimalIDs()` resolves transmitter codes to animals and attaches
#'     per-animal tagging dates and nominal delays; `matchDeployments()` matches detections to deployment
#'     windows and back-fills station names and coordinates. They key on different columns (transmitter vs
#'     receiver) and are \strong{order-independent}.}
#'   \item{\strong{4. Clean & bin} - \code{\link{filterDetections}}, \code{\link{correctPositions}},
#'     \code{\link{getTimeBins}}}{Where the data actually becomes analysis-ready: removing spurious
#'     detections, relocating positions that fall on land, and assigning regular time bins.}
#'   \item{\strong{5. Declare} - \code{\link{as_moby}}}{Records which column plays which role and attaches
#'     study metadata (`id.groups`, `epsg.code`, `land.shape`). Not really a stage - it is idempotent, so
#'     call it whenever you want to add metadata without losing what is already attached.}
#' }
#'
#' @section A typical run:
#' ```r
#' # 1 - read
#' det  <- importDetections(det_file, source = "vue")
#' tags <- importTags(tag_file, source = "vue")
#'
#' # 2 - enrich: animal identity, plus tagging dates and nominal delays
#' det <- assignAnimalIDs(det, tags)
#' det <- det[!is.na(det$ID), ]          # drop transmitters absent from 'tags'
#' det$ID <- droplevels(det$ID)
#'
#' # 3 - clean: note this returns a REPORT, not a data frame
#' qc  <- filterDetections(det, max.speed = 2, speed.unit = "m/s")
#' qc                                    # what was removed, and why
#' det <- qc$data
#'
#' # 4 - bin, and declare anything only you know
#' det$timebin <- getTimeBins(det$datetime, interval = "30 mins")
#' dataset <- as_moby(det, id.groups = groups, epsg.code = 32629)
#' ```
#'
#' Add the receiver-log branch when you need it - audit, then match:
#' ```r
#' deployments <- importDeployments(log_file, source = "vue")
#' checkDeployments(deployments, detections = det)     # report only; fix what it flags
#' det <- matchDeployments(det, deployments)
#' ```
#'
#' Already hold a tidy data frame? No moby object is required at any point:
#' ```r
#' qc <- filterDetections(df, tagging.dates = tagging_dates, max.speed = 2)
#' ```
#'
#' @section What happens if you skip a step:
#' \tabular{lll}{
#'   \strong{function} \tab \strong{status} \tab \strong{skip it and...} \cr
#'   `importDetections()`  \tab optional    \tab nothing breaks; supply canonical columns yourself \cr
#'   `importTags()`        \tab optional    \tab nothing breaks; you must supply tagging dates another way \cr
#'   `importDeployments()` \tab optional    \tab nothing breaks; only needed if you use the receiver log \cr
#'   `checkDeployments()`  \tab recommended \tab nothing breaks; you may match against a log with real errors \cr
#'   `matchDeployments()`  \tab optional    \tab coordinates stay as exported, unvalidated and ungapfilled \cr
#'   `assignAnimalIDs()`   \tab \strong{essential in practice} \tab `filterDetections()` \strong{errors} without tagging dates; the min_lag filter is unavailable \cr
#'   `filterDetections()`  \tab \strong{essential} \tab false detections, pre-tagging records and impossible speeds stay in \cr
#'   `correctPositions()`  \tab optional    \tab nothing breaks unless coordinates sit on land \cr
#'   `getTimeBins()`       \tab recommended \tab no per-bin analyses (COAs, wide tables, chronograms) \cr
#'   `as_moby()`           \tab recommended \tab no `id.groups`, projected CRS or land-aware distances \cr
#' }
#'
#' @section Ordering rules (there are only two):
#' \enumerate{
#'   \item \strong{Audit before you delete.} \code{\link{checkDeployments}} is read-only and idempotent, so
#'     it has no single correct slot - run it early to find problems, and again after corrections to confirm
#'     they landed. The one position that matters: run it \emph{before}
#'     `matchDeployments(drop.unmatched = TRUE)`, which removes the records that four of its five
#'     detection-side checks rely on. Under the default (`FALSE`) nothing is lost and order is free.
#'   \item \strong{Enrich before you filter.} \code{\link{filterDetections}} needs tagging dates, which
#'     \code{\link{assignAnimalIDs}} supplies. The two enrich steps themselves may run in either order.
#' }
#'
#' @section Three easy traps:
#' \itemize{
#'   \item \code{\link{filterDetections}} and \code{\link{correctPositions}} return **lists**, not data
#'     frames - take `$data`. Writing `det <- filterDetections(det)` leaves you holding a `mobyFilter`
#'     object and the next call fails.
#'   \item \code{\link{getTimeBins}} returns a **vector**; assign it to a column yourself.
#'   \item \code{\link{assignAnimalIDs}} leaves `NA` IDs (with a warning) for transmitters that are not in
#'     your tag table - typically other projects' tags. Remove them, or they travel into the analysis.
#' }
#'
#' @seealso \code{\link{moby_import_schema}} for the canonical column names; \code{\link{as_moby}} for the
#' data model; the package vignettes for worked examples.
#' @name moby_pipeline
NULL


#######################################################################################################
#######################################################################################################
#######################################################################################################
