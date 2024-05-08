#######################################################################################################
## Season shading #####################################################################################
#######################################################################################################

#' Return season boundaries for background shading.
#'
#' @description Returns a table containing season boundaries/limits, used to add background polygons
#' in \code{\link{plotAbacus}}.
#'
#' @param start.time A POSIXct object containing the earliest date used in the plot region.
#' @param end.time A POSIXct object containing the latest date used in the plot region.
#' @param interval Time-bins interval (in minutes).
#' @param color.pal Vector of 4 colors, one for each season
#' (in the following order: winter, spring, summer and autumn).
#' @param hemisphere Earth hemisphere for which to calculate seasons.
#' @seealso \code{\link{getSeason}}
#' @export

shadeSeasons <- function(start.time, end.time, interval, color.pal=c("white", "grey95", "grey80", "grey87"),
                         hemisphere="Northern") {

  color_pal <- data.frame("season"=c("winter", "spring", "summer", "autumn"),
                          "color"=color.pal)
  complete_seqs <- seq.POSIXt(from=start.time, to=end.time, by=interval*60)
  complete_seqs <- data.frame("timebin"=complete_seqs, "season"=getSeason(complete_seqs, hemisphere))
  complete_seqs$timebin <- as.character(complete_seqs$timebin)
  consec_seasons <- rle(as.character(complete_seqs$season))
  complete_seqs$season  <- paste0(complete_seqs$season, "_", rep(1:length(consec_seasons$lengths), consec_seasons$lengths))
  seasons_start <- aggregate(complete_seqs$timebin, by=list(complete_seqs$season), min, simplify=F)
  seasons_end <- aggregate(complete_seqs$timebin, by=list(complete_seqs$season), max, simplify=F)
  seasons_table <-  plyr::join(seasons_start, seasons_end, by="Group.1", type="left")
  colnames(seasons_table) <- c("season", "start", "end")
  seasons_table$season <- sub("\\_.*", "", seasons_table$season)
  date_format <- ifelse(any(nchar(seasons_table$start)>10), "%Y-%m-%d %H:%M:%S", "%Y-%m-%d")
  seasons_table$start <- as.POSIXct(unlist(seasons_table$start), date_format, tz="UTC")
  seasons_table$end <- as.POSIXct(unlist(seasons_table$end), date_format, tz="UTC")
  seasons_table <- seasons_table[order(seasons_table$start),]
  index <- 1:nrow(seasons_table)
  seasons_table$start[index!=1] <-  seasons_table$start[index!=1]  - (interval*60)/2
  seasons_table$end[index!=nrow(seasons_table)] <-  seasons_table$end[index!=nrow(seasons_table)] + (interval*60)/2
  season_table <- plyr::join(seasons_table, color_pal, by="season", type="left")
  return(season_table)
 }

#######################################################################################################
#######################################################################################################
#######################################################################################################
