#######################################################################################################
## Plot receiver deployment / operating-period timelines #############################################
#######################################################################################################

#' Plot receiver operating periods (deployment timeline)
#'
#' @description Draws a timeline of when each receiver (or station/site) was deployed and actively
#' listening over the study period - one row per unit, with a bar spanning each `[deploy, recover]`
#' interval and broken bars where a real coverage gap occurs. It is the visual companion to
#' \code{\link{checkDeployments}} (which reports metadata problems) and shares the visual language of
#' \code{\link{plotAbacus}} (calendar axis, labelled year band, group colouring). Because a gap in a
#' detection record can mean either "animal absent" or "receiver not listening", this figure supplies
#' the sampling-effort context needed to read the detection plots correctly.
#'
#' @details Rows are the unique values of `row.by` (a receiver, or a coarser unit such as a station or
#' site); a unit's deployments are merged into operating-period bars, with consecutive deployments
#' whose gap is `<= merge.gaps` days fused into one bar and longer gaps left as breaks. Rows can be
#' coloured and legended by a coarser grouping (`group.by`, e.g. array or region). An optional `events`
#' overlay draws neutral markers at points in time on the relevant rows (e.g. taggings, servicing
#' visits, downloads). The layout is device-stable (inch-anchored band, labels and legend), matching
#' the rest of the plotting family.
#'
#' @param deployments A receiver-deployment / station log (e.g. from \code{\link{importDeployments}}),
#' with at least the `row.by` and `deployment.deploy.col` columns (and usually `deployment.recover.col`).
#' @param row.by Name of the column defining the rows (the operating unit). Defaults to `"receiver"`;
#' set to `"station"` (or any column) to aggregate operating periods to a coarser unit.
#' @param group.by Optional column giving a coarser grouping used to colour the rows and their labels
#' (e.g. `"array"`, `"region"`). If NULL, all rows share `single.color`.
#' @param events Optional data frame of point events to overlay, with a column named as `row.by` (the
#' row key) and a date-time column (`time`, or the first POSIXct column). Study-neutral: taggings,
#' servicing, downloads, etc.
#' @template deploymentDateArgs
#' @param end Optional study-end date-time used to close still-active deployments. If NULL, the latest
#' date in the log is used.
#' @param merge.gaps Numeric. Consecutive deployments of a unit separated by `<=` this many days are
#' drawn as a single bar; longer gaps appear as breaks. Defaults to 1.
#' @param sort.by Row ordering: `"group"` (by `group.by`, then first deployment; default), `"start"`
#' (first deployment) or `"name"` (alphabetical).
#' @param color.pal Colours for the groups. If NULL, a colourblind-safe Okabe-Ito palette.
#' @param single.color Colour used when `group.by` is NULL. Defaults to "grey20".
#' @param bar.height Bar thickness as a fraction of the row spacing. Defaults to 0.45.
#' @param dividers Logical; draw a subtle separator line between groups (only when `group.by` is set
#' and `sort.by = "group"`, so groups are contiguous). Defaults to TRUE.
#' @param top.band A \code{\link[base]{strftime}} format for the labelled band above the panel (e.g.
#' "%Y" for years). FALSE to omit. Defaults to "%Y".
#' @param date.interval Controls the x-axis date labels. `"auto"` (default) chooses calendar-aware
#' breaks automatically. A positive integer switches to manual mode - a label every n-th
#' `date.format` unit (e.g. `date.format = "%b"`, `date.interval = 3` labels every third month) -
#' useful for very long or very short spans.
#' @param date.start Integer phase/offset for the first manual label (used with a numeric
#' `date.interval`). Defaults to 1.
#' @param date.format Optional \code{\link[base]{strftime}} format for the x-axis labels. If NULL,
#' calendar-aware labels are chosen automatically (and "%b" in manual mode).
#' @param background.color Panel background colour. Defaults to "grey96".
#' @param events.pch,events.color,events.cex,events.label Marker symbol, colour, size and legend label
#' for the `events` overlay. Defaults `pch = 8`, "grey15", 1, "event".
#' @param legend Logical; draw a group/event legend in the right margin. If NULL (default), on when
#' `group.by` is set or `events` are supplied.
#' @param main Plot title. Defaults to "Periods of operation".
#' @param cex Global expansion factor for all text. Defaults to 1.
#' @template deviceArgs
#'
#' @return Invisibly, a tidy per-row coverage table: `row`, `group`, `n_deployments`, `first`, `last`,
#' `active_days`, `gap_days` and `coverage` (active fraction of the first-to-last span).
#' @seealso \code{\link{checkDeployments}}, \code{\link{importDeployments}}, \code{\link{plotAbacus}}
#' @examples
#' # receiver operating-period timeline from the deployment log
#' plotDeployments(rays_deployments)
#' @export

plotDeployments <- function(deployments,
                            events = NULL,
                            deployment.deploy.col = "deploy",
                            deployment.recover.col = "recover",
                            group.by = NULL,
                            row.by = "receiver",
                            end = NULL,
                            merge.gaps = 1,
                            sort.by = c("group", "start", "name"),
                            color.pal = NULL,
                            single.color = "grey20",
                            bar.height = 0.45,
                            dividers = TRUE,
                            top.band = "%Y",
                            date.interval = "auto",
                            date.start = 1,
                            date.format = NULL,
                            background.color = "grey96",
                            events.pch = 8,
                            events.color = "grey15",
                            events.cex = 1,
                            events.label = "event",
                            legend = NULL,
                            main = "Periods of operation",
                            cex = 1,
                            file = NULL,
                            width = NULL,
                            height = NULL,
                            res = 300) {

  ##############################################################################
  # Checks + compute ###########################################################
  ##############################################################################

  sort.by <- match.arg(sort.by)
  d <- as.data.frame(deployments)
  errors <- c()
  if(!row.by %in% names(d)) errors <- c(errors, sprintf("'row.by' column ('%s') not found in 'deployments'.", row.by))
  if(!deployment.deploy.col %in% names(d)) errors <- c(errors, sprintf("'deployment.deploy.col' column ('%s') not found in 'deployments'.", deployment.deploy.col))
  if(!is.null(group.by) && !group.by %in% names(d)) errors <- c(errors, sprintf("'group.by' column ('%s') not found in 'deployments'.", group.by))
  if(length(errors) == 0 && !inherits(d[[deployment.deploy.col]], "POSIXct")) errors <- c(errors, "'deployment.deploy.col' must be a POSIXct date-time column.")
  if(length(errors) > 0) stop(paste(c("", paste0("- ", errors)), collapse = "\n"), call. = FALSE)
  if(!deployment.recover.col %in% names(d)) d[[deployment.recover.col]] <- as.POSIXct(NA)

  tz <- .dataTZ(d[[deployment.deploy.col]])
  if(is.null(legend)) legend <- !is.null(group.by) || !is.null(events)
  cex_axis <- 0.8 * cex; cex_lab <- 1.0 * cex; cex_legend <- 0.75 * cex; cex_band <- 0.75 * cex

  row_key <- as.character(d[[row.by]])
  start <- d[[deployment.deploy.col]]; stop_ <- d[[deployment.recover.col]]
  study_end <- if(!is.null(end)) as.POSIXct(end, tz = tz) else max(c(start, stop_), na.rm = TRUE)
  stop_[is.na(stop_)] <- study_end
  keep <- !is.na(start) & !is.na(row_key)
  if(!any(keep)) stop("No deployments with a valid row key and deploy date to plot.", call. = FALSE)
  row_key <- row_key[keep]; start <- start[keep]; stop_ <- stop_[keep]
  grp_vec <- if(!is.null(group.by)) as.character(d[[group.by]])[keep] else rep(NA_character_, length(row_key))
  row_group <- tapply(grp_vec, row_key, function(x){ x <- x[!is.na(x)]; if(length(x)) x[1] else NA_character_ })

  # merge each unit's intervals, folding gaps <= merge.gaps days; keep longer gaps as separate segments
  merge_row <- function(s, e){
    o <- order(s); s <- s[o]; e <- e[o]; segs <- list(); cs <- s[1]; ce <- e[1]
    for(i in seq_along(s)[-1]){
      if(as.numeric(difftime(s[i], ce, units = "days")) <= merge.gaps) ce <- max(ce, e[i])
      else { segs[[length(segs) + 1]] <- c(cs, ce); cs <- s[i]; ce <- e[i] }
    }
    segs[[length(segs) + 1]] <- c(cs, ce); segs
  }
  rows <- unique(row_key)
  segs_by_row <- lapply(rows, function(r) merge_row(start[row_key == r], stop_[row_key == r]))
  names(segs_by_row) <- rows

  first_by <- vapply(segs_by_row, function(s) min(vapply(s, `[`, numeric(1), 1)), numeric(1))
  ord <- switch(sort.by,
                group = order(as.character(row_group[rows]), first_by),
                start = order(first_by),
                name  = order(rows))
  rows <- rows[ord]; segs_by_row <- segs_by_row[rows]
  g <- as.character(row_group[rows]); n <- length(rows)

  groups <- unique(g[!is.na(g)])
  if(length(groups) == 0){ row_col <- rep(single.color, n); group_cols <- character(0) }
  else {
    if(is.null(color.pal)) color.pal <- .okabe_ito_pal(length(groups))
    group_cols <- stats::setNames(color.pal[seq_along(groups)], groups)
    row_col <- ifelse(is.na(g), single.color, group_cols[g])
  }

  # per-row coverage table (also the invisible return)
  as_t <- function(x) as.POSIXct(x, origin = "1970-01-01", tz = tz)
  coverage <- do.call(rbind, lapply(seq_len(n), function(i){
    s <- segs_by_row[[i]]
    active <- sum(vapply(s, function(z) as.numeric(difftime(as_t(z[2]), as_t(z[1]), units = "days")), numeric(1)))
    span <- as.numeric(difftime(as_t(max(vapply(s, `[`, numeric(1), 2))), as_t(min(vapply(s, `[`, numeric(1), 1))), units = "days"))
    data.frame(row = rows[i], group = g[i], n_deployments = length(s),
               first = as_t(min(vapply(s, `[`, numeric(1), 1))), last = as_t(max(vapply(s, `[`, numeric(1), 2))),
               active_days = round(active, 1), gap_days = round(span - active, 1),
               coverage = if(span > 0) round(active / span, 3) else 1, stringsAsFactors = FALSE)
  }))
  rownames(coverage) <- NULL

  # events overlay: resolve the row key and time columns
  ev_x <- ev_y <- numeric(0)
  if(!is.null(events)){
    ev <- as.data.frame(events)
    if(!row.by %in% names(ev)) stop(sprintf("'events' must contain a '%s' column (the row key).", row.by), call. = FALSE)
    tcol <- if("time" %in% names(ev)) "time" else names(ev)[vapply(ev, inherits, logical(1), "POSIXct")][1]
    if(is.na(tcol)) stop("'events' must contain a date-time column named 'time' (or a POSIXct column).", call. = FALSE)
    yi <- match(as.character(ev[[row.by]]), rows); ok <- !is.na(yi)
    ev_x <- as.numeric(ev[[tcol]])[ok]; ev_y <- yi[ok]
  }

  xr <- range(c(start, stop_))
  .printDeploymentsSummary(n_rows = n, row.by = row.by, groups = groups,
                           n_deployments = nrow(d[keep, , drop = FALSE]), xr = xr,
                           total_active = sum(coverage$active_days), n_events = length(ev_x))


  ##############################################################################
  # Device + inch-stable layout (ported from plotAbacus) #######################
  ##############################################################################

  if(!is.null(file)){
    .beginFileOutput(file, width, height, res,
                     w.rule = list(base = 9, slope = 0, n = 0, lo = 7, hi = 20),
                     h.rule = list(base = 1.6, slope = 0.28, n = n, lo = 3, hi = 40),
                     crowd.unit = "rows")
    on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  }
  original_par <- .savePar(); on.exit(.restorePar(original_par), add = TRUE)

  char_w_in <- par("cin")[1]; char_h_in <- par("cin")[2]
  band_gap_in <- 0.10; band_height_in <- 0.17
  bottom_in <- 0.75
  left_in <- max(strwidth(rows, units = "inches", cex = cex_axis)) + 0.30
  top_in <- 0.30 + (if(!isFALSE(top.band)) band_gap_in + band_height_in else 0) + (if(!is.null(main)) 0.32 else 0)
  if(legend){
    leg_labels <- c(if(length(groups)) groups, if(!is.null(events)) events.label)
    right_in <- max(nchar(leg_labels)) * char_w_in * cex_legend + 0.40
  } else right_in <- 0.30
  par(mai = c(bottom_in, left_in, top_in, right_in), mgp = c(2.5, 0.6, 0), xpd = TRUE)

  plot(NA, xlim = as.numeric(xr), ylim = c(n + 1, 0), axes = FALSE, xlab = "", ylab = "")
  usr <- par("usr"); x0 <- usr[1]; x1 <- usr[2]
  rect(x0, 0, x1, n + 1, col = background.color, border = NA)
  segments(x0, seq_len(n), x1, seq_len(n), lty = 3, lwd = 0.5, col = "grey55")   # row guides


  ##############################################################################
  # Draw: bars, events, labels, axis, band, title, legend ######################
  ##############################################################################

  # operating-period bars (broken across real gaps)
  for(i in seq_len(n)) for(s in segs_by_row[[i]])
    rect(s[1], i - bar.height / 2, s[2], i + bar.height / 2, col = row_col[i], border = NA)

  # subtle separator lines between contiguous groups (only meaningful when group-sorted)
  if(dividers && length(groups) > 1 && sort.by == "group" && n > 1){
    b <- which(g[-1] != g[-n])
    if(length(b)) segments(x0, b + 0.5, x1, b + 0.5, col = "grey45", lwd = 0.6)
  }

  if(length(ev_x)) points(ev_x, ev_y, pch = events.pch, col = events.color, cex = events.cex, lwd = 1.4, xpd = NA)

  # row labels, coloured by group, in the left margin
  for(i in seq_len(n)) mtext(rows[i], side = 2, at = i, las = 1, line = 0.3, col = row_col[i], cex = cex_axis)
  mtext(tools::toTitleCase(row.by), side = 2, line = left_in / char_w_in - 1.2, cex = cex_lab)   # y-title from row.by

  # calendar x-axis: "auto" = calendar-aware breaks; a numeric date.interval = every n-th label.
  # 'cdates' (daily grid over the span) is shared by the manual axis and the year band below.
  cdates <- seq.POSIXt(lubridate::floor_date(xr[1], "day"), lubridate::ceiling_date(xr[2], "day"), by = "day")
  if(identical(date.interval, "auto")){
    ax <- .prettyDateAxis(xr[1], xr[2], format = date.format, band.format = top.band)
    axis(1, at = as.numeric(ax$at), labels = ax$labels, cex.axis = cex_axis, tck = -0.02)
    if(length(ax$minor)) axis(1, at = as.numeric(ax$minor), labels = FALSE, tck = -0.012, lwd.ticks = 0.5)
  }else{
    if(!is.numeric(date.interval) || date.interval < 1) stop("'date.interval' must be \"auto\" or a positive integer.", call. = FALSE)
    fmt <- if(is.null(date.format)) "%b" else date.format
    lab <- strftime(cdates, fmt, tz = tz); rl <- rle(lab)
    key <- paste0(lab, "_", rep(seq_along(rl$lengths), rl$lengths)); uniq <- unique(key)
    first_idx <- vapply(uniq, function(u) min(which(key == u)), integer(1))
    disp <- uniq[seq(date.start, length(uniq), by = date.interval)]
    disp_idx <- vapply(disp, function(u) min(which(key == u)), integer(1))
    axis(1, at = as.numeric(cdates[disp_idx]), labels = sub("_.*", "", disp), cex.axis = cex_axis, tck = -0.02)
    axis(1, at = as.numeric(cdates[first_idx]), labels = FALSE, tck = -0.012, lwd.ticks = 0.5)
  }
  mtext("Date", side = 1, line = 2.0, cex = cex_lab)

  # black year band, inch-anchored above the panel (y = 0 is the panel top)
  box_top_in <- grconvertY(0, "user", "inches")
  if(!isFALSE(top.band)){
    r <- rle(strftime(cdates, top.band, tz = tz)); ends <- cumsum(r$lengths); starts <- c(1, utils::head(ends, -1) + 1)
    centers <- (starts + ends) / 2; wf <- r$lengths / length(cdates)
    y0m <- grconvertY(box_top_in + band_gap_in, "inches", "user")
    y1m <- grconvertY(box_top_in + band_gap_in + band_height_in, "inches", "user")
    keepl <- rep(TRUE, length(r$values))
    if(wf[1] < 0.06) keepl[1] <- FALSE
    if(wf[length(wf)] < 0.06) keepl[length(keepl)] <- FALSE
    rect(x0, y0m, x1, y1m, col = "black", border = "black", xpd = NA)
    if(length(ends) > 1) segments(as.numeric(cdates[ends[-length(ends)]]), y0m, y1 = y1m, col = "white", lwd = 1.5, xpd = NA)
    if(any(keepl)) text(as.numeric(cdates[round(centers[keepl])]), mean(c(y0m, y1m)),
                        labels = r$values[keepl], col = "white", cex = cex_band, font = 2, xpd = NA)
  }

  if(!is.null(main)){
    title_y_in <- box_top_in + (if(!isFALSE(top.band)) band_gap_in + band_height_in else 0) + 0.16
    text(mean(c(x0, x1)), grconvertY(title_y_in, "inches", "user"), labels = main, font = 2, cex = cex_lab * 1.25, xpd = NA)
  }

  # legend stacked in the right margin (event marker, then group swatches)
  if(legend && length(leg_labels)){
    x_leg <- grconvertX(grconvertX(x1, "user", "inches") + 0.12, "inches", "user")
    y_cur <- grconvertY(0, "user", "inches")
    place <- function(y_in, ...){
      co <- graphics::legend(x = x_leg, y = grconvertY(y_in, "inches", "user"), bty = "n", xpd = NA, cex = cex_legend, ...)
      y_in - (grconvertY(co$rect$top, "user", "inches") - grconvertY(co$rect$top - co$rect$h, "user", "inches")) - 0.10
    }
    if(!is.null(events)) y_cur <- place(y_cur, legend = events.label, pch = events.pch, col = events.color, pt.cex = 1.1, y.intersp = 1.4)
    if(length(groups)) place(y_cur, legend = groups, fill = group_cols[groups], border = NA, y.intersp = 1.4)
  }

  rect(x0, 0, x1, n + 1, col = NA, border = "black")   # panel box on top
  invisible(coverage)
}


#######################################################################################################
## Internal helper ####################################################################################

#' @keywords internal
#' @noRd
.printDeploymentsSummary <- function(n_rows, row.by, groups, n_deployments, xr, total_active, n_events){
  kv <- .kv
  .summaryOpen("Deployment timeline")
  kv("Rows", sprintf("%d (by %s)", n_rows, row.by))
  if(length(groups)) kv("Groups", sprintf("%d (%s)", length(groups), paste(groups, collapse = ", ")))
  kv("Deployments", n_deployments)
  kv("Span", sprintf("%s to %s", strftime(xr[1], "%Y-%m-%d"), strftime(xr[2], "%Y-%m-%d")))
  kv("Receiver-days", format(round(total_active), big.mark = ","))
  if(n_events > 0) kv("Events", n_events)
  .summaryClose()
}

#######################################################################################################
#######################################################################################################
