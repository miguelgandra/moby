################################################################################
## Continuous Wavelet Transform (CWT) Analysis & Plotting                     ##
################################################################################

#' Continuous wavelet transform (CWT)
#'
#' @description Function that analyzes and plots periodic patterns in a time series using a
#' Continuous Wavelet Transform (CWT) framework. CWT analysis provides an alternative method
#' to the Fast Fourier Transform (FFT) or other time-frequency decomposition techniques,
#' enabling the examination of periodicities over time. This function serves mostly as a
#' wrapper for the \code{\link[wavScalogram]{cwt_wst}} function from the `wavScalogram` package
#' (Bolos & Benitez, 2022).
#'
#' @inheritParams setDefaults
#' @param data A data frame containing binned animal detections or other time-based measurements.
#' The time series does not need to be regular; gaps between time bins are allowed.
#' If there are gaps between measurements (i.e., missing time bins), the handling method is
#' controlled by the \code{gap.handling} parameter.
#' @param variable Name of the column containing the numeric variable to be analyzed.
#' @param plot.title A string specifying the title of the plot. By default, this title is
#' automatically generated as "Wavelet Power Spectrum" followed by the `variable` name.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species).
#' @param wavelet.type A character string specifying the wavelet type to be used in the continuous wavelet transform.
#' This is passed to the `wname` argument in the \code{\link[wavScalogram]{cwt_wst}} function.
#' Possible values are: "MORLET", "DOG", "PAUL", "HAAR", or "HAAR2". The default is "MORLET".
#' @param gap.handling Method for handling missing time bins in the time series. Options are:
#' \itemize{
#'      \item "zero" (default): Fill gaps with zeros (appropriate for detection/count data)
#'      \item "mean": Fill gaps with the mean of the non-missing values for that individual
#'      \item "locf": Last observation carried forward (use the last valid observation)
#'      \item "interpolate": Linear interpolation between valid observations (recommended for environmental data)
#' }
#' @param power.scaling Controls how the signal power (Z-axis color) is scaled. Options are:
#' \itemize{
#'      \item "linear" (default): Raw power values. Large peaks may obscure smaller signals.
#'      \item "quantile": Quantile-based breaks. Maximizes contrast by distributing colors evenly across the distribution. Can appear "noisy".
#'      \item "log": Logarithmic scaling (log10). Good compromise: compresses large peaks and reveals lower power signals without maximizing noise.
#'      \item "sqrt": Square-root scaling. milder compression than log.
#' }
#' @param upper.value Optional numeric value. Threshold for winsorizing the raw power values (before scaling).
#' Values above this threshold are set to this value. Useful for removing extreme outliers that skew the color scale.
#' @param upper.quant Optional numeric value between 0 and 1. Percentile threshold for winsorizing the raw power values.
#' e.g., 0.99 caps values at the 99th percentile.
#' @param same.scale Forces same spectral scale (zlims) across all plots,
#' allowing for density comparison between individuals.
#' @param mask.coi Logical. If TRUE, regions outside the cone of influence (COI)
#' are excluded from color scale calculations and plotted in black. This prevents
#' edge effects from dominating the color scale. Defaults to FALSE.
#' @param period.range The range of period scales (y-axis limits) to be considered,
#' specified in the units defined by \code{time.unit}. Defaults to c(3, 48) in hours.
#' @param axis.periods Periods to include/highlight on the y-axis, specified in the units
#' defined by \code{time.unit}. Defaults to c(6,12,16,24,48) in hours.
#' @param time.unit Time unit for `period.range` and `axis.periods`.
#' Options are "mins", "hours" and "days". Defaults to "hours".
#' @param color.pal Color palette. Defaults to the Jet colormap.
#' @param min.days Discard individuals that were detected in less than x days.
#' @param detrend.method Method for detrending the time series prior to CWT analysis.
#' Options are:
#' \itemize{
#'   \item "none" (default): No detrending is applied.
#'   \item "diff": Applies first-order differencing (\code{\link[base]{diff}}) to remove
#'   linear trends. Note that this analyzes the *rate of change* of the signal.
#'   \item "loess": Fits a LOESS regression (using \code{\link[stats]{loess}}) and
#'   performs the CWT on the residuals. This is effective for removing complex,
#'   non-linear trends. The smoothness is controlled by \code{loess.span}.
#' }
#' @param loess.span The span (alpha) parameter for LOESS detrending, used only if
#' \code{detrend.method = "loess"}. Controls the degree of smoothing; larger values
#' fit a smoother, less "wiggly" trend line. Defaults to 0.75.
#' @param standardize Logical. If TRUE, the time series for each individual is
#' standardized to have a mean of 0 and a standard deviation of 1 (i.e.,
#' Z-scored using \code{\link[base]{scale}}) prior to CWT analysis.
#' This is applied *after* detrending. This can help reveal periodicities
#' in signals with low amplitude relative to their mean. Defaults to FALSE.
#' @param date.format Date-time format (as used in \code{\link[base]{strptime}}),
#' defining the x-axis labels. Defaults to month ("%d/%b").
#' @param date.interval Number defining the interval between each
#' displayed date (x-axis label). Defaults to 4.
#' @param date.start Integer defining the first displayed date (can be used in combination
#'  with 'date.interval" to better control the x-axis labels). Defaults to 1.
#' @param cex.main Determines the size of the title(s). Defaults to 1.2.
#' @param cex.lab Determines the size of the y-axis and y-axis labels. Defaults to 1.1.
#' @param cex.axis Determines the size of the text labels on the axes. Defaults to 1.
#' @param cex.legend Determines the size of the text labels on the color legend. Defaults to 0.9.
#' @param legend.xpos Relative position of left and right edge of color bar on first axis (0-1).
#' Defaults to c(0.90, 0.915).
#' @param legend.ypos Relative position of left and right edge of color bar on second axis (0-1).
#' Defaults to c(0.15, 0.85).
#' @param cols Number of columns to arrange plots in when multiple individuals are displayed.
#' Ignored if \code{par.args$mfrow} is provided, in which case the layout is controlled
#' entirely by the user via \code{par.args}.
#' @param par.args Optional named list of graphical parameters to pass to [graphics::par()].
#' This allows fine control over the multi-panel layout and plot spacing.
#' Any arguments supplied here will override the internal defaults used by the function
#' (such as \code{mfrow}, \code{mar}, \code{oma}, or \code{mgp}).
#' For example, to reduce margins and adjust outer spacing:
#' \preformatted{
#' par.args = list(
#'    mar = c(3, 4, 2, 1),
#'    oma = c(1, 1, 2, 1)
#' )
#' }
#' By default, an empty list (\code{list()}) is used, which applies the function's built-in layout settings.
#' @param cores Number of CPU cores to use for the computations. Defaults to 1, which
#' means no parallel computing (single core).  If set to a value greater than 1,
#' the function will use parallel computing to speed up calculations.
#' Run \code{parallel::detectCores()} to check the number of available cores.
#' @param ... Further arguments passed to the \code{\link[wavScalogram]{cwt_wst}} function.
#' (e.g. border_effects, waverad).
#'
#' @references
#' Bolos, V. J., & Benitez, R. (2022). wavScalogram: an R package with scalogram wavelet tools for time series analysis. The R Journal, 14(2), 164-185.
#'
#' @export

plotCWTs <- function(data,
                     variable,
                     plot.title = paste("Wavelet Power Spectrum -", tools::toTitleCase(variable)),
                     id.col = getDefaults("ID"),
                     timebin.col = getDefaults("timebin"),
                     id.groups = NULL,
                     wavelet.type = "MORLET",
                     gap.handling = "zero",
                     power.scaling = "linear",
                     upper.value = NULL,
                     upper.quant = NULL,
                     same.scale = FALSE,
                     mask.coi = FALSE,
                     period.range = c(3, 48),
                     axis.periods = c(6, 12, 16, 24, 48),
                     time.unit = "hours",
                     color.pal = NULL,
                     min.days = NULL,
                     detrend.method = "none",
                     loess.span = 0.75,
                     standardize = FALSE,
                     date.format = "%d/%b",
                     date.interval = 4,
                     date.start = 1,
                     cex.main = 1.2,
                     cex.lab = 1.1,
                     cex.axis = 1,
                     cex.legend = 0.9,
                     legend.xpos = c(0.895, 0.91),
                     legend.ypos = c(0.15, 0.85),
                     cols = 1,
                     par.args = list(),
                     cores = 1,
                     ...) {


  # ============================================================================
  # 1. Parameter Validation and Initial Checks
  # ============================================================================

  # Perform internal argument checks
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data

  # Validate dependencies and input types
  errors <- c()
  if(!requireNamespace("wavScalogram", quietly=TRUE)) errors <- c(errors, "The 'wavScalogram' package is required but is not installed. Please install it using install.packages('wavScalogram').")
  if(!inherits(data[[variable]], c("numeric", "integer"))) errors <- c(errors, "The signal variable must be of class numeric or integer.")
  if (!is.numeric(period.range) || length(period.range) != 2) errors <- c(errors, "`period.range` must be a numeric vector of length 2 (min, max).")
  if(!time.unit %in% c("mins", "hours", "days")) errors <- c(errors, "Invalid time.unit. Use 'mins', 'hours', or 'days'.")
  if(!gap.handling %in% c("zero", "mean", "locf", "interpolate")) errors <- c(errors, "Invalid gap.handling method.")
  if(!detrend.method %in% c("none", "diff", "loess")) errors <- c(errors, "Invalid detrend.method.")
  if(detrend.method == "loess" && (!is.numeric(loess.span) || loess.span <= 0)) errors <- c(errors, "`loess.span` must be a positive numeric value.")
  if(!is.logical(standardize)) errors <- c(errors, "`standardize` must be TRUE or FALSE.")

  # Power scaling validation
  if(!power.scaling %in% c("linear", "quantile", "log", "sqrt")) errors <- c(errors, "Invalid power.scaling. Use 'linear', 'quantile', 'log', or 'sqrt'.")

  # Threshold validation
  if(!is.null(upper.value) && (!is.numeric(upper.value) || upper.value < 0)) errors <- c(errors, "`upper.value` must be a positive numeric value.")
  if(!is.null(upper.quant) && (!is.numeric(upper.quant) || upper.quant <= 0 || upper.quant > 1)) errors <- c(errors, "`upper.quant` must be a numeric value between 0 and 1.")

  # Check for zoo package if interpolation is required
  if(gap.handling %in% c("locf", "interpolate") && !requireNamespace("zoo", quietly=TRUE)) {
    errors <- c(errors, paste0("The 'zoo' package is required for gap.handling='", gap.handling, "'. Please install it."))
  }

  if(length(errors) > 0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop(paste0("\n", paste0("- ", stop_message, collapse="\n")), call.=FALSE)
  }

  # Save current par settings to restore them on exit
  original_par <- par(no.readonly=TRUE)
  # Remove read-only parameters to prevent warnings on restoration
  original_par[c("pin", "cin", "cra", "csi", "cxy", "din", "page")] <- NULL
  on.exit(par(original_par), add = TRUE)

  .printConsole("Calculating Wavelet Periodograms")
  data[[id.col]] <- droplevels(data[[id.col]])

  # Convert period ranges to the base unit (minutes) for calculation
  if (time.unit == "mins") {
    period.range_calc <- period.range
    time_factor <- 1
    unit_abbrev <- "mins"
  } else if (time.unit == "hours") {
    period.range_calc <- period.range * 60
    time_factor <- 60
    unit_abbrev <- "h"
  } else if (time.unit == "days") {
    period.range_calc <- period.range * 1440
    time_factor <- 1440
    unit_abbrev <- "days"
  }

  # ============================================================================
  # 2. Data Preparation and Time Series Construction
  # ============================================================================

  # Robustly infer sampling interval (dt) from the smallest non-zero difference between timebins
  sorted_times <- sort(unique(data[[timebin.col]]))
  if (length(sorted_times) < 2) {
    stop("At least two unique time bins are required to determine the sampling interval.", call. = FALSE)
  }
  interval <- difftime(sorted_times, dplyr::lag(sorted_times), units = "min")
  interval <- as.numeric(min(interval[interval > 0], na.rm = TRUE))

  if (!is.finite(interval) || interval == 0) {
    stop("Could not determine a valid sampling interval (dt > 0). Check 'timebin.col'.", call. = FALSE)
  }
  cat(paste0("Inferred sampling interval (dt): ", interval, " minutes\n"))

  # Create a complete sequence of time bins to handle irregularities/gaps
  min_time <- min(data[[timebin.col]])
  max_time <- max(data[[timebin.col]])
  all_timebins <- seq(min_time, max_time, by = paste(interval, "min"))

  # Initialize synchronization table
  cwt_table <- data.frame(all_timebins)
  colnames(cwt_table) <- timebin.col

  # Pivot data: populate the table for each individual
  for(id in levels(data[[id.col]])) {

    # Extract data for current ID
    id_data <- data[data[[id.col]] == id, c(timebin.col, variable)]

    # Handle duplicates by averaging
    if(any(duplicated(id_data[[timebin.col]]))) {
      id_data <- aggregate(id_data[[variable]], by = list(id_data[[timebin.col]]), FUN = mean, na.rm = TRUE)
      colnames(id_data) <- c(timebin.col, variable)
    }

    # Align data to the complete time sequence
    values <- rep(NA_real_, length(all_timebins))
    matched_idx <- match(id_data[[timebin.col]], all_timebins)
    values[matched_idx] <- id_data[[variable]]

    # Mark pre-detection and post-detection periods as NaN (distinct from missing internal gaps)
    non_na_idx <- which(!is.na(values))
    if(length(non_na_idx) > 0) {
      first_idx <- non_na_idx[1]
      last_idx  <- non_na_idx[length(non_na_idx)]

      if(first_idx > 1) values[1:(first_idx-1)] <- NaN
      if(last_idx < length(values)) values[(last_idx+1):length(values)] <- NaN
    } else {
      values[] <- NaN
    }

    cwt_table[[id]] <- values
  }

  # Expand period range slightly to prevent edge truncation during CWT
  period.range_calc[1] <- period.range_calc[1] - 30

  # Filter individuals based on minimum data requirements (min.days)
  if(!is.null(min.days)){
    data$day <- strftime(data[[timebin.col]], "%Y-%m-%d", tz="UTC")
    days_detected <- by(data$day, data[[id.col]], function(x) length(unique(x)))
    selected_individuals <- which(as.numeric(days_detected) >= min.days)
    nindividuals <- length(selected_individuals)
    cat(paste(nindividuals, "individuals with >", min.days, "logged days\n"))
    cat(paste(nlevels(data[[id.col]]) - nindividuals, "individual(s) excluded\n"))
  } else {
    selected_individuals <- 1:nlevels(data[[id.col]])
    nindividuals <- length(selected_individuals)
  }

  # Extract time series for selected individuals and trim outer NaNs
  data_individual <- lapply(selected_individuals, function(i) {
    x <- cwt_table[, i+1]
    return(x[!is.nan(x)])
  })
  names(data_individual) <- levels(data[[id.col]])[selected_individuals]

  # Apply Gap Handling (crucial for valid CWT)
  cat(paste0("Handling gaps using method: '", gap.handling, "'\n"))
  data_individual <- lapply(data_individual, function(x) {
    if(gap.handling == "zero") {
      x[is.na(x)] <- 0
    } else if(gap.handling == "mean") {
      x[is.na(x)] <- mean(x, na.rm=TRUE)
    } else if(gap.handling == "locf") {
      x <- zoo::na.locf(x, na.rm=FALSE)
    } else if(gap.handling == "interpolate") {
      x <- zoo::na.approx(x, na.rm=FALSE)
    }
    return(x)
  })

  # Convert to ts objects
  data_ts <- lapply(data_individual, ts)


  # ============================================================================
  # 3. Layout Configuration
  # ============================================================================

  # Calculate grid layout based on groups or simple count
  if(!is.null(id.groups)){
    group_ids_selected <- lapply(id.groups, function(x) x[x %in% levels(data[[id.col]])[selected_individuals]])
    group_numbers <- lapply(group_ids_selected, length)
    group_rows <- lapply(group_numbers, function(x) ceiling(x/cols))
    rows <- do.call("sum", group_rows)
    group_plots <- lapply(group_rows, function(x) x*cols)

    # Map animals to plot indices
    animal_indexes <- mapply(function(nids, nplots) {
      if(nids < nplots) c(1:nids, rep(NA, nplots-nids)) else 1:nids
    }, nids=group_numbers, nplots=group_plots, SIMPLIFY=FALSE)

    for(i in 2:length(animal_indexes)){
      animal_indexes[[i]] <- animal_indexes[[i]] + max(animal_indexes[[i-1]], na.rm=TRUE)
    }
    animal_indexes <- unlist(animal_indexes, use.names=FALSE)
    background_pal <- grey.colors(length(id.groups), start=0.97, end=0.93)
  } else {
    rows <- ceiling(nindividuals/cols)
    nplots <- rows * cols
    animal_indexes <- if(nindividuals < nplots) c(1:nindividuals, rep(NA, nplots-nindividuals)) else 1:nindividuals
    background_col <- "gray96"
  }

  plot_layout <- matrix(animal_indexes, nrow=rows, ncol=cols, byrow=TRUE)
  plot_ids <- as.integer(apply(plot_layout, 1, function(x) x))
  plot_ids <- plot_ids[!is.na(plot_ids)]

  # Apply graphical parameters
  par.defaults <- list(mfrow = c(rows, cols), mar = c(4, 5, 4, 6), oma = c(2, 2, 3, 2), mgp = c(3, 0.8, 0))
  par.defaults[names(par.args)] <- par.args
  par(par.defaults)

  # Set default color palette if not provided
  if(is.null(color.pal)){
    color.pal <- .jet_pal(100)
    color.pal <- colorRampPalette(color.pal[5:100])(100)
  }


  # ============================================================================
  # 4. CWT Calculation
  # ============================================================================

  # Pre-processing: Detrending
  if (detrend.method == "diff") {
    cat("Detrending: Applying first-order differencing\n")
    data_ts <- lapply(data_ts, diff)
  } else if (detrend.method == "loess") {
    cat(paste0("Detrending: Applying LOESS smoothing (span = ", loess.span, ")\n"))
    for (id_name in names(data_ts)) {
      ts_data <- data_ts[[id_name]]
      if (length(ts_data) < 4) {
        warning(paste("Series", id_name, "too short for LOESS, skipping."), call. = FALSE)
        next
      }
      time_idx <- 1:length(ts_data)
      loess_fit <- try(stats::loess(ts_data ~ time_idx, span = loess.span), silent = TRUE)
      if (!inherits(loess_fit, "try-error")) {
        data_ts[[id_name]] <- stats::residuals(loess_fit)
      } else {
        warning(paste("LOESS failed for", id_name, "- using raw data."), call. = FALSE)
      }
    }
  }

  # Pre-processing: Standardization
  if(standardize) {
    cat("Standardizing: Applying Z-score scaling\n")
    data_ts <- lapply(data_ts, scale)
  }

  # Execute CWT (Parallel or Sequential)
  pb <- txtProgressBar(min=0, max=length(data_individual), initial=0, style=3)

  if(cores > 1) {
    cat(paste0("Starting parallel computation on ", cores, " cores\n"))
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    `%dopar%` <- foreach::`%dopar%`
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))

    cwts <- foreach::foreach(i=1:length(data_ts), .options.snow=opts, .packages="wavScalogram") %dopar% {
      wavScalogram::cwt_wst(data_ts[[i]], dt=interval, scales=c(period.range_calc, 20),
                            powerscales=TRUE, wname=wavelet.type, border_effects="BE",
                            makefigure=FALSE, energy_density=TRUE, figureperiod=TRUE, ...)
    }
  } else {
    cwts <- vector("list", length(data_ts))
    for(i in 1:length(data_ts)){
      cwts[[i]] <- wavScalogram::cwt_wst(data_ts[[i]], dt=interval, scales=c(period.range_calc, 20),
                                         powerscales=TRUE, wname=wavelet.type, border_effects="BE",
                                         makefigure=FALSE, energy_density=TRUE, figureperiod=TRUE, ...)
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)


  # ============================================================================
  # 5. Density and Scale Calculation
  # ============================================================================

  # Calculate power density values (Coefficient^2 / Scale)
  density_values <- vector("list", length(cwts))

  for(i in 1:length(cwts)) {
    cwt <- cwts[[i]]
    Z <- t(t(abs(cwt$coefs) ^ 2) / cwt$scales)

    # --- UPDATED: Apply Thresholds (Winsorizing) ---
    if(!is.null(upper.quant)){
      quant_limit <- quantile(Z, probs = upper.quant, na.rm=TRUE)
      Z[Z > quant_limit] <- quant_limit
    }
    if(!is.null(upper.value)){
      Z[Z > upper.value] <- upper.value
    }
    # -----------------------------------------------

    # --- UPDATED: Apply Power Transformations ---
    if(power.scaling == "log") {
      Z <- log10(Z + 1e-12) # Small offset to avoid log(0)
    } else if (power.scaling == "sqrt") {
      Z <- sqrt(Z)
    }
    # ------------------------------------------------

    # Optionally mask values outside the Cone of Influence (COI) for scale calculation
    if(mask.coi) {
      coi_mask <- matrix(TRUE, nrow = nrow(Z), ncol = ncol(Z))
      for(j in 1:nrow(Z)) {
        coi_scale <- cwt$coi_maxscale[j] * cwt$fourierfactor
        coi_mask[j, cwt$scales * cwt$fourierfactor > coi_scale] <- FALSE
      }
      density_values[[i]] <- Z[coi_mask]
    } else {
      density_values[[i]] <- as.vector(Z)
    }
  }

  # Determine global scale limits if same.scale is requested
  density_range <- range(unlist(density_values), na.rm=TRUE)
  global_breaks <- NULL

  if(same.scale && power.scaling == "quantile") {
    all_valid_z <- unlist(density_values)
    all_valid_z <- all_valid_z[is.finite(all_valid_z)]

    # Compute quantiles based strictly on positive signal to avoid background zero-inflation
    positive_z <- all_valid_z[all_valid_z > 0]

    if(length(positive_z) > 0) {
      global_breaks <- unique(quantile(positive_z, probs = seq(0, 1, length.out = length(color.pal) + 1)))
      global_breaks[1] <- min(all_valid_z, na.rm=TRUE) # Ensure coverage of 0
      global_breaks <- sort(unique(global_breaks))

      # Fallback if quantiles are degenerate
      if(length(global_breaks) != length(color.pal) + 1) {
        global_breaks <- seq(min(all_valid_z), max(all_valid_z), length.out = length(color.pal) + 1)
      }
    }
  }


  # ============================================================================
  # 6. Plot Generation
  # ============================================================================

  for (i in plot_ids) {

    if(is.na(i)) {
      plot.new()
      next
    }

    id <- selected_individuals[i]
    cwt <- cwts[[i]]
    Z <- t(t(abs(cwt$coefs) ^ 2) / cwt$scales)

    # --- UPDATED: Apply Thresholds (Winsorizing) ---
    if(!is.null(upper.quant)){
      quant_limit <- quantile(Z, probs = upper.quant, na.rm=TRUE)
      Z[Z > quant_limit] <- quant_limit
    }
    if(!is.null(upper.value)){
      Z[Z > upper.value] <- upper.value
    }
    # -----------------------------------------------

    # --- UPDATED: Apply Power Transformations ---
    if(power.scaling == "log") {
      Z <- log10(Z + 1e-12)
    } else if (power.scaling == "sqrt") {
      Z <- sqrt(Z)
    }
    # ------------------------------------------------

    # --- Determine Z-limits and Breaks ---
    if(same.scale) {
      zlim <- density_range
      current_breaks <- global_breaks
    } else {
      # Local scale calculation
      if(mask.coi) {
        coi_mask <- matrix(TRUE, nrow = nrow(Z), ncol = ncol(Z))
        for(j in 1:nrow(Z)) {
          coi_scale <- cwt$coi_maxscale[j] * cwt$fourierfactor
          coi_mask[j, cwt$scales * cwt$fourierfactor > coi_scale] <- FALSE
        }
        valid_Z <- Z[coi_mask]
      } else {
        valid_Z <- as.vector(Z)
      }

      zlim <- range(valid_Z, na.rm=TRUE)
      current_breaks <- NULL

      if(power.scaling == "quantile") {
        valid_Z <- valid_Z[is.finite(valid_Z)]
        positive_Z <- valid_Z[valid_Z > 0]

        if(length(positive_Z) > 0) {
          current_breaks <- unique(quantile(positive_Z, probs = seq(0, 1, length.out = length(color.pal) + 1)))
          current_breaks[1] <- min(zlim)
          current_breaks <- sort(unique(current_breaks))
        }
        # Fallback to linear if breaks are malformed
        if(!is.null(current_breaks) && length(current_breaks) != length(color.pal) + 1) {
          current_breaks <- seq(min(valid_Z), max(valid_Z), length.out = length(color.pal) + 1)
        }
      }
    }

    if(is.null(current_breaks)) {
      current_breaks <- seq(min(zlim), max(zlim), length.out = length(color.pal) + 1)
    }

    # --- Main Plot ---
    graphics::image(x=1:nrow(cwt$coefs), y=1:ncol(cwt$coefs), z=Z,
                    zlim = zlim,
                    breaks = current_breaks,
                    col=color.pal, axes=FALSE, useRaster=TRUE, main="", xlab="Date",
                    ylab = paste0("Period (", unit_abbrev, ")"),
                    frame.plot = TRUE, cex.lab = cex.lab, xaxs = "i", yaxs = "i")

    title(main=levels(data[[id.col]])[id], cex.main=cex.main, line=1)

    # --- X-Axis (Dates) ---
    all_dates <- cwt_table[,timebin.col][!is.nan(cwt_table[,id+1])]
    all_dates <- strftime(all_dates, date.format)

    # Identify unique date blocks for labeling
    consec_dates <- rle(all_dates)
    consec_dates <- paste0(all_dates, "_", rep(1:length(consec_dates$lengths), consec_dates$lengths))
    unique_dates <- unique(consec_dates)
    disp_dates <- unique_dates[seq(date.start, length(unique_dates), by=date.interval)]

    indexes <- unlist(lapply(unique_dates, function(x) min(which(consec_dates==x))))
    disp_indexes <- unlist(lapply(disp_dates, function(x) min(which(consec_dates==x))))
    disp_dates <- sub("\\_.*", "", disp_dates)

    axis(1, labels=disp_dates, at=disp_indexes, cex.axis=cex.axis)
    axis(1, labels=FALSE, at=indexes, tck=-0.02, lwd.ticks=0.5)

    # --- Y-Axis (Period) ---
    # Map log-scale periods to plot coordinates
    period_indexes <- .rescale(log2(axis.periods*time_factor),
                               from=log2(range(cwt$scales*cwt$fourierfactor)),
                               to=c(1, length(cwt$scales)))
    axis(2, at=period_indexes, labels=axis.periods, cex.axis=cex.axis, las=1)
    abline(h=period_indexes, lty=5, col="grey70", lwd=0.7)

    # --- Cone of Influence (COI) Masking ---
    # 1. Edge Effects COI (U-shaped)
    x_polygon <- 1:nrow(cwt$coefs)
    coi_y_values <- .rescale(log2(cwt$coi_maxscale * cwt$fourierfactor + 1e-20),
                             from = log2(range(cwt$scales * cwt$fourierfactor)),
                             to = c(1, length(cwt$scales)))

    y_bottom <- par("usr")[3]
    y_top    <- par("usr")[4]
    coi_y_clamped <- pmin(pmax(coi_y_values, y_bottom), y_top)

    mask_col <- if(mask.coi) "black" else adjustcolor("white", alpha.f = 0.5)

    polygon(c(x_polygon, rev(x_polygon)),
            c(coi_y_clamped, rep(y_top, length(x_polygon))),
            col = mask_col, border = TRUE)

    # 2. Minimum Scale COI (Bottom block)
    min_coi_log_value <- log2(cwt$coi_minscale * cwt$fourierfactor + 1e-20)
    min_coi_y_coord <- .rescale(min_coi_log_value,
                                from = log2(range(cwt$scales * cwt$fourierfactor)),
                                to = c(1, length(cwt$scales)))

    if (length(min_coi_y_coord) == 1 && is.numeric(min_coi_y_coord) && min_coi_y_coord > y_bottom) {
      rect(par("usr")[1], y_bottom, par("usr")[2], min_coi_y_coord, col = mask_col, border = TRUE)
    }

    # --- Color Legend ---
    # Draw legend only on the last plot of a set, or every plot if scales differ
    if (!same.scale || (same.scale && i == max(plot_ids))) {

      lab_scientific <- FALSE

      if(power.scaling == "quantile" && !is.null(current_breaks)){

        # ## Quantile Legend Logic ##
        leg_zlim <- c(0, 1)
        leg_zval <- seq(0, 1, length.out = 6)

        # Map ticks to break indices to find corresponding data values
        idx <- round(seq(1, length(current_breaks), length.out = 6))
        vals <- current_breaks[idx]
        leg_zlab <- vals # Initialize

        # Formatting
        sci_idx <- which(abs(vals) > 0 & (abs(vals) < 0.1 | abs(vals) >= 1000))
        norm_idx <- which(abs(vals) == 0 | (abs(vals) >= 0.1 & abs(vals) < 1000))

        if(length(sci_idx) > 0) leg_zlab[sci_idx] <- formatC(vals[sci_idx], format="e", digits=2)
        if(length(norm_idx) > 0) leg_zlab[norm_idx] <- formatC(vals[norm_idx], format="f", digits=2)

        leg_digit <- 3

      } else if (power.scaling == "log") {

        # ## UPDATED: Log Legend Logic ##
        # Use visually evenly spaced ticks based on the log-transformed data
        leg_zlim <- zlim
        leg_zval <- pretty(zlim)
        leg_zval <- leg_zval[leg_zval >= min(zlim) & leg_zval <= max(zlim)]

        # Calculate labels by inverse-transforming (anti-log) the tick positions
        # Note: zlim is already log10.
        vals <- 10^leg_zval
        leg_zlab <- vals

        # Apply the same strict formatting as Quantile to make large numbers readable
        sci_idx <- which(abs(vals) > 0 & (abs(vals) < 0.1 | abs(vals) >= 1000))
        norm_idx <- which(abs(vals) == 0 | (abs(vals) >= 0.1 & abs(vals) < 1000))

        if(length(sci_idx) > 0) leg_zlab[sci_idx] <- formatC(vals[sci_idx], format="e", digits=2)
        if(length(norm_idx) > 0) leg_zlab[norm_idx] <- formatC(vals[norm_idx], format="f", digits=2)

        leg_digit <- 3 # Passed but mostly overridden by zlab strings

      } else {
        # ## Linear / Sqrt Legend Logic ##
        # Standard mapping where visual position equals data value.

        leg_zlim <- zlim
        leg_zval <- pretty(zlim)
        leg_zval <- leg_zval[leg_zval >= min(zlim) & leg_zval <= max(zlim)]
        leg_zlab <- NULL # No override

        # Determine digits dynamically
        digits_vec <- .decimalPlaces(leg_zval)
        leg_digit <- if(all(is.na(digits_vec))) 2 else max(digits_vec, na.rm=TRUE)

        # Force scientific if precision requirements are too high
        if(leg_digit > 3) lab_scientific <- TRUE
      }

      .colorlegend(col = color.pal,
                   zlim = leg_zlim,
                   zval = leg_zval,
                   zlab = leg_zlab,
                   posx = legend.xpos,
                   posy = legend.ypos,
                   main = "",
                   main.cex = 1,
                   digit = leg_digit,
                   lab.scientific = lab_scientific,
                   cex = cex.legend,
                   xpd = NA)
    }
  }

  # Add Main Title
  mtext(text=plot.title, side=3, line=1.6, outer=T, cex=cex.main, font=2)

  # Add Group Legend if applicable
  if(!is.null(id.groups)){
    label_pos <- rev(unlist(lapply(group_rows, function(x) x/2)))
    for(i in 2:length(group_rows)){label_pos[i] <- label_pos[i] + rev(group_rows)[[i-1]]}
    label_pos <- grconvertY(label_pos/rows, "ndc", "user")
    text(x=grconvertX(0.01, "ndc", "user"), y=label_pos, labels=rev(names(id.groups)),
         srt=90, cex=cex.main, font=2, xpd=NA)
  }

}
