##########################################################################################################
## Create migrations table ###############################################################################
##########################################################################################################

#' Create migrations table
#'
#' @description Generates a table containing migration metrics such as the number of movements between sites,
#' the number of migrating individuals, and other related statistics. Optionally, visualizations of the
#' distribution of some of these metrics can also be produced (e.g., departure and arrival hour, departure and arrival month).
#'
#' @inheritParams setDefaults
#' @param data A data frame containing animal detections.
#' @param spatial.col Name of the variable/column containing the spatial information
#' used to calculate transitions. It can include receiver IDs in case all movements/transitions
#' are of interest, or contain (for example) location/habitat classes for broader-scale analyses.
#' @param id.metadata A data frame containing additional information about the tagged animals,
#' such as length, sex or transmitter type. It should contain a column indicating animal IDs,
#' as specified by the id.col argument. If the `von.bertalanffy` argument is set to TRUE, this data frame
#' must also include a "length" column, which specifies the length of each individual at the time of tagging.
#' @param id.groups Optional. A list specifying groups of animal IDs. This parameter allows for the estimation of metrics
#' independently for subsets of animals and facilitates visual aggregation of individuals belonging to the same category
#' (e.g., different species or life stages).
#' @param plot Boolean. If true a plot is generated containing all the stats together with temporal metrics.
#' If false, a summary data frame is returned.
#' @param plot.stats Indicates which variables to plot. Available options: "all", "depart.hour",
#' "arrival.hour", "depart.month", "arrival.month", "duration". Defaults to "all".
#' @param same.scale Standardizes scales across plots and id groups.
#' @param von.bertalanffy Optional. Predict lengths at the moment of departure using a Von Bertalanffy Growth model.
#' @param VBGF.params A list containing the required parameters to fit von Bertalanffy growth curves (VBGF), including parameters such as 'Linf', 'K', and 't0'.
#' See the \code{\link[TropFishR]{VBGF}} function for detailed information on these parameters. If the data includes multiple species,
#' provide separate parameter sets for each, corresponding to the order of `id.groups` (a list of lists, where each inner list
#' contains the parameters specific to each group or species).
#' @param cex.main Numeric. Determines the size of main titles (used only if `plot` is set to TRUE). Defaults to 1.1.
#' @param cex.axis Numeric. Determines the size of axis labels (used only if `plot` is set to TRUE). Defaults to 0.7.
#' @param cex.table Numeric. Determines the size of text in the summary table (used only if `plot` is set to TRUE). Defaults to 0.95.
#' @export


migrationsTable <- function(data, id.col=getDefaults("id"), datetime.col=getDefaults("datetime"),
                            spatial.col, tagging.dates=getDefaults("tagdates"),
                            id.metadata=NULL, id.groups=NULL, plot=TRUE, plot.stats="all",
                            same.scale=TRUE, von.bertalanffy=FALSE, VBGF.params=NULL,
                            cex.main=1.1, cex.axis=0.7, cex.table=0.95) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # perform argument checks and return reviewed parameters
  reviewed_params <- .validateArguments()
  data <- reviewed_params$data
  tagging.dates <- reviewed_params$tagging.dates


  # perform argument checks for von bertalanffy growth models
  if(von.bertalanffy==TRUE){
    if (!requireNamespace("TropFishR", quietly=TRUE)) stop("The 'TropFishR' package is required for this function but is not installed. Please install 'TropFishR' using install.packages('TropFishR') and try again.", call.=FALSE)
    if(is.null(id.metadata) | !c("length") %in% colnames(id.metadata)) stop("'id.metadata' with a 'length' column required to apply the Von Bertalanffy Growth Model", call.=FALSE)
    if(!id.col %in% colnames(id.metadata)) stop(paste("'", id.col, "' column required in 'id.metadata'", sep=""), call.=FALSE)
    if(is.null(VBGF.params)) stop("Please supply VBGF.params when von.bertalanffy is set to true", call.=FALSE)
    if(!inherits(VBGF.params[[1]], "list")) VBGF.params<-list(VBGF.params)
    if(is.null(id.groups) && length(VBGF.params)>1) stop("Multiple 'VBGF.params' supplied, but no 'id.groups' were defined. Please provide 'id.groups' or reduce 'VBGF.params' to a single element.", call.=FALSE)
    if(!is.null(id.groups) && length(id.groups)!=length(VBGF.params)){
      warning("The length of 'VBGF.params' does not match the length of 'id.groups'. The same parameters will be applied to all IDs", call.=FALSE)
      VBGF.params <- rep(VBGF.params, length(id.groups))
    }
    cat("Applying Von Bertalanffy Growth curve to predict lengths at departure times\n")
  }

  # check temporal stats
  available_stats <- c("depart.hour", "arrival.hour", "depart.month", "arrival.month", "duration")
  if(any(plot.stats=="all")) plot.stats <- available_stats
  if(any(!plot.stats %in% available_stats)) stop(paste("Invalid plot.stats argument. Please select one or more variables from the possible options:", paste(available_stats, collapse=", ")), call.=FALSE)

  # convert sites var to factor
  if(!inherits(data[,spatial.col], "factor")){
    warning(paste("Converting", spatial.col, "column to factor"), call.=FALSE)
    data[,spatial.col] <- as.factor(data[,spatial.col])
  }
  ordered_sites <- levels(data[,spatial.col])


  # set id groups and split data
  if(is.null(id.groups)) id.groups <- list(levels(data[,id.col]))
  data_group <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])
  nindividuals <- lapply(id.groups, length)


  # convert tagging dates to data frame
  tagging.dates <- data.frame("id"=levels(data[,id.col]), "tagging.dates"=tagging.dates)
  colnames(tagging.dates)[1] <- id.col


  ##############################################################################
  ## Calculate metrics #########################################################
  ##############################################################################

  final_table <- list()
  graphs_data <- list()

  for(g in 1:length(id.groups)){

    # calculate nº of transitions per individual
    data_individual <- split(data_group[[g]], f=data_group[[g]][,id.col], drop=T)
    transitions <- list()
    transition_times <- list()
    for(i in 1:length(data_individual)) {
      # count transitions #####################################
      sites_data <- data_individual[[i]]
      sites_data <- sites_data[order(sites_data[,datetime.col]),]
      sites_seq <- sites_data[, spatial.col]
      movements <- paste0(head(sites_seq,-1), " --> ", tail(sites_seq,-1))
      result <- table(movements)
      result <- result[order(match(names(result), unique(movements)))]
      site1 <- sub(" -->.*", "", names(result))
      site2 <- sub(".*--> ", "", names(result))
      stationary <- which(site1==site2)
      names(result)[stationary] <- site1[stationary]
      transitions[[i]] <- result
      # get indices of transitions  ##########################
      sites_rle <- rle(as.character(sites_data[,spatial.col]))
      transition_pos <- cumsum(sites_rle$lengths)
      transition_pos <- transition_pos[-length(transition_pos)]
      # retrieve datetimes of transitions
      consec_movements <- rle(movements)$values
      site1 <- sub(" -->.*", "", consec_movements)
      site2 <- sub(".*--> ", "", consec_movements)
      stationary <- which(site1==site2)
      if(length(transition_pos)>=1){
        transition_times[[i]] <- data.frame("transition"=consec_movements[-stationary],
                                            "departure"=sites_data[,datetime.col][transition_pos],
                                            "arrival"=sites_data[,datetime.col][transition_pos+1],
                                            "id"=rep(unique(sites_data[,id.col]), length(transition_pos)))
      }
    }

    # merge individual results
    migrations <- plyr::rbind.fill(lapply(transitions, function(x) as.data.frame.matrix(t(x))))
    rownames(migrations) <- names(data_individual)
    transition_times <- transition_times[!unlist(lapply(transition_times, is.null))]
    transition_times <- do.call("rbind", transition_times)
    colnames(transition_times)[ncol(transition_times)] <- id.col

    # calculate time since tagging
    if(!is.null(tagging.dates)){
      transition_times <- plyr::join(transition_times, tagging.dates, by=id.col, type="left")
      transition_times$days_between_depart_tag <- difftime(transition_times$departure, transition_times[,ncol(transition_times)], units="days")
      transition_times$days_between_depart_tag <- round(as.numeric(transition_times$days_between_depart_tag))
      transition_times <- transition_times[,-(ncol(transition_times)-1)]
    }

    # calculate temporal metrics
    transition_times$departure_hour <- strftime(transition_times$departure, "%H", tz="UTC")
    transition_times$arrival_hour <- strftime(transition_times$arrival, "%H", tz="UTC")
    transition_times$departure_month <- strftime(transition_times$departure, "%m", tz="UTC")
    transition_times$arrival_month <- strftime(transition_times$arrival, "%m", tz="UTC")
    transition_times$duration_h <- as.numeric(difftime(transition_times$arrival, transition_times$departure, units="h"))
    graphs_data[[g]] <- transition_times


    ############################################################################
    ## Create final table ######################################################
    ############################################################################

    # create table with nº movements per migration type
    movements <- colSums(migrations, na.rm=T)
    migrations_table <- reshape2::melt(movements)
    colnames(migrations_table) <- "Movements"
    migrations_table$Type <- rownames(migrations_table)
    rownames(migrations_table) <- NULL
    migrations_table <- migrations_table[,c(2,1)]

    # add nº individuals per migration type
    migrating_ids <- apply(migrations, 2, function(x) rownames(migrations)[which(!is.na(x))])
    migrations_table$Individuals <- reshape2::melt(unlist(lapply(migrating_ids, length)))$value
    nids_percent <- round(migrations_table$Individuals/nindividuals[[g]]*100)
    migrations_table$'Individuals (%)' <- paste0(nids_percent, "%")

    # calculate additional tag metrics per migration type
    if(!is.null(id.metadata)){
      migrating_ids_data <- lapply(migrating_ids, function(x) setNames(data.frame(x), id.col))
      migrating_ids_data <- lapply(migrating_ids_data, function(x) plyr::join(x, id.metadata, by=id.col, type="left"))
      col_classes <- reshape2::melt(lapply(id.metadata, class))
      colnames(col_classes) <- c("class", "column")
      col_classes <- col_classes[col_classes$column!=id.col,]
      numeric_cols <- col_classes$column[col_classes$class %in% c("numeric", "integer")]
      character_cols <- col_classes$column[col_classes$class %in% c("factor", "character")]
      for(n in numeric_cols){
        mean_values <- reshape2::melt(unlist(lapply(migrating_ids_data, function(x) mean(x[,n], na.rm=T))))$value
        se_values <- reshape2::melt(unlist(lapply(migrating_ids_data, function(x) plotrix::std.error(x[,n], na.rm=T))))$value
        digits <- max(.decimalPlaces(id.metadata[,n]), na.rm=T)+1
        mean_values <- data.frame("mean"=paste(sprintf(paste0("%.", digits, "f"), mean_values), "\u00b1", sprintf(paste0("%.", digits, "f"), se_values)))
        mean_values <- sapply(mean_values, function(x) gsub(" \u00b1 NA", "", x, fixed=T))
        colnames(mean_values) <- paste0("Mean ", tools::toTitleCase(n))
        migrations_table <- cbind(migrations_table, mean_values)
      }
      for(c in character_cols){
        count_values <- lapply(migrating_ids_data, function(x) table(x[,c]))
        count_values <- lapply(count_values, function(x) paste(x, names(x), collapse = " | "))
        na_values <- lapply(migrating_ids_data, function(x) length(which(is.na(x[,c]))))
        count_values <- reshape2::melt(unlist(mapply(function(a,b) paste(a, "|", b, "NA"), a=count_values, b=na_values)))$value
        count_values <- gsub(" | 0 NA", "", count_values, fixed=T)
        count_values <- data.frame(count_values)
        colnames(count_values) <- tools::toTitleCase(c)
        count_values[unlist(na_values)==migrations_table$individuals,] <- "NA"
        migrations_table <- cbind(migrations_table, count_values)
      }
    }

    # predict growth if required
    if(von.bertalanffy==T){
      VBGF_group <- VBGF.params[[g]]
      transition_times <- plyr::join(transition_times, id.metadata[,c(id.col, "length")], by=id.col, type="left")
      ages_at_tagging <- TropFishR::VBGF(param=VBGF_group, L=transition_times$length, na.rm=F)
      ages_at_departure <- ages_at_tagging + (transition_times$days_between_depart_tag/365)
      transition_times$lengths_at_departure <- TropFishR::VBGF(param=VBGF_group, t=ages_at_departure, na.rm=F)
      growth_table <- stats::aggregate(transition_times$lengths_at_departure, by=list(transition_times$transition), mean, na.rm=T)
      growth_table$x <- sprintf("%.1f",  growth_table$x)
      se_lengths_at_departure <- stats::aggregate(transition_times$lengths_at_departure, by=list(transition_times$transition), function(x) plotrix::std.error(x))
      se_lengths_at_departure$x <- sprintf("%.1f",  se_lengths_at_departure$x)
      growth_table$x <- paste( growth_table$x, "\u00b1", se_lengths_at_departure$x)
      growth_table$x <- gsub(" \u00b1 NA", "", growth_table$x, fixed=T)
      colnames(growth_table) <- c("Type", "Depart Length")
      migrations_table <- plyr::join(migrations_table, growth_table, by="Type", type="left")
    }

    # order rows by site
    migrations_table$site1 <- sub(" -->.*", "", migrations_table$Type)
    migrations_table$site2 <-  sub(".*--> ", "", migrations_table$Type)
    migrations_table$site1 <- factor(migrations_table$site1, levels=ordered_sites)
    migrations_table$site2 <- factor(migrations_table$site2, levels=ordered_sites)
    migrations_table$transitions <- grepl("-->", migrations_table$Type, fixed=T)
    migrations_table <- migrations_table[order(migrations_table$site1, migrations_table$transitions, migrations_table$site2 ),]
    migrations_table <- migrations_table[,-which(colnames(migrations_table) %in% c("site1", "site2", "transitions"))]


    # add mean transition duration
    durations <- stats::aggregate(transition_times$duration_h, by=list(transition_times$transition), mean)
    colnames(durations) <- c("Type", "Mean")
    durations$Error <- stats::aggregate(transition_times$duration_h, by=list(transition_times$transition), function(x) plotrix::std.error(x))$x
    if(mean(durations$Mean)>72){
      durations[,-1] <- apply(durations[,-1], 2, function(x) sprintf("%.1f", x/24))
      units <- "(d)"
    }else{
      durations[,-1] <- apply(durations[,-1], 2, function(x) sprintf("%.1f", x))
      units <- "(h)"
    }
    durations$Duration <- paste(durations$Mean, "\u00b1", durations$Error)
    durations$Duration <- gsub(" \u00b1 NA", "", durations$Duration)
    colnames(durations)[4] <- paste(colnames(durations)[4], units)
    durations <- durations[,c(1,4)]
    migrations_table <- plyr::join(migrations_table, durations, by="Type", type="left")


    # add title
    if(length(id.groups)>1){
      table_title <- migrations_table[0,]
      table_title[1,] <- ""
      table_title$Type <- names(id.groups)[g]
      migrations_table <- rbind(table_title, migrations_table)
    }

    # replace NAs
    migrations_table[is.na(migrations_table)] <- "-"

    # save results to list
    final_table[[g]] <- migrations_table
  }

  # return results
  final_table <- do.call("rbind", final_table)

  if(plot==FALSE) return(final_table)


  ##############################################################################
  ## Prepare graph variables ##################################################
  ##############################################################################

  # convert duration hours to days if necessary
  all_durations <- as.numeric(unlist(lapply(graphs_data, function(x) x$duration_h)))
  duration_unit <- "h"
  if(mean(all_durations)>48){
    duration_unit <- "days"
    graphs_data <- lapply(graphs_data, function(x) {x$duration_h <- x$duration_h/24; return(x)})
  }
  duration_range <- range(as.numeric(unlist(lapply(graphs_data, function(x) x$duration_h))))

  # get graphs max
  panel_max <- lapply(graphs_data, function(x) c(table(x$transition, x$departure_hour), table(x$transition, x$arrival_hour),
                                                 table(x$transition, x$departure_month), table(x$transition, x$arrival_month)))
  panel_max <- max(unlist(panel_max))



  ##############################################################################
  ## Prepare layout variables ##################################################
  ##############################################################################

  # aggregate cols
  final_table$Individuals <- paste0(final_table$Individuals, " (", final_table$`Individuals (%)`,")")
  final_table$Individuals[final_table$Movements==""] <- ""
  final_table <- final_table[,-which(colnames(final_table)=="Individuals (%)")]
  if(any(grepl("Duration", plot.stats, fixed=TRUE))) final_table <- final_table[,-which(grepl("Duration", colnames(final_table), fixed=TRUE))]

  # prepare layout
  mat_table <-  matrix(1:raster::ncell(final_table), nrow=nrow(final_table), byrow=TRUE)
  nplots <- length(plot.stats)
  n_cells <- prod(dim(final_table))
  mat_graphs <- matrix((n_cells+1):(n_cells+nplots*nrow(final_table)), nrow=nrow(final_table), byrow=T)

  # calculate left and bottom plots for axes placement
  stationary <- which(!grepl("-->", final_table$Type, fixed=TRUE))
  group_dividers <- which(final_table$Movements=="")
  if(length(group_dividers)>0) stationary <- stationary[-group_dividers]
  bottom_plots <- stationary - 1
  bottom_plots <- c(bottom_plots, nrow(final_table))

  # set layout
  total_cells <- raster::ncell(cbind(mat_table, mat_graphs))
  final_mat <- cbind(mat_table, rep(total_cells+1, nrow(mat_table)), mat_graphs)
  final_mat <- rbind(final_mat, matrix(data=max(final_mat)+1, ncol=ncol(final_mat), nrow=1))
  layout_widths <- c(3.2, rep(0.9, ncol(mat_table)-1), 0.5, rep(1.5, ncol(mat_graphs)))
  layout_heights <- rep(1, nrow(mat_table))
  layout_heights[group_dividers] <- 0.3
  layout_heights[nrow(layout_heights)] <- 0.05
  layout(mat=final_mat, widths=layout_widths, heights=layout_heights)

  # set plot pars
  cex.axis <- 0.7
  tick_length <- -0.05


  ##############################################################################
  ## Generate plot #############################################################
  ##############################################################################

  # set par
  par(mar=c(0,0,0,0), oma=c(3,1,6,1),  mgp=c(2.5, 0.5, 0))

  # set xpos
  xpos <- c(0.05, rep(0.5, ncol(final_table)-1))
  text_adj <- c(0, rep(0.5, ncol(final_table)-1))

  # set headers
  headers <- colnames(final_table)
  headers[1] <- ""

  #######################################################################
  # draw table ##########################################################

  for(r in 1:nrow(final_table)){
    for(c in 1:ncol(mat_table)) {
      plot.new()
      if(final_table[r,c] %in% names(id.groups)) {
        text(0.05, 0, final_table[r,c], cex=cex.main, font=2, adj=0, xpd=NA)
      }else{
        text(xpos[c], 0.5, final_table[r,c], cex=cex.table, adj=text_adj[c])
      }
      if(r==1) title(main=headers[c], line=1.25, font=2, cex.main=cex.main, xpd=NA)
    }
  }

  # save usr
  usr <- par("usr")


  #######################################################################
  ## add barplots #######################################################

  # set graph headers
  plot.stats <- tools::toTitleCase(gsub(".", " ", plot.stats, fixed=T))

  # cycle through each row
  current_group <- ifelse(length(id.groups)>1, 0, 1)
  for(i in 1:nrow(final_table)){

    # adjust par settings
    if(same.scale==TRUE) par(mar=c(0.35, 0.5, 0.35, 0.5))
    else par(mar=c(0.35, 0.6, 0.35, 0.6), mgp=c(2.5, 0.4, 0))
    transition <- final_table$Type[i]

    # if divider row, skip to next
    if(transition %in% names(id.groups)){
      current_group <- current_group+1
      for(n in 1:nplots){
        par(mar=c(0,0,0,0))
        plot.new()
        if(i==1) title(main=plot.stats[n], line=1.25, font=2, cex.main=cex.main, xpd=NA)
      }
      next
    }

    # reorder transitions and skip to next if intrahabitat movement
    group_data <- graphs_data[[current_group]]
    transitions_order <- unique(final_table$Type)
    transitions_order <- transitions_order[!transitions_order %in% names(id.groups)]
    group_data <- group_data[order(match(group_data$transition, transitions_order)),]
    if(!transition %in% group_data$transition){
      for(n in 1:nplots){
        par(mar=c(0,0,0,0))
        plot.new()
        if(i==1) title(main=plot.stats[n], line=1.25, font=2, cex.main=cex.main, xpd=NA)
      }
      next
    }

    # prepare plot stats
    depart_hour <- as.numeric(group_data$departure_hour[group_data$transition==transition])
    arrival_hour <- as.numeric(group_data$arrival_hour[group_data$transition==transition])
    depart_month <- as.numeric(group_data$departure_month[group_data$transition==transition])
    arrival_month <- as.numeric(group_data$arrival_month[group_data$transition==transition])
    duration <- as.numeric(group_data$duration_h[group_data$transition==transition])
    depart_hour <- table(factor(depart_hour, levels=0:23))
    arrival_hour <- table(factor(arrival_hour, levels=0:23))
    depart_month <- table(factor(depart_month, levels=1:12))
    arrival_month <- table(factor(arrival_month, levels=1:12))

    # set y-axis scale
    if(same.scale==T){
      graphs_max<-panel_max
    }else{
      graphs_max <- max(c(depart_hour, arrival_hour, depart_month, arrival_month))*1.05
      duration_range <- range(as.numeric(group_data$duration_h), na.rm=T)
    }
    graphs_axis <- pretty(c(0,graphs_max))
    graphs_axis <- unique(as.integer(graphs_axis))


    # set plot list
    plot_data_list <- list(
      "Depart Hour"=list(data=depart_hour, color=adjustcolor("skyblue", alpha.f=0.4), labels=paste0(0:23, "h")),
      "Arrival Hour"=list(data=arrival_hour, color=adjustcolor("orange", alpha.f=0.4), labels=paste0(0:23, "h")),
      "Depart Month"=list(data=depart_month, color=adjustcolor("skyblue", alpha.f=0.4), labels=month.abb),
      "Arrival Month"=list(data=arrival_month, color=adjustcolor("orange", alpha.f=0.4), labels=month.abb)
    )

    # draw graphs
    for(stat in plot.stats) {
      if(stat %in% names(plot_data_list)) {
        item <- plot_data_list[[stat]]
        graphics::barplot(item$data, axes=FALSE, col=FALSE, border=FALSE, names.arg=FALSE, ylim=c(0, graphs_max))
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NA)
        b <- graphics::barplot(item$data, axes=FALSE, col=item$color, names.arg=FALSE, ylim=c(0, graphs_max), add=TRUE)
        if(i==1) title(main=stat, line=1.25, font=2, cex.main=cex.main, xpd=NA)
        if(i %in% bottom_plots) axis(1, at=b, labels=item$labels, cex.axis=cex.axis, tck=tick_length, xpd=NA)
        if(same.scale==TRUE && plot.stats[1]==stat) axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=cex.axis, tck=tick_length, las=1)
        else if(same.scale==FALSE) axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=cex.axis, tck=tick_length, las=1)
        box()
      } else if(stat=="Duration") {
        boxplot(duration, horizontal=TRUE, axes=FALSE, ylim=duration_range)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NA)
        boxplot(duration, horizontal=TRUE, col=adjustcolor("purple", alpha.f=0.3), axes=FALSE, ylim=duration_range, outpch=16, add=TRUE)
        if(i==1) title(main=paste0("Duration (", duration_unit, ")"), line=1.25, font=2, cex.main=cex.main, xpd=NA)
        duration_labels <- pretty(duration_range, n=6)
        duration_labels <- duration_labels[duration_labels >= duration_range[1] & duration_labels <= duration_range[2]]
        if(i %in% bottom_plots) axis(1, at=duration_labels, labels=duration_labels, cex.axis=cex.axis, tck=tick_length, xpd=NA)
        if(same.scale==TRUE && plot.stats[1]=="Duration") axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=cex.axis, tck=tick_length, las=1)
        else if(same.scale==FALSE) axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=cex.axis, tck=tick_length, las=1)
        box()
      }
    }

  }

  # draw box
  box("inner", lwd=1.5)
  box("outer", lwd=1.5)

}


#######################################################################################################
#######################################################################################################
#######################################################################################################
