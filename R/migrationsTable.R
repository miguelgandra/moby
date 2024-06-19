##########################################################################################################
## Create migrations table ###############################################################################
##########################################################################################################

#' Create migrations table
#'
#' @description Creates a table containing  migration metrics
#' (nº movements between sites, nº migrating individuals, etc.).
#'
#' @param data A data frame containing binned animal detections.
#' @param sites.col Variable containing locations.
#' @param aggregate.by Variable containing locations.
#' @param tag.dates A data frame containing two columns, one containing the animal IDs
#' (column name defined by 'id.col') and the other the tagging dates in POSIXct format.
#' @param animals.info A data frame containing additional information about the tagged animals,
#' such as length, sex or transmitter type. If more than one row exists per animal,
#' variables are collapsed before being merged with the remaining stats.
#' @param id.groups Optional. A list containing ID groups, used to
#' visually aggregate animals belonging to the same class (e.g. different species or life stages).
#' @param id.col Name of the column containing animal IDs Defaults to 'ID'.
#' @param plot Boolean. If true a plot is generated containing all the stats together with temporal metrics.
#' If false, a summary data frame is returned.
#' @param plot.stats Indicates which variables to plot. Available options: "all", "none", "depart.hour",
#' "arrival.hour", "depart.month", "arrival.month", "duration". Defaults to "all".
#' @param same.scale Standardizes scales across plots and id groups.
#' @param von.bertalanffy Optional. Predict lengths at the moment of departure using a Von Bertalanffy Growth model.
#' @param VBGF.params A list containing the required parameters to fit the  Von Bertalanffy Growth model.
#' e.g.: 'Linf', 'K' and 't0'. Check \code{\link[TropFishR]{VBGF}} function for further details.
#' @export


migrationsTable <- function(data, sites.col, aggregate.by=NULL, tag.dates=NULL, animals.info=NULL,
                            id.groups=NULL, id.col="ID", plot=T, plot.stats="all", same.scale=T, von.bertalanffy=F, VBGF.params=NULL) {


  ##############################################################################
  ## Initial checks ############################################################
  ##############################################################################

  # check if data contains id.col
  if(!id.col %in% colnames(data)){
    stop("ID column not found. Please assign the correct column name with 'id.col'")
  }

  # check if data contains id.col
  if(!is.null(tag.dates) & !id.col %in% colnames(tag.dates)){
    stop("ID column not found in tag.dates. Please assign the correct column name with 'id.col'")
  }

  if(von.bertalanffy==T & (is.null(animals.info) | !c("length") %in% colnames(animals.info))){
    stop("Please supply 'animals.info' with a 'length' column in order to calculate Von Bertalanffy Growth Model")
  }

  if(von.bertalanffy==T){
    if(is.null(animals.info) | !c("length") %in% colnames(animals.info)){
      stop("'animals.info' with a 'length' column required to apply the Von Bertalanffy Growth Model")
    }
    if(is.null(tag.dates)){
      stop("'tag.dates' required to apply the Von Bertalanffy Growth Model")
    }
    if(is.null(VBGF.params)){
      stop("Please supply VBGF.params when von.bertalanffy is set to true")
    }
    cat("Applying Von Bertalanffy Growth curve to predict lengths at departure\n")
  }

  # convert sites var to factor
  if(class(data[,sites.col])!="factor"){
    cat(paste("Warning: converting", sites.col, "column to factor\n"))
    data[,sites.col] <- as.factor(data[,sites.col])
  }

  # set id groups
  if(is.null(id.groups)) {
    id.groups <- list(levels(data[,id.col]))
  } else {
    if(any(duplicated(unlist(id.groups)))) {stop("Repeated ID(s) in id.groups")}
    if(any(!unlist(id.groups) %in% levels(data[,id.col]))) {"Some of the ID(s) in id.groups don't match the IDs in the data"}
  }
  data_group <- lapply(id.groups, function(x) data[data[,id.col] %in% x,])

  # check aggregate.by groups
  if(!is.null(aggregate.by)){
    if(class(data[,aggregate.by])!="factor") {
      cat(paste("Warning: converting", aggregate.by, "column to factor\n"))
      data[,aggregate.by] <- as.factor(data[,aggregate.by])
    }
    sites_table <- unique(data[,c(sites.col, aggregate.by)])
    colnames(sites_table) <- c("site", "group")
    if(any(duplicated(sites_table$site))){
      stop("At least one of the supplied sites has > 1 matching group
           (within the 'aggregate.by' column)\n")
    }
  }

  available_stats <- c("depart.hour", "arrival.hour", "depart.month", "arrival.month", "duration")
  if(any(plot.stats=="none")){plot.stats<-F
  }else{
    if(any(plot.stats=="all")){plot.stats <- available_stats}
    if(any(!plot.stats %in% available_stats)){
      stop(paste("Invalid plot.stats argument. Please select one or more variables from the possible options:", paste(available_stats, collapse=", ")))
    }
  }

  # retrieve total nº of individuals
  if(!is.null(id.groups)){
    nindividuals <- lapply(id.groups, length)
  }else{
    nindividuals <- list(nlevels(data[,id.col]))
  }


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
      sites_data <- sites_data[order(sites_data$timebin),]
      sites_seq <- sites_data[, sites.col]
      movements <- paste0(head(sites_seq,-1), " --> ", tail(sites_seq,-1))
      result <- table(movements)
      result <- result[order(match(names(result), unique(movements)))]
      site1 <- sub(" -->.*", "", names(result))
      site2 <- sub(".*--> ", "", names(result))
      stationary <- which(site1==site2)
      names(result)[stationary] <- site1[stationary]
      if(!is.null(aggregate.by)){
        replacement_vector <- setNames(as.character(sites_table$group), sites_table$site)
        for (site in names(replacement_vector)) {
          new_names <- gsub(site, replacement_vector[[site]], names(result))
        }
        names(result) <- new_names
        result <- table(rep(names(result), result))
      }
      transitions[[i]] <- result
      # get datetimes of transitions  ##########################
      sites_rle <- rle(as.character(sites_data[,sites.col]))
      # get indices of transitions
      transition_pos <- cumsum(sites_rle$lengths)
      transition_pos <- transition_pos[-length(transition_pos)]
      # retrieve timebins of transitions
      if(length(transition_pos)>=1){
        transition_times[[i]] <- data.frame("transition"=unique(movements)[-stationary],
                                            "departure"=sites_data$timebin[transition_pos],
                                            "arrival"=sites_data$timebin[transition_pos+1],
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
    if(!is.null(tag.dates)){
      transition_times <- plyr::join(transition_times, tag.dates, by=id.col, type="left")
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


    # aggregate sites if required
    if(!is.null(aggregate.by)){
      replacement_vector <- setNames(as.character(sites_table$group), sites_table$site)
      transition_times$transition <- sapply(transition_times$transition, function(transition) {
        for (site in names(replacement_vector)) {
          transition <- gsub(site, replacement_vector[[site]], transition)
        }
        transition
      })
    }
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
    migrations_table$Individuals <- paste0(migrations_table$Individuals, " (",  nids_percent, "%)")

    # calculate additional tag metrics per migration type
    if(!is.null(animals.info)){
      migrating_ids_data <- lapply(migrating_ids, function(x) data.frame("ID"=x))
      migrating_ids_data <- lapply(migrating_ids_data, function(x) plyr::join(x, animals.info, by="ID", type="left"))
      col_classes <- reshape2::melt(lapply(animals.info, class))
      colnames(col_classes) <- c("class", "column")
      col_classes <- col_classes[col_classes$column!="ID",]
      numeric_cols <- col_classes$column[col_classes$class %in% c("numeric", "integer")]
      character_cols <- col_classes$column[col_classes$class %in% c("factor", "character")]
      for(n in numeric_cols){
        mean_values <- reshape2::melt(unlist(lapply(migrating_ids_data, function(x) mean(x[,n], na.rm=T))))$value
        se_values <- reshape2::melt(unlist(lapply(migrating_ids_data, function(x) plotrix::std.error(x[,n], na.rm=T))))$value
        digits <- max(.decimalPlaces(animals.info[,n]), na.rm=T)+1
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
      transition_times <- plyr::join(transition_times, animals.info[,c(id.col, "length")], by=id.col, type="left")
      ages_at_tagging <- TropFishR::VBGF(param=VBGF.params, L=transition_times$length, na.rm=F)
      ages_at_departure <- ages_at_tagging + (transition_times$days_between_depart_tag/365)
      transition_times$lengths_at_departure <- TropFishR::VBGF(param=VBGF.params, t=ages_at_departure, na.rm=F)
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
    if(is.null(aggregate.by)){
      ordered_sites <- levels(data[,sites.col])
    } else {
      ordered_sites <- levels(sites_table$group)
    }
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


  ##############################################################################
  ## Return table or generate plot #############################################
  ##############################################################################

  if(plot==F){
    return(final_table)

  } else{

  #prepare layout
  mat_table <- matrix(rep(1, ncell(final_table)), nrow=nrow(final_table))
  nplots <- length(plot.stats)
  mat_graphs <- matrix(2:(nplots*nrow(final_table)+1), nrow=nrow(final_table), byrow=T)
  layout(mat=cbind(mat_table, mat_graphs))

  # set variables
  table_rows <- nrow(final_table)
  table_cols <- ncol(final_table)
  boxlim <- table_rows*(0.8/13)
  axes_size <- 0.7
  tick_length <- -0.12

  par(mar=c(0,0,0,2), oma=c(3,1,6,1))
  plot(0,0, xlim=c(1, table_cols+2), ylim=c(1, table_rows), type="n", axes=F, main="", ylab="", xlab="")
  mtext("Transition Movements", side=3, outer=T, line=4, cex=0.75, font=2)
  header <- colnames(final_table)
  header[1] <- ""
  values <- unlist(final_table)
  values[is.na(values)] <- "-"
  if(length(id.groups)>1){
    group_names <- names(id.groups)
    group_indexes <- which(values %in% group_names)
    values[group_indexes] <- ""
  }
  first_col <- 1:nrow(final_table)
  ypos <- rep(table_rows:1, table_cols)
  if(table_cols==4){xpos <- rep(c(1, 3.5, 4.5, 5.5), each=table_rows)}
  if(table_cols>4){xpos <- rep(c(1, 3.5, 4.5, 5.5, 6.5), each=table_rows)}
  text(x=xpos[first_col], y=ypos[first_col], labels=values[first_col], cex=0.95, xpd=T, adj=0)
  text(x=xpos[-first_col], y=ypos[-first_col], labels=values[-first_col], cex=0.95, xpd=T, adj=0.5)
  if(length(id.groups)>1){text(x=xpos[group_indexes], y=ypos[group_indexes], labels=group_names, cex=1.05, xpd=T, adj=0, font=2)}
  mtext(text=header, at=unique(xpos), cex=0.7, line=0.7, font=2,  adj=0.5)
  usr <- par("usr")
  box()


  #######################################################################
  ## add barplots #######################################################

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


  # cycle through each row
  current_group <- ifelse(length(id.groups)>1, 0, 1)
  for(i in 1:nrow(final_table)){

    # adjust par settings
    par(mar=c(0.35, 0, 0.35, 1), mgp=c(2.5, 0.3, 0))
    transition <- final_table$Type[i]

    # if divider row, skip to next
    if(transition %in% names(id.groups)){
      current_group <- current_group+1
      for(n in 1:nplots){plot.new()}
      next
    }

    # Reorder transitions and skip to next if intrahabitat movement
    group_data <- graphs_data[[current_group]]
    transitions_order <- unique(final_table$Type)
    transitions_order <- transitions_order[!transitions_order %in% names(id.groups)]
    group_data <- group_data[order(match(group_data$transition, transitions_order)),]
    if(!transition %in% group_data$transition){for(n in 1:nplots){plot.new()}; next}

    # check if this is the last row to plot axes
    plot_axes <- F
    if(which(unique(group_data$transition)==transition)==length(unique(group_data$transition))){
      plot_axes <- T
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
    if(same.scale==T){
      graphs_max<-panel_max
    }else{
      graphs_max <- max(c(depart_hour, arrival_hour, depart_month, arrival_month))
      duration_range <- range(as.numeric(group_data$duration_h), na.rm=T)
    }
    graphs_axis <- pretty(c(0,graphs_max))
    graphs_axis <- unique(as.integer(graphs_axis))


    # depart hour
    if(any(plot.stats=="depart.hour")){
      graphics::barplot(depart_hour, axes=F, col=F, border=F, names.arg=F, ylim=c(0, graphs_max))
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NA)
      b <- graphics::barplot(depart_hour, axes=F, col=adjustcolor("skyblue", alpha.f=0.4),
                             names.arg=F, ylim=c(0, graphs_max), add=T)
      if(plot_axes) axis(1, at=b, labels=paste0(0:23,"h"), cex.axis=axes_size, tck=tick_length, pos=-0.5)
      axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=axes_size, tck=tick_length, las=1)
      box()
    }

    # arrival hour
    if(any(plot.stats=="arrival.hour")){
      graphics::barplot(arrival_hour, axes=F, col=F, border=F, names.arg=F, ylim=c(0, graphs_max))
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NA)
      b <- graphics::barplot(arrival_hour, axes=F, col=adjustcolor("orange", alpha.f=0.4),
                             names.arg=F, ylim=c(0, graphs_max), add=T)
      if(plot_axes) axis(1, at=b, labels=paste0(0:23,"h"),, cex.axis=axes_size, tck=tick_length, pos=-0.5)
      axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=axes_size, tck=tick_length, las=1)
      box()
    }

    # depart month
    if(any(plot.stats=="depart.month")){
      graphics::barplot(depart_month, axes=F, col=F, border=F, names.arg=F, ylim=c(0, graphs_max))
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NA)
      b <- graphics::barplot(depart_month, axes=F, col=adjustcolor("skyblue", alpha.f=0.4),
                             names.arg=F, ylim=c(0, graphs_max), add=T)
      if(plot_axes) axis(1, at=b, labels=month.abb, cex.axis=axes_size, tck=tick_length, pos=-0.5)
      axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=axes_size, tck=tick_length, las=1)
      box()
    }

    # arrival month
    if(any(plot.stats=="arrival.month")){
      graphics::barplot(arrival_month, axes=F, col=F, border=F, names.arg=F, ylim=c(0, graphs_max))
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NA)
      b <- graphics::barplot(arrival_month, axes=F, col=adjustcolor("orange", alpha.f=0.4),
                             names.arg=F, ylim=c(0, graphs_max), add=T)
      if(plot_axes) axis(1, at=b, labels=month.abb, cex.axis=axes_size, tck=tick_length, pos=-0.5)
      axis(2, at=graphs_axis, labels=graphs_axis, cex.axis=axes_size,tck=tick_length, las=1)
      box()
    }

    # duration
    if(any(plot.stats=="duration")){
      boxplot(duration, horizontal=T, axes=F, ylim=duration_range)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col="gray96", border=NA)
      boxplot(duration, horizontal=T, col=adjustcolor("purple", alpha.f=0.3), axes=F, ylim=duration_range, add=T)
      duration_labels <- pretty(duration_range, n=6)
      duration_labels <- duration_labels[duration_labels>=duration_range[1] & duration_labels<=duration_range[2]]
      if(plot_axes) axis(1, at=duration_labels, labels=duration_labels, cex.axis=axes_size, tck=tick_length, pos=0)
      box()
    }
  }

  # add graph headers
  plot.stats <- tools::toTitleCase(gsub(".", " ", plot.stats, fixed=T))
  if(any(plot.stats=="Duration")){
    plot.stats[plot.stats=="Duration"] <- paste0("Duration (", duration_unit, ")")
  }

  par(mfg=c(1,2), xpd=T)
  title(main=plot.stats[1], line=1.25, font=2, cex.main=1, xpd=NA)
  par(mfg=c(1,3))
  title(main=plot.stats[2], line=1.25, font=2, cex.main=1, xpd=NA)
  par(mfg=c(1,4))
  title(main=plot.stats[3], line=1.25, font=2, cex.main=1, xpd=NA)
  par(mfg=c(1,5))
  title(main=plot.stats[4], line=1.25, font=2, cex.main=1, xpd=NA)
  par(mfg=c(1,6))
  title(main=plot.stats[5], line=1.25, font=2, cex.main=1, xpd=NA)
  }

}



#######################################################################################################
#######################################################################################################
#######################################################################################################
