#######################################################################################################
## Plot summary table for pairwise overlap results ####################################################
#######################################################################################################

#' Generate a visual representation of pairwise overlap statistics.
#'
#' @description Plots a contingency table illustrating the pairwise distribution of overlap
#' scores/significance across individuals. IDs can be sorted by length or any other quantitative metric.
#'
#' @param random.results A list containing null model permutation results. Should be the output from the \code{\link{randomizeOverlap}} function.
#' @param type The type of metric to be plotted. Must be one of "significance" (default) or "mean overlap".
#' @param color.pal A vector of colors defining the palette for the table.
#' - If \code{type = "significance"}, the palette should include three colors
#' corresponding to positive ("+"), negative ("-"), and non-significant ("ns") results.
#' - For \code{type = "mean overlap"}, a gradient of at least 10 colors is recommended.
#' @param sort.by An optional vector to reorder IDs (e.g., based on their size or
#' another quantitative metric). Must have the same length as the IDs in \code{random.results}.
#' @param sort.by.title A string providing a title for the \code{sort.by} metric
#' (e.g., "Size (cm)"). If not supplied, the function attempts to extract the name
#' of the variable passed to \code{sort.by}.
#' @param full.scale Logical. If \code{TRUE} and \code{type = "mean overlap"},
#' scales overlap values from 0% to 100%. Defaults to \code{FALSE}.
#' @param discard.missing Logical. If \code{TRUE}, removes rows and columns with
#' no valid pairwise interactions. Defaults to \code{TRUE}.
#' @param table.cex Numeric. Character expansion factor for the table text.
#' Defaults to \code{0.7}.
#' @param legend.cex Numeric. Character expansion factor for the legend text.
#' Defaults to \code{0.7}.
#'
#' @seealso \code{\link{calculateOverlap}} for overlap computation and
#' \code{\link{randomizeOverlap}} for generating \code{random.results}.
#'
#' @export


plotOverlapTable <- function(random.results,
                             type = "significance",
                             color.pal = NULL,
                             sort.by = NULL,
                             sort.by.title = NULL,
                             full.scale = FALSE,
                             discard.missing = TRUE,
                             table.cex = 0.8,
                             legend.cex = 0.7) {


  ##############################################################################
  # Initial checks #############################################################
  ##############################################################################

  # validate parameters
  errors <- c()
  valid_ids <- TRUE
  if(!inherits(random.results, "list")) errors <- c(errors, "random.results must be a list, as expected from the output of the 'randomizeOverlaps' function.")
  if(!"metric" %in% names(attributes(random.results))) errors <- c(errors, "random.results format not recognized. Please make sure to use the output from the 'randomizeOverlaps' function.")
  else if(!is.null(sort.by)){
    if(length(sort.by) != length(attributes(random.results)$ids)) errors <- c(errors, "'sort.by' variable must have the same length as the number of IDs in random.results.")
    else if(!inherits(sort.by, c("numeric", "integer", "factor", "POSIXct"))) errors <- c(errors, "'sort.by' must contain sortable values (e.g., numeric, integer, factor, or POSIXct).")
  }
  if(!type %in% c("mean overlap", "significance")) errors <- c(errors, "Wrong type argument, please select either 'mean overlap', or 'significance'")
  if(!is.null(color.pal)){
    if(type=="significance"){
      if(length(color.pal)!=3) errors <- c(errors, "Color palette should contain 3 colors when type is set to 'significance'")
    }else if (type=="mean overlap"){
      if(length(color.pal)<10){warning("- Consider increasing the number of colors in color.pal", call. = FALSE)}
    }
  }
  # stop execution if there are errors, displaying all accumulated messages
  if(length(errors)>0){
    stop_message <- sapply(errors, function(x) paste(strwrap(x, width=getOption("width")), collapse="\n"))
    stop_message <- c("\n", paste0("- ", stop_message, collapse="\n"))
    stop(stop_message, call.=FALSE)
  }


  # attempt to retrieve sort.by.title if not explicitly supplied
  if(is.null(sort.by.title) && !is.null(sort.by)){
    # capture the name of the variable supplied to sort.by
    sort.by.expr <- substitute(sort.by)
    # check if the expression is a `$` call
    if(is.call(sort.by.expr) && as.character(sort.by.expr[[1]]) == "$") {
      # extract the second part of the `$` call (the variable name)
      sort.by.title <- as.character(sort.by.expr[[3]])
    } else {
      stop("- Please provide a 'sort.by.title'.", call. = FALSE)
    }
  }


  ##############################################################################
  # Prepare data ###############################################################
  ##############################################################################

  # extract attributes from random.results for further processing
  pairwise_results <- random.results$pairwise_results
  id.groups <- attributes(random.results)$id.groups
  ids <- as.character(attributes(random.results)$ids)
  nids <- length(ids)
  group_comparisons <- attributes(random.results)$group.comparisons

  # create a lookup table for id.groups if provided
  if(!is.null(id.groups)){
    id_lookup <- reshape2::melt(id.groups)
    colnames(id_lookup) <- c("id", "group")
  }

  # set a default color palette if none is supplied
  if(is.null(color.pal)) {
    if(type=="significance"){
      # significance type requires three colors
      color.pal <- c("lightseagreen", "darkred", "grey85")
      color.pal <- adjustcolor(color.pal, alpha.f=0.3)
    }else if (type=="mean overlap"){
      # mean overlap type requires a gradient
      color.pal <- .viridis_pal(100)
      color.pal <- c(adjustcolor(color.pal[1], alpha.f=0.65), adjustcolor(color.pal[2:100], alpha.f=0.7))
    }
  }


  ##############################################################################
  # Create contingency table ###################################################
  ##############################################################################

  # filter out unwanted comparisons based on group_comparisons setting
  if(!is.null(id.groups)){
    discard_rows <- c()
    for(r in 1:nrow(pairwise_results)){
      # map individuals to their respective groups
      group1 <- plyr::mapvalues(pairwise_results$id1[r], from=id_lookup$id, to=as.character(id_lookup$group), warn_missing=FALSE)
      group2 <- plyr::mapvalues(pairwise_results$id2[r], from=id_lookup$id, to=as.character(id_lookup$group), warn_missing=FALSE)
      # apply filtering rules for within-group or between-group comparisons
      if(group_comparisons=="within" & group1!=group2){discard_rows<-c(discard_rows, r)}
      if(group_comparisons=="between" & group1==group2){discard_rows<-c(discard_rows, r)}
    }
    # remove unwanted rows from pairwise_results
    if(length(discard_rows)>0){
      pairwise_results <- pairwise_results[-discard_rows,]
    }
  }

  # create contingency table based on the specified type of analysis
  if(type=="mean overlap"){
    # reshape pairwise results into a contingency table for overlap values
    contingency_table <- reshape2::dcast(pairwise_results, id1 ~ id2, value.var="overlap")
    rownames(contingency_table) <- contingency_table$id1
    contingency_table <- contingency_table[,-1]
    if(discard.missing){
      # discard rows and columns that are entirely NA
      contingency_table <- contingency_table[rowSums(is.na(contingency_table)) != ncol(contingency_table), ]
      contingency_table <- contingency_table[,colSums(is.na(contingency_table)) != nrow(contingency_table)]
    }
    # determine the range of overlap values for further scaling
    overlap_range <- range(pairwise_results, na.rm=TRUE)
  }else if(type=="significance"){
    # reshape pairwise results into a contingency table for significance values
    contingency_table <- reshape2::dcast(pairwise_results, id1 ~ id2, value.var="association")
    rownames(contingency_table) <- contingency_table$id1
    contingency_table <- contingency_table[,-1]
    # replace specific significance values with symbolic representations
    contingency_table[contingency_table=="positive"] <- "+"
    contingency_table[contingency_table=="negative"] <- "-"
    contingency_table[contingency_table=="non-significant"] <- "ns"
    if(discard.missing){
      # discard rows and columns that are entirely NA
      contingency_table <- contingency_table[rowSums(is.na(contingency_table)) != ncol(contingency_table), ]
      contingency_table <- contingency_table[,colSums(is.na(contingency_table)) != nrow(contingency_table)]
    }
  }

  # store row and column IDs for later use
  ids1 <- rownames(contingency_table)
  ids2 <- colnames(contingency_table)

  # add row and column IDs as the first row and column of the table
  contingency_table <- cbind(ids1, contingency_table)
  contingency_table <- rbind(c(NA, ids2), contingency_table)
  rownames(contingency_table) <- NULL
  colnames(contingency_table) <- NULL
  contingency_table[1,1] <- "IDs"

  # sort IDs if specified by the 'sort.by' parameter
  if(!is.null(sort.by)){
    # retrieve sorting variables for row and column IDs
    vars1 <- unlist(sapply(ids1, function(x) sort.by[which(ids==x)]))
    vars2 <- unlist(sapply(ids2, function(x) sort.by[which(ids==x)]))
    # add sorting variables to the table
    contingency_table <- cbind(contingency_table, c(NA, vars1))
    contingency_table <- rbind(contingency_table, c(NA, vars2, NA))
    # reorder the table based on the sorting variables
    cols <- ncol(contingency_table)
    rows <- nrow(contingency_table)
    contingency_table <- contingency_table[,c(1, cols, 2:(cols-1))]
    contingency_table <- contingency_table[c(1, rows, 2:(rows-1)),]
    contingency_table[-c(1:2),] <- contingency_table[-c(1:2),][order(vars1),]
    contingency_table[,-c(1:2)] <- contingency_table[,-c(1:2)][,order(vars2)]
    # set the title for the sorting column
    contingency_table[2,2] <- sort.by.title
    rownames(contingency_table) <- NULL
    colnames(contingency_table) <- NULL
    # format the values based on the decimal places required for sorting variables
    digits <- max(.decimalPlaces(c(vars1, vars2)), na.rm=TRUE)
    contingency_table[-c(1:2), 2] <- sprintf(paste0("%.", digits, "f"), as.numeric(unlist(contingency_table[-c(1:2), 2])))
    contingency_table[2, -c(1:2)] <- sprintf(paste0("%.", digits, "f"), as.numeric(unlist(contingency_table[2, -c(1:2)])))
  }

  # get the number of columns and rows in the contingency table
  n_cols <- ncol(contingency_table)
  n_rows <- nrow(contingency_table)

  # prepare values for visualization
  values <- expand.grid(1:n_rows, 1:n_cols)
  colnames(values) <- c("row", "col")
  values$value <- unlist(apply(values, 1, function(x) contingency_table[x[1], x[2]]))
  values$font <- 1
  values$font[values$row==1 | values$col==1] <- 2
  values$color <- NA
  # set cell colors based on significance or overlap type
  if(type=="significance"){
    values$color[values$value=="+"] <- color.pal[1]
    values$color[values$value=="-"] <- color.pal[2]
    values$color[values$value=="ns"] <- color.pal[3]
    values$color[is.na(values$value) & values$row>1 & values$col>1] <- color.pal[3]
  }else if (type=="mean overlap"){
    if(!is.null(sort.by)){overlap_cells <- which(values$row>2 & values$col>2)
    }else{overlap_cells <- which(values$row>1 & values$col>1)}
    overlaps <- as.numeric(values$value[overlap_cells])
    if(full.scale){
      # rescale the overlaps to fit the full color palette range
      overlaps <- round(.rescale(overlaps, from=c(0,100), to=c(1, length(color.pal)-1)))
    }else{
      # rescale the overlaps based on their actual range
      overlaps <- round(.rescale(overlaps, from=range(overlaps, na.rm=TRUE), to=c(1, length(color.pal)-1)))
    }
    # format the overlap values
    overlaps[is.na(overlaps)] <- 0
    values$color[overlap_cells] <- color.pal[overlaps+1]
    values$value[overlap_cells] <- sprintf("%.1f", as.numeric(values$value[overlap_cells]))
    values$value[values$value=="NA"] <- ""
  }

  # calculate the width of the table for dynamic positioning of the legend
  table_width <- n_cols * 1.1  # Adjust 1.1 based on your character width (might need tuning)

  # calculate inset dynamically: the distance from the table's right edge
  legend_inset_x <- -table_width * 0.05  # 5% of the table width
  legend_inset_y <- 0  # Keep vertical position as it is
  legend_position <- c(legend_inset_x, legend_inset_y)



  ##############################################################################
  # Plot the Table #############################################################
  ##############################################################################


  # set margins for the plot
  par(mar=c(1, 1, 1, 12))

  # create an empty plot with appropriate axis limits for the table dimensions
  plot(0, 0, xlim=c(0.5, n_cols+0.5), ylim=c(n_rows+0.5, 0.5), xaxs="i", yaxs="i",
       type="n", axes=FALSE, main="", ylab="", xlab="")

  # draw the table cells with appropriate colors
  rect(xleft=values$col-0.5 , xright=values$col+0.5 , ybottom=values$row-0.5,
       ytop=values$row+0.5, col=values$color, border=NA)
  text(x=values$col, y=values$row, labels=values$value, font=values$font, cex=table.cex, xpd=TRUE, adj=0.5)

  # draw gridlines if sorting is enabled
  if(!is.null(sort.by)){
    segments(x0=2.5, x1=ncol(contingency_table)+0.5, y0=2.5)
    segments(y0=2.5, y1=nrow(contingency_table)+0.5, x0=2.5)
  }else{
    segments(x0=1.5, x1=ncol(contingency_table)+0.5, y0=1.5)
    segments(y0=1.5, y1=nrow(contingency_table)+0.5, x0=1.5)
  }

  # add bounding box around the table
  segments(x0=0.5, x1=n_cols+0.5, y0=0.5, xpd=TRUE)
  segments(x0=0.5, x1=n_cols+0.5, y0=n_rows+0.5, xpd=TRUE)
  segments(y0=0.5, y1=n_rows+0.5, x0=0.5, xpd=TRUE)
  segments(y0=0.5, y1=n_rows+0.5, x0=n_cols+0.5, xpd=TRUE)

  # add legend based on the significance or overlap type
  if(type=="significance"){
    plus <- "overlap > random\n(joint resource utilization)"
    minus <- "overlap < random\n(spatiotemporal segregation)"
    ns <- "non-significant overlap"
    # calculate legend position
    legend_x <- par("usr")[2] + diff(par("usr")[1:2]) * 0.02
    legend_y <- mean(par("usr")[3:4])
    .legend(x=legend_x, y=legend_y, legend=c(plus, minus, ns), box.lwd=0.1,
            fill=color.pal, box.cex=c(1.5, 1), bty="n", cex=legend.cex, y.intersp=1.6, xpd=TRUE)
  }else if(type=="mean overlap"){
    if(full.scale) {
      overlap_range <- c(0,100)
      overlap_labs <- pretty(0:100, min.n=4)
    }else{
      overlap_labs <- pretty(overlap_range, min.n=4)
      overlap_labs <- overlap_labs[overlap_labs>=min(overlap_range) & overlap_labs<=max(overlap_range)]
    }
    overlap_digits <- max(.decimalPlaces(overlap_labs), na.rm=TRUE)
    .colorlegend(col=color.pal, zlim=overlap_range, zval=overlap_labs,
                            posx=c(0.88, 0.895), posy=c(0.2, 0.65), main="Overlap (%)",
                            main.cex=legend.cex+0.1, digit=overlap_digits, main.adj=0, cex=legend.cex)
  }

}



##################################################################################################
##################################################################################################
##################################################################################################
