% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMigrations.R
\name{plotMigrations}
\alias{plotMigrations}
\title{Plot migrations network}
\usage{
plotMigrations(
  data,
  id.col = getDefaults("ID"),
  lon.col = getDefaults("lon"),
  lat.col = getDefaults("lat"),
  datetime.col = getDefaults("datetime"),
  spatial.col,
  id.groups = NULL,
  id.metadata = NULL,
  land.shape = NULL,
  land.color = "gray50",
  epsg.code = getDefaults("epsg"),
  background.layer = NULL,
  background.pal = "#F3F7F7",
  title.cex = 1,
  color.nodes.by = "detection",
  nodes.color = c("darkblue", "black"),
  nodes.alpha = 0.8,
  nodes.size = c(0.05, 0.1),
  nodes.label.wrap = FALSE,
  nodes.label.cex = 0.5,
  nodes.label.color = "white",
  repel.nodes = TRUE,
  repel.buffer = 1.1,
  edge.type = "movements",
  edge.color = "darkblue",
  edge.curved = 0.5,
  edge.width = c(0.4, 3.5),
  edge.arrow.size = 0.5,
  edge.arrow.width = 1.5,
  edge.label.cex = 0.6,
  edge.label.color = "black",
  edge.label.font = 1,
  scale.km = NULL,
  scale.pos = "bottom",
  scale.inset = c(0, 0.05),
  scale.height = 1.5,
  scale.cex = 0.6,
  extent.factor = 1.1,
  cols = 1,
  ...
)
}
\arguments{
\item{data}{A data frame containing animal detections.}

\item{id.col}{Name of the column containing animal IDs. Defaults to 'ID'.}

\item{lon.col}{Name of the column containing longitude values (either projected or unprojected). Defaults to 'lon'.}

\item{lat.col}{Name of the column containing latitude values (either projected or unprojected). Defaults to 'lat'.}

\item{datetime.col}{Name of the column containing datetimes in POSIXct format. Defaults to 'datetime'.}

\item{spatial.col}{Name of the column in both \code{data} and \code{id.metadata} (if supplied)
that contains spatial information for calculating transitions. This can include receiver IDs
in case all movements/transitions are of interest, or contain (for example) location/habitat
classes for broader-scale analyses.}

\item{id.groups}{Optional. A list containing ID groups (e.g. different species or size-classes).
If supplied, a separate map will be plotted for each group.}

\item{id.metadata}{Optional. A data frame containing information on the tagging location of each animal.
Used to indicate the number of tagged individuals at each site.}

\item{land.shape}{Optional. A shape file containing coastlines.}

\item{land.color}{Optional. Color of land areas. Defaults to "gray50".}

\item{epsg.code}{An integer representing the EPSG code of the coordinate reference system (CRS). This parameter
serves two purposes: (1) to specify the projection system used when coordinates are already projected, or
(2) to define the target projection system when coordinates are supplied in geographic format. Note that
this must always refer to a \strong{projected CRS} (in meters).}

\item{background.layer}{Optional. A projected raster containing a variable to be displayed in the background
(e.g. bathymetry, temperature, etc.).}

\item{background.pal}{Color palette for the background layer. If a background layer is not supplied,
the first color is used as background. Defaults to "#F3F7F7".}

\item{title.cex}{Size of the plot title. Defaults to 1.}

\item{color.nodes.by}{Variable used to color nodes (locations). It can be set to any variable in the dataset.
Alternatively, it can be set to 'detections' (default) to distinguish between nodes with and without detections,
or to 'group' to assign different colors to different ID groups (one color per map). Defaults to "detection".}

\item{nodes.color}{Color(s) for the nodes. Defaults to c("darkblue", "black").}

\item{nodes.alpha}{Numeric. Transparency level for the nodes (0 to 1): 0=fully transparent, 1=fully opaque. Defaults to 0.8.}

\item{nodes.size}{A numeric vector of length 2. Represents the desired min and max vertex sizes relative
to the x-axis in terms of percentage (see details). See \code{\link[netdiffuseR]{rescale_vertex_igraph}}. Defaults to c(0.04, 0.08).}

\item{nodes.label.wrap}{Logical. If TRUE, splits node labels into multiple lines (to better fit within the nodes). Defaults to FALSE.}

\item{nodes.label.cex}{The font size for vertex labels. Defaults to 0.5.}

\item{nodes.label.color}{The font color for vertex labels. Defaults to white.}

\item{repel.nodes}{Logical. If TRUE, nodes are plotted using a repulsion algorithm to avoid overlap. Defaults to TRUE}

\item{repel.buffer}{Controls the amount of space between nodes if repel.nodes is set to TRUE. Defaults to 1.1.}

\item{edge.type}{A character string indicating the metric to be used to calculate the network edges.
It can be either "movements" to represent the number of movements between sites or "individuals"
to represent the number of individual animals moving between sites. Defaults to "movements".}

\item{edge.color}{Color(s) for the edges. Defaults to "darkblue".}

\item{edge.curved}{Specifies whether to draw curved edges, or not. This can be a logical or a numeric vector or scalar.
A numeric value specifies the curvature of the edge; zero curvature means straight edges, negative values mean the edge bends
clockwise, positive values the opposite. TRUE means curvature 0.5, FALSE means curvature zero. Defaults to 0.5.}

\item{edge.width}{A numeric vector of length 2. Represents the desired min and max
width/thickness of the edges. Defaults to c(0.4, 3.5).}

\item{edge.arrow.size}{The size of the arrows. Defaults to 0.5.}

\item{edge.arrow.width}{The width of the arrows. Defaults to 1.5.}

\item{edge.label.cex}{The font size for edge labels. Defaults to 0.6.}

\item{edge.label.color}{The color of the edge labels. Defaults to black.}

\item{edge.label.font}{The font for the edge labels. It is interpreted the same way
as the font graphical parameter: 1 is plain text, 2 is bold face, 3 is italic,
4 is bold and italic and 5 specifies the symbol font. Defaults to 1.}

\item{scale.km}{Distance covered by the scale bar, in kilometers. If NULL, it is
automatically defined as 20\% of the plot region.}

\item{scale.pos}{Position of the map scale, specified by keyword.
See \code{\link[grDevices]{xy.coords}}. Defaults to "bottom".}

\item{scale.inset}{Controls how far the scale bar is placed from the plot edges.
Can be specified as a single value (applied to both x and y directions) or as a two-element vector
(with x and y inset values specified separately).
\itemize{
\item \strong{Positive values} move the scale bar \emph{inward}, towards the center of the plot.
\item \strong{Negative values} move the scale bar \emph{outward}, potentially placing it outside the visible plot area.
}

Defaults to c(0, 0.05).}

\item{scale.height}{Controls the thickness of the scale bar. Defaults to 1.5.}

\item{scale.cex}{Size of the scale bar values. Defaults to 0.6.}

\item{extent.factor}{Numeric. Factor by which to adjust the extent of the plotting region,
defined based on the bounding box around animal positions/detections. A value of 1 keeps
the original bounding box, values greater than 1 increase the extent, and values less
than 1 decrease it. Defaults to 1.1 (10\% increase).}

\item{cols}{Number of columns in the plot panel layout (used in the 'mfrow' argument). Defaults to 1.}

\item{...}{Arguments passed to the \code{\link[igraph]{plot.igraph}} function.
See \code{\link[igraph]{igraph.plotting}} for the complete parameters list.}
}
\value{
A plot displaying the migrations network. Nodes correspond to different sites,
with coordinates reflecting their approximate location in geographic space.
Arrows illustrate movements between sites. The size of each node indicates the number of
individuals detected at that site, and labels on the arrows show either the number of
transition movements or the the number of individuals transiting between sites.
Optionally, the number of tagged individuals can be indicated in parentheses
within the nodes of the corresponding sites.
}
\description{
Plots a spatial network representing animal transitions/movements between receivers/sites.
If tagging information is supplied, the number of tagged individuals is indicated for each location/site, and
the release site is included in the calculated transitions.
}
\examples{
\dontrun{
pdf("./spatial-network.pdf", width=6, height=8)
plotMigrations(data=data_filtered,
               id.groups=species_ids,
               spatial.col="site",
               land.shape=coastline,
               edge.type="individuals",
               edge.color=c("darkblue", "darkred"),
               edge.width=c(0.2, 4),
               color.nodes.by="group",
               nodes.color=c("darkblue", "darkred"),
               nodes.size=c(0.05,0.1),
               nodes.label.wrap=T,
               repel.nodes=T)
dev.off()
}
}
