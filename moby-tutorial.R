###################################################################################################################
###################################################################################################################
## Miguel Gandra || CCMAR || m3gandra@gmail.com || Apr 2024 #######################################################
###################################################################################################################
###################################################################################################################

# moby: A R package for streamlined biotelemetry analyses and data visualization.

# Sample script to process and analyse passive acoustic telemetry data using the 'moby' package.
#
# ---> filter out spurious detections
# ---> assign a time bin to each detection, according with a given interval duration (e.g., 30 mins)
# ---> calculate COAS - Centres of activity
# ---> assign temporal classes (e.g., diel phase, season, reproductive status, etc.)
# ---> create summary table, including residency index
# ---> generate abacus plot (animal detections over time)
# ---> generate level plots (detections/individuals per hour and date)


#######################################################################################################
# Automatically install required libraries  ###########################################################
#######################################################################################################

if(!require(moby)){install.packages("moby"); library(moby)}
if(!require(lubridate)){install.packages("lubridate"); library(lubridate)}
if(!require(raster)){install.packages("raster"); library(raster)}
if(!require(sp)){install.packages("sp"); library(sp)}
if(!require(pals)){install.packages("pals"); library(pals)}


#######################################################################################################
# Load Data  ##########################################################################################
#######################################################################################################

# set epsg code for the study region
epsg_code <- CRS("+init=epsg:3063")

# import raw detections files
data_raw <- read.csv2("./data/acoustic_detections.csv", header=T)
selected_cols <- c("date_time", "receiver_id", "tag_id", "station_name", "deploy_latitude", "deploy_longitude", "network_project_code")
data_raw <- data_raw[,colnames(data_raw) %in% selected_cols]
colnames(data_raw) <- c("datetime", "receiver", "network", "transmitter", "station", "latitude", "longitude")
data_raw$datetime <- as.POSIXct(data_raw$datetime, "%Y-%m-%d %H:%M:%S", tz="UTC")
data_raw <- data_raw[order(data_raw$transmitter, data_raw$datetime),]

# import animals metadata
data_tags <- read.csv2("./data/tagged_animals.csv", header=T)
colnames(data_tags) <- c("ID", "length", "transmitter", "tagging_date", "dead_date", "tagging_location")
data_tags$tagging_date <- as.POSIXct(data_tags$tagging_date, "%Y-%m-%d", tz="UTC")

# import coastline shapefile
coastline <- raster::shapefile("./layers/coastline/GSHHS_shp/GSHHS_h_L1.shp")
eastern_atlantic <- matrix(c(-15, 8, 32, 55), ncol=2, byrow=T)
coastline <- crop(coastline, extent(eastern_atlantic))
coastline <- sp::spTransform(coastline, epsg_code)


#######################################################################################################
# Remove spurious detections ##########################################################################
#######################################################################################################

filter_results <- filterDetections(data=data_raw, tagging.dates=data_tags$tagging_date, min.detections=0,
                                   land.shape=coastline, hours.threshold=24, acoustic.range=800,
                                   max.speed=3, grid.resolution=1000, mov.directions=8)
data_filtered <- filter_results$data


#######################################################################################################
# Assign a time bin to each entry #####################################################################
#######################################################################################################

# specify time bins interval (in minutes)
interval <- "30 mins"

data_filtered$timebin <- as.POSIXct(NA, tz="UTC")
data_filtered$timebin <- lubridate::floor_date(data_filtered$datetime, interval)

# convert timebins to POSIXct
data_filtered$timebin <-as.POSIXct(data_filtered$timebin, format= "%Y-%m-%d %H:%M:%S", tz="UTC")

# reorder columns
n_cols <- ncol(data_filtered)
data_filtered <- data_filtered[,c(n_cols, 1:(n_cols-1))]

# save curated detections
write.csv2(data_filtered, file="./data/curated_detections.csv", row.names=F)


#######################################################################################################
# Restructure data and assign temporal variables ######################################################
#######################################################################################################

# create data matrix from binned data
data_table <- createWideTable(data=data_filtered, id.col="ID", start.dates=data_tags$tagging_date, value.col="station")

# define study region coordinates (used to retrieve sunrise and sunset times)
study_coords <- matrix(colMeans(coordinates_wsg84@coords), ncol=2)

# assign diel phase to each detection
data_table$timeofday <- getDielPhase(data_table$timebin, study_coords, phases=4)
data_filtered$timeofday <- getDielPhase(data_filtered$timebin, study_coords, phases=4)

# assign annual season to each detection
data_table$season <- getSeason(data_table$timebin)
data_filtered$season <- getSeason(data_filtered$timebin)

# assign reproductive season to each detection
data_table$reprod_season <- getReprodPeriod(data_table$timebin, spawning.start="Apr", spawning.end="Aug", format="%b")
data_filtered$reprod_season <- getReprodPeriod(data_filtered$timebin, spawning.start="Apr", spawning.end="Aug", format="%b")

# assign month and hour
data_table$month <- strftime(data_table$timebin, "%b", tz="UTC")
data_filtered$month <- strftime(data_filtered$timebin, "%b", tz="UTC")
month_ord <- order(match(unique(data_coas$month), month.abb))
data_table$month <- factor(data_table$month, levels=unique(data_table$month)[month_ord])
data_filtered$month <- factor(data_filtered$month, levels=unique(data_filtered$month)[month_ord])
data_table$hour <- strftime(data_table$timebin, "%H", tz="UTC")
data_filtered$hour <- strftime(data_filtered$timebin, "%H", tz="UTC")


#######################################################################################################
# Create summary table ################################################################################
#######################################################################################################

animal_info <- data_tags[,c("ID", "transmitter", "length", "sex", "size")]
colnames(animal_info) <- c("ID", "Transmitter", "Fork Length (cm)", "Sex", "Tag")
summary_table <- summaryTable(data=data_filtered, tagging.dates=data_tags$tagging_date, tags.info=animal_info, error.stat="se")

#save detection periods table
write.csv2(summary_table, "./summary_table.csv", row.names=F, fileEncoding="Windows-1252")


#######################################################################################################
# Generate abacus plot ################################################################################
#######################################################################################################

nstations <- nlevels(data_bins$receiver)
color_pal <- pals::ocean.haline(nstations)

pdf("./abacus_plot.pdf", height=7.5, width=10,  useDingbats=F)
plotAbacus(binned.data=data_filtered, tagging.dates=tagging_dates, color.by="receiver",
           color.pal=color_pal, discard.missing=T, season.shade=T, date.start=3,
           top.mural="%Y", date.format="%b", date.interval=6, pt.cex=1.25, highlight.isolated=T)
dev.off()


#######################################################################################################
# Plot level plots ####################################################################################
#######################################################################################################

pdf("./level_plot.pdf", height=8, width=16,  useDingbats=F)
par(mfrow=c(2,1))
plotLevels(data_filtered, variables=c("individuals", "co-occurrences"), style="points", date.format="%b/%y",
           color.pal=habitats_pal, date.interval=6, date.start=3, sunriset.coords=study_coords, lunar.info=F,
           color.by="receiver", polygons="season", background.col="grey98", grid=F,
           highlight.isolated=T, uniformize.scale=T, uniformize.dates=T, pt.cex=c(0.75,4))
dev.off()



###################################################################################################################
###################################################################################################################
###################################################################################################################
