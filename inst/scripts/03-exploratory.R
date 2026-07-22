## -----------------------------------------------------------------------------
#| message: false
library(moby)

# the cleaned, analysis-ready checkpoint (a mobyData object: see ?rays)
data(rays)

# the two-species grouping is stored in the object's metadata
id_groups <- mobyMeta(rays)$id.groups
names(id_groups)


## -----------------------------------------------------------------------------
str(id_groups)


## -----------------------------------------------------------------------------
#| fig-height: 6
plotAbacus(rays, id.groups = id_groups)


## -----------------------------------------------------------------------------
#| fig-height: 6
plotStationStats(rays, type = c("detections", "individuals"), id.groups = id_groups)


## -----------------------------------------------------------------------------
# representative coordinates of the array (centre of the MPA)
mpa_coords <- matrix(c(mean(rays$lon), mean(rays$lat)), nrow = 1)

rays$diel <- getDielPhase(rays$datetime, coords = mpa_coords, phases = 2)
table(rays$diel)


## -----------------------------------------------------------------------------
#| fig-height: 6
# chronograms operate on time-binned detections; bin to the hour first
rays$timebin <- getTimeBins(rays$datetime, interval = "60 mins")

# one panel per species (split.by), points coloured by receiver
plotChronogram(rays, coords = mpa_coords, split.by = "species",
               color.by = "station")


## -----------------------------------------------------------------------------
# hourly detection counts per individual: a long ID x time-bin x count table
binned <- aggregate(list(detections = rep(1L, nrow(rays))),
                    by = list(ID = rays$ID, timebin = rays$timebin), FUN = sum)


## -----------------------------------------------------------------------------
#| fig-height: 6.5
plotPeriodogram(binned, id.col = "ID", timebin.col = "timebin", id.groups = id_groups)


## -----------------------------------------------------------------------------
#| fig-height: 6
plotScalogram(binned[binned$ID %in% id_groups[[1]], ],
              variable = "detections", id.col = "ID", timebin.col = "timebin", ncol = 2)


## -----------------------------------------------------------------------------
monitoring_end <- as.POSIXct("2023-07-05", tz = "UTC")

residency <- calculateResidency(rays, last.monitoring.date = monitoring_end)
head(residency)


## -----------------------------------------------------------------------------
summaryTable(rays, last.monitoring.date = monitoring_end, id.groups = id_groups)

