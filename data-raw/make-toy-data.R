###############################################################################################
## Generate the toy teaching dataset shipped with moby ########################################
###############################################################################################
##
## Two demersal elasmobranchs (Raja clavata, Dasyatis pastinaca) tracked with an acoustic
## array inside a small coastal MPA. The data are FULLY SYNTHETIC - coordinates, dates and
## detection patterns are simulated, so there are no sensitivity/ethics concerns with the
## locations of real protected species. The simulation is only meant to be ecologically
## plausible enough to exercise every step of the moby workflow in the tutorials.
##
## Run with:  source("data-raw/make-toy-data.R")
## (re-generates the objects in data/ and the raw CSVs in inst/extdata/).
###############################################################################################

set.seed(20240601)

## ---- array (6 receivers in a small bay; SYNTHETIC coordinates) ----------------------------
stations <- data.frame(
  station   = sprintf("ST%02d", 1:6),
  lon       = -9.00 + c(-0.020, -0.008,  0.004,  0.015,  0.028,  0.010),
  lat       =  38.45 + c( 0.004,  0.014,  0.006,  0.016,  0.002, -0.008),
  receiver  = sprintf("VR2W-%04d", 1001:1006),
  stringsAsFactors = FALSE
)

## ---- tagged animals (4 per species) -------------------------------------------------------
species_levels <- c("Raja clavata", "Dasyatis pastinaca")
animals <- data.frame(
  ID       = c(sprintf("R%02d", 1:4), sprintf("D%02d", 1:4)),
  species  = rep(species_levels, each = 4),
  stringsAsFactors = FALSE
)
# tagging spread over the first two weeks; each animal released at a "home" station
animals$tagging_date <- as.POSIXct("2023-04-01", tz = "UTC") +
  sample(0:13, nrow(animals), replace = TRUE) * 86400
animals$tagging_station <- stations$station[sample(seq_len(nrow(stations)), nrow(animals), replace = TRUE)]
animals$transmitter <- sprintf("A69-1602-%05d", 30001:(30000 + nrow(animals)))

## ---- simulate detections ------------------------------------------------------------------
# each animal has a preferred subset of stations (Raja = central, Dasyatis = peripheral),
# visits them in presence "bouts", with a diel modulation (more nocturnal detections).
study_end <- as.POSIXct("2023-06-30 23:00", tz = "UTC")

simulate_individual <- function(id, sp, tag_date, home) {
  # station preference weights
  centrality <- abs(seq_len(6) - 3.5)               # ST03/ST04 most central
  pref <- if (sp == "Raja clavata") 1 / (centrality + 0.5) else (centrality + 0.5)
  pref[stations$station == home] <- max(pref) * 1.6  # bias toward release site
  pref <- pref / sum(pref)

  n_bouts <- sample(25:45, 1)
  bout_starts <- sort(tag_date + runif(n_bouts, 0, as.numeric(difftime(study_end, tag_date, units = "secs"))))
  recs <- list()
  for (b in seq_along(bout_starts)) {
    st  <- sample(stations$station, 1, prob = pref)
    n   <- rpois(1, 8) + 1                            # detections within the bout
    t0  <- bout_starts[b]
    times <- t0 + cumsum(rexp(n, rate = 1 / 1200))    # ~20 min spacing
    times <- times[times <= study_end]
    if (!length(times)) next
    # diel thinning: keep night detections more often (hours 20-06)
    hr <- as.integer(format(times, "%H", tz = "UTC"))
    keep <- runif(length(times)) < ifelse(hr >= 20 | hr < 6, 0.95, 0.55)
    times <- times[keep]
    if (!length(times)) next
    recs[[b]] <- data.frame(ID = id, datetime = times, station = st, stringsAsFactors = FALSE)
  }
  do.call(rbind, recs)
}

det <- do.call(rbind, Map(simulate_individual,
                          animals$ID, animals$species, animals$tagging_date, animals$tagging_station))
det <- det[order(det$ID, det$datetime), ]
rownames(det) <- NULL

# attach station coordinates + species + transmitter
det <- merge(det, stations[, c("station", "lon", "lat", "receiver")], by = "station", all.x = TRUE)
det <- merge(det, animals[, c("ID", "species", "transmitter")], by = "ID", all.x = TRUE)
det <- det[order(det$ID, det$datetime), c("ID", "datetime", "station", "lon", "lat", "species", "receiver", "transmitter")]
rownames(det) <- NULL

###############################################################################################
## Assemble the exported objects ##############################################################
###############################################################################################

## raw detections (as a plain data.frame, mimicking a harmonised import) ----------------------
rays_detections <- det

## tag metadata -------------------------------------------------------------------------------
rays_tags <- animals[, c("ID", "transmitter", "species", "tagging_date", "tagging_station")]
rownames(rays_tags) <- NULL

## receiver-deployment log, in the RAW European Tracking Network (ETN) export format ----------
## This is the layout users most commonly receive; importDeployments(source = "etn") harmonises it
## to moby's canonical schema. It is written to inst/extdata/ so the import tutorial can demonstrate
## that step on a realistic file (all values SYNTHETIC).
n_st <- nrow(stations); na_chr <- rep(NA_character_, n_st); na_num <- rep(NA_real_, n_st); na_lgl <- rep(NA, n_st)
rays_deployments_raw <- data.frame(
  deployment_id              = 5000L + seq_len(n_st),
  receiver_id                = stations$receiver,
  acoustic_project_code      = "RAYS-MPA",
  station_name               = stations$station,
  station_description        = paste("Station", stations$station),
  station_manager            = na_chr,
  deploy_date_time           = "2023-03-25 09:00:00",
  deploy_latitude            = stations$lat,
  deploy_longitude           = stations$lon,
  intended_latitude          = na_num,
  intended_longitude         = na_num,
  mooring_type               = "bottom-mooring",
  bottom_depth               = na_num,
  riser_length               = na_num,
  deploy_depth               = c(12, 18, 15, 22, 9, 14),
  battery_installation_date  = "2023-03-25",
  battery_estimated_end_date = "2024-03-25",
  activation_date_time       = na_chr,
  recover_date_time          = "2023-07-05 15:00:00",
  recover_latitude           = stations$lat,
  recover_longitude          = stations$lon,
  download_date_time         = "2023-07-05 15:00:00",
  download_file_name         = sprintf("%s_20230705_1.vrl", gsub("-", "_", stations$receiver)),
  valid_data_until_date_time = na_chr,
  sync_date_time             = na_chr,
  time_drift                 = na_num,
  # instrument-config / acoustic-release fields: present in the ETN export layout but empty for
  # these plain VR2W bottom-mooring receivers (matches how a real export fills them - all NA).
  ar_battery_installation_date  = na_lgl,
  ar_confirm                    = na_lgl,
  transmit_profile              = na_lgl,
  transmit_power_output         = na_lgl,
  log_temperature_stats_period  = na_lgl,
  log_temperature_sample_period = na_lgl,
  log_tilt_sample_period        = na_lgl,
  log_noise_stats_period        = na_lgl,
  log_noise_sample_period       = na_lgl,
  log_depth_stats_period        = na_lgl,
  log_depth_sample_period       = na_lgl,
  comments                   = na_chr,
  stringsAsFactors = FALSE
)

## cleaned, analysis-ready checkpoint (a mobyData) --------------------------------------------
# built with moby so later tutorials can load a ready object and stay independent of T01/T02.
pkgload::load_all(".", quiet = TRUE)

# harmonise the raw ETN deployment log to moby's canonical schema (exactly what the import tutorial does)
rays_deployments <- importDeployments(rays_deployments_raw, source = "etn")

id_groups <- split(rays_tags$ID, rays_tags$species)[species_levels]   # named by species

rays <- as_moby(
  rays_detections,
  id.col        = "ID",
  datetime.col  = "datetime",
  station.col   = "station",
  lon.col       = "lon",
  lat.col       = "lat",
  epsg.code     = 32629,   # UTM 29N (coords ~ -9 lon, off Portugal) -> projected analyses (KUDs, distances)
  tagging.dates = setNames(rays_tags$tagging_date, rays_tags$ID),
  id.groups     = id_groups
)
rays$timebin <- getTimeBins(rays$datetime, interval = "1 hour")   # adds the canonical 'timebin' column

###############################################################################################
## Save ########################################################################################
###############################################################################################

save_obj <- function(obj, name) {
  assign(name, obj)
  save(list = name, file = file.path("data", paste0(name, ".rda")), compress = "xz")
}
save_obj(rays_detections, "rays_detections")
save_obj(rays_tags, "rays_tags")
save_obj(rays_deployments, "rays_deployments")
save_obj(rays, "rays")

# raw CSVs for the import tutorial:
#  - detections in a 'generic' source format (non-canonical column names, hand-mapped)
csv <- rays_detections
names(csv) <- c("animal_id", "timestamp", "station_name", "deploy_longitude",
                "deploy_latitude", "scientific_name", "receiver_id", "transmitter")
write.csv(csv, "inst/extdata/rays_detections.csv", row.names = FALSE)
#  - the receiver-deployment log in raw ETN export format (harmonised via source = "etn")
write.csv(rays_deployments_raw, "inst/extdata/rays_deployments.csv", row.names = FALSE)

cat("Toy dataset generated:\n",
    " detections:", nrow(rays_detections), "rows |",
    nlevels(factor(rays_detections$ID)), "animals |",
    length(unique(rays_detections$station)), "stations\n")
