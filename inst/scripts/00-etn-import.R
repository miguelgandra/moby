## -----------------------------------------------------------------------------
# run once
pak::pak("inbo/etn")
install.packages("usethis")


## -----------------------------------------------------------------------------
usethis::edit_r_environ()


## -----------------------------------------------------------------------------
library(etn)
etn_projects <- etn::get_animal_projects()
View(etn_projects)


## -----------------------------------------------------------------------------
id_metadata <- etn::get_animals(animal_project_code = "Inforbiomares")

id_deployments <- etn::get_acoustic_deployments(acoustic_project_code = "Inforbiomares")

all_detections <- etn::get_acoustic_detections(
  animal_project_code = "Inforbiomares",
  scientific_name     = c("Dasyatis pastinaca", "Raja clavata")
)
all_detections <- as.data.frame(all_detections)

# record the download time for reproducibility
attr(all_detections, "download.date") <- Sys.time()


## -----------------------------------------------------------------------------
write.csv2(all_detections, "INFORBIOMARES_detections.csv", row.names = FALSE)

