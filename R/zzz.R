
mobyEnv <- NULL

# this will run automatically when package is loaded
.onLoad <- function(libname, pkgname){

  # create new environment to store needed variables
  mobyEnv <<- new.env()

  # store default arguments
  mobyEnv$defaults <- list(
    id.col = "ID",
    datetime.col = "datetime",
    timebin.col = "timebin",
    station.col = "station",
    lon.col = "lon",
    lat.col = "lat",
    epsg.code = NULL,
    tagging.dates = NULL
  )
}


# this will run automatically when package is attached
.onAttach <- function(libname, pkgname){

  version <- utils::packageVersion("moby")
  version_line <- sprintf("** Version: %-47s **", version)

  packageStartupMessage("
  **************************************************************
  *******************    Welcome to moby    ********************
  **************************************************************
  ** Biotelemetry analyses and data visualization             **
  ", version_line, "
  ** To get started, check out the documentation:             **
  ** - help(package='moby') for a list of available functions **
  ** Please cite this package if you use it in publications:  **
  ** - citation('moby') for details                           **
  **************************************************************
  **************************************************************
  ")
}
