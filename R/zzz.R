
mobyEnv <- NULL

# this will run automatically when package is loaded
.onLoad <- function(libname, pkgname){

  version <- utils::packageVersion("moby")
  version_line <- sprintf("** Version: %-47s **", version)

  packageStartupMessage("
  **************************************************************
  *******************    Welcome to moby    ********************
  **************************************************************
  ** Biotelemetry analyses and data visualization             **
  ", version_line, "
  ** To get started, check out the documentation:             **
  ** - vignette('moby') for a comprehensive guide             **
  ** - help(package='moby') for a list of available functions **
  ** Please cite this package if you use it in publications:  **
  ** - citation('moby') for details                           **
  **************************************************************
  **************************************************************
  ")

  # create new environment to store needed variables
  mobyEnv <<- new.env()

  # store default arguments
  mobyEnv$defaults <- list(
    id.col = "ID",
    datetime.col = "datetime",
    timebin.col = "timebin",
    tagging.dates = NULL,
    tag.durations = NULL
  )

}
