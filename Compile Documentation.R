library(devtools)
library(roxygen2)

devtools::document()
devtools::check()


devtools::use_vignette("my-vignette")
