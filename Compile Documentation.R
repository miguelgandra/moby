library(devtools)
library(roxygen2)

devtools::document()
rcheck <- devtools::check()


usethis::use_gpl3_license()
#devtools::use_vignette("my-vignette")

#use_readme_md(open = rlang::is_interactive())



remotes::install_github("eddelbuettel/dang")
dang::checkPackageAsciiCode(dir = ".")


