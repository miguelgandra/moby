library(devtools)
library(roxygen2)

devtools::document()
devtools::check()


devtools::use_vignette("my-vignette")


use_readme_md(open = rlang::is_interactive())



remotes::install_github("eddelbuettel/dang")
dang::checkPackageAsciiCode(dir = ".")
