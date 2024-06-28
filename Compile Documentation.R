library(devtools)
library(roxygen2)

devtools::document()
rcheck <- devtools::check()

cat(rcheck$warnings[1])
cat(rcheck$notes)

usethis::use_gpl3_license()
usethis::use_vignette("moby-vignette")
#devtools::build_vignettes()

#use_readme_md(open = rlang::is_interactive())

#usethis::use_github_action_check_standard()

remotes::install_github("eddelbuettel/dang")
dang::checkPackageAsciiCode(dir = ".")

stringi::stri_escape_unicode("º")
stringi::stri_escape_unicode("±")


#usethis::use_build_ignore("Commands")
#usethis::use_build_ignore("Compile Documentation.R")
#usethis::use_build_ignore("moby-tutorial.R")
