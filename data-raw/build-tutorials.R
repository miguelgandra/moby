###############################################################################################
## Maintainer helper: generate runnable .R tutorial scripts from the .qmd articles ###########
###############################################################################################
##
## Single source of truth = the Quarto articles in vignettes/articles/. This script extracts
## their R code into plain, classroom-friendly scripts in inst/scripts/, so the two never
## drift apart. Re-run after editing any article.
##
##   source("data-raw/build-tutorials.R")
##
## Requires a recent {knitr} (handles the `#|` chunk options used by Quarto).
###############################################################################################

src_dir <- "vignettes/articles"
out_dir <- "inst/scripts"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

qmd <- sort(list.files(src_dir, pattern = "\\.qmd$", full.names = TRUE))

for (f in qmd) {
  out <- file.path(out_dir, sub("\\.qmd$", ".R", basename(f)))
  knitr::purl(f, output = out, quiet = TRUE, documentation = 1L)
  cat("wrote", out, "\n")
}
