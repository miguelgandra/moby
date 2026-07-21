# moby tutorial scripts

Runnable R scripts for the moby tutorial series, intended for line-by-line use in RStudio
(e.g. during teaching sessions). They mirror the rendered tutorials on the package website:

<https://miguelgandra.github.io/moby/>

Each script is **generated** from the corresponding Quarto article in `vignettes/articles/`
(via `data-raw/build-tutorials.R`) — edit the `.qmd`, not the `.R`.

| Script | Topic |
|--------|-------|
| `00-etn-import.R`           | Getting real data from ETN (optional appendix) |
| `01-import.R`               | Raw detections → a `mobyData` object |
| `02-quality-control.R`      | QC, deployment matching & filtering |
| `03-exploratory.R`          | Exploratory analysis, residency, summaries |
| `04-space-use.R`            | Movement metrics & home ranges (KUDs) |
| `05-movement-networks.R`    | Movement (spatial) networks |
| `06-association-networks.R` | Association (social) networks |

All tutorials run on the bundled toy dataset (`data(rays)`), so they work offline with no
credentials. See `?rays` for details.
