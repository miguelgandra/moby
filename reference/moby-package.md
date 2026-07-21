# moby: Streamlined Acoustic Telemetry Analyses and Data Visualization

moby provides a comprehensive suite of functions and methods tailored
for processing, visualizing, and interpreting passive acoustic telemetry
data, with a primary focus on marine environments. It offers tools
spanning data cleaning, temporal and spatial classification, movement
and home-range analysis, residency metrics, spatiotemporal overlap and
social-network analysis, and a wide range of publication-ready plots.

Most functions accept user-defined column names (e.g. `id.col`,
`datetime.col`, `station.col`). To avoid repeating these on every call,
wrap your detection data once with
[`as_moby`](https://miguelgandra.github.io/moby/reference/as_moby.md):
the resulting `mobyData` object carries the column mapping, coordinate
reference system, tagging dates and ID groups as metadata, which
downstream functions read automatically. Plain data frames remain fully
supported.

## Getting started

See [`help(package = "moby")`](https://rdrr.io/pkg/moby/man) for the
full list of functions, and the package vignette
(`browseVignettes("moby")`) for an introductory workflow.

## See also

Useful links:

- <https://github.com/miguelgandra/moby>

- <https://miguelgandra.github.io/moby/>

- Report bugs at <https://github.com/miguelgandra/moby/issues>

## Author

**Maintainer**: Miguel Gandra <m3gandra@gmail.com>
([ORCID](https://orcid.org/0000-0003-1506-5141))

Authors:

- Miguel Gandra <m3gandra@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-1506-5141))
