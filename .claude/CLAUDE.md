# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**phytorisk** is an R package for modeling the spatially explicit invasion risk of *Phytophthora cinnamomi* (Pc). It implements four dispersal mechanisms and produces raster outputs representing infection risk surfaces.

## Common Commands

All development tasks use `devtools` and roxygen2. Run these in an R session:

```r
devtools::document()          # Regenerate NAMESPACE and man/ from roxygen2
devtools::build_readme()      # Rebuild README.md from README.Rmd
devtools::check()             # Full R CMD check locally
devtools::load_all()          # Load package for interactive testing
devtools::install()           # Install package locally
```

## Architecture

### Core Mechanism Functions (`R/`)

Four exported functions implement dispersal mechanisms, each accepting spatial inputs (rasters via `terra`, vectors via `sf`) and returning a normalized risk raster (0–1):

| Function | File | Mechanism |
|---|---|---|
| `mec_rootcontact()` | `mec-rootcontact.R` | Root-to-root transmission via 3×3 spatial window |
| `mec_soilwater()` | `mec-soilwater.R` | Soil water transport using flow direction/accumulation (requires `flowdem`) |
| `mec_surfacewater()` | `mec-surfacewater.R` | Surface water spread using TWI (Topographic Water Index) |
| `mec_zoospread()` | `mec-zoospread.R` | Animal movement simulation with trajectory analysis |

Two ensemble functions combine mechanism outputs:
- `phytorisk_ensemble()` — weighted combination of pre-computed mechanism rasters
- `phytorisk_ensemble_raw()` — one-step wrapper that calls all mechanisms then ensembles

### Internal Helpers (`R/utils-not-exported.R`)

Key internal algorithms (not exported):
- `find_connected()` / `is_connected()` — flood-fill connectivity for root contact spread
- `extract_neighbours()` / `isolate_pixels()` — 8-neighbor pixel operations for flow analysis
- `twi()` — Topographic Water Index: `log(1/slope)`
- `Animal()`, `move_towards_food()`, `stay_within_area()` — animal movement simulation objects/functions for `mec_zoospread`

Input validation is centralized in `R/utils-assert.R` with spatial assertion helpers.

### Roxygen Templates (`man/roxygen/templates/`)

Shared parameter documentation is defined as roxygen2 templates (not inline): `aoi.R`, `poi.R`, `dem.R`, `treecover.R`, `quiet.R`, `th.R`, `buffer.R`. Use `@template <name>` in roxygen blocks rather than duplicating parameter descriptions.

### Sample Data (`inst/spatial/`)

The package ships sample rasters and vectors for examples:
- `dem.tiff` / `dem_light.tiff` — digital elevation models (light version for fast examples)
- `trees.tiff` / `trees_light.tiff` — tree cover rasters
- `poi.geojson` — points of interest (infection sources)
- `tejera.geojson` — study area polygon (AOI)

Access via `system.file("spatial", "<file>", package = "phytorisk")`.

## Key Dependencies

- **terra** — all raster operations
- **sf** — vector spatial data
- **cli** — progress bars and user messages
- **flowdem** (GitHub: `KennethTM/flowdem`) — required for `mec_soilwater()`; install via `pak::pak("KennethTM/flowdem")`

## Development Notes

- `internal/build.R` contains the full CRAN submission workflow checklist.
- `internal/scrap.R` holds experimental code — not part of the package.
- Documentation site is built via pkgdown and deployed automatically to GitHub Pages on push to `main` (see `.github/workflows/pkgdown.yaml`).
- Commit style: semantic prefixes (`feat:`, `fix:`, `docs:`, `refactor:`).
