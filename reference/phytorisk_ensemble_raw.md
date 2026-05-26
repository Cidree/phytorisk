# Raw Phytophthora ensemble risk

Calculates Phytophthora cinnamomi's ensemble risk from raw data

## Usage

``` r
phytorisk_ensemble_raw(
  aoi,
  poi,
  dem,
  treecover,
  weights = "equal",
  th = 100,
  buffer = 50,
  include_zoospread = FALSE,
  append_mec = FALSE,
  ...,
  quiet = FALSE
)
```

## Arguments

- aoi:

  A `sf` polygon representing the area of interest. Used to mask the
  tree cover

- poi:

  A single-point `sf` object denoting the point of interest to run the
  simulations

- dem:

  A single-band `SpatRaster` with a digital elevation model

- treecover:

  A single-band `SpatRaster` where 1 represents host trees, and 0
  represents background area

- weights:

  weights of the ensemble model. The default uses the same weights for
  each model. The argument accepts a numeric vector with the
  corresponding weights. See
  [phytorisk_ensemble](https://cidree.github.io/phytorisk/reference/phytorisk_ensemble.md)

- th:

  Threshold of flow accumulation to delineate streams

- buffer:

  A buffer in meters to extend the spread in every direction

- include_zoospread:

  logical. Whether to include the optional module of
  [mec_zoospread](https://cidree.github.io/phytorisk/reference/mec_zoospread.md)

- append_mec:

  Logical. Whether to append to the results of each individual module to
  the output SpatRaster

- ...:

  arguments passed to
  [mec_zoospread](https://cidree.github.io/phytorisk/reference/mec_zoospread.md)

- quiet:

  A logical value. If `TRUE`, suppresses any informational messages.
  Defaults to `FALSE`.

## Value

A `SpatRaster`

## Examples

``` r
# \donttest{
## load packages
library(phytorisk)
library(sf)
library(terra)

## load data
poi_sf <- st_read(
  system.file("spatial/poi.geojson", package = "phytorisk"),
  quiet = TRUE
)
dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))
aoi_sf <- st_read(
  system.file("spatial/tejera.geojson", package = "phytorisk"),
  quiet = TRUE
)

## calculate ensemble risk, returning individual mechanisms
risk_equal_sr <- phytorisk_ensemble_raw(
  aoi = aoi_sf,
  poi = poi_sf,
  dem = dem_sr,
  treecover = trees_sr,
  append_mec = TRUE
)
#> 
#> ── Mec Ii - Spread in soil water ───────────────────────────────────────────────
#> ℹ Filling DEM...
#> ✔ DEM filled [26ms]
#> 
#> ℹ Filling basins...
#> ✔ Basins filled [41ms]
#> 
#> ℹ Removing depressions...
#> ✔ Depressions removed [28ms]
#> 
#> ℹ Filling depressions...
#> ✔ Depressions filled [18ms]
#> 
#> ℹ Getting flow directions...
#> ✔ Flow directions [18ms]
#> 
#> ℹ Calculating flow accumulation...
#> ✔ Flow accumulation calculated [17ms]
#> 
#> ℹ Delineating streams...
#> ✔ Streams delineated [28ms]
#> 
#> ℹ Determining the wet front
#> ✔ Wet front determined [763ms]
#> 
#> 
#> ── Mec Iii - Root-to-root contact ──────────────────────────────────────────────
#> ℹ Preparing tree data...
#> ✔ Tree data prepared [52ms]
#> 
#> ℹ Finding root-to-root contact...
#> ✔ Finished [13.5s]
#> 
#> 
#> ── Mec II - Spread in surface water ────────────────────────────────────────────
#> ℹ Calculating natural drainage network...
#> ✔ Natural drainage network calculated [54ms]
#> 
#> ℹ Identifying surface water close to foci...
#> ✔ Surface water close to foci identified [57ms]
#> 
#> ℹ Finding connected pixels...
#> ✔ Finished [833ms]
#> 

## visualize results
plot(risk_equal_sr)

# }
```
