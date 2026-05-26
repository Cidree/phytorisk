# Transmission through root-to-root contact

Simulates diffusive inoculum spread through root contact between
neighbouring trees

## Usage

``` r
mec_rootcontact(treecover, aoi, poi, quiet = FALSE)
```

## Arguments

- treecover:

  A single-band `SpatRaster` where 1 represents host trees, and 0
  represents background area

- aoi:

  A `sf` polygon representing the area of interest. Used to mask the
  tree cover

- poi:

  A single-point `sf` object denoting the point of interest to run the
  simulations

- quiet:

  A logical value. If `TRUE`, suppresses any informational messages.
  Defaults to `FALSE`.

## Value

A `SpatRaster`

## Details

This function models inoculum movement due to root-to-root contact based
on tree continuity. It calculates vegetation continuity by counting
connected pixels within a 3×3 window. If continuity exists, a new raster
is created with connected pixels marked as 1. The function then
identifies all pixels connected to the foci within the study area. If no
direct pixel-to-pixel contact between trees exists, the infection risk
is considered zero.

## References

Cardillo, E., Abad, E., Meyer, S., 2021. Iberian oak decline caused by
Phytophthora cinnamomi: A spatiotemporal analysis incorporating the
effect of host heterogeneities at landscape scale. For. Pathol. 51,
e12667. [doi:10.1111/efp.12667](https://doi.org/10.1111/efp.12667)

Cardillo, E., Acedo, A., Abad, E., 2018. Topographic effects on
dispersal patterns of Phytophthora cinnamomi at a stand scale in a
Spanish heathland. PloS One 13, e0195060.

## Examples

``` r
# \donttest{
## load packages
library(phytorisk)
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(terra)
#> terra 1.9.27

## load data
study_area_sf <- st_read(
  system.file("spatial/tejera.geojson", package = "phytorisk"),
  quiet = TRUE
)

poi_sf <- st_read(
  system.file("spatial/poi.geojson", package = "phytorisk"),
  quiet = TRUE
)

trees_sr <- rast(system.file("spatial/trees_light.tiff", package = "phytorisk"))

## simulate mechanism
mec_rootcontact_sr <- mec_rootcontact(trees_sr, study_area_sf, poi_sf)
#> ℹ Preparing tree data...
#> ✔ Tree data prepared [45ms]
#> 
#> ℹ Finding root-to-root contact...
#> ✔ Finished [13.4s]
#> 
# }
```
