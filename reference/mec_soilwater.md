# Transport through soil water

Simulates the spread of *Phytophthora cinnamomi* inoculum within the
soil matrix through the soil water

## Usage

``` r
mec_soilwater(dem, poi, th = 100, quiet = FALSE)
```

## Arguments

- dem:

  A single-band `SpatRaster` with a digital elevation model

- poi:

  A single-point `sf` object denoting the point of interest to run the
  simulations

- th:

  Threshold of flow accumulation to delineate streams

- quiet:

  A logical value. If `TRUE`, suppresses any informational messages.
  Defaults to `FALSE`.

## Value

A `SpatRaster`

## Details

This function models inoculum movement in soil, considering flow
direction and altitude under wet and intermediate moisture conditions.
It identifies pixels where flow direction matches or is adjacent to the
foci, then processes connected pixels based on altitude differences. The
result is a binary raster showing the spatial pattern of inoculum
dispersal, constrained by topography and moisture.

## References

Ristaino, J., Gumpertz, M., 2000. New Frontiers in the Study of
Dispersal and Spatial Analysis of Epidemics Caused by Species in the
Genus Phytophthora. Annu. Rev. Phytopathol. 38, 541–576.
[doi:10.1146/annurev.phyto.38.1.541](https://doi.org/10.1146/annurev.phyto.38.1.541)

Vannini, A., Natili, G., Anselmi, N., Montaghi, A., Vettraino, A.M.,
2010. Distribution and gradient analysis of Ink disease in chestnut
forests. For. Pathol. 40, 73–86.
[doi:10.1111/j.1439-0329.2009.00609.x](https://doi.org/10.1111/j.1439-0329.2009.00609.x)

Vannini, A., Natili, G., Thomidis, T., Belli, C., Morales-Rodriguez, C.,
2021. Anthropogenic and landscape features are associated with ink
disease impact in Central Italy. For. Pathol. 51, e12722.
[doi:10.1111/efp.12722](https://doi.org/10.1111/efp.12722)

## Examples

``` r
## load packages
library(phytorisk)
library(sf)
library(terra)

## load data
dem_sr <- rast(system.file("spatial/dem_light.tiff", package = "phytorisk"))
poi_sf <- st_read(
  system.file("spatial/poi.geojson", package = "phytorisk"),
  quiet = TRUE
)

## simulate mechanism
mec_soilwater_sr <- mec_soilwater(dem_sr, poi_sf)
#> Loading required namespace: flowdem
#> ℹ Filling DEM...
#> ✔ DEM filled [59ms]
#> 
#> ℹ Filling basins...
#> ✔ Basins filled [36ms]
#> 
#> ℹ Removing depressions...
#> ✔ Depressions removed [29ms]
#> 
#> ℹ Filling depressions...
#> ✔ Depressions filled [30ms]
#> 
#> ℹ Getting flow directions...
#> ✔ Flow directions [30ms]
#> 
#> ℹ Calculating flow accumulation...
#> ✔ Flow accumulation calculated [20ms]
#> 
#> ℹ Delineating streams...
#> ✔ Streams delineated [20ms]
#> 
#> ℹ Determining the wet front
#> ✔ Wet front determined [930ms]
#> 
```
