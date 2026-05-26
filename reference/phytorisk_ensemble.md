# Phytophthora ensemble risk

Calculates Phytophthora cinnamomi's ensemble risk from mechanisms

## Usage

``` r
phytorisk_ensemble(
  mec_soilwater,
  mec_rootcontact,
  mec_surfacewater,
  mec_zoospread = NULL,
  weights = "equal"
)
```

## Arguments

- mec_soilwater:

  The result of
  [mec_soilwater](https://cidree.github.io/phytorisk/reference/mec_soilwater.md)

- mec_rootcontact:

  The result of
  [mec_rootcontact](https://cidree.github.io/phytorisk/reference/mec_rootcontact.md)

- mec_surfacewater:

  The result of
  [mec_surfacewater](https://cidree.github.io/phytorisk/reference/mec_surfacewater.md)

- mec_zoospread:

  Optional module. The result of
  [mec_zoospread](https://cidree.github.io/phytorisk/reference/mec_zoospread.md)

- weights:

  Weights of the ensemble model. The default uses the same weights for
  each model. The argument accepts a numeric vector with the
  corresponding weights. See *Details*

## Value

A `SpatRaster`

## Details

This function calculates the ensemble risk as the weighted sum of the
individual risk models as follows:

\\Risk = w_1 \cdot Mec\_{soilwater} + w_2 \cdot Mec\_{rootcontact} + w_3
\cdot Mec\_{mec_surfacewater} + w_4 \cdot Mec\_{zoospread}\\.

Being \\w_i = 0.25\\ when
[mec_zoospread](https://cidree.github.io/phytorisk/reference/mec_zoospread.md)
is not null, and \\w_i = \frac{1}{3}\\ when
[mec_zoospread](https://cidree.github.io/phytorisk/reference/mec_zoospread.md)
is null. Additionally, the user can specify a numeric vector of length 3
or 4 specifying the desired weights based on expert criteria.

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

## first, calculate the individual mechanisms
mec_soilwater_sr <- mec_soilwater(dem_sr, poi_sf)
#> ℹ Filling DEM...
#> ✔ DEM filled [16ms]
#> 
#> ℹ Filling basins...
#> ✔ Basins filled [22ms]
#> 
#> ℹ Removing depressions...
#> ✔ Depressions removed [15ms]
#> 
#> ℹ Filling depressions...
#> ✔ Depressions filled [17ms]
#> 
#> ℹ Getting flow directions...
#> ✔ Flow directions [16ms]
#> 
#> ℹ Calculating flow accumulation...
#> ✔ Flow accumulation calculated [22ms]
#> 
#> ℹ Delineating streams...
#> ✔ Streams delineated [18ms]
#> 
#> ℹ Determining the wet front
#> ✔ Wet front determined [498ms]
#> 
mec_surface_lst   <- mec_surfacewater(dem_sr, mec_soilwater_sr, poi_sf)
#> ℹ Calculating natural drainage network...
#> ✔ Natural drainage network calculated [35ms]
#> 
#> ℹ Identifying surface water close to foci...
#> ✔ Surface water close to foci identified [37ms]
#> 
#> ℹ Finding connected pixels...
#> ✔ Finished [548ms]
#> 
mec_rootcontact_sr <- mec_rootcontact(trees_sr, aoi_sf, poi_sf)
#> ℹ Preparing tree data...
#> ✔ Tree data prepared [108ms]
#> 
#> ℹ Finding root-to-root contact...
#> ✔ Finished [13.6s]
#> 

## calculate ensemble risk using equal weights
risk_equal_sr <- phytorisk_ensemble(
 mec_soilwater = mec_soilwater_sr,
 mec_rootcontact = mec_rootcontact_sr,
 mec_surfacewater = mec_surface_lst
)

## assign more weight to root-to-root contact
risk_weighted_sr <- phytorisk_ensemble(
 mec_soilwater = mec_soilwater_sr,
 mec_rootcontact = mec_rootcontact_sr,
 mec_surfacewater = mec_surface_lst,
 weights = c(.3, .4, .3)
)
# }
```
