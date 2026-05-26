# Inoculum spread in surface water

Simulates the spread of the inoculum dispersal through surface water by
integrating hydrological flow paths and wetness conditions

## Usage

``` r
mec_surfacewater(dem, mec_soilwater, poi, buffer = 50, quiet = FALSE)
```

## Arguments

- dem:

  A single-band `SpatRaster` with a digital elevation model

- mec_soilwater:

  The result of
  [mec_soilwater](https://cidree.github.io/phytorisk/reference/mec_soilwater.md)

- poi:

  A single-point `sf` object denoting the point of interest to run the
  simulations

- buffer:

  A buffer in meters to extend the spread in every direction

- quiet:

  A logical value. If `TRUE`, suppresses any informational messages.
  Defaults to `FALSE`.

## Value

A list with the surface flow risk (SpatRaster), and the surface water
(sf)

## Details

This function models inoculum movement in downslope streams using the
natural drainage network and wetlands. It creates a binary raster for
the closest stream to the foci and selects associated downslope streams.
The topographic water index (TWI) is used as a proxy for soil moisture
to identify wetlands, creating a binary raster for pixels above the 90th
percentile of TWI. The wetland raster is then combined with the stream
raster to define potential areas for inoculum spread.

## References

Li, A.Y., Williams, N., Fenwick, S.G., Hardy, G.E.St.J., Adams, P.J.,
2014b. Potential for dissemination of Phytophthora cinnamomi by feral
pigs via ingestion of infected plant material. Biol. Invasions 16,
765–774.
[doi:10.1007/s10530-013-0535-7](https://doi.org/10.1007/s10530-013-0535-7)

Ruiz-Gómez, F.J., Pérez-de-Luque, A., Navarro-Cerrillo, R.M., 2019. The
involvement of Phytophthora root rot and drought stress in holm oak
decline: from ecophysiology to microbiome influence. Curr. For. Rep. 5,
251–266.

## Examples

``` r
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

## simulate mechanism
mec_soilwater_sr <- mec_soilwater(dem_sr, poi_sf)
#> ℹ Filling DEM...
#> ✔ DEM filled [18ms]
#> 
#> ℹ Filling basins...
#> ✔ Basins filled [24ms]
#> 
#> ℹ Removing depressions...
#> ✔ Depressions removed [17ms]
#> 
#> ℹ Filling depressions...
#> ✔ Depressions filled [31ms]
#> 
#> ℹ Getting flow directions...
#> ✔ Flow directions [31ms]
#> 
#> ℹ Calculating flow accumulation...
#> ✔ Flow accumulation calculated [41ms]
#> 
#> ℹ Delineating streams...
#> ✔ Streams delineated [34ms]
#> 
#> ℹ Determining the wet front
#> ✔ Wet front determined [848ms]
#> 
mec_surface_sr <- mec_surfacewater(dem_sr, mec_soilwater_sr, poi_sf)
#> ℹ Calculating natural drainage network...
#> ✔ Natural drainage network calculated [95ms]
#> 
#> ℹ Identifying surface water close to foci...
#> ✔ Surface water close to foci identified [62ms]
#> 
#> ℹ Finding connected pixels...
#> ✔ Finished [928ms]
#> 
```
