# Dispersion by animals

Optional Module. Simulates dispersion by domestic and wild animal
movement

## Usage

``` r
mec_zoospread(
  aoi,
  poi,
  mec_surfacewater,
  n_animals = 5,
  n_steps = 100,
  pixel_size = 1,
  n_iter = 10,
  dist = 5,
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

- mec_surfacewater:

  The result of
  [mec_surfacewater](https://cidree.github.io/phytorisk/reference/mec_surfacewater.md)

- n_animals:

  Number of simulated animals

- n_steps:

  Number of steps of each animal for the simulations

- pixel_size:

  Size of the movement in pixels

- n_iter:

  Number of random iterations for each animal

- dist:

  Filter trajectories less than the specified meters away of the POI
  (susceptible to inoculum)

- quiet:

  A logical value. If `TRUE`, suppresses any informational messages.
  Defaults to `FALSE`.

## Value

A `SpatRaster`

## Details

This function models animal movement as a series of straight-line steps
towards resource sources, while avoiding exclusion areas and staying
within defined boundaries. It randomly assigns initial animal positions
and calculates movement direction towards resources. The direction is
normalized and randomized using a rotation matrix. Animal trajectories
are stored in vector format, and those crossing the foci are identified
for further analysis. Simulation parameters, such as the number of
animals and steps, are set at the start.

## References

Kliejunas, J.T., Ko, W.H., 1976. Dispersal of Phytophthora cinnamomi on
the island of Hawaii. Phytopathology 66, 457–460.

Li, A.Y., Williams, N., Fenwick, S.G., Hardy, G.E.St.J., Adams, P.J.,
2014b. Potential for dissemination of Phytophthora cinnamomi by feral
pigs via ingestion of infected plant material. Biol. Invasions 16,
765–774.
[doi:10.1007/s10530-013-0535-7](https://doi.org/10.1007/s10530-013-0535-7)

Cardillo, E., Acedo, A., Abad, E., 2018. Topographic effects on
dispersal patterns of Phytophthora cinnamomi at a stand scale in a
Spanish heathland. PloS One 13, e0195060.

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

## first, calculate the soil water and surface water dispersal mechanisms
mec_soilwater_sr <- mec_soilwater(dem_sr, poi_sf)
#> ℹ Filling DEM...
#> ✔ DEM filled [25ms]
#> 
#> ℹ Filling basins...
#> ✔ Basins filled [25ms]
#> 
#> ℹ Removing depressions...
#> ✔ Depressions removed [18ms]
#> 
#> ℹ Filling depressions...
#> ✔ Depressions filled [19ms]
#> 
#> ℹ Getting flow directions...
#> ✔ Flow directions [32ms]
#> 
#> ℹ Calculating flow accumulation...
#> ✔ Flow accumulation calculated [30ms]
#> 
#> ℹ Delineating streams...
#> ✔ Streams delineated [34ms]
#> 
#> ℹ Determining the wet front
#> ✔ Wet front determined [732ms]
#> 
mec_surface_sr   <- mec_surfacewater(dem_sr, mec_soilwater_sr, poi_sf)
#> ℹ Calculating natural drainage network...
#> ✔ Natural drainage network calculated [69ms]
#> 
#> ℹ Identifying surface water close to foci...
#> ✔ Surface water close to foci identified [96ms]
#> 
#> ℹ Finding connected pixels...
#> ✔ Finished [843ms]
#> 

## calculate the spread by animals (dummy example)
mec_zoospread_sr <- mec_zoospread(
  aoi = aoi_sf,
  poi = poi_sf,
  mec_surface = mec_surface_sr,
  n_animals = 5,
  n_steps = 5,
  pixel_size = 1,
  n_iter = 2,
  dist = 5
)
#> 
#> ── Starting animal movement simulation ──
#> 
#> ✔ Simulation completed. 2 trajectories generated.
#> ℹ Preparing results
#> ! There are no trajectories within 5 meters of the inoculum
#> ℹ Preparing results
#> ✔ Success [26ms]
#> 
# }
```
