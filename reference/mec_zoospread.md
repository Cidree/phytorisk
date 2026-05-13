# Dispersion by animals

Optional Module. Simulates dispersion by domestic and wild animal
movement

## Usage

``` r
mec_zoospread(
  aoi,
  poi,
  mec_surface,
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

- mec_surface:

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
## load packages
# TODO
```
