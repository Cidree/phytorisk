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

- ...:

  arguments passed to
  [mec_zoospread](https://cidree.github.io/phytorisk/reference/mec_zoospread.md)

- quiet:

  A logical value. If `TRUE`, suppresses any informational messages.
  Defaults to `FALSE`.

## Value

A `SpatRaster`

## Details

\#TODO
