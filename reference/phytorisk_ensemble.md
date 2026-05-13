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
## TODO
```
