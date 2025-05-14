
<!-- README.md is generated from README.Rmd. Please edit that file -->

# phytorisk

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/phytorisk)](https://CRAN.R-project.org/package=phytorisk)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<!-- badges: end -->

`phytorisk` is an R package developed to quantify the risk of invasion
of the oomycete *Phytophthora cinnamomi* Rands (Pc hereafter). Pc is a
pathogen that infects plant’s roots, and in many cases, kills the host.
The dispersal of the pathogen depends on multiple factors, which are
depicted in the next figure:

<figure id="fig-dispersal">
<img src="man/figures/Fig1.jpg" width="825"
alt="Potential mechanisms of dispersal of Pc" />
<figcaption aria-hidden="true">Potential mechanisms of dispersal of
Pc</figcaption>
</figure>

1)  Movement of inoculum in the field due to diffusive root-to-root
    contact `mec_rootcontact()`

2)  Inoculum movement to roots in water particles within soil
    `mec_soilwater()`

3)  Dispersion of inoculum by domestic and wild animal movement
    `mec_zoospread()`

4)  Inoculum spread in surface water `mec_surfaceflow()`

## Installation

You can install the last stable version of `phytorisk` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("Cidree/phytorisk")
```

## 
