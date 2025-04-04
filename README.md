
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ALGR

<!-- badges: start -->

[![R-CMD-check](https://github.com/pogoyoly/ALGR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pogoyoly/ALGR/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/pogoyoly/ALGR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/pogoyoly/ALGR?branch=main)
<!-- badges: end -->

ALGR is an agricultural landcover generator that is designed to ease the
integration of artificial agricultural landcover maps in the work flow
of ecological modellers. The package is designed to generate landscapes
using different algorithms, and then store the information in an output
file that includes both a raster layer, and a list containing all
information of the fields and farmers. This package allows to generate
landscapes in a systematic reproducable way while controlling for
different variables.

## Installation

You can install the development version of ALGR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pogoyoly/ALGR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ALGR)
## basic example code
r<-generate_pn(500,500,1,1,3,0.005,TRUE, "land_percentage", percetange = 70)
output<-establish_pac(potential_space= r,
                      cell_size=1,
                      includsion_value = 1,
                      mean_field_size = 1000,
                      sd_field_size = 500,
                      distribution = "norm",
                      mean_shape_index = 4,
                      sd_shape_index = 1,
                      percent = 80)


return_by_field(output, method = 1)
```

<img src="man/figures/README-example-1.png" width="100%" />
