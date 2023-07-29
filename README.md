
# heatwave3 <img src="logo.png" width=200 align="right" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/robwschlegel/heatwave3/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/robwschlegel/heatwave3/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/robwschlegel/heatwave3/branch/main/graph/badge.svg)](https://app.codecov.io/gh/robwschlegel/heatwave3?branch=main)
<!-- badges: end -->

Welcome to heatwave3. The goal of this package is to improve on the serial application of the Hobday et al. 2016, 2018 marine heatwave (MHW) definition on time series data. Traditionally we have had to detect MHW within one time series at a time. For large gridded datasets 

## Installation

You can install the development version of heatwave3 from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("robwschlegel/heatwave3")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(heatwave3)

# Set a file pathway
my_nc <- "~/Documents/ncdf/data.nc"
detect3(file_in = my_nc, file_out = "~/Documents/ncdf/res.nc")
```

## Code of Conduct

Please note that the heatwave3 project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

