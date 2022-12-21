
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MS2extract

<!-- badges: start -->
<!-- badges: end -->

The goal of MS2extract is to create a tool to import MS2 data of known
standards or material and targeted extract the MS2 expectra in order to
create an in-house MS2 library.

## Installation

You can install the development version of MS2extract from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("CooperstoneLab/MS2extract")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MS2extract)
## calculating ppm range

chlorogenic_acid_pos <- 355.1023
ppm_error = 10
# Calculate ranges
ppm_baund <- ppm_range(mz = chlorogenic_acid_pos, ppm = ppm_error)
```
