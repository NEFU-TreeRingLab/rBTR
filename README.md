<!-- README.md is generated from README.Rmd. Please edit that file -->

# rBTR

<!-- badges: start -->

<!-- badges: end -->

**This is the Beta version of BTR model**

rBTR is an Rpackage of the Broadleaved Tree-Ring (BTR) model. 



## Installation

You can install the development version of BTR model like so:

```r
# devtools::install_github("kdoodk/rBTR")
```
## Note

Now rBTR can utilize multi-core parallel computing. 

```r
rBTR::btr_parallel( clims, param ,age, Cores = 12 ) 
```
