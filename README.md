** Other Language Versions: [简体中文](README_zh_cn.md)
# rBTR

`rBTR` is an R package implementing the Broadleaved Tree-Ring (BTR) model. For more details, please refer to the following paper:

> Zhao, B. *et al.* A process-based model of climate-driven xylogenesis and tree-ring formation in broad-leaved trees (BTR). *Tree Physiol.* **44**, tpae127 (2024).[DOI:[10.1093/treephys/tpae127](http://dx.doi.org/10.1093/treephys/tpae127) ]

## Installation

You can install the `rBTR` package using:
``` r
# Install directly from GitHub
# devtools::install_github("NEFU-TreeRingLab/rBTR")
```
---
## How to Use
#### Running the Model

The `rBTR` model includes two main functions:
1. `btr()`:Runs the model in a single thread.
2. `btr_parallel()`:Supports parallel computation for faster simulations.
Both functions work the same way, but `btr()` is better for short simulations, while `btr_parallel()` is ideal for long simulations.
```r
# Running the BTR model
# Single-threaded version
# btr(clim = Clims, parameters = BPparam, age = BPage, writeRes = F, Named = "TestData")

# Parallel version
# btr_parallel(clim = Clims, parameters = BPparam, age = BPage, writeRes = F, Named = "TestData")

# Notes:
# `clim`   - Input meteorological data (data.frame)
# `paramters` - Model parameters (data.frame)
# `age`    - Age trend data (data.frame)
# `writeRes` - Boolean. If TRUE, the results are automatically saved in a folder named  "res_(running time)" in the working directory.
# `Named`  - String. If `writeRes` is TRUE, results are saved in a folder named after this string.
```

- Model parameterization methods can be consulted in the R package `BTRtools` [Link:BTRtools](https://github.com/NEFU-TreeRingLab/BTRtools)

####  Preparing Input Data

- Meteorological Data

The input meteorological data should follow this format:
| Year | Month | Day  | TEM   | VPD   | soilM |
| ---- | ----- | ---- | ----- | ----- | ----- |
| 2020 | 1     | 1    | -23.7 | 0.029 | 0.31  |

***VPD*** (Vapor Pressure Deficit) can be replaced with ***RH***  (Relative Humidity).
***soilM*** (Soil Moisture) can be estimated using ***PRE*** (precipitation data). However, the soil moisture estimation module hasn’t been rigorously tested yet.

Before running the model, the input data needs to be pre-processed. Specify the  ***lat*** (site latitude) and ***rootd*** (root zone depth) during this step.
```r
# Example: Pre-processing meteorological data
# Load sample meteorological data
# data(Climdata)

# Pre-process the data
# Clims <- Compute_clim(Climdata, lat = 47.13, rootd = 1000)

# Example with alternative meteorological data
# data(Climdata.2)
# Clims <- Compute_clim(Climdata.2, lat = 47.13, rootd = 1000)
```

The pre-processed input data will look like this:
| Year | DOY  | TEM    | VPD   | soilM | gE    | Ls    | dL_i   | rootd |
| ---- | ---- | ------ | ----- | ----- | ----- | ----- | ------ | ----- |
| 2020 | 1    | --23.7 | 0.029 | 0.31  | 0.532 | 0.697 | 0.0144 | 1000  |

***Ls***: Daylight hours (can be replaced with field-measured values).
***rootd***: Root zone depth (field-measured values can improve simulation accuracy).
Example data is included in the package.

#### Age Trend

The model doesn’t simulate age trends directly, so you must provide this as input. The data format is as follows:
| Year | age  | Tage      | Lage      |
| ---- | ---- | --------- | --------- |
| 2000 | 41   | 0.1718823 | 0.9606754 |

***Tage*** is the trend of cambial activity;

***Lage*** is the trend of cell enlargement.



Example data is included in the package:

```r
# Example: Load age trend data
# data(BPage)
```
#### Model Parameters

The model requires a set of input parameters, which are summarized in a data.frame. A default parameter table is included in the package and can be exported to ***CSV*** or ***Excel*** for modification.
```r
# Example: Model parameters
# data(BPparam)

# Save the parameter table locally for editing
# write.csv(BPparam, "BPparam.csv")  # Save as CSV
# openxlsx::write.xlsx(BPparam, "BPparam.xlsx")  # Save as Excel
```
