<!-- README.md is generated from README.Rmd. Please edit that file -->

# rBTR

<!-- badges: start -->

<!-- badges: end -->

**Rpackage for the BTR model is still under development.**

We expect an update of the example section of the readme file by June 8th.



rBTR is an Rpackage of the Broadleaved Tree-Ring (BTR) model. 



## Installation

You can install the development version of BTR model like so:

```r
# devtools::install_github("kdoodk/rBTR")
```

## Example

This is a basic example which shows you how to solve a common problem:

#### Simulation

```r
library(rBTR)

## computer daylength, soil moisture and VPD
Climate_data <- Compute_clim( climdata = LS_clim , parameters = Clim_param , syear = NA, eyear = NA )
# or
Climate_data <- Compute_clim( climdata = LS_climdata , parameters = Clim_param , syear = NA, eyear = NA )

## Simulate tree growth
Res <- btr( clim = Climate_data, parameters = BP_param, syear = NA, eyear = NA, intraannual = F)

```

#### Required data frame organization

1. ###### Gorwth parameters list
   
   Example growth parameters list. Key cols are required columns in the data frame, they  are called in the BTR model.
   Use `data(BP_param)` or` data(FM_param)` view sample parameter data frame.
   
   ![](man/FIgs/readme_df_Gparam.png)

2. ###### Climate parameters list

   Example climate parameters list. 
   Use `data(clim_param)` view sample parameter data frame.
   If you only want to calculate the day length, set the parameter ***latitude***. The remaining parameters are used to           calculate soil moisture using the **cpc-leaky bucket** model.

   <img src = "man/FIgs/readme_df_Cparam.png" , width = "80%" / >
   ![](man/FIgs/readme_df_Cparam.png)

3. ###### Climate data
   `Computer_clim()`'function could use ***latitude*** to calculate ***Li (daylength)*** and ***gE (relative rate of cell growth driven by daylength)***, ***MAT (mean air temperature)*** and ***PRE (precipitation)*** to calculate ***soilM (soil moisture)***, and ***RH (Relative Humidity)*** and ***MAT*** to calculate ***VPD (vapor pressure deficit)***.
   Example climate data input.
   Use `data(LS_clim)` view sample climate data.
    
   <img src = "man/FIgs/readme_df_readme_df_clim1.png" , width = "80%" / >
   ![](man/FIgs/readme_df_readme_df_clim1.png)

   If you have more reliable soil moisture or VPD data for the sample site, `Computer_clim()` will only calculate the missing parts.
   Use `data(LS_climdata)` view sample climate data.

   <img src = "man/FIgs/readme_df_readme_df_clim2.png" , width = "80%" / >
   ![](man/FIgs/readme_df_readme_df_clim2.png)