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

## Example

This is a basic example which shows you how to solve a common problem:

#### Simulation

```r
library(rBTR)

## computer daylength, gE, soil moisture, VPD, etc.
Climate_data <- Compute_clim( climdata = LS_clim , parameters = Clim_param , syear = NA, eyear = NA )
# or computer daylength, gE , etc.
Climate_data <- Compute_clim( climdata = LS_climdata , parameters = Clim_param , syear = NA, eyear = NA )

## Simulate tree growth
Res <- btr( clim = Climate_data, parameters = BP_param, syear = NA, eyear = NA, intraannual = F)

```

#### Result

```r
## The main output are 4 tables
names(Res)
#> [1] "annaulRing"      "xylem_trait"     "microclim"       "dailyParameters"

## "annaulRing" Summarized annual results
## Year(A),  RingWidth(mm) , maxVA(μm^2)
Res$annaulRing
#> # A tibble: 11 × 6
#>     Year CellLayer CellNumber RingWidth VesselNumber maxVA
#>    <dbl>     <dbl>      <dbl>     <dbl>        <dbl> <dbl>
#>  1  2000       532     28380       4.34          380  995.
#>  2  2001       476     25393.      3.85          340  993.
#>  3  2002       500     26673.      4.11          357  995.
#>  4  2003       482     25712.      3.98          344  997.
#>  ......

## "xylem_trait" Reault of cell anatomy traits in each layer
# Year(A), cell_L: cell layer, CA(μm^2): cell area, CRD(μm): cell radial diameter, 
# CTD(μm): cell tangential diameter, CV(μm^2): cell lumen area, WA(μm^2): cell wall area,
# LWA(μm^2): lignified wall area, WT(μm): cell wall thickness,
# EDOY: Day of cell start enlargement,
# TDOY: Day of cell wall start thickening,
# DDOY: Day of cell mature,
# NoV； Number of vessels in this layer
# Raddist(μm): Radial distance of this layer from the previous ring
# 'XXX_v': trait of vessel cell.
colnames(Res$xylem_trait)
#>  [1] "Year"    "cell_L"  "CA"      "CRD"     "CTD"     "CV"      "DDOY"    "EDOY"    "LWA"     "TDOY"    "WA"      "WT"     
#> [13] "CA_v"    "CRD_v"   "CTD_v"   "CV_v"    "DDOY_v"  "EDOY_v"  "LWA_v"   "NoV_v"   "TDOY_v"  "VN_v"    "WA_v"    "WT_v"   
#> [25] "Raddist"

##  "microclim" input climate data
Res$microclim 
#> Year Month Day   TEM soilM        VPD DOY        gE        Ls rootd       dL_i gT gM gV    L_i.fiber   L_i.vessel
#> 1  2000     1   1 -26.1     0 0.01590113   1 0.2226430 0.6947605     0 0.01305703  0  0  1 7.395351e-09 6.627039e-06
#> 2  2000     1   2 -18.8     0 0.03317144   2 0.2238368 0.6959656     0 0.01446143  0  0  1 3.003902e-09 3.837076e-06
#> 3  2000     1   3 -19.4     0 0.03018919   3 0.2251499 0.6972894     0 0.01588563  0  0  1 1.152818e-09 2.146492e-06
#> 4  2000     1   4 -27.0     0 0.01596437   4 0.2265824 0.6987306     0 0.01729450  0  0  1 4.270088e-10 1.175217e-06
#> 5  2000     1   5 -28.0     0 0.01515112   5 0.2281342 0.7002878     0 0.01868669  0  0  1 1.527123e-10 6.298890e-07
#> 6  2000     1   6 -18.1     0 0.02788169   6 0.2298052 0.7019596     0 0.02006089  0  0  1 5.276033e-11 3.306052e-07
#> ......

## "dailyParameters" Daily values of dynamic parameters, to help check that the model parameters are working properly.
Res$dailyParameters
#> .....
```



#### Required data frame organization

1. ###### Gorwth parameters list
   
   Example growth parameters list.<br>
   Key cols are required columns in the data frame, they  are called in the BTR model.<br>
   Use `data(BP_param)` or` data(FM_param)` view sample parameter data frame.
   ![readme_df_Gparam.png](./man/Figs/readme_df_Gparam.png)

2. ###### Climate parameters list
   
   Example climate parameters list. <br>
   Use `data(clim_param)` view sample parameter data frame.<br>
   If you only want to calculate the day length, set the parameter ***latitude***. <br>
   The remaining parameters are used to calculate soil moisture using the **cpc-leaky bucket** model.<br>
   Use `data(clim_param)` view sample parameter data frame.
   ![readme_df_Cparam.png](./man/Figs/readme_df_Cparam.png)

3. ###### Climate data
   
   `Computer_clim()` function is used to generate the clim input data for the BTR model.<br>
   `Computer_clim()`function could use ***latitude*** to calculate ***Li (daylength)*** and ***gE (relative rate of cell growth driven by daylength)***,<br>
   ***MAT (mean air temperature)*** and ***PRE (precipitation)*** to calculate ***soilM (soil moisture)***, <br>
   ***RH (Relative Humidity)*** and ***MAT*** to calculate ***VPD (vapor pressure deficit)***.<br>
   Example climate data input.<br>
   Use `data(LS_clim)` view sample climate data.
   ![readme_df_clim1.png](./man/Figs/readme_df_clim1.png)
   If you have more reliable ***soilM (soil moisture)*** or ***VPD*** data for the sample site, `Computer_clim()` will only calculate the missing parts.<br>
   Use `data(LS_climdata)` view sample climate data2.
   ![ readme_df_clim2.png ](man/Figs/readme_df_clim2.png)


