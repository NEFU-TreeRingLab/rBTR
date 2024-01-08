### gM 与gT计算 -------------------------------------

#' Computer gM & gT & gP
#'
#' @param clim climates about Temperature & Soil Moisture
#' @param growth_Param parameters of climates model
#'
#' @importFrom dplyr select mutate
#' @importFrom tidyr spread
#'
#' @return gM & gT

Compute_gR <- function( clim , growth_Param){

  growthParam <-growth_Param |> dplyr::select(c("parameter","values")) |>
    tidyr::spread(key = parameter, value = values)

  ### gT

  clim <- dplyr::mutate(clim, gT = NA , gM = NA, gV = NA, L_i.fiber =NA , L_i.vessel = NA)

  clim$gT[ clim$TEM < growthParam$T1 & is.na(clim$gT) ] <- 0

  clim$gT[ clim$TEM > growthParam$T4 & is.na(clim$gT) ] <- 0

  clim$gT[ clim$TEM < growthParam$T2 & is.na(clim$gT) ] <-
    (clim$TEM[ clim$TEM < growthParam$T2 & is.na(clim$gT) ] - growthParam$T1) / (growthParam$T2-growthParam$T1)

  clim$gT[ clim$TEM > growthParam$T3 & is.na(clim$gT) ] <-
    (growthParam$T4 - clim$TEM[ clim$TEM > growthParam$T3 & is.na(clim$gT) ]) / (growthParam$T4-growthParam$T3)

  clim$gT[ clim$TEM <= growthParam$T3 & is.na(clim$gT) ] <- 1

  ### gM

  clim$gM[ clim$soilM < growthParam$M1 & is.na(clim$gM) ] <- 0

  clim$gM[ clim$soilM > growthParam$M4 & is.na(clim$gM) ] <- 0

  clim$gM[ clim$soilM < growthParam$M2 & is.na(clim$gM) ] <-
    ( clim$soilM[ clim$soilM < growthParam$M2 & is.na(clim$gM) ] - growthParam$M1 ) / (growthParam$M2-growthParam$M1)

  clim$gM[ clim$soilM > growthParam$M3 & is.na(clim$gM) ] <-
    (growthParam$M4 - clim$soilM[ clim$soilM > growthParam$M3 & is.na(clim$gM) ]) / (growthParam$M4-growthParam$M3)

  clim$gM[ clim$soilM <= growthParam$M3 & is.na(clim$gM) ] <- 1

  ### gT


  clim$gV[ clim$VPD < growthParam$VPD1 & is.na(clim$gV) ] <- 0

  clim$gV[ clim$VPD > growthParam$VPD4 & is.na(clim$gV) ] <- 0

  clim$gV[ clim$VPD < growthParam$VPD2 & is.na(clim$gV) ] <-
    (clim$VPD[ clim$VPD < growthParam$VPD2 & is.na(clim$gV) ] - growthParam$VPD1) / (growthParam$VPD2-growthParam$VPD1)

  clim$gV[ clim$VPD > growthParam$VPD3 & is.na(clim$gV) ] <-
    (growthParam$VPD4 - clim$VPD[ clim$VPD > growthParam$VPD3 & is.na(clim$gV) ]) / (growthParam$VPD4-growthParam$VPD3)

  clim$gV[ clim$VPD <= growthParam$VPD3 & is.na(clim$gV) ] <- 1

  ## L_i.fiber , gL_vessel

  clim$L_i.fiber  <- growthParam$a.fiber * exp(
    -exp(  growthParam$b.fiber - growthParam$c.fiber * clim$dL_i ) )

  clim$L_i.vessel <- growthParam$a.vessel * exp(
    -exp(  growthParam$b.vessel - growthParam$c.vessel * clim$dL_i ) )



  return(clim)

}  ## Compute_gR end --------------------------------------

### gM 与gT计算 -------------------------------------

#' Computer gM & gT & gP
#'
#' @param clim climates about Temperature & Soil Moisture
#' @param growth_Param parameters of climates model
#'
#' @importFrom dplyr select mutate
#' @importFrom tidyr spread
#'
#' @return gM & gT
Compute_gR2 <- function( clim , growth_Param){

  growthParam <-growth_Param |> dplyr::select(c("parameter","values")) |>
    tidyr::spread(key = parameter, value = values)

  ### gT

  clim <- dplyr::mutate(clim, gT = NA , gM = NA, gV = NA, L_i.fiber =NA , L_i.vessel = NA)

  Ta <- clim$TEM + 273.15
  R <- 8.314

  clim$gT <- ( Ta * exp( - growthParam$deltaH_A_Da / ( R * Ta) ) /
                 (1 + exp( growthParam$deltaS_D / R - growthParam$deltaH_D /(R*Ta) ) ) ) |>
    nor( Zeros = T)
  clim$gT[clim$TEM <= growthParam$T1| clim$TEM >= growthParam$T4 ] <- 0



  ### gM

  clim$gM[ clim$soilM < growthParam$M1 & is.na(clim$gM) ] <- 0

  clim$gM[ clim$soilM > growthParam$M4 & is.na(clim$gM) ] <- 0

  clim$gM[ clim$soilM < growthParam$M2 & is.na(clim$gM) ] <-
    (clim$soilM[ clim$soilM < growthParam$M2 & is.na(clim$gM) ] - growthParam$M1) / (growthParam$M2-growthParam$M1)

  clim$gM[ clim$soilM > growthParam$M3 & is.na(clim$gM) ] <-
    (growthParam$M4 - clim$soilM[ clim$soilM > growthParam$M3 & is.na(clim$gM) ]) / (growthParam$M4-growthParam$M3)

  clim$gM[ clim$soilM <= growthParam$M3 & is.na(clim$gM) ] <- 1

  ### gT


  clim$gV[ clim$VPD < growthParam$VPD1 & is.na(clim$gV) ] <- 0

  clim$gV[ clim$VPD > growthParam$VPD4 & is.na(clim$gV) ] <- 0

  clim$gV[ clim$VPD < growthParam$VPD2 & is.na(clim$gV) ] <-
    (clim$VPD[ clim$VPD < growthParam$VPD2 & is.na(clim$gV) ] - growthParam$VPD1) / (growthParam$VPD2-growthParam$VPD1)

  clim$gV[ clim$VPD > growthParam$VPD3 & is.na(clim$gV) ] <-
    (growthParam$VPD4 - clim$VPD[ clim$VPD > growthParam$VPD3 & is.na(clim$gV) ]) / (growthParam$VPD4-growthParam$VPD3)

  clim$gV[ clim$VPD <= growthParam$VPD3 & is.na(clim$gV) ] <- 1

  ## L_i.fiber , gL_vessel

  clim$L_i.fiber  <- growthParam$a.fiber * exp(
    -exp(  growthParam$b.fiber - growthParam$c.fiber * clim$dL_i ) )

  clim$L_i.vessel <- growthParam$a.vessel * exp(
    -exp(  growthParam$b.vessel - growthParam$c.vessel * clim$dL_i ) )



  return(clim)

}  #
