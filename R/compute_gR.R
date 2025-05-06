### gM 与gT计算 -------------------------------------
#' Computer gM & gT & gP
#'
#' @export
#'
#' @param clim climates about Temperature & Soil Moisture
#' @param growth_Param parameters of climates model
#' @param gTmethod formula of gT
#'
#' @importFrom dplyr select mutate
#' @importFrom tidyr spread
#'
#' @return gM & gT
compute_gR <- function( clim , growth_Param, gTmethod ){

  ## 计算gT
  if (gTmethod == 'VS') {
    clim <- clim |>
      dplyr::mutate(gT = dplyr::case_when( TEM < growth_Param$T1 | TEM >= growth_Param$T4 ~ 0,
                                           TEM >= growth_Param$T2 & TEM <= growth_Param$T3 ~ 1,
                                           TEM < growth_Param$T2 & TEM >= growth_Param$T1 ~
                                             ( TEM - growth_Param$T1 ) / (growth_Param$T2-growth_Param$T1),
                                           TEM < growth_Param$T4 & TEM > growth_Param$T3 ~
                                             ( growth_Param$T4 - TEM ) / (growth_Param$T4-growth_Param$T3) )
      )

  } else {

    clim$gT <- ( Ta * exp( - growthParam$deltaH_A_Da / ( 8.314 * (clim$TEM + 273.15)) ) /
                   (1 + exp( growthParam$deltaS_D / 8.314 - growthParam$deltaH_D /( 8.314 * (clim$TEM + 273.15)) ) ) ) |>
      nor( Zeros = T)
    clim$gT[clim$TEM <= growthParam$T1| clim$TEM >= growthParam$T4 ] <- 0
  }

  ## 计算 gM,gV,L_i.fiber,L_i.vessel
  clim <- clim |>
    dplyr::mutate(gM = dplyr::case_when( soilM < growth_Param$M1 | soilM >= growth_Param$M4 ~ 0,
                                         soilM >= growth_Param$M2 & soilM <= growth_Param$M3 ~ 1,
                                         soilM < growth_Param$M2 & soilM >= growth_Param$M1 ~
                                             ( soilM - growth_Param$M1 ) / (growth_Param$M2-growth_Param$M1),
                                         soilM < growth_Param$M4 & soilM > growth_Param$M3 ~
                                             ( growth_Param$M4 - soilM ) / (growth_Param$M4-growth_Param$M3)
                                        ),
                  gV = dplyr::case_when( VPD < growth_Param$VPD1 | VPD >= growth_Param$VPD4 ~ 0,
                                         VPD >= growth_Param$VPD2 & VPD <= growth_Param$VPD3 ~ 1,
                                         VPD < growth_Param$VPD2 & VPD >= growth_Param$VPD1 ~
                                             ( VPD - growth_Param$VPD1 ) / (growth_Param$VPD2-growth_Param$VPD1),
                                         VPD < growth_Param$VPD4 & VPD > growth_Param$VPD3 ~
                                             ( growth_Param$VPD4 - VPD ) / (growth_Param$M4-growth_Param$VPD3)
                                        ),
                  L_i.fiber  = growth_Param$a.fiber * exp(
                    -exp(  growth_Param$b.fiber - growth_Param$c.fiber * dL_i ) ),
                  L_i.vessel = growth_Param$a.vessel * exp(
                    -exp(  growth_Param$b.vessel - growth_Param$c.vessel * dL_i ) )
      )

  ## 计算有效积温和气候生长季长度
  ## 生长季开始::积温 > 阈值 & 第一次连续5天温度 > T1
  ## 生长季结束:: 最高温度日期以后 & 积温 > 阈值 & 第一次连续5天温度 < T1

  clim <- clim |>
    dplyr::group_by(Year) |>
    dplyr::arrange(Year,DOY) |> dplyr::mutate( aT  = TEM - temThreshold,
                                               aT = dplyr::case_when(
                                                 aT <0 ~ 0,
                                                 aT >= 0 ~ aT
                                               ),
                                               aaT = cumsum(aT),
                                               tn = c(dplyr::case_when(aT == 0 ~ 1 ,aT != 0 ~ 0  ) |>
                                                        embed(5) |> rowSums() , 5,5,5,5 ) ,
                                               HotDay = which.max(TEM),
                                               Fs1 = which(aaT >= fixparam.divi$AAT & tn == 0    )[1],
                                               Fs2 = which( aaT >= fixparam.divi$AAT & tn == 5 & DOY >= HotDay  )[1],
                                               GS = dplyr::case_when( DOY >= Fs1 & DOY <= Fs2  ~ 'GR',
                                                                      .default = "UN")
    ) |> dplyr::select(-tn,-Fs1,-Fs2,-HotDay)



}
