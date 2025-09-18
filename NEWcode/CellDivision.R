#' cell_division function
#'
#' @param Doy_GR Today growth rate data
#' @param fixparam.divi Fix parameters
#' @param fixparam.growth.fiber Fix parameters of fiber cells
#' @param fixparam.growth.vessel Fix parameters
#' @param dynparam.growth.t Dynamic parameters
#' @param cells cells anatomy
#' @param vessels vessels anatomy
#' @param CZgR if climate limited Cambial cell growth
#'
#' @return today cells layer number & vessels number
#'
## #' @export
#'
#' @importFrom dplyr bind_rows
#'
#'

### 函数需要计算当日的生长速率，输出细胞和导管的数量 checkDiv 和 checkDef 在输入木质部后重置
cell_division <- function( Doy_GR , Res_cambial){   ## Fixp_cambi,, TA

  dt_AllCambial <- Res_cambial ## 形成层细胞生长 ## |> data.table::as.data.table()
  Doy_GR$deltaD <- ifelse( dt_AllCambial$sumDiv == 0,1, Doy_GR$deltaD ) ## error catch
  ###
  # 形成层扩大（生长活动）
  # dt_AllCambial$CA <- (1 + Doy_GR$v_cz) * dt_AllCambial$CA
  #
  # # 判断是否分裂并更新相关列
  # dt_AllCambial$checkDiv <- ifelse(dt_AllCambial$CA >= 2 * dt_AllCambial$OriCA, 1, 0)
  # dt_AllCambial$sumDiv <- ifelse(dt_AllCambial$CA >= 2 * dt_AllCambial$OriCA,
  #                                dt_AllCambial$sumDiv + 1,
  #                                dt_AllCambial$sumDiv)
  #
  # # 生长活动后更新
  # dt_AllCambial$CA <- ifelse(dt_AllCambial$checkDiv == 1,
  #                            dt_AllCambial$OriCA,
  #                            dt_AllCambial$CA)
  # dt_AllCambial$checkDef <- ifelse(dt_AllCambial$dtDef >= Doy_GR$deltaD,
  #                                  dt_AllCambial$dtDef %/% Doy_GR$deltaD,
  #                                  dt_AllCambial$checkDef)
  #
  # # 分裂后更新
  # dt_AllCambial$sumDef <- ifelse(dt_AllCambial$checkDef != 0,
  #                                dt_AllCambial$sumDef + dt_AllCambial$checkDef,
  #                                dt_AllCambial$sumDef)
  # dt_AllCambial$dtDef <- ifelse(dt_AllCambial$checkDef != 0,
  #                               dt_AllCambial$dtDef %% Doy_GR$deltaD,
  #                               dt_AllCambial$dtDef)

  dt_AllCambial[
    , `:=`(CA = (1 + Doy_GR$v_cz) * CA,
           DOY = Doy_GR$DOY )
    ## 形成层扩大（生长活动）
  ][
    , `:=`(checkDiv = data.table::fifelse( CA >= 2*OriCA, 1 , 0 ),## 判断是否分裂
           dtDef = data.table::fifelse( CA >= 2*OriCA, dtDef + 1 , dtDef ),## 分裂则分化计数+1
           sumDiv = data.table::fifelse( CA >= 2*OriCA, sumDiv + 1, sumDiv ) ## 汇总层数。
           )
  ][ ## 生长活动后
    , `:=`(
      CA = data.table::fifelse( checkDiv == 1,OriCA,CA ), ## 如果分裂则重置形成层细胞大小
      checkDef = data.table::fifelse( dtDef >= Doy_GR$deltaD, dtDef %/% Doy_GR$deltaD, checkDef ) )## 检查是否分裂

  ][## 分裂后
    , `:=`(
      sumDef = data.table::fifelse( checkDef != 0 ,sumDef + checkDef , sumDef  ),
      dtDef = data.table::fifelse( checkDef != 0 ,dtDef %% Doy_GR$deltaD, dtDef ) )
  ]


  dt_AllCambial

  return( dt_AllCambial )


} ## cell_division end
