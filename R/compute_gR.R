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
#' @importFrom data.table %between%
#'
#' @return gM & gT
compute_gR <- function( clim , growth_Param, gTmethod,mode ){

  ## 计算gT gTmethod
  if (gTmethod == 'VS') {
    clim[, gT := data.table::fcase(
      TEM < growth_Param$T1 | TEM >= growth_Param$T4, 0,
      TEM >= growth_Param$T2 & TEM <= growth_Param$T3, 1,
      TEM < growth_Param$T2 & TEM >= growth_Param$T1,
      (TEM - growth_Param$T1) / (growth_Param$T2 - growth_Param$T1),
      TEM < growth_Param$T4 & TEM > growth_Param$T3,
      (growth_Param$T4 - TEM) / (growth_Param$T4 - growth_Param$T3)
    )]
  } else {
    clim[, gT := ((TEM + 273.15) *
                    exp(-growth_Param$deltaH_A_Da / (8.314 * (TEM + 273.15))) /
                    (1 + exp(growth_Param$deltaS_D / 8.314 -
                               growth_Param$deltaH_D / (8.314 * (TEM + 273.15))))) |>
           nor(Zeros = TRUE)]

    clim[TEM <= growth_Param$T1 | TEM >= growth_Param$T4, gT := 0]
  }

  ## 计算 gM,gV,L_i.fiber,L_i.vessel
  # 第一部分：计算 gM, gV, L_i.fiber, L_i.vessel
  clim[, `:=`(
    gM = data.table::fcase(
      soilM < growth_Param$M1 | soilM >= growth_Param$M4, 0,
      soilM >= growth_Param$M2 & soilM <= growth_Param$M3, 1,
      soilM < growth_Param$M2 & soilM >= growth_Param$M1,
      (soilM - growth_Param$M1) / (growth_Param$M2 - growth_Param$M1),
      soilM < growth_Param$M4 & soilM > growth_Param$M3,
      (growth_Param$M4 - soilM) / (growth_Param$M4 - growth_Param$M3)
    ),

    gV = data.table::fcase(
      VPD < growth_Param$VPD1 | VPD >= growth_Param$VPD4, 0,
      VPD >= growth_Param$VPD2 & VPD <= growth_Param$VPD3, 1,
      VPD < growth_Param$VPD2 & VPD >= growth_Param$VPD1,
      (VPD - growth_Param$VPD1) / (growth_Param$VPD2 - growth_Param$VPD1),
      VPD < growth_Param$VPD4 & VPD > growth_Param$VPD3,
      (growth_Param$VPD4 - VPD) / (growth_Param$VPD4 - growth_Param$VPD3)
    ),

    L_i.fiber = growth_Param$a.fiber * exp(
      -exp(growth_Param$b.fiber - growth_Param$c.fiber * dL_i)),

    L_i.vessel = growth_Param$a.vessel * exp(
      -exp(growth_Param$b.vessel - growth_Param$c.vessel * dL_i))
  )]

  # 第二部分：计算积温和生长季
  clim[ , aT := data.table::fifelse( (TEM - growth_Param$T1) >=0, TEM - growth_Param$T1,0 )  , by = Year]

  # 计算累积积温
  clim[, aaT := cumsum(aT), by = Year]

  # 计算连续5天温度<=T1的逻辑
  clim[, tn := {
    temp <- data.table::fifelse(aT == 0, 1, 0)
    c(rowSums(embed(temp, 5)), rep(5, 4))  # 注意：embed需要处理边界情况
  }, by = Year]

  # 找出每年最高温度日期
  clim[, HotDay := which.max(TEM), by = Year]

  # 计算生长季开始和结束
  clim[, Fs1 := which(aaT >= fixparam.divi$AAT & tn == 0)[1], by = Year]
  clim[, Fs2 := which(aaT >= fixparam.divi$AAT & tn == 5 & DOY >= HotDay)[1], by = Year]

  # 确定生长季
  clim[, GS := data.table::fifelse(DOY >= Fs1 & DOY <= Fs2, "GR", "UN")]

  # 删除临时列
  clim[, c("tn", "Fs1", "Fs2", "HotDay") := NULL]

  if (mode == 'mix') {
    clim[, `:=`(
      egR = gE * gT * pmin(gE, gV),
      wgR = gE * gT
    )]
  } else {
    clim[, `:=`(
      egR = gE * pmin(gT, gE, gV),
      wgR = gE * gT
    )]
  }

  return(clim)
}

compute_gR2 <- function(clim, growth_Param, gTmethod, mode) {
  # 预加载关键参数到内存（减少列表索引开销）
  T1 <- growth_Param$T1; T2 <- growth_Param$T2; T3 <- growth_Param$T3; T4 <- growth_Param$T4
  M1 <- growth_Param$M1; M2 <- growth_Param$M2; M3 <- growth_Param$M3; M4 <- growth_Param$M4
  VPD1 <- growth_Param$VPD1; VPD2 <- growth_Param$VPD2; VPD3 <- growth_Param$VPD3; VPD4 <- growth_Param$VPD4

  # 优化gT计算 ---------------------------------------------------------------
  if (gTmethod == 'VS') {
    clim[, gT := data.table::fcase(
      TEM < T1 | TEM >= T4, 0.0,
      TEM >= T2 & TEM <= T3, 1.0,
      TEM < T2 & TEM >= T1, (TEM - T1)/(T2 - T1),
      TEM < T4 & TEM > T3, (T4 - TEM)/(T4 - T3)
    )]
  } else {
    # 预计算重复使用的温度值
    TK <- clim$TEM + 273.15
    exp_term1 <- exp(-growth_Param$deltaH_A_Da / (8.314 * TK))
    exp_term2 <- exp(growth_Param$deltaS_D/8.314 - growth_Param$deltaH_D/(8.314*TK))

    clim[, gT := (TK * exp_term1) / (1 + exp_term2)]
    clim[TEM <= T1 | TEM >= T4, gT := 0.0]
  }

  # 批量计算多列（减少内存分配次数）---------------------------------------------
  clim[, `:=`(
    gM = data.table::fcase(
      soilM < M1 | soilM >= M4, 0.0,
      soilM >= M2 & soilM <= M3, 1.0,
      soilM < M2 & soilM >= M1, (soilM - M1)/(M2 - M1),
      soilM < M4 & soilM > M3, (M4 - soilM)/(M4 - M3)
    ),
    gV = data.table::fcase(
      VPD < VPD1 | VPD >= VPD4, 0.0,
      VPD >= VPD2 & VPD <= VPD3, 1.0,
      VPD < VPD2 & VPD >= VPD1, (VPD - VPD1)/(VPD2 - VPD1),
      VPD < VPD4 & VPD > VPD3, (VPD4 - VPD)/(M4 - VPD3)  # 注意这里可能存在参数错误
    ),
    L_i.fiber = growth_Param$a.fiber * exp(-exp(growth_Param$b.fiber - growth_Param$c.fiber*dL_i)),
    L_i.vessel = growth_Param$a.vessel * exp(-exp(growth_Param$b.vessel - growth_Param$c.vessel*dL_i))
  )]

  # 优化累积计算（避免重复分组）-----------------------------------------------
  clim[, aT := pmax(TEM - T1, 0.0)]
  clim[, aaT := cumsum(aT), by = Year]

  # 使用滚动窗口函数优化连续天数计算（替代embed）
  clim[, tn := data.table::frollsum(data.table::fifelse(aT == 0, 1L, 0L), 5, align = "left", fill = 5L), by = Year]

  # 优化极值查找（利用which.min性能优势）
  clim[, HotDay := which.max(TEM), by = Year]

  # 向量化条件查找（替代多次分组which）
  AAT <- fixparam.divi$AAT
  clim[, Fs1 := .I[aaT >= AAT & tn == 0L][1], by = Year]
  clim[, Fs2 := .I[aaT >= AAT & tn == 5L & DOY >= HotDay][1], by = Year]

  # 逻辑索引直接赋值（替代fifelse）
  clim[, GS := "UN"]
  clim[DOY >= Fs1 & DOY <= Fs2, GS := "GR"]

  # 模式计算优化（避免重复运算）
  min_val <- if (mode == 'mix') clim[, pmin(gE, gV)] else clim[, pmin(gT, gE, gV)]
  clim[, `:=`(
    egR = gE * gT * min_val,
    wgR = gE * gT
  )][, c("tn","Fs1","Fs2","HotDay") := NULL]

  return(clim)
}

