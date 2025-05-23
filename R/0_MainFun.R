##########################
##### Main function1 #####
##########################

#' Main function calculate tree-ring anatomic
#'
#' @param clim  climate data must have Year，DOY，Temp ，Prec
#' @param parameters model parameter date ，we ues a excel data to save it.
#' @param syear modeling start year
#' @param eyear modeling end year
#' @param Cores 核心数
#' @param writeRes auto write res_data, logical value, if 'Ture', output result list as xlsx in file.
#' @param intraannual output daily result, logical, if 'Ture', output daily result in file.
#' @param gTmethod method of calculate gT ,'VS' and 'Jonhson'
#'
#' @param division dealtDivision 'fix' is constant, 'limit' is gM/gV & RCTA control.
#' @param testLim use Lage limit mw & ml , Lage in age param_data, regulation cell enlargement by CAmax param .
#' @param CZgR growth limit : climate-Cambial Cell, cell wall-gM, cell wall-gV, NewgRtype; 0 for False, 1 for True
#' @param Pbar Whether to show progress bar
#' @param testMod test Mode, logical, if 'False', Lage&Dage use ratio
#' @param Dcase How to calculate Division Limit, including: "min","mean","multiply".
#'
#' @return  a excel data
#'
#' @importFrom dplyr filter select select left_join bind_rows summarise
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom data.table rbindlist fwrite as.data.table
#' @importFrom openxlsx write.xlsx
#' @importFrom stats na.omit cor.test setNames
#' @importFrom purrr map2 map
#' @importFrom parabar start_backend par_lapply
#'
#' @export
#'
btr <- function(  clim, parameters, age, syear = NA, eyear = NA ,
                  writeRes = F, intraannual = F, gTmethod = "Jonhson" , division = "limit" ,
                  testLim = F,
                  CZgR = c(1,0,0,1 ) , Pbar = F , testMod = F ,Dcase = "min", Named = NULL ) { ## functions start  ring_width,

  writeRes[intraannual == T] <- T

  if (writeRes == T) {
    redir <- paste0(  Named,  "res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))

    dir.create(path = eval(redir) )
  }
  ## 提取分裂模型参数：

  fixparam.divi <- parameters[ parameters$Module == "CambialActivity" ,
                               c("Parameter","Values") ]   |>
    tidyr::spread( key = 'Parameter', value = 'Values')

  fixparam.growth.fiber <- parameters[parameters$Module == "FiberGrowth" ,
                                      c("Parameter","Values")]  |>
    tidyr::spread( key = 'Parameter', value = 'Values')

  fixparam.growth.vessel <- parameters[parameters$Module == "VesselGrowth" ,
                                       c("Parameter","Values") ] |>
    tidyr::spread( key = 'Parameter', value = 'Values')
  ## 提取微气候模型参数
  parameters$Values <- as.numeric(parameters$Values)
  growth_Param  <- parameters[parameters$Module == "GrowthRate",] ## 生长速率阈值参数

  ## error-catching
  if (is.na(syear) ) {syear = max( min(clim$Year), min(age$Year)  )   } else {
    if( syear < max( min(clim$Year), min(age$Year)  )   ){ stop( paste("syear is wrong:", syear,"ClimY:", min(clim$Year),"AgeY:" ,min(age$Year) ) ) }
  }

  if (is.na(eyear)) {eyear = min(max(clim$Year) , max(age$Year))      } else {
    if( eyear > min(max(clim$Year) , max(age$Year)) ){ stop( paste(  "eyear is wrong:",eyear, "ClimY:",max(clim$Year), "AgeY:",max(age$Year)  ) ) }
  }

  ## 提取对应年份数据
  clim <- clim[clim$Year %in% c(syear:eyear), ]
  age <- age[age$Year %in% c(syear:eyear), ]

  temThreshold <- dplyr::filter(growth_Param, Parameter == "T1")$Values

  if ( all( gTmethod == "VS" ) ) {
    microclim <- Compute_gR(clim , growth_Param) |>
      dplyr::left_join(age,by="Year")  # |>
    # data.table::as.data.table() ##VS
  } else { microclim <- Compute_gR2(clim , growth_Param) |>
    dplyr::left_join(age,by="Year") # |>
  # data.table::as.data.table()
  } ## 'Jonhson'

  ## 计算有效积温和气候生长季长度
  ## 生长季开始::积温 > 阈值 & 第一次连续5天温度 > T1
  ## 生长季结束:: 最高温度日期以后 & 积温 > 阈值 & 第一次连续5天温度 < T1

  microclim <- microclim |> dplyr::group_by(Year) |>
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

  ##### dynparam 未来动参在模型内指定，移除出参数表
  # dynparam.growth.0 <- parameters[parameters$modul == "division" & parameters$paramtype == "dynamic" ,
  #                                 c("parameter","values") ] |>
  #   tidyr::spread( key = 'parameter', value = 'values')
  dynparam.growth.0 <- data.frame( Parameter = c('L_i.fiber','L_i.vessel','dCA_cz','SumCL','SumVL',
                                                 'SumV','v_c.fiber','v_w.fiber','v_l.fiber',
                                                 'v_c.vessel','v_w.vessel','v_l.vessel',
                                                 'deltaVN','Vcz','grwothSeason','Age','T_age','czgR','egR', 'wgR' ),
                                   Values = 0 ) |>
    tidyr::spread( key = 'Parameter', value = 'Values')

  fixparam.growth.origin <- parameters[parameters$ParamType == "FixedParam",
                                       c("Parameter","Module", "Values")]

  dynparam.growth.0$Age <- age$age[ age$Year == syear] ## error catch
  ## 年循环计算开始前
  ## RCTA line 计算 RCTA 平衡曲线
  RCTA <- fixparam.divi$maxRCTA  *
    nor( microclim$L_i.vessel[microclim$Year == syear ] *-1 )
  ## 20240803 使用L+2 替换 L+1 ：this is L116
  # RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <- 99
  RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <- fixparam.divi$maxRCTA * fixparam.divi$RCTADivT
  ## __end ----

  ## 设置各类初始值： 纤维细胞和导管初始值

  Cambial <- matrix(NA,ncol = 12, nrow = 1,
                    dimnames = list(c("1"),c('Year','cell_L','CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY' )) ) |>
    as.data.frame()

  Cambial$CV  <- parameters$Values[ parameters$Parameter == 'CV' ]
  Cambial$WT  <- parameters$Values[ parameters$Parameter == 'WT' ]
  Cambial$CTD <- parameters$Values[ parameters$Parameter == 'CTD' ]
  Cambial$CRD <- Cambial$CV/(Cambial$CTD - 2*Cambial$WT   ) + 2*Cambial$WT
  Cambial$CA  <- Cambial$CTD*Cambial$CRD
  Cambial$WA  <- Cambial$CA-Cambial$CV

  cells <- Cambial
  cells[is.na(cells)] <- 0

  vessels <- Cambial |> dplyr::mutate( NoV = NA , VN = NA   )
  vessels[is.na(vessels)] <- 0


  ## 年循环 #####
  if (syear ==  eyear ) {
    pb <- utils::txtProgressBar(title = "PB",min = syear, max = (eyear+1), style = 3, file = stderr())
  } else {
     pb <- utils::txtProgressBar(title = "PB",min = syear, max = eyear, style = 3, file = stderr())
  }


  yRes <- list()


  for (years in syear:eyear) { ## year cycle start  年循环计算 准备改并行  # years = syear

    yRes[[ as.character(years) ]] <- year_growth( x = years,
                                                  microclim ,testMod,testLim, intraannual, writeRes , division,Dcase, CZgR,
                                                  cells, vessels,
                                                  dynparam.growth.0,RCTA,
                                                  fixparam.divi,fixparam.growth.fiber, fixparam.growth.origin, fixparam.growth.vessel
    )

    utils::setTxtProgressBar(pb = pb, value = years)

  }  ## year cycle end 年生长计算结束------------------

  close(con = pb)


  ## 合并所有结果;
  Cells <- list()
  IntraAnnual <- list()
  # InterAnnual <- list()
  dailyParam <- list()

  for (ys in 1:length(yRes)) {
    Cells[[ as.character(ys)  ]] <- yRes[[ys]]$summaryYears
    IntraAnnual[[ as.character(ys) ]] <- yRes[[ys]]$AnnualGrowth
    dailyParam[[ as.character(ys) ]] <- yRes[[ys]]$dailyParameters
  }

  AnnualGrowth <- IntraAnnual |>
    data.table::rbindlist() |>
    dplyr::rename(VesselFraction = RCTA)

  AnnaulRing <- AnnualGrowth[,1:16] |>
    dplyr::group_by(Year) |>
    dplyr::mutate( StartDoy = min(DOY),.before = DOY )|>
    dplyr::filter(DOY == max(DOY)) |>
    dplyr::rename(EndDoy = DOY)

  summaryYears <- data.table::rbindlist(Cells)

  dailyParameters <- data.table::rbindlist(dailyParam)

  Outputs <- list(annaulRing = AnnaulRing,
                  xylem_trait = summaryYears , IntraAnnualGrowth = AnnualGrowth,
                  microclim = microclim, dailyParameters =  dailyParameters ) |>
    dev.outputs()

  resp <- list(parameters = parameters,age = age)

  if (writeRes == T) {
    openxlsx::write.xlsx(Outputs,paste0( redir, "/Outputs.xlsx" ) )
    openxlsx::write.xlsx(resp,file = paste0(redir,"/Parameters.xlsx" ))
    print( redir)
  }

  return(Outputs)

} ##  TRXXX2 funtion end -----------------------------------------












#########################@
##### Main function2 #####
#########################@

#' Main function calculate tree-ring anatomic
#'
#' @param clim  climate data must have Year，DOY，Temp ，Prec
#' @param parameters model parameter date ，we ues a excel data to save it.
#' @param syear modeling start year
#' @param eyear modeling end year
#' @param Cores 核心数
#' @param writeRes auto write res_data, logical value, if 'Ture', output result list as xlsx in file.
#' @param intraannual output daily result, logical, if 'Ture', output daily result in file.
#' @param gTmethod method of calculate gT ,'VS' and 'Jonhson'
#' @param division dealtDivision 'fix' is constant, 'limit' is gM/gV & RCTA control.
#' @param testLim use Lage limit mw & ml , Lage in age param_data, regulation cell enlargement by CAmax param .
#' @param CZgR growth limit : climate-Cambial Cell, cell wall-gM, cell wall-gV, NewgRtype; 0 for False, 1 for True
#' @param Pbar Whether to show progress bar
#' @param testMod test Mode, logical, if 'False', Lage&Dage use ratio
#' @param Dcase How to calculate Division Limit, including: "min","mean","multiply".
#'
#' @return  a excel data
#'
#' @importFrom dplyr filter select select left_join bind_rows summarise
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom data.table rbindlist fwrite as.data.table
#' @importFrom openxlsx write.xlsx
#' @importFrom stats na.omit
#' @importFrom purrr map2 map
#' @importFrom parabar start_backend par_lapply set_option stop_backend
## #' @importFrom . year_growth
#' @export
#'
btr_parallel <- function(  clim, parameters, age, syear = NA, eyear = NA ,Cores = 5,
                           writeRes = F,intraannual = F, gTmethod = "Jonhson" , division = "limit" ,
                           testLim = F,
                           CZgR = c(1,0,0,1 ) , Pbar = T , testMod = F ,Dcase = "min", Named = NULL ) { ## functions start  ring_width,
  ## 检查变量名


  ##

  ##
  writeRes[intraannual == T] <- T

  if (writeRes == T) {
    redir <- paste0(  Named,  "res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))

    dir.create(path = eval(redir) )
  }
  ## 提取分裂模型参数：

  fixparam.divi <- parameters[ parameters$Module == "CambialActivity" ,
                               c("Parameter","Values") ]   |>
    tidyr::spread( key = 'Parameter', value = 'Values')

  fixparam.growth.fiber <- parameters[parameters$Module == "FiberGrowth" ,
                                      c("Parameter","Values")]  |>
    tidyr::spread( key = 'Parameter', value = 'Values')

  fixparam.growth.vessel <- parameters[parameters$Module == "VesselGrowth" ,
                                       c("Parameter","Values") ] |>
    tidyr::spread( key = 'Parameter', value = 'Values')
  ## 提取微气候模型参数
  parameters$Values <- as.numeric(parameters$Values)
  growth_Param  <- parameters[parameters$Module == "GrowthRate",] ## 生长速率阈值参数

  ## error-catching
  if (is.na(syear) ) {syear = max( min(clim$Year), min(age$Year)  )   } else {
    if( syear < max( min(clim$Year), min(age$Year)  )   ){ stop( paste("syear is wrong:", syear,"ClimY:", min(clim$Year),"AgeY:" ,min(age$Year) ) ) }
  }

  if (is.na(eyear)) {eyear = min(max(clim$Year) , max(age$Year))      } else {
    if( eyear > min(max(clim$Year) , max(age$Year)) ){ stop( paste(  "eyear is wrong:",eyear, "ClimY:",max(clim$Year), "AgeY:",max(age$Year)  ) ) }
  }

  ## 提取对应年份数据
  clim <- clim[clim$Year %in% c(syear:eyear), ]
  age <- age[age$Year %in% c(syear:eyear), ]

  temThreshold <- dplyr::filter(growth_Param, Parameter == "T1")$Values

  if ( all( gTmethod == "VS" ) ) {
    microclim <- Compute_gR(clim , growth_Param) |>
      dplyr::left_join(age,by="Year")  # |>
    # data.table::as.data.table() ##VS
  } else { microclim <- Compute_gR2(clim , growth_Param) |>
    dplyr::left_join(age,by="Year") # |>
  # data.table::as.data.table()
  } ## 'Jonhson'

  ## 计算有效积温和气候生长季长度
  ## 生长季开始::积温 > 阈值 & 第一次连续5天温度 > T1
  ## 生长季结束:: 最高温度日期以后 & 积温 > 阈值 & 第一次连续5天温度 < T1

  microclim <- microclim |> dplyr::group_by(Year) |>
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

  ##### dynparam 未来动参在模型内指定，移除出参数表
  # dynparam.growth.0 <- parameters[parameters$modul == "division" & parameters$paramtype == "dynamic" ,
  #                                 c("parameter","values") ] |>
  #   tidyr::spread( key = 'parameter', value = 'values')
  dynparam.growth.0 <- data.frame( Parameter = c('L_i.fiber','L_i.vessel','dCA_cz','SumCL','SumVL',
                                                 'SumV','v_c.fiber','v_w.fiber','v_l.fiber',
                                                 'v_c.vessel','v_w.vessel','v_l.vessel',
                                                 'deltaVN','Vcz','grwothSeason','Age','T_age','czgR','egR', 'wgR' ),
                                   Values = 0 ) |>
    tidyr::spread( key = 'Parameter', value = 'Values')

  fixparam.growth.origin <- parameters[parameters$ParamType == "FixedParam",
                                       c("Parameter","Module", "Values")]

  dynparam.growth.0$Age <- age$age[ age$Year == syear] ## error catch
  ## 年循环计算开始前
  ## RCTA line 计算 RCTA 平衡曲线
  RCTA <- fixparam.divi$maxRCTA  *
    nor( microclim$L_i.vessel[microclim$Year == syear ] *-1 )
  ## 20240803 使用L+2 替换 L+1 ：this is L116
  # RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <- 99
  RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <- fixparam.divi$maxRCTA * fixparam.divi$RCTADivT
  ## __end ----

  ## 设置各类初始值： 纤维细胞和导管初始值

  Cambial <- matrix(NA,ncol = 12, nrow = 1,
                    dimnames = list(c("1"),c('Year','cell_L','CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY' )) ) |>
    as.data.frame()

  Cambial$CV  <- parameters$Values[ parameters$Parameter == 'CV' ]
  Cambial$WT  <- parameters$Values[ parameters$Parameter == 'WT' ]
  Cambial$CTD <- parameters$Values[ parameters$Parameter == 'CTD' ]
  Cambial$CRD <- Cambial$CV/(Cambial$CTD - 2*Cambial$WT   ) + 2*Cambial$WT
  Cambial$CA  <- Cambial$CTD*Cambial$CRD
  Cambial$WA  <- Cambial$CA-Cambial$CV

  cells <- Cambial
  cells[is.na(cells)] <- 0

  vessels <- Cambial |> dplyr::mutate( NoV = NA , VN = NA   )
  vessels[is.na(vessels)] <- 0

  ## 年循 #####

  if (Pbar == T) {
    parabar::set_option("progress_track", T)
  } else {parabar::set_option("progress_track", T)}

  # ifelse( Pbar == T , parabar::set_option("progress_track", T), parabar::set_option("progress_track", FALSE))

  backend <- parabar::start_backend(cores = Cores, cluster_type = "psock", backend_type = "async")

  yRes <-
    parabar::par_lapply(backend, x = syear:eyear,
                        fun = rBTR:::year_growth ,
                        microclim ,testMod,testLim, intraannual, writeRes ,division,Dcase, CZgR,
                        cells, vessels,
                        dynparam.growth.0,RCTA,
                        fixparam.divi,fixparam.growth.fiber, fixparam.growth.origin, fixparam.growth.vessel  )

  parabar::stop_backend(backend)

  ## 合并所有结果;
  Cells <- list()
  IntraAnnual <- list()
  # InterAnnual <- list()
  dailyParam <- list()

  for (ys in 1:length(yRes)) {
    Cells[[ as.character(ys)  ]] <- yRes[[ys]]$summaryYears
    IntraAnnual[[ as.character(ys) ]] <- yRes[[ys]]$AnnualGrowth
    dailyParam[[ as.character(ys) ]] <- yRes[[ys]]$dailyParameters
  }

  AnnualGrowth <- IntraAnnual |>
    data.table::rbindlist() |>
    dplyr::rename(VesselFraction = RCTA)

  AnnaulRing <- AnnualGrowth[,1:16] |>
    dplyr::group_by(Year) |>
    dplyr::mutate( StartDoy = min(DOY),.before = DOY )|>
    dplyr::filter(DOY == max(DOY)) |>
    dplyr::rename(EndDoy = DOY)

  summaryYears <- data.table::rbindlist(Cells)

  dailyParameters <- data.table::rbindlist(dailyParam)

  Outputs <- list(annaulRing = AnnaulRing,
                  xylem_trait = summaryYears , IntraAnnualGrowth = AnnualGrowth,
                  microclim = microclim, dailyParameters =  dailyParameters ) |>
    dev.outputs()

  resp <- list(parameters = parameters,age = age)

  if (writeRes == T) {
    openxlsx::write.xlsx(Outputs,paste0( redir, "/Outputs.xlsx" ) )
    openxlsx::write.xlsx(resp,file = paste0(redir,"/Parameters.xlsx" ))
    print( redir)
  }

  return(Outputs)

} ##  TRXXX2 funtion end -----------------------------------------
