#### BTR模型主程序 ####
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
#' @param mode 气象因子组合模式 "min","multiply".
#' @param Dcase How to calculate Division Limit, including: "min","mean","multiply".
#' @param Named Named save document
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
btr <- function(  clim, parameters, age, syear = NA, eyear = NA, Cores = 1,
                  writeRes = F, intraannual = F, gTmethod = "Jonhson", mode = 'mix',
                  Dcase = "min",Tradeoff = 'Ks', Named = NULL ) { ## functions start  ring_width,
  #### 运行前检查 ####
  # writeRes intraannual 任意为 TRUE 时 建立文件夹保存数据
  if ( any( writeRes,intraannual ) ) {
    redir <- paste0(  Named,  "res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))
    dir.create(path = eval(redir) )
  }

  ## error-catching
  syear[is.na(syear)] <- max( min(clim$Year), min(age$Year) )
  eyear[is.na(eyear)] <- min( max(clim$Year), max(age$Year) )


  if ( all( c( syear,eyear ) %in% c( max( min(clim$Year), min(age$Year)  ): min(max(clim$Year) , max(age$Year)) )) == F ) {
     stop( paste( "The simulated time span exceeds the range permitted by the data, please check syear and eyear!" ) )
  }

  #### 模型计算 ####

  # 拆分参数表
  #### system.time({ }) ####
  parameters$Values <- as.numeric(parameters$Values)

  fixparam.divi <- parameters[ parameters$Module == "CambialActivity" ,
                               c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame()

  fixparam.growth.fiber <- parameters[ parameters$Module == "FiberGrowth" ,
                                       c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame()

  fixparam.growth.vessel <- parameters[ parameters$Module == "VesselGrowth" ,
                                        c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame()

  ## 提取微气候模型参数
  growth_Param <- parameters[ parameters$Module == "GrowthRate" ,
                              c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame()

  # age <- age[age$Year %in% c(syear:eyear), ]

  age <- data.table::as.data.table(age[age$Year %in% c(syear:eyear), ])
  data.table::setkey(age, Year)

  temThreshold <- growth_Param$T1

  # 提取对应年份数据并计算微气候


  # system.time({
  clim1 <- data.table::as.data.table(clim)
  clim2 <- compute_gR(clim1[Year %in% c(syear:eyear),] , growth_Param, gTmethod, mode )
  # })

  #  system.time({
  #   clim4 <- data.table::as.data.table(clim)
  #   clim3 <- compute_gR2(clim = clim4[Year %in% c(syear:eyear),] , growth_Param, gTmethod, mode )
  # })

  # 建立表 动态参数、初始值
  params <- c('dCA_cz','v_c.fiber','v_w.fiber','v_l.fiber',
              'v_c.vessel','v_w.vessel','v_l.vessel','dVN_cz','Vcz')

  # 单步完成：创建模板 + 合并 + 添加年龄列
  dynparam.growth.0 <- clim[GS != 'UN', list('Year', 'DOY')][,
                                                      c(params) := 0L ][,  # 直接添加参数列并初始化为0
                                                                        Age := age['syear', age, on = 'Year']  # 使用快速合并提取年龄
                                                      ]


  # Create initial data.frame as data.table
  dynparam.growth.0 <- data.table::data.table(
    Parameter = c('dCA_cz','v_c.fiber','v_w.fiber','v_l.fiber',
                  'v_c.vessel','v_w.vessel','v_l.vessel',
                  'dVN_cz','Vcz'),
    Values = 0
  )

  # Spread the data (using dcast instead of tidyr::spread)
  dynparam.growth.0 <- data.table::dcast(dynparam.growth.0,
                             . ~ Parameter,
                             value.var = "Values")[, . := NULL]  # Remove the dummy column

  # Bind with the clim data (using cbind instead of bind_cols)
  dynparam.growth.0 <- cbind(dynparam.growth.0,
                             clim[GS != 'UN', .(Year, DOY)])

  # Add age column
  dynparam.growth.0[, Age := age$age[age$Year == syear]]

  ## 年循环计算开始前
  ## 计算 权衡的年内平衡曲线
  TradeLine <- fixparam.divi$maxRCTA  *
    nor( microclim$L_i.vessel[microclim$Year == syear ] *-1 )

  ## 20240803 使用L+2 替换 L+1 ：this is L116
  # RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <- 99
  RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <-
    fixparam.divi$maxRCTA * fixparam.divi$RCTADivT
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

} ## btr end ------------






