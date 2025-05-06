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
#' @param mode test Mode, logical, if 'False', Lage&Dage use ratio
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
                  Dcase = "min", Named = NULL ) { ## functions start  ring_width,
  #### 运行前检查 ####
  # writeRes intraannual 任意为 TRUE 时 建立文件夹保存数据
  if ( any( writeRes,intraannual ) ) {
    redir <- paste0(  Named,  "res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))
    dir.create(path = eval(redir) )
  }

  ## error-catching
  syear[is.na(syear)] <- max( min(clim$Year), min(age$Year) )
  eyear[is.na(eyear)] <- min( max(clim$Year), max(age$Year) )

  if ( all( c( syear,eyear ) %in% c( max( min(clim$Year), min(age$Year)  ): min(max(clim$Year) , max(age$Year))   )) == F ) {
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

  ## 提取对应年份数据
  clim <- clim[clim$Year %in% c(syear:eyear), ]
  age <- age[age$Year %in% c(syear:eyear), ]

  temThreshold <- growth_Param$T1

  microclim <- compute_gR(clim , growth_Param, gTmethod )



} ## btr end ------------






