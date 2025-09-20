#### BTR模型主程序 ####
# R包文档通常包括以下几个部分：
#
# 标题和描述：简要描述函数或数据集的作用。
# 参数说明：详细描述每个参数的含义和类型。
# 返回值：描述函数的返回值。
# 示例代码：提供使用函数的示例代码。
# 参考：列出相关的参考资料或函数。
# 内部参数/变量使用 “.”开头驼峰命名，尽量少用缩写,内部函数使用mypkg_开头
# 外部参数/变量蛇形命名"."间隔，外部函数蛇形命名“_”间隔,动-名词方式命名
#
#


#' Simulate tree-ring anatomic by BTR model
#'
#' Main function of Broad-leaved tree-ring model
#'
#' @param clim  气象数据 ，使用computer_clim 的运行结果作为输入
#' @param parameters 模型运行参数，使用data.frame储存
#' @param age 年龄趋势输入
#' @param start.year 模拟起始年份
#' @param end.year 模拟结束年份
#' @param write.result 布尔值 默认为F， 当为TURE时将模拟结果保存为excel文件
#' @param gT.method 温度的相对生长率计算方式(梯形函数或者jonhson的酶促反应方程)
#' @param division 模型对细胞分裂的计算方式，默认为 fix 固定分裂间期，还可以选min或者mix
#' @param Named 当write.result 为TURE时保存的文件名
#' set.layers 模拟生长的最大细胞层数，用于优化计算
#'
#' @return  a excel data
#'
#' @import data.table
#' @importFrom dplyr filter select select left_join bind_rows summarise
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom data.table rbindlist fwrite as.data.table fifelse
#' @importFrom openxlsx write.xlsx
#' @importFrom stats na.omit cor.test setNames
#' @importFrom purrr map2 map
#'
#' @export
#'


### BTR模型主程序 尽量使用base包和向量计算 以优化运行速度
#' 以下为说明文件完成开发后替换
#' 内部参数/变量使用 “.”开头驼峰命名，尽量少用缩写,内部函数使用mypkg_开头
#' 外部参数/变量蛇形命名"."间隔，外部函数蛇形命名“_”间隔,动-名词方式命名
#'
#'

run_btr <- function( clim, parameters, age, start.year = NA,end.year = NA,
                     write.result = FALSE, gT.method = "Jonhson",
                     division = 'fix', Named = NULL ,...
                    ){ ### set.layers = 300
  ## check input data
  # write.result or save.daily 为 TRUE 时 建立文件夹保存数据
  if ( write.result ) {
    .redir <- paste0(  Named,  "res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))
    dir.create(path = eval(.redir) )
  }

  ## error-catching
  start.year[is.na(start.year)] <- max( min(clim$Year), min(age$Year) )
  end.year[is.na(end.year)] <- min( max(clim$Year), max(age$Year) )
  if ( all( c( start.year,end.year ) %in% c( max( min(clim$Year), min(age$Year)  ): min(max(clim$Year) , max(age$Year)) )) == F ) {
    stop( paste( "The simulated time span exceeds the range permitted by the data, please check start.year and end.year!" ) )
  }

  #### 模型计算 ####
  # 拆分参数表,避免同名参数引用错误
  parameters$Values <- as.numeric(parameters$Values) ## 确保为数值数据

  .fixParamDivi <- parameters[ parameters$Module == "CambialActivity" ,
                               c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame() ## 形成层分裂相关参数

  # 提取 Module 为 "CambialActivity" 的行，并选择 Parameter 和 Values 列
  .fixParamDivi <- parameters[parameters$Module == "CambialActivity", c("Parameter", "Values")]
  row.names(.fixParamDivi) <- .fixParamDivi$Parameter# 将 Parameter 列设为行名
  .fixParamDivi$Parameter <- NULL # 删除 Parameter 列（因为已作为行名）
  .fixParamDivi <- as.data.frame(t(.fixParamDivi)) # 转置数据框并转换为 data.frame（此时行变为列）

  .fixParamGrowthFiber <- parameters[parameters$Module == "FiberGrowth", c("Parameter", "Values")]
  row.names(.fixParamGrowthFiber) <- .fixParamGrowthFiber$Parameter
  .fixParamGrowthFiber$Parameter <- NULL
  .fixParamGrowthFiber <- as.data.frame(t(.fixParamGrowthFiber))

  .fixParamGrowthVessel <- parameters[parameters$Module == "VesselGrowth", c("Parameter", "Values")]
  row.names(.fixParamGrowthVessel) <- .fixParamGrowthVessel$Parameter
  .fixParamGrowthVessel$Parameter <- NULL
  .fixParamGrowthVessel <- as.data.frame(t(.fixParamGrowthVessel))

  ## 提取微气候模型参数
  .growthParam <- parameters[parameters$Module == "GrowthRate", c("Parameter", "Values")]
  row.names(.growthParam) <- .growthParam$Parameter
  .growthParam$Parameter <- NULL
  .growthParam <- as.data.frame(t(.growthParam))

  ## 整理年龄趋势
  AgeTrend <- age[age$Year %in% c(syear:eyear), ] |> data.table::as.data.table()
  data.table::setkey(AgeTrend, Year)
  ## 生理 0℃
  temThreshold <- growth_Param$T1

  # 提取对应年份数据并计算微气候
  # system.time({
  microclim <- clim[clim$Year %in% c(syear:eyear),] |>
    data.table::as.data.table() |>
    compute_gR(growth_Param, gTmethod, mode )






} ## end run_btr
