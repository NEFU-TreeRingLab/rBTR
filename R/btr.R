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
#' @param clim  climate data must have Year，DOY，Temp ，Prec
#' @param parameters model parameter date ，we ues a excel data to save it.
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

run_btr <- function( clim, paramters, age, start.year = NA,end.year = NA,
                     write.result = FALSE, save.daily = FALSE,
                     gR.method = "Jonhson",division = 'fix', Named = NULL ,
                     set.layers = 300   ){
  ## check input data
  # write.result or save.daily 为 TRUE 时 建立文件夹保存数据
  if ( any( write.result, save.daily ) ) {
    redir <- paste0(  Named,  "res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))
    dir.create(path = eval(redir) )
  }

  ## error-catching
  syear[is.na(syear)] <- max( min(clim$Year), min(age$Year) )
  eyear[is.na(eyear)] <- min( max(clim$Year), max(age$Year) )
  if ( all( c( syear,eyear ) %in% c( max( min(clim$Year), min(age$Year)  ): min(max(clim$Year) , max(age$Year)) )) == F ) {
    stop( paste( "The simulated time span exceeds the range permitted by the data, please check syear and eyear!" ) )
  }






} ## end run_btr
