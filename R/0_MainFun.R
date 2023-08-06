#########################
##### Main function #####
#########################

#' Main function calculate tree-ring anatomic
#'
#' @param clim  climate data must have Year，DOY，Temp ，Prec
#' @param parameters model parameter date ，we ues a excel data to save it.
#' @param syear modeling start year
#' @param eyear modeling end year
#' @param intraannual output daily result
#'
#' @return  a excel data
#'
#' @importFrom dplyr filter select select left_join bind_rows
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom data.table rbindlist fwrite
#' @importFrom openxlsx write.xlsx
#' @importFrom stats na.omit
#'
#'
#' @export
#'

btr <- function(  clim, parameters, syear = NA, eyear = NA ,intraannual = F){ ## functions start  ring_width,

  redir <- paste0("res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))


  dir.create(path = eval(redir) )

  ## 提取微气候模型参数

  parameters$values <- as.numeric(parameters$values)

  growth_Param  <- dplyr::filter(parameters, modul == "gR"  ) ## 生长速率阈值参数

  ## error-catching
  if (is.na(syear) ) {syear = min(clim$Year)} else {
    if( syear < min(clim$Year) ){ print("syear is wrong") }
  }

  if (is.na(eyear)) {eyear = max(clim$Year)} else {
    if( eyear > max(clim$Year) ){ print("eyear is wrong") }
  }
  ## 提取对应年份数据
  clim <- dplyr::filter( clim, Year >= syear & Year <= eyear)

  microclim <- Compute_gR(clim , growth_Param)

  ##### 细胞大小计算 ——————

  ## 提取分裂模型参数：

  fixparam.divi <- dplyr::filter(parameters ,  modul == "division" & paramtype == "fixed") %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = 'parameter', value = 'values')

  fixparam.growth.fiber <- dplyr::filter(parameters ,  modul == "growthC" & paramtype == "fixed") %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = 'parameter', value = 'values')

  fixparam.growth.vessel <- dplyr::filter(parameters ,  modul == "growthV" & paramtype == "fixed") %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = 'parameter', value = 'values')

  dynparam.growth.0 <- dplyr::filter(parameters ,  modul == "division" & paramtype == "dynamic") %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = 'parameter', value = 'values')

  ## 年循环计算

  cells <- dplyr::filter(parameters ,  grepl("C0", modul ) ) %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = 'parameter', value = 'values')
  cells[is.na(cells)] <- 0

  vessels <- dplyr::filter(parameters ,  grepl("V0", modul ) ) %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = 'parameter', value = 'values')
  vessels[is.na(vessels)] <- 0

  temThreshold <- dplyr::filter(growth_Param, parameter == "T1")$values

  ##test use
  dailyParameters = logical()
  # tt1 = matrix(NA,ncol = 15,nrow = 366)
  # colnames(tt1) <-c("years","Today",colnames(dynparam.growth.t))
  ##

  summaryYears <- list()


  for (years in syear:eyear) { ## year cycle start  年循环计算





    xylogenesis <- list(cells[-1,],vessels[-1,])
    names(xylogenesis) <- c("cells","vessels")

    gR.year <- dplyr::filter(microclim , Year == years)

    ### 进度条
    print(  paste0("Modeling ", years," / ", "Final ", eyear )    )
    pb <- utils::txtProgressBar(title = "PB",min = 1, max = nrow(gR.year), style = 3, file = stderr())

    ## 年初重置部分参数
    dynparam.growth.t <- dynparam.growth.0
    dynparam.growth.t$Age <- dynparam.growth.t$Age + years - syear
    # dynparam.growth.t$Dd = 0
    # dynparam.growth.t$SumV = 0
    # dynparam.growth.t$SumCL = 1
    # dynparam.growth.t$SumNV = 1

    # j  = 1  ## 正在生长的细胞层数计数
    # k  = 1  ## 正在生长的导管层数计数
    # dcl = 0  ## 完成生长的细胞层数计数
    # dvl = 0  ## 完成生长的导管层数计数
    interimDailyParam <- matrix(NA,ncol = ncol(dynparam.growth.t)+2,nrow = 366)
    colnames(interimDailyParam) <- c("years","Today",colnames(dynparam.growth.t))
    ## 重置结束
    actAccumulatedTem <- 0 ## 活动积温 active accumulated temperature
    # TA <- 0

    summaryDaily <- list()

    for ( Today in 1:nrow(gR.year) ) { ## daily sycle start
      # for (Today in 1:130) {


      # print("")
      # print(c(years,Today))
      # print(dynparam.growth.t)

      # cat( c(years, Today,"\n",as.character(dynparam.growth.t)) , file = "E:\\文档-正在使用\\论文内容\\202203\\code\\11.txt",
      #      sep = " ", fill = T, labels = paste0("{", 1:10, "}:"),append = T)

      clim.today <- gR.year[Today,]

      aatToday <- (clim.today$TEM - temThreshold)
      aatToday[aatToday<0] <- 0
      actAccumulatedTem <- actAccumulatedTem + aatToday
      # print( paste( Today, actAccumulatedTem ,clim.today$rootd )  )

      # dynparam.growth.t <- L_balance(clim.today, dynparam.growth.t ,fixparam.light)

      if (actAccumulatedTem >= fixparam.divi$AAT & clim.today$rootd > 0 ) { ## if start 可以生长时计算细胞分裂分化

          # print( paste( "START",Today    ) )

        newCells <-  cell_division(clim.today, fixparam.divi, fixparam.growth.fiber, fixparam.growth.vessel, dynparam.growth.t,
                             cells, vessels) ### , TVA,TA

        newCells[[1]]  <- stats::na.omit(newCells[[1]])
        newCells[[2]]  <- stats::na.omit(newCells[[2]])
        # TC <- TC + newCells[[4]][1]


        dynparam.growth.t <- newCells[[3]]
        # tt2 <- cbind(years,Today,dynparam.growth.t)
        interimDailyParam[Today,] <- as.numeric(c(years,Today,dynparam.growth.t))
        ##

        if ( nrow(newCells[[1]]) != 0  ) {
          xylogenesis[["cells"]] <- dplyr::bind_rows(xylogenesis[["cells"]],newCells[[1]])
          if (nrow(newCells[[2]])!= 0 ) {
            xylogenesis[["vessels"]] <- dplyr::bind_rows(xylogenesis[["vessels"]],newCells[[2]])
          }
        } ## newCells != 0

        ## 细胞生长部分 #####

        # j <- dcl
        # k <- dvl
        layer.max <- nrow(xylogenesis[["cells"]])
        #
        # nrow(xylogenesis[["cells"]]$DDOY = 0 )
        growthing.cells <- dplyr::filter( xylogenesis[["cells"]] , DDOY == 0  )
        if ( layer.max != 0 & nrow( growthing.cells ) != 0 ) {



          growthing.cells <- cells_growth(cell = growthing.cells, CorV = "C" , clim.today, layer.max,
                            fixparam.growth.fiber, fixparam.growth.vessel, dynparam.growth.t )
          # growthing.cells$DDOY = clim.today$DOY


          rowOfCells <- growthing.cells$cell_L

          xylogenesis[["cells"]][rowOfCells,] <- growthing.cells

        } ## cell end-------------

        layer.max <- nrow(xylogenesis[["vessels"]])
        growthing.cells <- dplyr::filter( xylogenesis[["vessels"]], DDOY == 0  )

        if( layer.max != 0 & nrow(growthing.cells) != 0  ){  ## 有导管时导管生长



          growthing.cells <- cells_growth(cell = growthing.cells, CorV = "V" , clim.today, layer.max,
                            fixparam.growth.fiber,fixparam.growth.vessel, dynparam.growth.t)

          rowOfCells <- growthing.cells$NoV

          xylogenesis[["vessels"]][rowOfCells,] <- growthing.cells

        } ## 有导管时导管生长 end ---------------------


        todayVessels <- xylogenesis[["vessels"]]
        todayFibers <- xylogenesis[["cells"]]


        colnames(todayVessels) <- paste("V", sep = "" , colnames(todayVessels) )

        dailyResult <- dplyr::left_join(todayFibers, todayVessels, by = c("Year" = "VYear", "cell_L" = "Vcell_L") )

        # dailyResult
        #%>% mutate(doy = "Today", .after = "Year")

        #dailyResult$doy[ dailyResult$doy == "Today"] <- Today

        summaryDaily[paste0(Today)] <- list(dailyResult)




      } ## if (actAccumulatedTem >= fixparam.divi$AT) end  可生长日期判定结束 -----------------

      utils::setTxtProgressBar(pb = pb, value = Today)
      # todayVessels <- xylogenesis[["vessels"]]
      # todayFibers <- xylogenesis[["cells"]]
      #
      #
      # colnames(todayVessels) <- paste("V", sep = "" , colnames(todayVessels) )
      #
      # dailyResult = left_join(todayFibers, todayVessels, by = c("Year" = "VYear", "cell_L" = "Vcell_L") )
      #
      # # dailyResult
      # #%>% mutate(doy = "Today", .after = "Year")
      #
      # #dailyResult$doy[ dailyResult$doy == "Today"] <- Today
      #
      # finald[paste0(Today)] <- list(dailyResult)
      #
      # Today  = Today +1
      # Today

    } ## today sycle end  每日生长计算结束 ------------------

    ## 这里需要对年结束进行汇总。

    # colnames(xylogenesis[["vessels"]]) <- paste("V", sep = "" , colnames(xylogenesis[["vessels"]]) )
    #
    # final = left_join(xylogenesis[["cells"]], xylogenesis[["vessels"]], by = c("Year" = "VYear", "cell_L" = "Vcell_L") )
    #
    # orderd = c("Year","cell_L",
    #            "CA" , "CRD","CTD", "CV" ,"CWT", "LWA", "DH", "MORK", "TDOY", "EDOY", "DDOY",
    #            "VVN","VCA","VCRD","VCTD","VCV","VCWT","VLWA","VDH","VMORK","VEDOY","VTDOY","VDDOY")
    #

    dailyParameters <- rbind(dailyParameters,interimDailyParam)
    if (intraannual == T) {
    data.table::rbindlist(summaryDaily,use.names = T,fill = T,idcol = "doy") %>% dplyr::select(c("Year","doy","cell_L",everything())) %>%
        data.table::fwrite(paste0(redir,"/", as.character(years), ".csv" )   )
    }


    summaryYears[[as.character(years)]]  <- dailyResult %>%  dplyr::select(c("Year","cell_L",everything()))

    close(con = pb)


  }  ## year cycle end 年生长计算结束------------------

  ## 合并所有结果

  summaryYears <- data.table::rbindlist(summaryYears)

  parameters$values[ parameters$parameter == "width"]

  annaulRing <-
    summaryYears %>% filter(Year >= syear) %>% group_by(Year) %>%summarise(
    CellLayer = max(cell_L,na.rm = T),
    CellNumber = mean( parameters$values[ parameters$parameter == "width"]  / CTD * CellLayer ,na.rm = T) +
      max(VNoV,na.rm = T) ,
    RingWidth = (mean( parameters$values[ parameters$parameter == "width"] / CTD ,na.rm = T) *
             sum(CA,na.rm = T) + sum(VCA,na.rm = T)) /
             parameters$values[ parameters$parameter == "width"] / 1000,
    VesselNumber = max(VNoV,na.rm = T )

  )


  Outputs <- list(annaulRing, summaryYears , microclim, as.data.frame( stats::na.omit(dailyParameters) ) )

  names(Outputs) <- c('annaulRing','xylem_trait', 'microclim', 'dailyParameters' )


  openxlsx::write.xlsx(Outputs,paste0( redir, "/Outputs.xlsx" ) )

  openxlsx::write.xlsx(parameters,file = paste0(redir,"/Parameters.xlsx" ))


  print( redir)

  return(Outputs)

} ##  TRXXX2 funtion end -----------------------------------------
