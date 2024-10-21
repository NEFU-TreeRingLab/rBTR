#############################
#### year Groeth modual ####
#############################


### return ：生长完成的细胞，日生长汇总表， 日动参表 ，选择write，日细胞表

## interFun of Daily growth
#' yearGrowth
#' @param years 年 输入
#' @param microclim climate data
#' @param testMod locic testMode
#' @param testLim loigc Li
#' @param intraannual logic output intraannual cells
#' @param writeRes logic write results
#' @param division 分裂
#' @param Dcase 忘了1
#' @param CZgR 忘了2
#'
#' @param cells fiber cells t0
#' @param vessels vessel t0
#' @param dynparam.growth.0 dynamic param 0
#' @param dynparam.growth.t dynmaic param t
#'
#' @param fixparam.divi Fix parameters of cell division
#' @param RCTA 标准RCTA
#' @param fixparam.growth.origin Fis param
#' @param fixparam.growth.fiber Fix parameters of fiber cells
#' @param fixparam.growth.vessel Fix parameters
#'
#' @importFrom dplyr filter select select left_join bind_rows summarise mutate
#' @importFrom magrittr %>%
#' @importFrom data.table rbindlist fwrite as.data.table
#' @importFrom openxlsx write.xlsx
#' @importFrom stats na.omit setNames cor.test
#' @importFrom purrr map2 map
## #' @importFrom . cell_division daily_grwoth
#'
#' @return today cells anatomy
#'
## #' @export

year_growth <- function( x, microclim ,testMod,testLim,intraannual, writeRes , division,Dcase, CZgR,
                         cells, vessels,
                         dynparam.growth.0,RCTA,
                         fixparam.divi,fixparam.growth.fiber, fixparam.growth.origin, fixparam.growth.vessel  ){
  years <- x
  ## 数据初始化
  RCTAt <- 0 ### RCTA for dD
  gR.year <- microclim[microclim$Year == years,]


  ## check max VCA
  if ( !is.na( unique(gR.year$Lage))  ) {
    ## VCA - age

    if( testMod == T ){

      if (testLim == T ) { # 先不改，不要用
        fixparam.growth.vessel$ml <- gR.year$Lage[1] /fixparam.growth.vessel$CAmax * fixparam.growth.vessel$ml
        fixparam.growth.vessel$mw <- gR.year$Lage[1] /fixparam.growth.vessel$CAmax * fixparam.growth.vessel$mw
      }

      fixparam.growth.vessel$CAmax <- gR.year$Lage[1]
      fixparam.divi$deltaD         <- gR.year$Dage[1]

    }else{
      fixparam.growth.vessel$CAmax <- gR.year$Lage[1] *
        fixparam.growth.origin$Values[fixparam.growth.origin$Parameter == "CAmax" & fixparam.growth.origin$Module =="VesselGrowth" ]

      ## Dage use  ####
      fixparam.divi$deltaD <- gR.year$Lage[1]^fixparam.divi$LAtoDiv * # gR.year$Dage[1] *
        fixparam.growth.origin$Values[fixparam.growth.origin$Parameter == "deltaD" & fixparam.growth.origin$Module =="CambialActivity" ]

      if (testLim == T ) { # 先不改，不要用
        fixparam.growth.vessel$ml <- unique(gR.year$Lage) * fixparam.growth.vessel$ml
        fixparam.growth.vessel$mw <- unique(gR.year$Lage) * fixparam.growth.vessel$mw
      } ##

    }
    ## diff - age
    # fixparam.divi$deltaD <- unique(gR.year$Dage) * fixparam.divi$deltaD
  }

  # fixparam.growth.vessel$WAmax <- age$Lage[age$Year == years ] * fixparam.growth.vessel$WAmax

  ## 年初重置部分参数
  dynparam.growth.t <- dynparam.growth.0
  # dynparam.growth.t$Age <- dynparam.growth.t$Age
  growth.day <- 1

  ## 每年的空白生长表
  dailyCells <- list(
    dailyFiber = matrix(data = 0, ncol = ncol(cells), nrow = 400) |>
      as.data.frame() |>
      dplyr::rename( !!!setNames(  paste0("V", seq(1,12,1)), colnames(cells))  ) |>
      dplyr::mutate( Year = years ,cell_L = c(1:400) ),
    dailyVessels = matrix(data = 0, ncol = ncol(vessels), nrow = 400) |>
      as.data.frame() |>
      dplyr::rename( !!!setNames(  paste0("V", seq(1,14,1)), colnames(vessels))  ) |>
      dplyr::mutate( Year = years ,cell_L = c(1:400) )
  )
  ## 生长当年的每日生长记录
  summaryDaily <- purrr::map( 1 : 366 ,
                              function(x  ){ c <- data.frame( seq(1: 400) , rep( x,400  ),matrix(nrow =400,ncol = 27   ) )
                              colnames(c) <- c( colnames(cells), paste0("V", colnames(vessels[,-1:-2])) ,"VAs","CAs","Dh","Kh","Raddist"   )
                              return( c ) }   )

  ## 生成每日活动变量记录表
  dailyParameters <- data.frame(Year = NA,
                                DOY = c(1:366)  ,
                                matrix(NA , nrow = 366, ncol= 20   ))

  colnames(dailyParameters) <- c('years',	'Today',	'Age',	'czgR',	'dCA_cz'	,'deltaVN',	'egR','wgR',	'grwothSeason',
                                 'L_i.fiber',	'L_i.vessel',	'SumCL',	'SumV',	'SumVL',	'T_age',
                                 'v_c.fiber',	'v_c.vessel',	'v_l.fiber',	'v_l.vessel',	'v_w.fiber',	'v_w.vessel',	'Vcz')
  dailyParameters$years <- years
  # ## 生成年内细胞表
  # summaryYears <- matrix( nrow = 400 ,ncol = 29  ) |> data.frame() |>
  #   dplyr::rename( !!!setNames( paste0("X", c(1:29)) ,c( colnames(cells), paste0("V", colnames(vessels[,-1:-2])) ,
  #                                                        "VAs" , "CAs" ,"Dh" ,"Kh","Raddist"   ))  )
  # summaryYears$cell_L <- c(1:400)
  #
  ## 生成年内细胞汇总表
  AnnualGrowth <- matrix(nrow = 366, ncol= 22   ) |> as.data.frame()  ## 16 - 22

  colnames(AnnualGrowth ) <- c( 'Year' ,'DOY' , 'RingArea' ,'RingWidth' ,'CellLayer' ,
                                'MeanVesselLumenArea' , 'MaxVesselLumenArea', 'VesselNumber' ,'CellNumber' , 'VesselTotalLumenArea',
                                'VesselDensity', 'RCTA', 'MeanDh' , 'MeanKh' ,'Ks','dD', ## 1:16 old part
                                "fiberInE", "fiberInT","fiberMature", "rwInE", "rwInT","rwMature" ) ## +6 new part
  AnnualGrowth$DOY  <- c(1:366)

  AnnualGrowth$Year <- years

  ## Inner Year ######
  for ( Today in 1: nrow(gR.year)) { ## daily sycle start
    ## 每日 clim
    clim.today <- gR.year[Today,]

    if (clim.today$aaT >= fixparam.divi$AAT & clim.today$rootd > 0 ) { ## if start 可以生长时计算细胞分裂分化

      # print( paste( "START",Today    ) )

      ## Division period start
      if (division == 'fix') {
        deltaD_T <- fixparam.divi$deltaD
      } else {
        waterDivLim  <- fixparam.divi$a1 * min(clim.today$gM ,clim.today$gV )

        rctaDivLim <- fixparam.divi$a2 *RCTAt/RCTA[ Today ]
        # rctaDivLim[ RCTA[ Today ]  == 99 ] <- waterDivLim ## ERROR CATCH mabey don't use it :20240803 修正后没有99

        Div_limit <-  switch( Dcase,
                              "min" = min(waterDivLim,rctaDivLim),
                              "mean" = mean( c(waterDivLim,rctaDivLim)),
                              "multiply" = waterDivLim*rctaDivLim )

        deltaD_T <- fixparam.divi$deltaD * ( 1 +  fixparam.divi$Div_alpha-
                                               exp( log(1/fixparam.divi$Div_alpha) * - Div_limit )  )
      } ## if div fix end -----
      ### Division  period end

      newCells <- rBTR:::cell_division(clim.today = clim.today, fixparam.divi = fixparam.divi,
                                         fixparam.growth.fiber = fixparam.growth.fiber,
                                         fixparam.growth.vessel = fixparam.growth.vessel,
                                         dynparam.growth.t = dynparam.growth.t,
                                         cells = cells, vessels = vessels,CZgR = CZgR,deltaD_T = deltaD_T)

      dailyCells <- rBTR:::daily_grwoth( newCell = newCells[[1]], newVessel = newCells[[2]], vesselsNum = newCells[[3]],
                                           dailyCells , cells, vessels,
                                           clim.today,
                                           fixparam.growth.fiber, fixparam.growth.vessel ,
                                           dynparam.growth.t =  newCells[[4]])

      dailyParameters[dailyParameters$Today == Today,c(-1,-2) ] <- dynparam.growth.t <- newCells[[4]]


      #### 汇总日数据 #####
      if( any( dailyCells$dailyFiber$DDOY != 0)  ){ ### 有细胞成熟后计算？

        growth.day <- as.numeric(Today)

        summaryDaily[[ Today ]][,1:12] <- dailyCells$dailyFiber
        summaryDaily[[ Today ]][,13:24] <- dailyCells$dailyVessels[,-1:-2]
        summaryDaily[[ Today ]] <- filter( summaryDaily[[ Today ]], CA != 0  )

        summaryDaily[[Today]] <-
          dplyr::mutate(summaryDaily[[Today]], VAs=VCA *VVN,CAs =  fixparam.growth.origin$Values[ fixparam.growth.origin$Parameter == "TreeRingWidth" ] / CTD * CA )

        summaryDaily[[Today]] <-
          dplyr::mutate( summaryDaily[[Today]] ,Dh = (VCRD-2*VWT), Kh = (10^-24 * pi * 998.2)/(128*1.002*10^-9) *(VCRD-2*VWT)^4  ,
                         Raddist =  round( (cumsum( VAs + CAs ) - VAs - CAs)/fixparam.growth.origin$Values[ fixparam.growth.origin$Parameter == "TreeRingWidth" ] ,3 ) )

        # |> mutate(doy = "Today", .after = "Year")

        daily.t <- summaryDaily[[Today]]

        daily.t$VCV[daily.t$VDDOY == 0 ] <- 0 ## 不成熟的导管面积为0 （计算Ks和RCTA时不计入）

        AnnualGrowth[ AnnualGrowth$DOY == Today, 3:15 ] <-  daily.t   |> ## 不筛选成熟细胞
          dplyr::summarise(
            RingArea = (mean( fixparam.growth.origin$Values[ fixparam.growth.origin$Parameter == "TreeRingWidth"] / CTD ,na.rm = T) *
                          sum(CA,na.rm = T) + sum( VCA *VVN ,na.rm = T)) / 10^6  ,
            RingWidth = RingArea / fixparam.growth.origin$Values[ fixparam.growth.origin$Parameter == "TreeRingWidth"] * 1000,
            CellLayer = max(cell_L,na.rm = T),
            MeanVesselLumenArea = mean( VCV+0 ,na.rm = T ),
            MaxVesselLumenArea = max( VCV+0 ,na.rm = T  ),
            VesselNumber = sum(VVN,na.rm = T ),
            CellNumber = mean( fixparam.growth.origin$Values[ fixparam.growth.origin$Parameter == "TreeRingWidth"]  / CTD * CellLayer ,na.rm = T) + max(VNoV,na.rm = T) ,
            VesselTotalLumenArea =  sum( VCV *VVN ,na.rm = T)/10^6,
            VesselDensity = VesselNumber / RingArea,
            RCTA = VesselTotalLumenArea / RingArea,
            MeanDh = sum( (VCRD-2*VWT) ^5 ,na.rm = T ) / sum((VCRD-2*VWT)^4 ,na.rm = T )*10^-6, ##
            # MeanKh = (10^-24 * pi * 998.2)/(128*1.002*10^-9) *sum( (VCRD-2*VWT)^4  ,na.rm = T  ) , ##
            MeanKh = ( pi * 998.21)/(128*1.002*10^-9) * sum( (  (VCRD-2*VWT)*10^-6   )^4  ,na.rm = T  )   * VesselDensity *10^3 , ##
            Ks = MeanKh/ (fixparam.growth.origin$Values[ fixparam.growth.origin$Parameter == "TreeRingWidth"] * RingWidth),
          ) ## end summarise -----

        AnnualGrowth[ AnnualGrowth$DOY == Today, 16 ] <- deltaD_T  ## dD

        AnnualGrowth[ AnnualGrowth$DOY == Today, 22 ] <-
          max( daily.t$Raddist[ daily.t$DDOY != 0 ]) ##RW mature fiber
        AnnualGrowth[   AnnualGrowth$DOY == Today, 21 ] <-
          max( daily.t$Raddist[ daily.t$TDOY != 0 ]) - max( daily.t$Raddist[ daily.t$DDOY != 0 ]) ##RW Thickness fiber
        AnnualGrowth[   AnnualGrowth$DOY == Today, 20 ] <-
          AnnualGrowth[   AnnualGrowth$DOY == Today, 4 ] -
          max( daily.t$Raddist[ daily.t$TDOY != 0 ])/1000 ##RW Enlargement fiber

        AnnualGrowth[   AnnualGrowth$DOY == Today, 19 ] <-
          max( daily.t$cell_L[ daily.t$DDOY != 0 ] ) ##cell_L mature fiber
        AnnualGrowth[   AnnualGrowth$DOY == Today, 18 ] <-
          max( daily.t$cell_L[ daily.t$TDOY != 0 ]) - max( daily.t$cell_L[ daily.t$DDOY != 0 ]) ##cell_L Thickness fiber
        AnnualGrowth[   AnnualGrowth$DOY == Today, 17 ] <-
          AnnualGrowth[   AnnualGrowth$DOY == Today, 5 ] -
          max( daily.t$cell_L[ daily.t$TDOY != 0 ]) ##cell_L Enlargement fiber

        # is.na(  AnnualGrowth[   AnnualGrowth$DOY == Today, ])

        RCTAt <- AnnualGrowth[   AnnualGrowth$DOY == Today, 3:15 ]$RCTA
      } ## if nrow( dailyResult ) end----

    } ## if (actAccumulatedTem >= fixparam.divi$AT) end  可生长日期判定结束 -----------------

  } ## today sycle end  每日生长计算结束 ------------------

  # AnnualGrowth <- AnnualGrowth[ !is.na(AnnualGrowth$RingArea ), ]

  ## 这里需要对年结束进行汇总。
  if (intraannual == T & writeRes == T) {
    data.table::rbindlist(summaryDaily,use.names = T,fill = T,idcol = "doy") |>
      dplyr::filter(CA != 0) |>
      ##filter(  !is.na( EDOY)) |>
      dplyr::select(c("Year","doy","cell_L",everything())) |>
      data.table::fwrite(paste0(redir,"/", as.character(years), ".csv" )   )
  } ## if intraannual & writeRes  end -------

  #   summaryYears <- summaryDaily[[growth.day]]
  dailyParameters$years <- years

  summaryYears <- summaryDaily[[growth.day]] |> dplyr::mutate(VNoV = cumsum(VVN))
  summaryYears[ summaryYears == 0  ] <- NA
  summaryYears$Raddist[ is.na(summaryYears$Raddist )  ] <- 0 ###ERROR CATCH
  # summaryYears$VNoV[summaryYears$VCA == 0 ] <- 0

  res <- list(
    summaryYears = summaryYears, ## 年轮细胞表
    AnnualGrowth = AnnualGrowth[ !is.na(AnnualGrowth$RingArea ), ] , ## 年内生长汇总
    dailyParameters = dailyParameters[ !is.na( dailyParameters$czgR),  ]
  )
  return( res )
}




