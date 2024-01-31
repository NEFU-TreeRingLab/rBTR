#########################
##### Main function #####
#########################

#' Main function calculate tree-ring anatomic
#'
#' @param clim  climate data must have Year，DOY，Temp ，Prec
#' @param parameters model parameter date ，we ues a excel data to save it.
#' @param age age trend of simulate ring
#' @param syear modeling start year
#' @param eyear modeling end year
#' @param writeRes auto write res_data, logical value, if 'Ture', output result list as xlsx in file.
#' @param intraannual output daily result, logical, if 'Ture', output daily result in file.
#' @param gTmethod method of calculate gT ,'VS' and 'Jonhson et al.'
#'
#' delete this part
#' @param division dealtDivision 'fix' is constant, 'limit' is gM/gV & RCTA control. ## cancel this part
#' @param testLim use Lage limit mw & ml , Lage in age param_data, regulation cell enlargement by CAmax param .#'
#' @param CZgR growth limit : climate-Cambial Cell, cell wall-gM, cell wall-gV, NewgRtype; 0 for False, 1 for True#'
#' @param Dcase How to calculate Division Limit, including: "min","mean","multiply".

#' @param Pbar Whether to show progress bar
#' @param testMod test Mode, logical, if 'False', Lage&Dage use ratio

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
#' @importFrom purrr map map2
#'
#'
#' @export
#'


btr <- function(  clim, parameters, age, syear = NA, eyear = NA ,
                  writeRes = T,intraannual = F, gTmethod = "Jonhson" ,
                  division = "limit" ,
                  testLim = F,
                  CZgR = c(1,0,0,1 ) ,
                  Pbar = F ,
                  testMod = F ,
                  Dcase = "min",
                  Named = NULL ) { ## functions start  ring_width,

  ##
  writeRes[intraannual == T] <- T

  if (writeRes == T) {
    redir <- paste0(  Named,  "res_",format(Sys.time(), "%Y%m%d_%H-%M-%OS"))

    dir.create(path = eval(redir) )
  }

   ## 提取微气候模型参数

  parameters$values <- as.numeric(parameters$values)
  growth_Param  <- parameters[parameters$modul == "gR",] ## 生长速率阈值参数

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

  if ( all( gTmethod == "VS" ) ) {
    microclim <- Compute_gR(clim , growth_Param) |>
      dplyr::left_join(age,by="Year")  |>
      data.table::as.data.table() ##VS
  } else { microclim <- Compute_gR2(clim , growth_Param) |>
    dplyr::left_join(age,by="Year")  |>
    data.table::as.data.table() } ## 'Jonhson'


  ##### 细胞大小计算 ——————

  ## 提取分裂模型参数：

  fixparam.divi <- parameters[ parameters$modul == "division" & parameters$paramtype == "fixed",
                               c("parameter","values") ]   |>
    tidyr::spread( key = 'parameter', value = 'values')

  fixparam.growth.fiber <- parameters[parameters$modul == "growthC" & parameters$paramtype == "fixed" ,
                                      c("parameter","values")]  |>
    tidyr::spread( key = 'parameter', value = 'values')

  fixparam.growth.vessel <- parameters[parameters$modul == "growthV" & parameters$paramtype == "fixed",
                                       c("parameter","values") ] |>
    tidyr::spread( key = 'parameter', value = 'values')

  dynparam.growth.0 <- parameters[parameters$modul == "division" & parameters$paramtype == "dynamic" ,
                                  c("parameter","values") ] |>
    tidyr::spread( key = 'parameter', value = 'values')

  fixparam.growth.origin <- parameters[parameters$paramtype == "fixed",
                                       c("parameter","modul", "values")]

  #### VA & Differentiation
  MaxVesselArea <- fixparam.growth.vessel$CAmax
  MaxdeltaD <- fixparam.divi$deltaD
  #### VA & Differentiation ---

  ## 年循环计算
  ## RCTA line

  RCTA <- fixparam.divi$maxRCTA  *
    nor( microclim$L_i.vessel[microclim$Year == syear ] *-1 )
  RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <- 99
  ## __end ----

  cells <- dplyr::filter(parameters ,  grepl("C0", modul ) ) |>
    dplyr::select( c("parameter","values") ) |>
    tidyr::spread( key = 'parameter', value = 'values') |>
    dplyr::select( "cell_L","Year", everything())
  cells[is.na(cells)] <- 0

  vessels <- dplyr::filter(parameters ,  grepl("V0", modul ) ) |>
    dplyr::select( c("parameter","values") ) |>
    tidyr::spread( key = 'parameter', value = 'values')|>
    dplyr::select( "cell_L","Year", everything())
  vessels[is.na(vessels)] <- 0

  temThreshold <- dplyr::filter(growth_Param, parameter == "T1")$values

  ##test use
  # dailyParameters = logical()
  dailyParameters <- data.frame(rep( syear:eyear,each =366  ),
                                rep( 1:366,time = (eyear - syear +1)  ),
                                matrix(nrow = (eyear - syear +1) *366, ncol= 19   ))

  colnames(dailyParameters) <- c('years',	'Today',	'Age',	'czgR',	'dCA_cz'	,'deltaVN',	'egR',	'grwothSeason',
                                 'L_i.fiber',	'L_i.vessel',	'SumCL',	'SumV',	'SumVL',	'T_age',
                                 'v_c.fiber',	'v_c.vessel',	'v_l.fiber',	'v_l.vessel',	'v_w.fiber',	'v_w.vessel',	'Vcz')

  # prow <- ceiling( 10/log( 2, (age$Tage * fixparam.divi$va_cz + 1) ))*
  #                microclim[ gT > 0,.(x = max(DOY) - min(DOY)) ,by = Year ][,x]  ###
  prow <- rep( 366, (eyear - syear+1) )


  summaryYears <- purrr::map2( seq(syear,eyear,1) , prow ,
                               function( x ,prow){ c <- data.frame( seq(1: prow) , rep( x,prow  ),matrix(nrow =prow,ncol = 27   ) )
                               colnames(c) <- c( colnames(cells), paste0("V", colnames(vessels[,-1:-2])) ,
                                                 "VAs" , "CAs" ,"Dh" ,"Kh","Raddist"   )
                               return( c ) }  )

  names( summaryYears ) <- as.character( syear:eyear )

  ## summary data.table
  AnnualGrowth <- matrix(nrow = (eyear - syear +1) *366, ncol= 22   ) |> as.data.frame()  ## 16 - 22

  colnames(AnnualGrowth ) <- c( 'Year' ,'DOY' , 'RingArea' ,'RingWidth' ,'CellLayer' ,
                                'MeanVesselLumenArea' , 'MaxesselLumenArea', 'VesselNumber' ,'CellNumber' , 'VesselTotalLumenArea',
                                'VesselDensity', 'RCTA', 'MeanDh' , 'MeanKh' ,'Ks','dD', ## 1:16 old part
                                "fiberInE", "fiberInT","fiberMature", "rwInE", "rwInT","rwMature" ) ## +6 new part


  AnnualGrowth$Year <- rep( syear:eyear,each =366  )
  AnnualGrowth$DOY  <- rep( 1:366,time = (eyear - syear +1)  )


  for (years in syear:eyear) { ## year cycle start  年循环计算


    RCTAt <- 0 ### RCTA for dD
    gR.year <- microclim[microclim$Year == years,]

    Xy.fiber  <- data.frame( seq(1: prow[years-syear+1] ) , rep( years,prow[years-syear+1]  ) , matrix( ncol = 10 , nrow = prow[years-syear+1] )   )
    Xy.vessel <- data.frame( seq(1: prow[years-syear+1]) , rep( years,prow[years-syear+1]  ) , matrix( ncol = 12 , nrow = prow[years-syear+1] )   )



    colnames(Xy.fiber) <- colnames(cells)
    colnames(Xy.vessel) <- colnames(vessels)

    ## check max VCA
    if ( is.na( unique(gR.year$Lage)) == F ) {
      ## VCA - age

      if( testMod == T ){

        fixparam.growth.vessel$CAmax <- gR.year$Lage[1]
        fixparam.divi$deltaD         <- gR.year$Dage[1]

      }else{

        fixparam.growth.vessel$CAmax <- MaxVesselArea * gR.year$Lage[1]
        fixparam.divi$deltaD <- MaxdeltaD * gR.year$Dage[1]

        ## old type
        # fixparam.growth.vessel$CAmax <- gR.year$Lage[1] *
          # fixparam.growth.origin$values[fixparam.growth.origin$parameter == "CAmax" & fixparam.growth.origin$modul == "growthV" ]
        # fixparam.divi$deltaD <- gR.year$Dage[1] *
          # fixparam.growth.origin$values[fixparam.growth.origin$parameter == "deltaD" & fixparam.growth.origin$modul =="division" ]
      } ## test mod end ---

    }

    ### 进度条
    if (Pbar == T ) {
      print(  paste0("Modeling ", years," / ", "Final ", eyear )    )
      pb <- utils::txtProgressBar(title = "PB",min = 1, max = nrow(gR.year), style = 3, file = stderr())
    }


    ## 年初重置部分参数
    dynparam.growth.t <- dynparam.growth.0
    dynparam.growth.t$Age <- dynparam.growth.t$Age + years - syear
    growth.day <- 1

    ## 重置结束
    actAccumulatedTem <- 0 ## 活动积温 active accumulated temperature
    # TA <- 0
    summaryDaily <- purrr::map( 1 :366 ,
                                function(x , prow ){ c <- data.frame( seq(1: prow) , rep( x,prow  ),matrix(nrow =prow,ncol = 27   ) )
                                colnames(c) <- c( colnames(cells), paste0("V", colnames(vessels[,-1:-2])) ,"VAs","CAs","Dh","Kh","Raddist"   )
                                return( c ) } ,prow[years - syear+1]  )

    for ( Today in 1:nrow(gR.year) ) { ## daily sycle start

      clim.today <- gR.year[Today,]

      aatToday <- (clim.today$TEM - temThreshold)
      aatToday[aatToday<0] <- 0
      actAccumulatedTem <- actAccumulatedTem + aatToday

      if (actAccumulatedTem >= fixparam.divi$AAT & clim.today$rootd > 0 ) { ## if start 可以生长时计算细胞分裂分化

        ## Division period start
        climateDivLim  <- fixparam.divi$a1 * min(clim.today$gM ,clim.today$gV )

        rctaDivLim <- fixparam.divi$a2 *RCTAt / RCTA[ Today ]
        rctaDivLim[ RCTA[ Today ]  == 99 ] <- climateDivLim ## ERROR CATCH

        Div_limit <- min(climateDivLim,rctaDivLim)

        deltaD_T <- fixparam.divi$deltaD * ( 1 +  fixparam.divi$Div_alpha-
                                               exp( log(1/fixparam.divi$Div_alpha) * - Div_limit )  )
        ### Division  period end

        newCells <-  cell_division(clim.today = clim.today, fixparam.divi = fixparam.divi,
                                   fixparam.growth.fiber = fixparam.growth.fiber,
                                   fixparam.growth.vessel = fixparam.growth.vessel,
                                   dynparam.growth.t = dynparam.growth.t,
                                   cells = cells, vessels = vessels,CZgR = CZgR,deltaD_T = deltaD_T)
        ### , TVA,TA

        newCells[[1]]  <- stats::na.omit(newCells[[1]])
        newCells[[2]]  <- stats::na.omit(newCells[[2]])

        dynparam.growth.t <- newCells[[3]]

        dailyParameters[ dailyParameters$years == years & dailyParameters$Today == Today , ] <-  c(years,Today,dynparam.growth.t)
        ##

        Xy.fiber[ Xy.fiber$cell_L %in% newCells[[1]]$cell_L,  ] <- newCells[[1]]

        Xy.vessel[ Xy.vessel$cell_L %in% newCells[[2]]$cell_L,  ] <- newCells[[2]]

        ## 细胞生长部分 #####


        layer.max.f <- nrow( Xy.fiber[ !is.na(Xy.fiber$EDOY),  ]  )

        growthing.fiber <- Xy.fiber[ !is.na(Xy.fiber$EDOY)& Xy.fiber$DDOY == 0  ,  ]
        if ( layer.max.f != 0 & nrow( growthing.fiber ) != 0 ) {



          growthing.fiber <- cells_growth(cell = growthing.fiber, CorV = "C" , clim.today, layer.max.f,
                                          fixparam.growth.fiber, fixparam.growth.vessel, dynparam.growth.t )

          Xy.fiber[ Xy.fiber$cell_L %in% growthing.fiber$cell_L, ] <- growthing.fiber

        } ## cell end-------------

        layer.max.v <- nrow(Xy.vessel[ !is.na(Xy.vessel$EDOY),  ])
        growthing.vessel <- Xy.vessel[ !is.na(Xy.vessel$EDOY)& Xy.vessel$DDOY == 0  ,  ]

        if( layer.max.v != 0 & nrow(growthing.vessel) != 0  ){  ## 有导管时导管生长



          growthing.vessel <- cells_growth(cell = growthing.vessel, CorV = "V" , clim.today, layer.max.v,
                                           fixparam.growth.fiber,fixparam.growth.vessel, dynparam.growth.t)

          Xy.vessel[Xy.vessel$cell_L %in% growthing.vessel$cell_L ,] <- growthing.vessel

        } ## 有导管时导管生长 end ---------------------

        if(  layer.max.v+layer.max.f != 0 ){
          growth.day <- as.numeric(Today)
          summaryDaily[[ Today ]][,1:12] <- Xy.fiber
          summaryDaily[[Today]][,13:24] <- Xy.vessel[,-1:-2]
          summaryDaily[[Today]] <-
            summaryDaily[[Today]][ !is.na( summaryDaily[[Today]]$EDOY ),]

          summaryDaily[[Today]] <-
            dplyr::mutate(summaryDaily[[Today]], VAs=VCA *VVN,CAs = parameters$values[ parameters$parameter == "width"]  / CTD * CA   )

          summaryDaily[[Today]]$VAs[ is.na(summaryDaily[[Today]]$VAs) ] <-  0

          summaryDaily[[Today]] <-
            dplyr::mutate( summaryDaily[[Today]] ,Dh = (VCRD-2*VWT) ,  Kh = (10^-24 * pi * 998.2)/(128*1.002*10^-9) *(VCRD-2*VWT)^4  ,
                           Raddist =  round(  (cumsum( VAs +  CAs ) - VAs -  CAs)/parameters$values[ parameters$parameter == "width"]   ,3  )    )

          daily.t <- summaryDaily[[Today]]

          daily.t$VCV[daily.t$VDDOY == 0 ] <- 0
          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 3:15 ] <-  daily.t   |> ## 筛选成熟细胞
            dplyr::summarise(
              RingArea = (mean( parameters$values[ parameters$parameter == "width"] / CTD ,na.rm = T) *
                            sum(CA,na.rm = T) + sum( VCA *VVN ,na.rm = T)) / 10^6  ,
              RingWidth = RingArea / parameters$values[ parameters$parameter == "width"] * 1000,
              CellLayer = max(cell_L,na.rm = T),
              MeanVesselLumenArea = mean( VCV+0 ,na.rm = T ),
              MaxVesselLumenArea = max( VCV+0 ,na.rm = T  ),
              VesselNumber = sum(VVN,na.rm = T ),
              CellNumber = mean( parameters$values[ parameters$parameter == "width"]  / CTD * CellLayer ,na.rm = T) + max(VNoV,na.rm = T) ,
              VesselTotalLumenArea =  sum( VCV *VVN ,na.rm = T)/10^6,
              VesselDensity = VesselNumber / RingArea,
              RCTA = VesselTotalLumenArea / RingArea,
              MeanDh = sum( (VCRD-2*VWT) ^5 ,na.rm = T ) / sum((VCRD-2*VWT)^4 ,na.rm = T )*10^-6, ##
              # MeanKh = (10^-24 * pi * 998.2)/(128*1.002*10^-9) *sum( (VCRD-2*VWT)^4  ,na.rm = T  ) , ##
              MeanKh = ( pi * 998.21)/(128*1.002*10^-9) * sum( (  (VCRD-2*VWT)*10^-6   )^4  ,na.rm = T  )   * VesselDensity *10^3 , ##
              Ks = MeanKh/ (parameters$values[ parameters$parameter == "width"] * RingWidth),
            ) ## end summarise -----

          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 16 ] <- deltaD_T

          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 22 ] <-
            max( daily.t$Raddist[ daily.t$DDOY != 0 ]) ##RW mature fiber
          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 21 ] <-
            max( daily.t$Raddist[ daily.t$TDOY != 0 ]) - max( daily.t$Raddist[ daily.t$DDOY != 0 ]) ##RW Thickness fiber
          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 20 ] <-
            AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 4 ] -
            max( daily.t$Raddist[ daily.t$TDOY != 0 ])/1000 ##RW Enlargement fiber

          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 19 ] <-
            max( daily.t$cell_L[ daily.t$DDOY != 0 ]) ##cell_L mature fiber
          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 18 ] <-
            max( daily.t$cell_L[ daily.t$TDOY != 0 ]) - max( daily.t$cell_L[ daily.t$DDOY != 0 ]) ##cell_L Thickness fiber
          AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 17 ] <-
            AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 5 ] -
            max( daily.t$cell_L[ daily.t$TDOY != 0 ]) ##cell_L Enlargement fiber

          RCTAt <- AnnualGrowth[ AnnualGrowth$Year == years & AnnualGrowth$DOY == Today, 3:15 ]$RCTA
        } ## if nrow( dailyResult ) end----

      } ## if (actAccumulatedTem >= fixparam.divi$AT) end  可生长日期判定结束 -----------------

      # } ## Today sycle  end ------

      ## 进度条结束

      if (Pbar == T ) {
        utils::setTxtProgressBar(pb = pb, value = Today)
      }


    } ## today sycle end  每日生长计算结束 ------------------

    ## 这里需要对年结束进行汇总。

    if (intraannual == T & writeRes == T) {
      data.table::rbindlist(summaryDaily,use.names = T,fill = T,idcol = "doy") |>filter(  !is.na( EDOY)) |>
        dplyr::select(c("Year","doy","cell_L",everything())) |>
        data.table::fwrite(paste0(redir,"/", as.character(years), ".csv" )   )
    } ## if intraannual & writeRes  end -------

    summaryYears[[as.character(years)]] <- summaryDaily[[growth.day]]

    if (Pbar == T ) {
      close(con = pb)
    } ## if Pbar end --------



  }  ## year cycle end 年生长计算结束------------------

  ## 合并所有结果;

  AnnualGrowth <- na.omit( AnnualGrowth ) #%>% as.data.table()


  AnnaulRing <- AnnualGrowth[,c(1:16)] |>
    dplyr::group_by(Year) |>
    dplyr::mutate( StartDoy = min(DOY),.before = DOY )|>
    dplyr::filter(DOY == max(DOY)) |>
    dplyr::rename(EndDoy = DOY)

  summaryYears <- data.table::rbindlist(summaryYears) |>filter( !is.na(EDOY) )

  Outputs <- list(AnnaulRing, summaryYears , AnnualGrowth, microclim, as.data.frame( stats::na.omit(dailyParameters) ) )

  names(Outputs) <- c('annaulRing','xylem_trait','AnnualGrowth', 'microclim', 'dailyParameters' )




  resp <- list(parameters,age)
  names(resp) <- c( 'parameters','age')

  if (writeRes == T) {
    openxlsx::write.xlsx(Outputs,paste0( redir, "/Outputs.xlsx" ) )
    openxlsx::write.xlsx(resp,file = paste0(redir,"/Parameters.xlsx" ))
    print( redir)

  }





  return(Outputs)

} ##  TRXXX2 funtion end -----------------------------------------
