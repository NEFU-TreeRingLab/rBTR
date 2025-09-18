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
#' @param division dealtDivision 'fix' is constant, 'limit' is gM/gV & RCTA control.
#' @param mode 气象因子组合模式 "min","multiply".
#' @param Dcase How to calculate Division Limit, including: "min","mean","multiply".
#' @param Named Named save document
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
btr <- function(  clim, parameters, age, syear = NA, eyear = NA, Cores = 1,
                  writeRes = F, intraannual = F, gTmethod = "Jonhson", mode = 'mix',
                  Dcase = "min",division = 'fix', Named = NULL ,setlayers = 300 ) { ## functions start  ring_width,
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
  parameters$Values <- as.numeric(parameters$Values) ## 确保为数值数据

  fixparam.divi <- parameters[ parameters$Module == "CambialActivity" ,
                               c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame() ## 形成层分裂相关参数

  fixparam.growth.fiber <- parameters[ parameters$Module == "FiberGrowth" ,
                                       c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame() ## 纤维细胞生长相关参数

  fixparam.growth.vessel <- parameters[ parameters$Module == "VesselGrowth" ,
                                        c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame() ## 导管细胞生长相关参数

  ## 提取微气候模型参数
  growth_Param <- parameters[ parameters$Module == "GrowthRate" ,
                              c("Parameter","Values") ]   |>
    tibble::remove_rownames()|>
    tibble::column_to_rownames('Parameter') |>
    t() |>
    as.data.frame()

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


  # })

  #  system.time({
  #   clim4 <- data.table::as.data.table(clim)
  #   clim3 <- compute_gR2(clim = clim4[Year %in% c(syear:eyear),] , growth_Param, gTmethod, mode )
  # })

  ## 建立形成层细胞表取代动态参数表
  # 建立表 动态参数、初始值
  params <- c('dCA_cz','v_c.fiber','v_w.fiber','v_l.fiber',
              'v_c.vessel','v_w.vessel','v_l.vessel','dVN_cz','Vcz')

  # 单步完成：创建模板 + 合并 + 添加年龄列
  GrowthRate <- microclim[, c('Year', 'DOY','GS','gM','gV','L_i.fiber','L_i.vessel','egR','wgR')][,
                           c(params) := 0L ][  ## 使用快速合并提取年龄
                           AgeTrend,, on = 'Year'][  ## error catch
                           ,grwothSeason := data.table::fifelse(  (egR <= 0.01 & DOY > 200) | GS == 'UN' ,1,0 )][ # grwothSeason 0:生长季; 1:非生长季
                           , `:=`( egR = data.table::fifelse(grwothSeason == 1 , 0 , egR),# 生长季快结束时细胞不扩大
                                       wgR = data.table::fifelse( egR == 0 & grwothSeason ==1 & GS != 'UN' , max(0.01,wgR), wgR )) ][ # 非生长季wgR增加，快速木质化  & 计算速率
                           , `:=`( v_c.fiber = egR * fixparam.growth.fiber$va_c.fiber,
                                   v_w.fiber = wgR * fixparam.growth.fiber$va_w.fiber * ( 1 + L_i.fiber ),
                                   v_l.fiber = wgR * fixparam.growth.fiber$va_l.fiber * ( 1 + L_i.fiber ),
                                   v_c.vessel = egR * fixparam.growth.vessel$va_c.vessel,
                                   v_w.vessel = wgR * fixparam.growth.vessel$va_w.vessel * ( 1 + L_i.vessel ),
                                   v_l.vessel = wgR * fixparam.growth.vessel$va_l.vessel * ( 1 + L_i.vessel ),
                                   Vcz = egR * fixparam.divi$va_cz )]   #

   ## 检查并计算 形成层扩大速率的年龄趋势
   if ( is.null(GrowthRate$Tage)    ) {
     GrowthRate$Tage <- fixparam.divi$alpha_age *
       exp(fixparam.divi$beta_age * GrowthRate$age)
   }
   ## 形成层扩大速率计算
   GrowthRate[,czgR := fixparam.divi$alpha_cz * exp( fixparam.divi$beta_cz * ( egR ) ) ][ # 计算czgR
     , czgR := data.table::fifelse( czgR > egR, egR, czgR  ) ][ # 形成层活动速率不超过总体速率
     , v_cz := data.table::fifelse( grwothSeason == 1, 0, ( fixparam.divi$va_cz * Tage )* czgR ) ] # 计算形成层扩大率

   ##
   ## 细胞生长上限趋势 和 导管分化趋势校准
   if ( all(!is.na( unique(GrowthRate$Lage))) ) {     ## VCA - age
     GrowthRate[ ,`:=`(CAmax = Lage * fixparam.growth.vessel$CAmax,
                       deltaD = Lage^fixparam.divi$LAtoDiv * fixparam.divi$deltaD ) ]
   }

  ## 年循环计算开始前
  ## 计算 权衡的年内平衡曲线
  TradeLine <- fixparam.divi$maxRCTA  *
    nor( microclim$L_i.vessel[microclim$Year == syear ] *-1 )
  ## 20240803 使用L+2 替换 L+1 ：this is L116
  # RCTA[ RCTA <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <- 99
  TradeLine[ TradeLine <=  fixparam.divi$maxRCTA * fixparam.divi$RCTADivT  ] <-
    fixparam.divi$maxRCTA * fixparam.divi$RCTADivT # error catch
  ## __end ----

  ## 设置各类初始值： 纤维细胞和导管初始值

  Cambial <- matrix(NA,ncol = 10, nrow = 1,
                    dimnames = list(c("1"),c('CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY' )) ) |>
    as.data.frame()

  Cambial$CV  <- parameters$Values[ parameters$Parameter == 'CV' ]
  Cambial$WT  <- parameters$Values[ parameters$Parameter == 'WT' ]
  Cambial$CTD <- parameters$Values[ parameters$Parameter == 'CTD' ]
  Cambial$CRD <- Cambial$CV/(Cambial$CTD - 2*Cambial$WT   ) + 2*Cambial$WT
  Cambial$CA  <- Cambial$CTD*Cambial$CRD
  Cambial$WA  <- Cambial$CA-Cambial$CV

  cells <- Cambial |> as.data.table()
  cells[is.na(cells)] <- 0

  # vessels <- Cambial |> dplyr::rename_with(~ paste0("V", .), dplyr::everything()) |>
  #   dplyr::mutate( NoV = NA , VN = NA   )
  # vessels[is.na(vessels)] <- 0
  vessels <- matrix(0,ncol = 12, nrow = 1,
                    dimnames = list(c("1"),c('VCA','VCV','VWA','VLWA','VWT','VCRD','VCTD','VEDOY','VTDOY','VDDOY','NoV','VN' )) ) |>
    as.data.frame()
  # 提取生长季开始最早一年
  s_doy <- min(GrowthRate$DOY[GrowthRate$GS == 'GR']  )
  # 提取生长季结束时间
  GSend <- GrowthRate[GS == 'GR', max(DOY), by = Year]+1

  ## 生成空list N年，每年预留400层细胞 先生成年份*400行的表格，依次填充年份与层序，再按年拆分为list
  Daily_expanded <- expand.grid(Year = eyear:syear, DOY = s_doy:max(GSend$V1),cell_L = c(1:setlayers))
  AllXylem <- cbind( Daily_expanded,cbind(cells,vessels)[rep(1, times = nrow(Daily_expanded) ), ] ) |>
    data.table::as.data.table()
  AllXylem <- split(AllXylem,AllXylem$DOY)
  # AllCells <- cbind( Daily_expanded, cells[rep(1, times = nrow(Daily_expanded) ), ] )   |>
  #   as.data.table()
  # AllVessels <- cbind( Daily_expanded, vessels[rep(1, times = nrow(Daily_expanded) ), ] )   |>
  #   as.data.table()
  AllCambial <- cbind( expand.grid(Year = eyear:syear, DOY = s_doy:max(GSend$V1)) , CA = cells$CA   )  |>
    as.data.table()
  AllCambial <- AllCambial[,`:=`(OriCA = CA,dtDef = 0 ,
                                 checkDiv = 0, checkDef=0,sumDiv = 0, sumDef = 0,RCTAt = 0 )  ]     # [,c('Year','CA','OriCA','dtDiv')]
  ## 生成用于计算的表格
  Cambials <- AllCambial[ DOY == s_doy,][,DOY := 0]
  Fibers <- cbind( expand.grid(Year = eyear:syear, cell_L = c(1:setlayers)), cells[rep(1, times = length(eyear:syear)* setlayers), ] ) |>
    data.table::as.data.table()
  Vessels <- cbind( expand.grid(Year = eyear:syear, cell_L = c(1:setlayers)), vessels[rep(1, times = length(eyear:syear)* setlayers ), ] ) |>
    data.table::as.data.table()
  ### 开始循环，由于各年间生长独立，因此对各年使用向量化计算
  pb <- utils::txtProgressBar(title = "PB",min = s_doy, max = max(GSend$V1) , style = 3, file = stderr())
  for (i_doy in s_doy:max(GSend$V1)) { ## i_doy= s_doy i_doy= 170 ##

    Doy_GR <- GrowthRate[ DOY == i_doy,   ]

    ## Division period start
    if (division != 'fix') {
      waterDivLim  <- fixparam.divi$a1 * pmin(Doy_GR$gM ,Doy_GR$gV )

      rctaDivLim <-  data.table::fifelse( Cambials$RCTAt ==0,0, fixparam.divi$a2 * TradeLine[ i_doy ]/Cambials$RCTAt )
      # rctaDivLim[ RCTA[ Today ]  == 99 ] <- waterDivLim ## ERROR CATCH mabey don't use it :20240803 修正后没有99

      Div_limit <-  switch( Dcase,
                            "min" = pmin(waterDivLim,rctaDivLim),
                            "mean" = pmean( c(waterDivLim,rctaDivLim)),
                            "multiply" = waterDivLim*rctaDivLim )

      Doy_GR[, deltaD:= data.table::fifelse( Div_limit == 0,
                                             deltaD , deltaD * ( 1 +  fixparam.divi$Div_alpha -
                                                            exp( log(1/fixparam.divi$Div_alpha) * - Div_limit ) ) )  ]

    } ## if div fix end -----

    ##
    # Res_cambial <- Cambials
    for (i in 1:10) {
      #### 分裂进入木质部
      Cambials <- cell_division(Doy_GR , Cambials ) ## 输入气候因子，参数和形成层大小，输出分裂出多少个细胞和多少个导管

      ### 查看形成的细胞数量,输入木质部表中
      if (any(Cambials$checkDiv != 0 )) {
        Fibers[ Cambials[checkDiv ==1,], EDOY := i.DOY, on = .( Year = Year, cell_L = sumDiv  )  ]
        if (any(Cambials$checkDef != 0 ) ) {
          Vessels[ Cambials[checkDef ==1& checkDef != 0 ,],
                   `:=`(VEDOY = i.DOY,NoV = checkDef, VN = sumDef,
                        VCA = cells$CA, VCV  = cells$CV   , VWA = cells$WA, VLWA = cells$LWA,
                        VWT = cells$WT, VCRD = cells$CRD, VCTD = cells$CTD ) ,
                   on = .( Year = Year, cell_L = sumDiv  )  ]
        }
        Cambials[,`:=`( checkDiv = 0, checkDef = 0 )] ## 重置计数
      }

      # #### 细胞生长部分
      if (any( Fibers$EDOY != 0 & Fibers$DDOY ==0 )) {
        Fibers[  cells_growth( Doy_GR,cell = Fibers[EDOY != 0 & DDOY ==0,], param = fixparam.growth.fiber,CorV = "C",i_doy ),
                 c('CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY') :=
                   mget(paste0("i.", c('CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY'))),
                 on = .(Year,cell_L)   ]


      } ## if fiber growth

      if (any( Vessels$VEDOY != 0 & Vessels$VDDOY ==0 )) {
        Vessels[  cells_growth( Doy_GR,cell = Vessels[VEDOY != 0 & VDDOY ==0,], param = fixparam.growth.vessel, CorV = "V",i_doy ),
                  c('VCA','VCV','VWA','VLWA','VWT','VCRD','VCTD','VEDOY','VTDOY','VDDOY','NoV','VN') :=
                   mget(paste0("i.", c('VCA','VCV','VWA','VLWA','VWT','VCRD','VCTD','VEDOY','VTDOY','VDDOY','NoV','VN'))),
                 on = .(Year,cell_L)   ]
      }




    } ## For i end -------

    ### 汇总每日形成层和木质部特征
    AllCambial[Cambials,
                  c('CA','OriCA','dtDef','sumDiv','sumDef','RCTAt') :=
                 mget(paste0("i.", c('CA','OriCA','dtDef','sumDiv','sumDef','RCTAt'))),
               on = .(Year,DOY)    ]
    AllXylem[[as.character(i_doy)]][  Fibers,
               c('CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY') :=
               mget(paste0("i.", c('CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY'))),
               on = .(Year,cell_L)   ]
    AllXylem[[as.character(i_doy)]][  Vessels,
               c('VCA','VCV','VWA','VLWA','VWT','VCRD','VCTD','VEDOY','VTDOY','VDDOY','NoV','VN') :=
                 mget(paste0("i.", c('VCA','VCV','VWA','VLWA','VWT','VCRD','VCTD','VEDOY','VTDOY','VDDOY','NoV','VN'))),
               on = .(Year,cell_L)   ]

    ## 计算RCTA ### 和 Ks
    if ( any(Vessels$VDDOY > 0 )   ) { ## cRCTAt
      R_daily <- Fibers[,.(Af= sum(CA)* parameters$Values[parameters$Parameter == 'Twidth']/ mean(CTD)),by='Year'
                        ][ Vessels[ VDDOY > 0 , .(VA = sum(VCV * NoV ), Av = sum(VCA * NoV )), by='Year'], on='Year'
                        ][ ,RCTAt := VA/(Av+Af),by = 'Year' ]

      Cambials[R_daily, RCTAt:= i.RCTAt,on = 'Year' ]
    } ## end cRCTAt
    utils::setTxtProgressBar(pb = pb, value = i_doy)
  } ## end for doy
  close(con = pb)
  #### 分类汇总全数据并输出
  AllXylem_summary <- data.table::rbindlist(AllXylem)[EDOY>0,][ ## 清洗表格
    ## 汇总计算
    , `:=`(
      RAcl = CA * parameters$Values[parameters$Parameter == 'Twidth']/CTD+VCA * NoV,
      RWcl = (CA * parameters$Values[parameters$Parameter == 'Twidth']/CTD+VCA * NoV)/
        parameters$Values[parameters$Parameter == 'Twidth'],
      FV = ceiling(parameters$Values[parameters$Parameter == 'Twidth']/CTD),
      Dh = VCRD-2*VWT,
      VAcl = VCA * NoV,
      khf = (10^-24 * pi * 998.2)/(128*1.002*10^-9) *(CRD-2*WT)^4,
      Kh = (10^-24 * pi * 998.2)/(128*1.002*10^-9) *(VCRD-2*VWT)^4
    )
  ]

  ## 年终生长量
  AnnualGrowth<- AllXylem_summary[DOY == i_doy,][,DOY := NULL][ ## 清洗表格
    ## 汇总计算
    , .(
      RingArea = sum(RAcl,na.rm = T)/10^6,
      RingWidth = sum(RAcl,na.rm = T)/parameters$Values[parameters$Parameter == 'Twidth'],
      CellLayer = max(cell_L,na.rm = T),
      MeanVesselLumenArea = mean(VCV ,na.rm = T),
      MaxVesselLumenArea = max(VCV ,na.rm = T ),
      VesselNumber = max(VN,na.rm = T ),
      CellNumber = sum(FV) + max(VN,na.rm = T ),
      VesselTotalLumenArea =  sum( VAcl ,na.rm = T)/10^6,
      MeanDh = sum( Dh ^5 ,na.rm = T ) / sum(Dh^4 ,na.rm = T )*10^-6,
      MeanKh = ( pi * 998.21)/(128*1.002*10^-9) * sum( ( Dh*10^-6 )^4 ,na.rm = T ) * max(VN,na.rm = T )/sum(RAcl,na.rm = T) *10^3 , ##
      Ks = ( pi * 998.21)/(128*1.002*10^-9) * sum( ( Dh*10^-6 )^4 ,na.rm = T ) * max(VN,na.rm = T )/sum(RAcl,na.rm = T) *10^3/sum(RAcl,na.rm = T)
    ),by="Year"
  ][  ## 再次汇总计算
    , `:=`(
      VesselDensity = VesselNumber / RingArea,
      RCTA = VesselTotalLumenArea / RingArea
    ),by = "Year"
  ]

  xylem_trait <- AllXylem_summary[DOY == i_doy,][,DOY := NULL][AnnualGrowth[,c('Year','RingWidth','VesselNumber')],on= 'Year'][
    ## 汇总计算
    , `:=`(
      Raddist = cumsum( RWcl ),
      RRaddist = cumsum( RWcl )/RingWidth*100,
      RLayer = cell_L / max(cell_L) *100,
      RVN = VN / max(VN)*100
    ),by="Year"
  ]

  InterAnnualGrowth <- AllXylem_summary[
    ,`:=`(Type = data.table::fcase(
      TDOY == 0 & DDOY == 0 , "E",
      TDOY != 0 & DDOY == 0 , "T",
      DDOY != 0, "M"
    )  )
  ][ ## 清洗表格
    ## 汇总计算
    , .(
      RingArea = sum(RAcl,na.rm = T)/10^6,
      RingWidth = sum(RAcl,na.rm = T)/parameters$Values[parameters$Parameter == 'Twidth'],
      CellLayer = max(cell_L,na.rm = T),
      MeanVesselLumenArea = mean(VCV ,na.rm = T),
      MaxVesselLumenArea = max(VCV ,na.rm = T ),
      VesselNumber = max(VN,na.rm = T ),
      CellNumber = sum(FV) + max(VN,na.rm = T ),
      VesselTotalLumenArea =  sum( VAcl ,na.rm = T)/10^6,
      MeanDh = sum( Dh ^5 ,na.rm = T ) / sum(Dh^4 ,na.rm = T )*10^-6,
      MeanKh = ( pi * 998.21)/(128*1.002*10^-9) * sum( ( Dh*10^-6 )^4 ,na.rm = T ) * max(VN,na.rm = T )/sum(RAcl,na.rm = T) *10^3 , ##
      Ks = ( pi * 998.21)/(128*1.002*10^-9) * sum( ( Dh*10^-6 )^4 ,na.rm = T ) * max(VN,na.rm = T )/sum(RAcl,na.rm = T) *10^3/sum(RAcl,na.rm = T),
      RWE = sum(RAcl[Type == 'E'],na.rm = T)/parameters$Values[parameters$Parameter == 'Twidth'],
      RWT = sum(RAcl[Type == 'T' ],na.rm = T)/parameters$Values[parameters$Parameter == 'Twidth'],
      RWM = sum(RAcl[Type == 'M'],na.rm = T)/parameters$Values[parameters$Parameter == 'Twidth']
      ),by= c("Year",'DOY')
  ][  ## 再次汇总计算
    , `:=`(
      RWTM = RWM + RWT,
      VesselDensity = VesselNumber / RingArea,
      RCTA = VesselTotalLumenArea / RingArea
    ),by = c("Year",'DOY')
  ]

  Res <- list(AllXylem_summary=AllXylem_summary,
              AnnualGrowth =AnnualGrowth ,
              xylem_trait =xylem_trait ,
              InterAnnualGrowth =InterAnnualGrowth ,
              GrowthRate = GrowthRate )



return(Res)


} ## btr end ------------





