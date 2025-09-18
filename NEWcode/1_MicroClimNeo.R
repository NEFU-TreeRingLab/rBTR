##### MicroClimate #####
#' MicroClimate
#' @param climdata Climsdata
#' @param parameters Parameters
#' @param syear Start year
#' @param eyear End year
#' @param lat Latitude of sample sites
#' @param rootd Root zone depth
#' @return gE , SoilM soil moisture, Ls, potEv et.al. climate factor
#' SoilM soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
#' potEv Potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm).
#'
#' @importFrom dplyr filter select mutate
#' @importFrom tidyr spread
#'
#' @export
#'
#'
Compute_clim <- function(climdata , syear = NA, eyear = NA , lat = NA, rootd = 1000,... ){

  ## Column Name
  colnames(climdata) <- toupper(colnames(climdata) )
  if( any(colnames(climdata) %in% 'SOILM' )  ){
    climdata <- dplyr::rename(climdata, Year = YEAR, Month = MONTH, Day =DAY ,soilM = SOILM  )
  } else {
    climdata <- dplyr::rename(climdata, Year = YEAR, Month = MONTH, Day =DAY  )
  }


  #### error-catching ####

  lat[ lat>90 | lat < -90   ] <- NA

  if (is.na(syear) ) {syear <- min(climdata$Year) } else {
    if( syear > max(climdata$Year) ){
      stop("syear is wrong")}
  }

  if (is.na(eyear)) {eyear <- max(climdata$Year)} else {
    if( eyear < min(climdata$Year) ){
      stop("syear is wrong")
    }
  }

  ## 提取对应年份数据
  climdata <- dplyr::filter( climdata,Year %in% c( syear:eyear) )


  #### 检查气象数据，看是否需要通过计算补齐 ####
  dtColname <- colnames( climdata )

  ## check date col
  ifelse( all(c('Year', 'DOY') %in% dtColname ),C1 <- 2 ,
          ifelse(all( c('Year', 'Month', 'Day') %in% dtColname ), C1 <- 1, C1 <- 0 )     )
  ## check Climfactor
  ctem <- cVPD <- csoilM <- cLs <- crootd <- 0

  cVPD[ 'RH' %in% dtColname ] <- 1
  csoilM[ 'PRE' %in% dtColname ] <- 1
  cLs[ !is.na(lat)  ] <- 1
  crootd[ !is.na(rootd)  ] <- 1

  ctem[ 'TEM' %in% dtColname ] <- 2
  cVPD[ 'VPD' %in% dtColname ] <- 2
  csoilM[ 'soilM' %in% dtColname ] <- 2
  cLs[ 'Ls' %in% dtColname & !is.na(lat)   ] <- 2
  crootd[ 'crootd' %in% dtColname ] <- 2

  if( any(c(C1,ctem,cVPD,csoilM,cLs,crootd) == 0 )  ) {stop("Some data missing")}

  #### 补齐数据 ####
  if (C1 == 1 ) {
    climdata <- dplyr::mutate(climdata, DOY = as.POSIXlt(ISOdate(Year ,Month, Day))$yday+1,.after = 'Year' )
  }

  if (cVPD ==1 ) {
    climdata <- dplyr::mutate(climdata,VPD = 0.61078*exp(17.27*TEM/(TEM+237.3))*(1-RH))
    climdata$VPD[  climdata$TEM <= 0 ] <- 0
  }

  if (cLs == 1) {
    climdata <- Compute_dayLength(climdata , lat) |>
      dplyr::mutate( gE = L/max(L), Ls = L , dL_i = c( 0 , diff( L ) ) *12    )
  }

  if ( cLs == 2 ){
    climdata <- Compute_dayLength(climdata , lat) |>
      dplyr::mutate( gE = Ls/max(L),  dL_i = c( 0 , diff( L ) ) *12    )
  }

  if ( crootd == 1 ) {
    climdata <- Compute_rootd( climdata, rootd ,... ) ##
  }

  if ( csoilM == 1 ) {
    climdata <- Computer_soliM( climdata,... ) ##
  }
  return(climdata )
} ##Compute_clim end -----------------------

#### 光周期 ####
#' Compute_daylengthfactor
#'
#' @param climdata data
#' @param lat Latitude
#' @importFrom dplyr group_by reframe ungroup left_join
#'
#' @return day length and standardization day length
Compute_dayLength <- function(climdata,lat ){

  #单位转换与常数设置
  lats  <-  lat * pi / 180# change to radians 角度改为弧度制
  I_0  <-  118.109#太阳常数MJ/s

    #日序公转角度
    # theta_0  <-  2 * pi * (doy - 1) / DAL[ly]

    dt <- climdata |> dplyr::group_by(Year) |>
      dplyr::reframe(DOY = DOY , theta_0 = 2 * pi * (DOY - 1) /max(DOY)  ) |>
      dplyr::ungroup() |> data.frame()##

    #地球轨道偏心率订正因子
    dt$rho_2  <-  1.000110 + 0.034221 * cos(dt$theta_0) + 0.001280 * sin(dt$theta_0) +
      0.000719 * cos(2 * dt$theta_0) + 0.000077 * sin(2 * dt$theta_0)
    #赤纬/太阳偏角(180/pi)*
    dt$delta  <-  (0.006918-0.399912 * cos(dt$theta_0) + 0.070257 * sin(dt$theta_0) -
                  0.006758 * cos(2 * dt$theta_0) + 0.000970 * sin(2 * dt$theta_0)-
                  0.002697 * cos(3 * dt$theta_0) + 0.000148 * sin(3 * dt$theta_0))
    #时角
    dt$omega_0  <-  acos(-tan(lats) * tan(dt$delta))
    #最大日照时数
    dt$N  <-  24 / pi * dt$omega_0
    dt$L  <-  dt$N / 12


    climdata <- dplyr::left_join( climdata, dt[ ,c('Year', "DOY",'L')] )

    # #逐日天文辐射量
    # dt$s_0  <-  ((I_0) / pi) * dt$rho_2 * (dt$omega_0 * sin(lat) * sin(dt$delta) +
    #                                    cos(lat) * cos(dt$delta) * sin(dt$omega_0))

    # #日照百分率
    # if (is.na(fixparam.Mclim$n) == F ) {
    #   dynparam.Mclim$s  <-  n / N
    # }

    ## 无效部分在算gE时会被约掉
    # 总模型 A-P model
    # Ogelman model不同辐射模型的比较及其对APSIM模型模拟效果的影响
    # _毛洋洋[37]https://doi.org/10.1016/j.agrformet.2006.02.001
    # Q = s_0 *(a+b*s+c*(s^2))
    # dt$Q  <-  dplyr::case_when( SRM == "AP" ~ dt$s_0 *( paramSRM['a'] + paramSRM['b'] * paramSRM['s']),
    #                             SRM == "Oge" ~ dt$s_0 *( paramSRM['a'] + paramSRM['b'] * paramSRM['s'] + paramSRM['c']*paramSRM['s']^2) )
  return(climdata)

} ## 光周期 end -----


#### 活动土层深度  --------------------------------------------------
#' Active soil depth
#'
#' @param climdata climsdata
#' @param rootd0 root zone depth
#' @param paramSoilMelt Soil melt parameters, Where p1 is Soil thawing coefficient 1, p2 is Soil thawing coefficient 2.
#' @param Tm Temperature of soil melt
#'
#' @importFrom dplyr group_by mutate case_when select
#' @return Active soil depth

Compute_rootd  <-  function(climdata, rootd0, paramSoilMelt = c( p1 = NA , p2 = NA ), Tm = NA   ){

  ####
  ezRootd <- 0
  ezRootd[ !is.na(Tm) ] <- 1
  ezRootd[ all( !is.na( paramSoilMelt ) ) ] <- 2

  climdata <- climdata |> dplyr::group_by(Year) |>
    dplyr::mutate(rootd = 0,
                  winterIndx = dplyr::case_when(TEM > 0~1, TEM<=0 ~ 0),
                  winIdx = 0,Winter = NA
                  )

  climdata$winterIndx <- c( rep( 0,4  ),  sapply(5:length(climdata$winterIndx), function(i) sum(climdata$winterIndx[(i - 5 + 1):i]))  )
  climdata$winIdx[ climdata$winterIndx == 5 ] <- 1

  Soil_melt_switch <- ifelse( climdata$TEM[1] >0, T, F )
  winter_switch <- ifelse( climdata$TEM[1] <= 0, T, F )

  for (days_i in 2: nrow(climdata)) {

    ### 计算季节与土壤状态
    if ( winter_switch ) {
      ## 当前为冬季时
      cuntT <- 0
      winter_switch[climdata$winIdx[ days_i ] == 1 ] <- F ## 当日达到连续5天均温0上时关闭冬季
      climdata$Winter[ days_i ] <- "Ungrowth"
      climdata$Winter[ days_i ][ winter_switch==F  ] <- "Growth"
    }else{
      ## 打开冬季的条件
      winter_switch[climdata$winterIndx[ days_i ] == 0 ] <- T ## 当日达到连续5天均温0上时关闭冬季
      climdata$Winter[ days_i ] <- "Growth"
      climdata$Winter[ days_i ][ winter_switch==T ] <- "Ungrowth"

      if ( ezRootd == 1 ) {
        cuntT[ Soil_melt_switch == F ] <- cuntT +  max(climdata$TEM[ days_i ],0)
        Soil_melt_switch[ cuntT >= Tm  ] <- T
      }
    }
    ### 计算土壤解冻深度



    if( ezRootd == 2 ){
      ### 判断土壤状态
      if ( climdata$TEM[days_i] > 0 & winter_switch == F  ) {
        climdata$rootd[days_i]  <- min( (  climdata$rootd[days_i-1] + climdata$TEM[days_i] * paramSoilMelt['p1'] * exp( -paramSoilMelt['p2'] * climdata$rootd[days_i - 1])), ##VSM土壤解冻方程
                                     rootd0  ) ### 土壤开始融化
      }else{ climdata$rootd[days_i]  <-  0 } ### 土壤融化
    }
    if( ezRootd == 1 ){
      ### 判断土壤状态
      if (Soil_melt_switch ) {
        climdata$rootd[days_i]  <- rootd0 ### 土壤融化开关为 TRUE 时土壤开始融化
      }else{ climdata$rootd[days_i]  <-  0 } ### 土壤融化开关为 FALSE 时土壤开始融化
    } ##
    if( ezRootd == 0 ){ climdata$rootd[days_i][ winter_switch == F ]  <- rootd0 }


  }
  climdata<- dplyr::select(climdata, -winterIndx, -winIdx, -Winter )
  climdata$rootd <-  round(climdata$rootd)
  return(climdata)
} ## 活动土层 end --------------------

#### 雨雪折算  -----------------------------------------------------------
#' Computer snow to precipitation
#'
#' @param climdata climsdata
#'
#' @return climsdata

Compute_P <-  function(climdata  ){  ## 雨雪模块 start
  P_snow <- 0
  ##降雪过程
  for (t in 1:nrow(climdata)) {###日循环
    P_snow[ climdata$TEM[t] <= 0  ] <- climdata$PRE[t] + P_snow
    climdata$PRE[t][ climdata$TEM[t] <= 0  ] <- 0
    climdata$PRE[t][ climdata$TEM[t] > 0 & P_snow != 0 ] <- climdata$PRE[t] + P_snow
    P_snow[ climdata$TEM[t] > 0 & P_snow != 0  ] <- 0
  }
  return(climdata)
} ## 雨雪模块 end ------------------------------------------------------------


#### Using Precipitation to Estimate Soil Moisture ####
#### cpc-leaky bucket
#' cpc-leaky bucket
#'
#' The main formula is based on VSLiteR cpc-leaky bucket model (submonthly version) and has been modified to fit the daily value data.
#' VSLiteR:[https://github.com/suztolwinskiward]
#'
#' @param climdata Climsdata
#' @param paramSoilM The parameters related to estimation of soil moisture, where 'M0' is the initial soil moisture, 'dp' is the rainfall infiltration step, 'alph' is the scalar runoff , 'Mmax' is the maximum soil water holding capacity, 'Mmin' is the minimum soil moisture, 'mu.th' is the subsurface runoff parameter, 'm.th' is the surface runoff parameter.
#'
#'
#' @return gE , SoilM soil moisture, Ls, potEv et.al. climate factor
#' SoilM soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
#' potEv Potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm).
#'
Computer_soliM <- function(climdata,
                           paramSoilM = c(M0 = 0.35, dp = 2, alph = 0.0013 , Mmax = 0.6, Mmin = 0.042 ,
                                          mu.th = 1.09, m.th = 4.1 ) ){

  climdata <- Compute_P(climdata)
  summaryMicroClim <- logical()

  for ( y in unique(climdata$Year) ) {  ## 按年拆分数据循环
    ydata <- climdata[ climdata$Year == y, ]

    soilM <- potEv <- gT <- gM <- matrix(NA, nrow(ydata), 1)
    deltaWaterFlow <- CG <- CR1 <- CR2 <- matrix(0, nrow(ydata), 1)

    istar <- (ydata$TEM / 5) ^ 1.514 #1.514 = 53/35
    istar[ydata$TEM < 0] <- NA
    I <- mean(istar,na.rm=T)
    ap <- (6.75e-7) * I ^ 3 - (7.71e-5) * I ^ 2 + (1.79e-2) * I + 0.49

    # attach(ydata$)

    if(paramSoilM['M0'] < 0 ){ paramSoilM['M0'] <- 0.3 } #核对并修正初始值在有效范围

    for (day in 1: nrow(ydata)) {

      if (ydata$rootd[day] > 0) {  ##土壤解冻后在进行计算

        if (ydata$TEM[day] <= 0){Ep <- 0}
        if (ydata$TEM[day] > 0 && ydata$TEM[day] < 26.5){Ep <- 16 * ydata$Ls[day] * (10 * ydata$TEM[day] / I) ^ ap}
        if (ydata$TEM[day] >= 26.5){Ep <-  -415.85 + 32.25*ydata$TEM[day] - 0.43* ydata$TEM[day]^2}
        Ep  <-  Ep/30.42

        potEv[day]  <-  Ep

        #参数初始化：
        #dp = 2.0 # mm of precip per increment 每次计算步长递增降雨量
        nstep <- floor(ydata$PRE[day]/paramSoilM['dp']) +1 # number of sub-monthly substeps 步进计算长度
        Pinc <- ydata$PRE[day]/nstep # precip per substep每个子步骤
        alphinc <- paramSoilM['alph'] /nstep # runoff rate per substep time interval每个子步时间间隔的径流率
        Epinc <- Ep/nstep # potential evapotrans per substep.每个子步骤潜在的蒸发量。

        sm0 <- paramSoilM['M0']##初始土壤湿度= M0 ##首年第一天为直接赋值M0                ***

        for(istep in 1:nstep){  ## 步长计算湿度
          #计算蒸散：潜在蒸散Ep*(土壤湿度/最大持水量) 注：不同林分使用同样潜在蒸散量是否合适
          Etrans <- Epinc*sm0*ydata$rootd[day]/(paramSoilM['Mmax']*ydata$rootd[day])

          #计算下渗G：μ*α/(1+μ)*土壤湿度 ##实测排水速度为0.138mm/day
          G  <-  ( paramSoilM['mu.th']*alphinc/(1+paramSoilM['mu.th'])*sm0*ydata$rootd[day] )# /30

          #计算径流R，地面径流+地下径流
          #地表径流=降雨*（土壤湿度/最大持水量）^m系数
          #地下径流=alpha系数/（1+μ系数）*土壤湿度
          R1  <-  ( Pinc*(sm0*ydata$rootd[day]/(paramSoilM['Mmax']*ydata$rootd[day]))^paramSoilM['m.th'] )# /30# +
          R2  <-  ( (alphinc/(1+paramSoilM['mu.th']))*sm0*ydata$rootd[day] ) #/30

          # 总模型：土壤湿度=降雨（Pinc）- 蒸散（Etrans）-径流（R）-下渗（G）
          dWdt <- Pinc - Etrans - R1 - R2 - G #总模型
          sm1 <- sm0 + dWdt/ydata$rootd[day] #根深湿度
          ##
          deltaWaterFlow[day]  <-  deltaWaterFlow[day] + dWdt
          CG[day]  <-  CG[day] + G   ##
          CR1[day]  <-  CR1[day] + R1
          CR2[day]  <-  CR2[day] + R2
          ##

          ##初步erro catch
          sm0 <- max(sm1,paramSoilM['Mmin'])
          sm0 <- min(sm0,paramSoilM['Mmax'])

        } ## 湿度步进计算结束

        soilM[day]  <-  sm0##不知道是否有效
        # error-catching:错误捕捉
        if (soilM[day] <= paramSoilM['Mmin']){soilM[day]  <-  paramSoilM['Mmin'];}
        if (soilM[day] >= paramSoilM['Mmax']){soilM[day]  <-  paramSoilM['Mmax'];}
        if (is.na(soilM[day])==1){soilM[day]  <-  paramSoilM['Mmin'];}

        paramSoilM['M0']  <-  soilM[day]

      }else{soilM[day]  <-  0;potEv[day]  <-  0} ##rootd[t]判断结束

    } ## 年内
    microclim <- data.frame(soilM, potEv ,CG,CR1,CR2,deltaWaterFlow   )  ## ,gT,gM
    summaryMicroClim <- rbind(summaryMicroClim,microclim)
    # detach(ydata)
  }  ##  年循环 out
  climdata <- cbind(climdata,summaryMicroClim) #|> select(c("Year", "Month","Day","DOY","gE","gT","gM","soilM","TEM","Ls","dailyPrec"))
  return(climdata)
} ## soliM end --------







