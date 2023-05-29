##### 微气候 #####

#### cpc-leaky bucket
#' cpc-leaky bucket
#'
#' @param climdata climsdata
#' @param parameters parameters
#' @param syear start year
#' @param eyear end year
#'
#' @return gE , gT ,gM , Ls et.al. climate factor
#'
#' @importFrom dplyr filter select mutate
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#'
#' @export

Compute_clim = function(climdata , parameters, syear = NA, eyear = NA ){

  #### form MainFun
  parameters$values <- as.numeric(parameters$values)

  fixparam.Mclim <- dplyr::filter(parameters ,  modul == "microclim" & paramtype == "fixed") %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = parameter, value = values)

  dynparam.Mclim <- dplyr::filter(parameters ,  modul == "microclim" & paramtype == "dynamic") %>%
    dplyr::select( c("parameter","values") ) %>% tidyr::spread( key = parameter, value = values)

  ## error-catching
  if (is.na(syear) ) {syear = min(climdata$Year) } else {
    if( syear < min(climdata$Year) ){
      stop("syear is wrong")}
  }

  if (is.na(eyear)) {eyear = max(climdata$Year)} else {
    if( eyear > max(climdata$Year) ){
      stop("syear is wrong")
    }
  }
  ## 提取对应年份数据
  climdata <- dplyr::filter( climdata, Year >= syear & Year <= eyear)

  summaryMicroClim <- vector()

  iyear <- syear:eyear

  # 合并数据
  rootd <- Compute_rootd(climdata, fixparam.Mclim, dynparam.Mclim)
  dailyPrec <- Compute_P(climdata, fixparam.Mclim, dynparam.Mclim,rootd)

  gE <- Compute_gE(climdata, fixparam.Mclim, dynparam.Mclim,syear,eyear)##[,c(-1,-2)] ##syear, eyear,

  data2 <- data.frame(gE, dailyPrec, rootd)

  rm("dailyPrec" , "rootd" , "gE")

  for (y in syear:eyear) {  ## 按年拆分数据循环
    ydata <- dplyr::filter(data2, Year == y )

    soilM <- potEv <- gT <- gM <- matrix(NA, nrow(ydata), 1)
    deltaWaterFlow <- CG <- CR1 <- CR2 <- matrix(0, nrow(ydata), 1)

    # attach(ydata$)

    if(dynparam.Mclim$M0 < 0 ){ dynparam.Mclim$M0 = 300 / fixparam.Mclim$rootd0} #核对并修正初始值

    for (day in 1: nrow(ydata)) {

      if (ydata$rootd[day] > 0) {  ##土壤解冻后在进行计算

        istar <- (ydata$TEM / 5) ^ 1.514 #1.514 = 53/35
        istar[ydata$TEM < 0] <- NA
        I <- mean(istar,na.rm=T)
        ap <- (6.75e-7) * I ^ 3 - (7.71e-5) * I ^ 2 + (1.79e-2) * I + 0.49

        if (ydata$TEM[day] <= 0){Ep = 0}
        if (ydata$TEM[day] > 0 && ydata$TEM[day] < 26.5){Ep <- 16 * ydata$Ls[day] * (10 * ydata$TEM[day] / I) ^ ap}
        if (ydata$TEM[day] >= 26.5){Ep <-  -415.85 + 32.25*ydata$TEM[day] - 0.43* ydata$TEM[day]^2}
        Ep  <-  Ep/30.5

        potEv[day]  <-  Ep

        #参数初始化：
        #dp = 2.0 # mm of precip per increment 每次计算步长递增降雨量
        nstep <- floor(ydata$dailyPrec[day]/fixparam.Mclim$dp) +1 # number of sub-monthly substeps 步进计算长度
        Pinc <- ydata$dailyPrec[day]/nstep # precip per substep每个子步骤
        alphinc <- fixparam.Mclim$alph /nstep # runoff rate per substep time interval每个子步时间间隔的径流率
        Epinc <- Ep/nstep # potential evapotrans per substep.每个子步骤潜在的蒸发量。

        sm0 <- dynparam.Mclim$M0##初始土壤湿度= M0 ##首年第一天为直接赋值M0                ***

        for(istep in 1:nstep){  ## 步长计算湿度
          #计算蒸散：潜在蒸散Ep*(土壤湿度/最大持水量) 注：不同林分使用同样潜在蒸散量是否合适
          Etrans <- Epinc*sm0*ydata$rootd[day]/(fixparam.Mclim$Mmax*ydata$rootd[day])

          #计算下渗G：μ*α/(1+μ)*土壤湿度 ##实测排水速度为0.138mm/day
          G  <-  ( fixparam.Mclim$mu.th*alphinc/(1+fixparam.Mclim$mu.th)*sm0*ydata$rootd[day] )# /30

          #计算径流R，地面径流+地下径流
          #地表径流=降雨*（土壤湿度/最大持水量）^m系数
          #地下径流=alpha系数/（1+μ系数）*土壤湿度
          R1  <-  ( Pinc*(sm0*ydata$rootd[day]/(fixparam.Mclim$Mmax*ydata$rootd[day]))^fixparam.Mclim$m.th )# /30# +
          R2  <-  ( (alphinc/(1+fixparam.Mclim$mu.th))*sm0*ydata$rootd[day] ) #/30

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
          sm0 <- max(sm1,fixparam.Mclim$Mmin)
          sm0 <- min(sm0,fixparam.Mclim$Mmax)

        } ## 湿度步进计算结束

        soilM[day]  <-  sm0##不知道是否有效
        # error-catching:错误捕捉
        if (soilM[day] <= fixparam.Mclim$Mmin){soilM[day]  <-  fixparam.Mclim$Mmin;}
        if (soilM[day] >= fixparam.Mclim$Mmax){soilM[day]  <-  fixparam.Mclim$Mmax;}
        if (is.na(soilM[day])==1){soilM[day]  <-  fixparam.Mclim$Mmin;}

        dynparam.Mclim$M0  <-  soilM[day]


      }else{soilM[day]  <-  0;potEv[day]  <-  0} ##rootd[t]判断结束


    } ## 年内


    microclim <- data.frame(soilM, potEv ,CG,CR1,CR2,deltaWaterFlow   )  ## ,gT,gM
    summaryMicroClim <- rbind(summaryMicroClim,microclim)
    # detach(ydata)
  }  ##  年循环 out

  mclim <- cbind(data2,summaryMicroClim) #%>% select(c("Year", "Month","Day","DOY","gE","gT","gM","soilM","TEM","Ls","dailyPrec"))

  # detach(fixparam.Mclim)

  dL_i <- c( 0 , diff( mclim$Ls ) )

  Microclim <- cbind( mclim , dL_i*12 )

  Microclim <- dplyr::mutate(Microclim,VPD = 0.61078*exp(17.27*TEM/(TEM+237.3))*(1-RH))


  return(Microclim)

} ##Compute_clim end -----------------------


### 标准太阳辐射计算 ----------------------------------------------------------
# Compute_daylengthfactor <-function(phi, a, b, n, s){
#' Compute_daylengthfactor
#'
#' @param fixparam.Mclim Fix parameters
#' @param dynparam.Mclim Dynamic parameters
#'
#' @return day length and standardization day length

compute_daylengthfactor <-function( fixparam.Mclim, dynparam.Mclim ){  ## daylength start

  gE  <-  matrix(NaN, nrow =  366, ncol =  4)
  colnames(gE)  <-  c('gE365', 'gE366', 'L365', 'L366')

  #单位转换与常数设置
  lat  <-  fixparam.Mclim$phi * pi / 180# change to radians 角度改为弧度制
  I_0  <-  118.109#太阳常数MJ/s
  DAL  <-  c(365, 366)
  for (ly in 1:2) {
    doy  <-  (1:DAL[ly])
    #日序公转角度
    theta_0  <-  2 * pi * (doy - 1) / DAL[ly]
    #地球轨道偏心率订正因子
    rho_2  <-  1.000110 + 0.034221 * cos(theta_0) + 0.001280 * sin(theta_0) +
      0.000719 * cos(2 * theta_0) + 0.000077 * sin(2 * theta_0)
    #赤纬/太阳偏角(180/pi)*
    delta  <-  (0.006918-0.399912 * cos(theta_0) + 0.070257 * sin(theta_0) -
               0.006758 * cos(2 * theta_0) + 0.000970 * sin(2 * theta_0)-
               0.002697 * cos(3 * theta_0) + 0.000148 * sin(3 * theta_0))
    #时角
    omega_0  <-  acos(-tan(lat) * tan(delta))
    #最大日照时数
    N  <-  24 / pi * omega_0
    L  <-  N / 12
    #逐日天文辐射量
    s_0  <-  ((I_0) / pi) * rho_2 * (omega_0 * sin(lat) * sin(delta) +
                                    cos(lat) * cos(delta) * sin(omega_0))
    #日照百分率
    if (is.na(fixparam.Mclim$n) == F ) {
      dynparam.Mclim$s  <-  n / N
    }

    #总模型 A-P model
    Q  <-  s_0 *(fixparam.Mclim$a + fixparam.Mclim$b * dynparam.Mclim$s)
    # Ogelman model不同辐射模型的比较及其对APSIM模型模拟效果的影响
    # _毛洋洋[37]https://doi.org/10.1016/j.agrformet.2006.02.001
    # Q = s_0 *(a+b*s+c*(s^2))
    ndl  <-  Q / max(Q)
    #gE
    for (t in doy){
      gE[t, ly] <-  ndl[t]
      gE[t, ly + 2]  <-  L[t]
    }

  }

  return(gE)

} ## daylength end -----------------



#### 太阳辐射因子 gE ------------------------------------------------------
#' Sun radiation factor gE
#'
#' @param climdata climsdata
#' @param fixparam.Mclim Fix parameters
#' @param dynparam.Mclim Dynamic parameters
#' @param syear start year
#' @param eyear end year
#'
#' @importFrom dplyr rename left_join
#' @importFrom data.table rbindlist
#' @importFrom tibble rownames_to_column
#'
#' @return Sun radiation factor gE

Compute_gE <- function(climdata, fixparam.Mclim, dynparam.Mclim, syear,eyear ){  ## gE start

  basedata  <-  Compute_daylengthfactor( fixparam.Mclim, dynparam.Mclim ) %>% data.frame() %>% tibble::rownames_to_column("DOY")
  basedata$DOY <- as.numeric(basedata$DOY )
  d365  <-  basedata[,c(1,2,4)] %>% dplyr::rename(gE = `gE365`, Ls = `L365`)
  d366  <-  basedata[,c(1,3,5)] %>% dplyr::rename(gE = `gE366`, Ls = `L366`)

  dclim <- base::split(climdata,climdata$Year)

  for (i in 1:length(dclim)) {
    if ( nrow(dclim[[i]]) == 365   ) {
      dclim[[i]] <- dplyr::left_join(dclim[[i]], d365)
    } else {
      dclim[[i]] <- dplyr::left_join(dclim[[i]], d366)
    }

  }

  gE <- data.table::rbindlist(dclim)

  return(gE)
} ## Compute_gE end -----------------------------

#### 活动土层深度  --------------------------------------------------
#' Active soil depth
#'
#' @param climdata climsdata
#' @param fixparam.Mclim Fix parameters
#' @param dynparam.Mclim Dynamic parameters

#'
#' @return Active soil depth

Compute_rootd  <-  function(climdata, fixparam.Mclim, dynparam.Mclim ){
  ####
  d_T  <-  climdata$TEM
  len  <-  length(d_T)
  dep  <-  matrix(NA, len, 1)##生成数据集
  ####
  for (t in 2:len) {
  # for (t in 1:len) {
    ###判断土壤状态
    if (dynparam.Mclim$Soil_melt_switch == 1 ) {
      dep[t]  <-  fixparam.Mclim$rootd0

      ### 逐渐溶解？


    }else{
      if (d_T[t] >= 0) {
        dep[t]  <-  (dep[t-1] + d_T[t] * fixparam.Mclim$p11 * exp( -fixparam.Mclim$p12 * dep[t - 1]))##VSM土壤解冻方程
        if (dep[t] < fixparam.Mclim$rootd0 ){dep[t]  <-  0 }
        if (dep[t] >= fixparam.Mclim$rootd0) { dep[t]  <-  fixparam.Mclim$rootd0 }# error catching


      }else{
        dep[t]  <-  0
      }

    }

    ###计算第二天土壤状态
    if (d_T[t] < 0) {
      dynparam.Mclim$winter  <-  dynparam.Mclim$winter + 1
      if (dynparam.Mclim$winter == 10) { dynparam.Mclim$Soil_melt_switch  <-  0; dynparam.Mclim$winter  <-  0; dynparam.Mclim$Tsum  <-  0 }
    }else{
      dynparam.Mclim$Tsum  <-  dynparam.Mclim$Tsum + d_T[t]
      if (dynparam.Mclim$Tsum > fixparam.Mclim$Tm) { dynparam.Mclim$Soil_melt_switch  <-  1 ; dynparam.Mclim$winter  <-  0}
    }

  }

  rootd  <-  round(dep)
  return(rootd)
} ## 活动土层 end --------------------



#### 雨雪折算  -----------------------------------------------------------
#' Computer snow to precipitation
#'
#' @param climdata climsdata
#' @param fixparam.Mclim Fix parameters
#' @param dynparam.Mclim Dynamic parameters
#' @param rootd  Active soil depth
#'
#' @return precipitation

Compute_P <-  function(climdata, fixparam.Mclim, dynparam.Mclim ,rootd){  ## 雨雪模块 start
  ##climdata格式为：year	month	day	DOY	TEM	PRE
  len  <-  length(climdata$Year)
  d_T <-  climdata$TEM
  d_Pre <-  climdata$PRE
  ####
  P_all <-  matrix(0, len, 1) ##土壤补水量
  P_sum  <-  0
  ##降雪过程
  for (t in 1:len) {###日循环
    if (rootd[t] <= 0   ) {
      P_all[t]  <-  0
      P_sum  <-  P_sum + d_Pre[t]
    } else {
      P_all[t]  <-  d_Pre[t] + P_sum
      P_sum  <-  0
    }


  }

  P_input  <-  P_all
  return(P_input)
} ## 雨雪模块 end ------------------------------------------------------------




