#' normalization data
#'
#' @param dat a numeric data.
#' @param TPH a logical value, if TPH is TRUE, calculate the data as 0-100, otherwise calculate as 0-1.
#' @param ROU a logical value, if TPH is TRUE, data will be maintained to two decimal places.
#' @param min1 a logical value,if TPH is TRUE, data minimum value is not zero.
#' @param Zeros a logical value,if TPH is TRUE, raw data minimum value set to 0.
#'
#' @return normalization data.
#' @export

nor <- function( dat,TPH = F ,ROU = F,min1 = F, Zeros = F ){
  a <- 0
  if (min1 == T) {
    a <- 0.01
  }

  dat <- as.numeric(dat)

  mind <- min(dat,na.rm = T)
  mind[Zeros == T]  <- 0

  nordat <- (dat-mind +a)/(max(dat,na.rm =T) - mind +a)

  if (TPH == T) {
    nordat <- nordat * 100
  }
  if (ROU == T) {
    nordat <- round(nordat,2)
  }

  return(nordat)
} ## nor end -------


#' Model test
#'
#' @param y_actual Observes data
#' @param y_predicted Simulates data
#' @param method cor.test method
#'
#' @return test result
#' @export

mod.test<-function(y_actual,y_predicted,method = "pearson", ... ){
  avr_y_actual <- mean(y_actual)

  ## SSR
  ss_regression <- sum( (y_predicted - avr_y_actual)^2)

  ## SSE
  ss_residuals <- sum( (y_actual - y_predicted)^2 )

  ## SST
  ss_total <- sum( (y_actual - avr_y_actual)^2 )


  ####
  ss_total2 <- ss_regression + ss_residuals
  ss_total3 <-  ss_total + ss_residuals

  rsquare <-  1- ss_residuals / ss_total

  #return(rsquare)#当模型偏差过大，rsquare很小时，不采用rsquare统计
  n1<-length(y_actual)
  n2<-length(y_predicted)#要求n1 == n2
  meansquare <- ss_residuals/(n1-2)
  #参考王辰勇译《线性回归分析导论》12页
  #return(meansquare)#MS残
  COR <- cor.test(y_actual, y_predicted, method = method,... )$estimate
  pvalue <- cor.test(y_actual, y_predicted, method = method ,...)$p.value

  RMSE<-(ss_residuals/n1)^0.5
  NRMSD<-RMSE/avr_y_actual


  MAE <- mean( abs(y_predicted- y_actual)  )
  MAPE <- mean( abs( (y_predicted- y_actual)/y_actual )   ) *100
  # res <- data.frame(rsquare,rsquare1,rsquare2,rsquare3,meansquare,RMSE,NRMSD,COR,MAPE  )
  res <- data.frame(rsquare,meansquare, COR,pvalue,RMSE,NRMSD, MAE, MAPE )

  return(res)
}

#' Develop Outputs
#'
#' @param Result rBTR Outputs list
#'
#' @return DevResult result
#' @export
#'

dev.outputs <- function(BtrResList = Outputs){

  BtrResList[['xylem_trait']] <- BtrResList[['xylem_trait']] |>
    dplyr::left_join( BtrResList[['annaulRing']] |> dplyr::select(Year,RingWidth,CellLayer,VesselNumber) ) |>
    dplyr::mutate( rrFiber =  (Raddist+0.5*CRD) /1000/RingWidth*100,
                   rrVessel =  (Raddist+0.5*VCRD) /1000/RingWidth*100 ,
                   rrLayer = cell_L/CellLayer*100 ,
                   rrVN = VNoV/VesselNumber*100,100  ) |>
    dplyr::select(-RingWidth,-CellLayer,-VesselNumber )
  BtrResList$xylem_trait$rrFiber[ BtrResList$xylem_trait$rrFiber >100 ] <- 100
  BtrResList$xylem_trait$rrVessel[ BtrResList$xylem_trait$rrVessel >100 ] <- 100

  return(BtrResList)

}

#' Model test Cells
#'
#' @param devOutputs Simulates data
#' @param OBSdata Observes data From ROXAS cell sheet
#' @param cor.test time windows (e.g. c(1960:2020) ).
#' @param celltype 'vessel' or 'fiber'
#' @param CellParam cell parameter e.g. "LA","WT","WA"
#' @param method cor.test method
#' @importFrom mgcv gam predict.gam
#' @importFrom dplyr rename select filter mutate
#'
#' @return test result
#' @export

mod.test.cells <- function( devOutputs=BtrResList , OBSdata ,
                            celltype = 'vessel',CellParam = "LA", Years , method = "pearson", k = 5, ... ){

  # head(Simdata)


  MtRes <- matrix(ncol = 11, nrow = length(Years) ) |> as.data.frame()
  # MtRes <- logical()

  if ( celltype == 'vessel') {
    Simdata <- devOutputs[["xylem_trait"]] |>
      dplyr::select(Year,cell_L,VCA, VCV,VWA, VWT, rrVessel, rrVN) |>
      dplyr::rename( CA = VCA, LA = VCV ,WA = VWA , CWTall =VWT, RR = rrVessel, RN = rrVN ) |>
      na.omit()
  } else {
    Simdata <- devOutputs[["xylem_trait"]] |>
      dplyr::select(Year,cell_L,CA,LA= CV,WA,CWTall= WT, rrFiber , rrLayer ) |>
      dplyr::rename( RR = rrFiber , RN = rrLayer  )
  }


  i=1
  for ( corTimes in Years) {

    GAMdata <- OBSdata |> dplyr::filter(Year == corTimes ) |>
      dplyr::select( Year,eval( CellParam),RRadDistR ) |> na.omit()
    colnames( GAMdata ) <- c('Year', 'y' , 'RR')

    Resdata <- Simdata |> dplyr::filter(Year == corTimes ) |>
      dplyr::select( Year,eval( CellParam),RR ) |> na.omit()
    colnames( Resdata ) <- c('Year', 'y' , 'RR')

    if (nrow(GAMdata ) != 0) {
      GAM1 <- mgcv::gam( y ~ s(RR , k = k )  , data =  GAMdata  )

      Resdata$trend <- mgcv::predict.gam(object = GAM1 , newdata = Resdata )

      res  <- mod.test( y_actual = Resdata$trend, y_predicted = Resdata$y ,method = method,...  ) |>
        dplyr::mutate( Year = corTimes,celltype = celltype, CellParam=CellParam,  .before =1 )

      MtRes[i ,] <- res
      i=i+1
    } else {
      warning(  paste("Missing Year:", corTimes ))
    } ## end if else nrow(GAMdata )




  } ## end for corTimes



  MtRes <- na.omit(MtRes)
  colnames(MtRes) <- colnames(res)

  MtRes <- rbind(MtRes , c( 'all',  celltype,  CellParam , apply(MtRes[4:11],2,mean)))


  return(MtRes)

} ## end func





