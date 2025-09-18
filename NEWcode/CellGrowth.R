#' cells_growth function
#'
#' @param cell one of growthing cell layer
#' @param CorV cells or vessels
#' @param clim.today Today climate data
#' @param layer.max max layer in grwoth area ## Del
#' @param fixparam.growth.fiber Fix parameters
#' @param fixparam.growth.vessel Fix parameters
#' @param dynparam.growth.t Dynamic parameters
#'
#' @return today cells anatomy
#'
## #' @export

cells_growth <- function( Doy_GR, cell , CorV = "C" ,  param, i_doy  ){ ## function start

  if(CorV == "C"){
    cell <- Doy_GR[,.(Year, vc = v_c.fiber,vw = v_w.fiber, vl = v_l.fiber,CAmax)][cell,on = 'Year'   ][,CAmax := param$CAmax ]
  }else{
    data.table::setnames(cell,c('VCA','VCV','VWA','VLWA','VWT','VCRD','VCTD','VEDOY','VTDOY','VDDOY'),
                         c('CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY'))
    cell <- Doy_GR[,.(Year, vc = v_c.vessel,vw = v_w.vessel, vl = v_l.vessel,CAmax)][cell,on = 'Year'   ]
  }

  ## CAt max 细胞壁弹性允许的最大CA
    CPeri <- cell$WA/param$WTmin

    ifelse(CorV == "C",
           CAtmax <- cell$CRD *(0.5*CPeri-cell$CRD) ,
           CAtmax <- (0.25*CPeri)^2 )

    limitWTa <- cell$WT/param$WTa
    limitWTa[ limitWTa>1 ] <- 1

    ## CA 细胞面积

    dCA_dt <-  cell$vc * cell$CA * ( 1 - cell$CA/cell$CAmax ) * ( 1- limitWTa )
    dCA_dt[dCA_dt < 0 ] <- 0 ## error catching

    cell$TDOY[dCA_dt == 0 & cell$TDOY == 0 & cell$CA > 0 ] <- i_doy

    dCA_dt[ cell$TDOY != 0  ] <- 0 ## error catching

    ### WA 细胞壁面积

    dWA_dt <- cell$vw * ( 1 - cell$WA/param$WAmax ) * ( 1 - 1/( 1 + (cell$CA-cell$WA)/param$mw )^param$sw  )

    dWA_dt[dWA_dt < 0 ] <- 0 ## error catching
    dWA_dt[cell$WT >= param$WTmax ] <- 0 ## error catching


    ### 细胞壁生长大于细胞总膨胀时扩大结束 dCA <= 0
    cell$TDOY[ (dCA_dt - dWA_dt ) <=0 & cell$TDOY == 0 & cell$CA > 0 ] <- i_doy
    cell$TDOY[ (dCA_dt - dWA_dt ) > 0 & cell$TDOY != 0 ] <- 0 ## 可逆的

    ### LWA 细胞壁木质化量
    dLWA_dt <- cell$vl *  ( 1- 1/( 1 + (cell$CA - cell$WA )/param$ml )^param$sl  )

    dLWA_dt[dLWA_dt > 0 &  dLWA_dt < 0.05 ] <- 0.05 ## error catching
    dLWA_dt[dLWA_dt < 0 ] <- 0 ## error catching
    dLWA_dt[cell$LWA >= cell$WA ] <- 0 ## error catching

    ###  summary
    ##
    dCW <- cell$CA - cell$WA

    cell$CA <- cell$CA +round(dCA_dt,3 )

    cell$CA[cell$CA >= CAtmax ] <- CAtmax[cell$CA >= CAtmax] ## error catching


    cell$WA <- cell$WA + round(dWA_dt,3)

    cell$WA[cell$WA > param$WAmax ] <- param$WAmax ## error catching
    cell$WA[cell$WA > cell$CA] <- cell$CA[cell$WA > cell$CA] ## error catching

    cell$LWA <-  cell$LWA + round(dLWA_dt,3)

    cell$LWA[cell$LWA > cell$WA ] <- cell$WA[cell$LWA > cell$WA ]

    ## DDOY & TODY

    cell$DDOY[ cell$LWA >= cell$WA   & cell$DDOY == 0 ] <- i_doy

    cell$TDOY[ cell$WA  >= param$WAmax & cell$TDOY == 0] <- i_doy

    cell$TDOY[ cell$WT  >= param$WTmax & cell$TDOY == 0] <- i_doy

    # cell$DDOY[cell$LWA < cell$WA & cell$DDOY == 0 & dLWA_dt == 0 ] = clim.today$DOY+9000

    ## WT & CV & CRD 径向细胞大小

    cell$CTD[CorV == "V"] <- floor(cell$CA ^ 0.5)
    cell$CRD <- cell$CA/cell$CTD##


    cell$WT <- ( 2*(cell$CTD + cell$CRD) - (4*(cell$CTD + cell$CRD)^2 - 16 *cell$WA)^0.5) / 8

    cell$CV <- cell$CA - cell$WA

    cell[ is.na(cell) ] <- 0

    if(CorV == "V") {
      data.table::setnames(cell,c('CA','CV','WA','LWA','WT','CRD','CTD','EDOY','TDOY','DDOY'),
                           c('VCA','VCV','VWA','VLWA','VWT','VCRD','VCTD','VEDOY','VTDOY','VDDOY')
                           )
    }
    cell[,c('vc','vw','vl','CAmax'):=NULL]
  return(cell)

} ### FUNCTION END -------------------------------
