#' cells_growth function
#'
#' @param cell one of growthing cell layer
#' @param CorV cells or vessels
#' @param clim.today Today climate data
#' @param layer.max max layer in grwoth area
#' @param fixparam.growth.fiber Fix parameters
#' @param fixparam.growth.vessel Fix parameters
#' @param dynparam.growth.t Dynamic parameters
#'
#' @return today cells anatomy

cells_growth <- function( cell = GC, CorV = "C" , clim.today, layer.max,
                          fixparam.growth.fiber,fixparam.growth.vessel, dynparam.growth.t ){ ## function start

  if(CorV == "C"){

    FGP <- fixparam.growth.fiber

    vc <- dynparam.growth.t$v_c.fiber

    vw <- dynparam.growth.t$v_w.fiber

    vl <- dynparam.growth.t$v_l.fiber

  }else{

    FGP <- fixparam.growth.vessel

    vc <- dynparam.growth.t$v_c.vessel

    vw <- dynparam.growth.t$v_w.vessel

    vl <- dynparam.growth.t$v_l.vessel

  }

  Death <- rep(1,nrow(cell))

  for(i in 1:10){
    ### limit factor

    ## CAt max
    CPeri <- cell$WA/FGP$WTmin

    ifelse(CorV == "C",
           CAtmax <- cell$CRD *(0.5*CPeri-cell$CRD) ,
           CAtmax <- (0.25*CPeri)^2 )

    limitWTa <- cell$WT/FGP$WTa
    limitWTa[ limitWTa>1 ] <- 1

    ## Death

    Death[cell$DDOY != 0] <- 0


    ## CA 细胞面积

    dCA_dt <- vc * cell$CA * ( 1 - cell$CA/FGP$CAmax ) * ( 1- limitWTa ) * Death

    dCA_dt[dCA_dt < 0 ] <- 0 ## error catching

    cell$TDOY[dCA_dt == 0 & cell$TDOY == 0 ] <- clim.today$DOY

    dCA_dt[ cell$TDOY != 0 ] <- 0

    ### WA 细胞壁面积

    dWA_dt <- vw * ( 1 - cell$WA/FGP$WAmax ) * ( 1 - 1/( 1 + (cell$CA-cell$WA)/FGP$mw )^FGP$sw  ) * Death

    dWA_dt[dWA_dt < 0 ] <- 0 ## error catching
    dWA_dt[cell$WT >= FGP$WTmax ] <- 0 ## error catching


    ### LWA 细胞壁木质化量
    dLWA_dt <- vl *  ( 1- 1/( 1 + (cell$CA - cell$WA )/FGP$ml )^FGP$sl  )*Death

    dLWA_dt[dLWA_dt > 0 &  dLWA_dt < 0.05 ] <- 0.05 ## error catching
    dLWA_dt[dLWA_dt < 0 ] <- 0 ## error catching
    dLWA_dt[cell$LWA >= cell$WA ] <- 0 ## error catching


    ###  summary
    ##
    dCW <- cell$CA - cell$WA

    cell$CA <- cell$CA + dCA_dt

    cell$CA[cell$CA >= CAtmax] <- CAtmax[cell$CA >= CAtmax] ## error catching


    cell$WA <- cell$WA + dWA_dt

    cell$WA[cell$WA > FGP$WAmax ] <- FGP$WAmax ## error catching
    cell$WA[cell$WA > cell$CA] <- cell$CA[cell$WA > cell$CA] ## error catching

    cell$LWA <-  cell$LWA + dLWA_dt

    cell$LWA[cell$LWA > cell$WA] <- cell$WA[cell$LWA > cell$WA]

    ## DDOY & TODY

    cell$DDOY[cell$LWA >= cell$WA & cell$DDOY == 0 ] <- clim.today$DOY

    cell$TDOY[ cell$WA >= FGP$WAmax & cell$TDOY == 0 ] <- clim.today$DOY

    cell$TDOY[ cell$WT >= FGP$WTmax & cell$TDOY == 0 ] <- clim.today$DOY


    ## WT & CV & CRD 径向细胞大小

    cell$CTD[CorV == "V"] <- floor(cell$CA ^ 0.5)
    cell$CRD <- cell$CA/cell$CTD##

    cell$WT <- ( 2*(cell$CTD + cell$CRD) - (4*(cell$CTD + cell$CRD)^2 - 16 *cell$WA)^0.5) / 8

    cell$CV <- cell$CA - cell$WA

    if ( all(cell$DDOY != 0 ) ) {
      break
    }

  } ## for 1:10 end --------

  return(cell)

} ### FUNCTION END -------------------------------
