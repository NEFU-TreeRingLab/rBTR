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

cells_growth <- function( cell , CorV = "C" , clim.today, # layer.max,
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

  cell$EDOY[cell$CA != 0 &cell$EDOY == 0 ] <- clim.today$DOY

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

    Death[cell$DDOY != 0 ] <- 0


    ## CA 细胞面积

    dCA_dt <- vc * cell$CA * ( 1 - cell$CA/FGP$CAmax ) * ( 1- limitWTa ) * Death
    # dCA_dt[dCA_dt != 0  ] <- round(dCA_dt[dCA_dt != 0  ] , 3)

    dCA_dt[dCA_dt < 0 ] <- 0 ## error catching

    cell$TDOY[dCA_dt == 0 & cell$TDOY == 0 & cell$CA > 0 ] <- clim.today$DOY

    dCA_dt[ cell$TDOY != 0  ] <- 0 ## error catching

    ### WA 细胞壁面积

    dWA_dt <- vw * ( 1 - cell$WA/FGP$WAmax ) * ( 1 - 1/( 1 + (cell$CA-cell$WA)/FGP$mw )^FGP$sw  ) * Death
    # dWA_dt[dWA_dt != 0  ] <- round(dWA_dt[dWA_dt != 0  ] , 3)

    dWA_dt[dWA_dt < 0 ] <- 0 ## error catching
    dWA_dt[cell$WT >= FGP$WTmax ] <- 0 ## error catching
    # dWA_dt[cell$DDOY != 0 ] <- 0 ## error catching


    ### LWA 细胞壁木质化量
    dLWA_dt <- vl *  ( 1- 1/( 1 + (cell$CA - cell$WA )/FGP$ml )^FGP$sl  )*Death
    # dLWA_dt[dLWA_dt != 0  ] <- round(dLWA_dt[dLWA_dt != 0  ] , 3)

    dLWA_dt[dLWA_dt > 0 &  dLWA_dt < 0.05 ] <- 0.05 ## error catching
    dLWA_dt[dLWA_dt < 0 ] <- 0 ## error catching
    dLWA_dt[cell$LWA >= cell$WA ] <- 0 ## error catching
    # dLWA_dt[cell$DDOY != 0 ] <- 0 ## error catching


    ###  summary
    ##
    dCW <- cell$CA - cell$WA

    # li =  ( 1 - 1/( 1 + (cell$CA - cell$WA )/FGP$mw ) ^ FGP$sw )
    # jl =  ( 1 - 1/( 1 + (cell$CA - cell$WA )/FGP$ml ) ^ FGP$sl )
    ##

    cell$CA <- cell$CA +round(dCA_dt,3 )

    cell$CA[cell$CA >= CAtmax & cell$CA > 0   ] <- CAtmax[cell$CA >= CAtmax & cell$CA > 0] ## error catching


    cell$WA <- cell$WA + round(dWA_dt,3)

    cell$WA[cell$WA > FGP$WAmax & cell$CA > 0 ] <- FGP$WAmax ## error catching
    cell$WA[cell$WA > cell$CA & cell$CA > 0] <- cell$CA[cell$WA > cell$CA & cell$CA > 0] ## error catching

    cell$LWA <-  cell$LWA + round(dLWA_dt,3)

    cell$LWA[cell$LWA > cell$WA & cell$CA > 0] <- cell$WA[cell$LWA > cell$WA & cell$CA > 0]

    ## DDOY & TODY

    # cell$TDOY[dCA_dt == 0 & cell$TDOY == 0 ] = clim.today$DOY

    cell$DDOY[cell$LWA >= cell$WA & cell$DDOY == 0  & cell$CA > 0  ] <- clim.today$DOY

    cell$TDOY[ cell$WA >= FGP$WAmax & cell$TDOY == 0  & cell$CA > 0 ] <- clim.today$DOY

    cell$TDOY[ cell$WT >= FGP$WTmax & cell$TDOY == 0  & cell$CA > 0 ] <- clim.today$DOY

    # cell$DDOY[cell$LWA < cell$WA & cell$DDOY == 0 & dLWA_dt == 0 ] = clim.today$DOY+9000

    ## WT & CV & CRD 径向细胞大小

    cell$CTD[CorV == "V"] <- floor(cell$CA ^ 0.5)
    cell$CRD <- cell$CA/cell$CTD##


    cell$WT <- ( 2*(cell$CTD + cell$CRD) - (4*(cell$CTD + cell$CRD)^2 - 16 *cell$WA)^0.5) / 8

    cell$CV <- cell$CA - cell$WA

    cell[ is.na(cell) ] <- 0

    # if (nrow(cell) != 0 ) {
    #   GD = data.frame(clim.today$Year,clim.today$DOY,cell$cell_L,i, CorV, dCW,vc , vw, vl ,
    #                   cell$CA, cell$WA,cell$LWA, dCA_dt,jw, dWA_dt,jl,dLWA_dt,cell$WT)
    #   for(j in 1: nrow(GD)){
    #     cat( c( as.character(GD[j,] ),"\n") ,
    #          file = paste0(redir,"/GD.csv") ,
    #          sep = ",", fill = F,append = T)
    #   }
    # }

    if ( all(cell$DDOY[cell$CA != 0] != 0 ) ) {
      break
    }

  } ## for 1:10 end --------

  return(cell)

} ### FUNCTION END -------------------------------
