
#' cell_division function
#'
#' @param clim.today Today climate data
#' @param fixparam.divi Fix parameters
#' @param fixparam.growth.fiber Fix parameters of fiber cells
#' @param fixparam.growth.vessel Fix parameters
#' @param dynparam.growth.t Dynamic parameters
#' @param cells cells anatomy
#' @param vessels vessels anatomy
#'
#' @return today cells layer number & vessels number
#'
#' @importFrom dplyr bind_rows
#'

cell_division <- function( clim.today,
                           fixparam.divi,fixparam.growth.fiber,fixparam.growth.vessel,
                           dynparam.growth.t, cells, vessels, deltaD_T){   ## Fixp_cambi,, TA

  egR <- clim.today$gE* clim.today$gT* min(  clim.today$gM, clim.today$gV)  ## Vessels growth rate
  wgR <- clim.today$gE* clim.today$gT    ##       ## fiber  growth rate

  ## error catch
  dynparam.growth.t$grwothSeason[egR <= 0.05 & clim.today$DOY > 200 ] <- 1 ##
  egR[dynparam.growth.t$grwothSeason == 1 ] <- 0.001 ## 生长季快结束时细胞不扩大

  ## if age influence wgR
  ## growth velocity
  dynparam.growth.t$v_c.fiber  <- fixparam.divi$va_c.fiber  * egR # ( wgR^(1 + dynparam.growth.t$L_just) )

  dynparam.growth.t$v_w.fiber  <- fixparam.divi$va_w.fiber  * wgR * (1+clim.today$L_i.fiber) # ( wgR^(1 - dynparam.growth.t$L_just) )

  dynparam.growth.t$v_l.fiber  <- fixparam.divi$va_l.fiber * wgR * (1+clim.today$L_i.fiber) # ( wgR^(1 - dynparam.growth.t$L_just) )


  dynparam.growth.t$v_c.vessel <- fixparam.divi$va_c.vessel * egR # ( wgR^(1 + dynparam.growth.t$L_just2) )

  dynparam.growth.t$v_w.vessel <- fixparam.divi$va_w.vessel * wgR * (1+clim.today$L_i.vessel) # ( wgR^(1 - dynparam.growth.t$L_just2) )

  dynparam.growth.t$v_l.vessel <- fixparam.divi$va_l.vessel * wgR * (1+clim.today$L_i.vessel) # ( wgR^(1 - dynparam.growth.t$L_just2) )




  ### cambium grwoth from a bind model （）

  if ( is.null(clim.today$Tage)    ) {
    dynparam.growth.t$T_age <- fixparam.divi$alpha_age *
      exp(fixparam.divi$beta_age * dynparam.growth.t$Age)
  } else {
    dynparam.growth.t$T_age <- clim.today$Tage
  }

  dynparam.growth.t$czgR <- fixparam.divi$alpha_cz *
    exp(fixparam.divi$beta_cz * ( egR )  ) ## new

  dynparam.growth.t$czgR[ dynparam.growth.t$czgR > egR ] <- egR ## error catch

  # CCMAX <- ( Fixp_divi$DCC  +  Dynp_growth$GR_age   )   *Dynp_growth$GR_cz   ##  * GRT

  ifelse( CZgR[1] == 1, cambial.cell.gR <- egR,  cambial.cell.gR <- 1)

  ## not good , Var in old tree is too small
  # v_cz <- ( fixparam.divi$va_cz * cambial.cell.gR + dynparam.growth.t$T_age )* dynparam.growth.t$czgR    ## new
  v_cz <- ( fixparam.divi$va_cz * dynparam.growth.t$T_age )* dynparam.growth.t$czgR    ## 3.0ver

  dynparam.growth.t$egR <- egR ## check 以后修改

  v_cz[dynparam.growth.t$grwothSeason == 1 ] <- 0 ## error catch 非生长季不生长

  dynparam.growth.t$dCA_cz[dynparam.growth.t$dCA_cz == 0 ] <- cells$CA

  dynparam.growth.t$Vcz <- v_cz

  Layers <- 0
  # CA_cz = cells$CA + dynparam.growth.t$dCA_cz
  CA_cz <- dynparam.growth.t$dCA_cz

  for (i in 1:10) {

    CA_cz <- (1+ v_cz) * CA_cz

    if ( CA_cz >= (2*cells$CA) ) {
      Layers <- Layers + 1
      CA_cz <- cells$CA
    }
  }

  dynparam.growth.t$dCA_cz <- CA_cz

  Ct <- cells[-1,] # 生成空白矩阵
  Vt <- vessels[-1,]
  Nv <- 0

  if (Layers > 0 ) { ## 判读是否分裂

    for (i in 1 : Layers ) { ## 分裂便计算分化
      ## 1，细胞层数+1
      dynparam.growth.t$SumCL <- dynparam.growth.t$SumCL +1
      cells$Year <- clim.today$Year
      cells$cell_L <- dynparam.growth.t$SumCL
      cells$EDOY <- clim.today$DOY

      Ct <- dplyr::bind_rows(Ct , cells)
      ## 1 end


      if (deltaD_T != 0) { #### 分化
        ## 2, 计算分化
        # 分裂间期（层数计）+1

        dynparam.growth.t$deltaVN <- dynparam.growth.t$deltaVN + 1

        # 计算该间期下分裂的层数
        Nv <- dynparam.growth.t$deltaVN %/% deltaD_T

        # 计算新间期
        dynparam.growth.t$deltaVN <- dynparam.growth.t$deltaVN %% deltaD_T

        # 每年首次分裂产生至少1个导管

        Nv[dynparam.growth.t$SumV == 0 & dynparam.growth.t$v_c.vessel != 0 & Nv == 0] <- 1

        dynparam.growth.t$SumV <- dynparam.growth.t$SumV + Nv ## 总导管数增加

        if (Nv != 0) { ## 判断分化导管数量

          vessels$Year <- clim.today$Year
          vessels$cell_L <- cells$cell_L
          vessels$EDOY <- clim.today$DOY
          vessels$VN <- Nv
          # dynparam.growth.t$SumV <- dynparam.growth.t$SumV + Nv
          vessels$NoV <- dynparam.growth.t$SumV
          Vt <- dplyr::bind_rows(Vt , vessels)

        }else{
          Vt[i,] <- NA
        }
        ## 导管生长判断结束

      }#### end if 分化



    } ## 每层生长结束 end for i -------
    ###

    Today_G <- list(Ct,Vt,dynparam.growth.t)

  }else{
    Ct[1,] <- NA
    Vt[1,] <- NA
    Today_G <- list(Ct,Vt,dynparam.growth.t)
  } ## 层生长结束


  # detach(fixparam.divi)

  return( Today_G )

} ## cell_division end ----------------
