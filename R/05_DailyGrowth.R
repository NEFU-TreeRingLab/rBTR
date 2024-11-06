#############################
#### daily Groeth modual ####
#############################

## interFun of Daily growth
#' interFun of Daily growth
#' @param newCells 形成层活动结果
#' @param newVessel 形成层活动结果
#' @param vesselsNum 形成层活动结果
#' @param dailyCells 细胞们
#' @param cells cells anatomy
#' @param vessels vessels anatomy
#' @param clim.today Today climate data
#' @param fixparam.growth.fiber Fix parameters of fiber cells
#' @param fixparam.growth.vessel Fix parameters
#' @param dynparam.growth.t dynparam.growth
#'
#'
#' @return today cells anatomy
#'
## #' @importFrom . cells_growth
#'
# #' @export





daily_grwoth <- function( newCell, newVessel, vesselsNum,
                          dailyCells , cells, vessels,
                          clim.today,
                          fixparam.growth.fiber, fixparam.growth.vessel,
                          dynparam.growth.t){

  ## 加入产生的细胞
  ifelse( all(newCell != 0 ),
          dailyCells$dailyFiber[ dailyCells$dailyFiber$cell_L %in% as.numeric( newCell ),c(3:12) ] <- cells[,c(3:12)], NA )
  ifelse( all(newVessel != 0 ), {
    dailyCells$dailyVessels[ dailyCells$dailyVessels$cell_L %in% as.numeric( newVessel ), c(3:14)] <- vessels[,c(3:14)]
    dailyCells$dailyVessels$VN[ dailyCells$dailyVessels$cell_L %in% as.numeric( newVessel ) ] <- vesselsNum } , NA )

  ## cell grwoth
  ifelse( any(dailyCells$dailyFiber$CA != 0 ) , dailyCells$dailyFiber <-
            rBTR:::cells_growth(cell = dailyCells$dailyFiber , CorV = "C" , clim.today,
                                  fixparam.growth.fiber, fixparam.growth.vessel,dynparam.growth.t  )  , NA   )
  ## fiber growth
  ifelse( any(dailyCells$dailyVessels$CA != 0 ) , dailyCells$dailyVessels <-
            rBTR:::cells_growth(cell = dailyCells$dailyVessels , CorV = "V" , clim.today,
                                  fixparam.growth.fiber, fixparam.growth.vessel,dynparam.growth.t  )   , NA)

  return(dailyCells)

}

