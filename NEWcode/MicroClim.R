#### 基于向量和data.table包计算微气候模块 ####
#### 计算气候 ####

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




}## Compute_clim end ----
