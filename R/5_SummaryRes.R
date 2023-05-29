#' normalization data
#'
#' @param dat a numeric data.
#' @param TPH a logical value, if TPH is TRUE, calculate the data as 0-100, otherwise calculate as 0-1.
#' @param ROU a logical value, if TPH is TRUE, data will be maintained to two decimal places.
#' @param min1 a logical value,if TPH is TRUE, data minimum value is not zero.
#' @param Zeros a logical value,if TPH is TRUE, raw data minimum value set to 0.
#'
#' @return normalization data.
#'

nor <- function( dat,TPH = F ,ROU = F,min1 = F, Zeros = F ){
  a = 0
  if (min1 == T) {
    a = 0.01
  }

  dat <- as.numeric(dat)

  mind <- min(dat,na.rm =T)
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
