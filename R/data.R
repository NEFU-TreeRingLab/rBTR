##' Cell anatomical characteristics for \emph{Betula platyphylla} in Liangshui
##'
##'
##' This dataset gives the fiber cell lumen area, fiber cell wall thickness
##' and vessel cell lumen area for \emph{Betula platyphylla}
##' at each relative ring width in each year.
##'
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##'
##'
##' @docType data
##' @keywords datasets
##' @name BP_cells
##' @usage data(BP_cells)
##' @format A \code{data.frame} containing five columns with year(Year),
##' relative ring width(RRingWidth), lumen area(LA) ,cell wall thickness(CWT)
##' and cell type(CellType).
##'
'BP_cells'

##' Cell anatomical characteristics for \emph{Fraxinus mandshurica} in Liangshui
##'
##'
##' This dataset gives the fiber cell lumen area, fiber cell wall thickness
##' and vessel cell lumen area for \emph{Fraxinus mandshurica}
##' at each relative ring width in each year.
##'
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##'
##'
##' @docType data
##' @keywords datasets
##' @name FM_cells
##' @usage data(FM_cells)
##' @format A \code{data.frame} containing five columns with year(Year),
##' relative ring width(RRingWidth), lumen area(LA) ,cell wall thickness(CWT)
##' and cell type(CellType).
##'
'FM_cells'

##' Tree ring width for \emph{Betula platyphylla} in Liangshui
##'
##'
##' This dataset gives the mean tree ring width for \emph{Betula platyphylla}
##' in each year.
##'
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##'
##'
##' @docType data
##' @keywords datasets
##' @name BP_ring
##' @usage data(BP_ring)
##' @format A \code{data.frame} containing four columns with year(Year),
##' mean tree ring width(RingWidth), Standard Deviation(SD) and sample depth(N).
##'
##'
'BP_ring'

##' Tree ring width for \emph{Fraxinus mandshurica} in Liangshui
##'
##'
##' This dataset gives the mean tree ring width for \emph{Fraxinus mandshurica}
##' in each year.
##'
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##'
##'
##' @docType data
##' @keywords datasets
##' @name FM_ring
##' @usage data(FM_ring)
##' @format A \code{data.frame} containing four columns with year(Year),
##' mean tree ring width(RingWidth), Standard Deviation(SD) and sample depth(N).
'FM_ring'


##' BTR model parameters for \emph{Betula platyphylla}
##'
##'
##' This dataset gives the BTR model parameters list for simulate
##' \emph{Betula platyphylla} xylogenesis in Liangshui.
##'
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##'
##'
##' @docType data
##' @keywords datasets
##' @name BP_param
##' @usage data(BP_param)
##' @format A \code{data.frame} containing 7 columns :
##' @format **in_formals** : is the parameter name in Manuscript,
##' @format **parameter** : is the parameter name in Rcode,
##' @format **description** : is the description of parameter,
##' @format **note** : is the object for which the parameter is used,
##' @format **paramtype** : is the parameter type,
##' @format **modul** : is the module running the parameter,
##' @format **values** : is the parameter value.
##'
##'
##'
'BP_param'


##' BTR model parameters for \emph{Fraxinus mandshurica}
##'
##'
##' This dataset gives the BTR model parameters list for simulate
##' \emph{Fraxinus mandshurica} xylogenesis in Liangshui.
##'
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##'
##'
##' @docType data
##' @keywords datasets
##' @name FM_param
##' @usage data(FM_param)
##' @format A \code{data.frame} containing 7 columns :
##' @format **in_formals** : is the parameter name in Manuscript,
##' @format **parameter** : is the parameter name in Rcode,
##' @format **description** : is the description of parameter,
##' @format **note** : is the object for which the parameter is used,
##' @format **paramtype** : is the parameter type,
##' @format **modul** : is the module running the parameter,
##' @format **values** : is the parameter value.

##'
'FM_param'

##' BTR model parameters for microclim module
##'
##'
##' This dataset gives the BTR model parameters list for simulate microclim in Liangshui.
##' Original code from "Run CPC Leaky Bucket model (submonthly version)"
##' We have modified some of the parameters to allow the model to estimate data on a daily basis.
##'
##' @source Modifications of Suz Tolwinski-Ward's monthly time-step code by Nick Graham in 2011.
#' Ported to R by SETW in 2015. Implementation of CPC Leaky Bucket model as described in
#' Huang et al., 'Analysis of Model-Calculated Soil Moisture over the United States
#' (1931-1993) and Applications to Long-Range Temperature Forecasts,' J. Clim. (1995)
##'
##'
##' @docType data
##' @keywords datasets
##' @name Clim_param
##' @usage data(Clim_param)
##' @format A \code{data.frame} containing 6 columns :
##' @format **parameter** : is the parameter name in Rcode,
##' @format **description** : is the description of parameter,
##' @format **note** : is the note of parameter,
##' @format **paramtype** : is the parameter type,
##' @format **modul** : is the module running the parameter,
##' @format **values** : is the parameter value.

##'
'Clim_param'


##' Daily Mean Temperature, Total Precipitation and Relative Humidity for
##' Liangshui,China
##'
##' This dataset gives the daily mean temperature, total
##' precipitation and relative humidity at Liangshui, Yichun, Heilongjiang, China.
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##' @docType data
##' @keywords datasets
##' @name LS_clim
##' @usage data(LS_clim)
##' @format A \code{data.frame} containing 7 columns with year (Year),
##' month (Month),day (Day), Day of year (DOY), precipitation (PRE),
##' temperature (TEM) and relative humidity (RH).
##'
'LS_clim'

##' Daily Mean Temperature, soil moisture and VPD for
##' Liangshui,China
##'
##' This dataset gives the daily mean temperature, soil Moisture and
##' VPD at Liangshui, Yichun, Heilongjiang, China.
##' @source A process-based model of xylogenesis and tree-ring formation of
##' broad-leaved trees (BTR) driven by climate. (not published yet)
##' @docType data
##' @keywords datasets
##' @name LS_climdata
##' @usage data(LS_climdata)
##' @format A \code{data.frame} containing 6 columns with year (Year),
##' month (Month),day (Day), soil moisture (soilM),
##' temperature (TEM) and VPD (VPD).
##'
'LS_climdata'
