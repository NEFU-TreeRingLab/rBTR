##' BTR model parameters for \emph{Betula platyphylla}
##'
##'
##' This dataset gives the BTR model parameters list for simulate
##' \emph{Betula platyphylla} xylogenesis in Liangshui.
##'
##' @source A process-based model of climate-driven xylogenesis and tree-ring formation in broad-leaved trees (BTR).[DOI:10.1093/treephys/tpae127]
##'
##'
##' @docType data
##' @keywords datasets
##' @name BPparam
##' @usage data(BPparam)
##' @format A \code{data.frame} containing 7 columns :
##' @format **Parameter** : is the parameter name in Rcode,
##' @format **ParamType** : is the parameter type,
##' @format **Module** : is the module running the parameter,
##' @format **Note** : is the object for which the parameter is used,
##' @format **Description** : is the description of parameter,##'
##' @format **Values** : is the parameter value.
##' @format **Unit** : is the unit parameter value.
##'
'BPparam'

##' BTR model age input for \emph{Betula platyphylla}
##'
##'
##' This dataset gives the BTR model age trend input for simulate
##' \emph{Betula platyphylla} xylogenesis in Liangshui.
##'
##' @source A process-based model of climate-driven xylogenesis and tree-ring formation in broad-leaved trees (BTR).[DOI:10.1093/treephys/tpae127]
##'
##' @docType data
##' @keywords datasets
##' @name BPage
##' @usage data(BPage)
##' @format A \code{data.frame} containing 4 columns :
##' @format **Year** : is the year of growth,
##' @format **age** : is the age of trees corresponding to the year,
##' @format **Tage** : is the trend of cambial activaty,
##' @format **Lage** : is the trend of cell enlargement.
##'
'BPage'

##' Meteorological data, input data for the BTR model from Liangshui,China
##'
##' Model primary inputs include date, temperature, VPD soil moisture at Liangshui, Yichun, Heilongjiang, China.
##' @source A process-based model of climate-driven xylogenesis and tree-ring formation in broad-leaved trees (BTR).[DOI:10.1093/treephys/tpae127]
##'
##' @docType data
##' @keywords datasets
##' @name Climdata
##' @usage data(Climdata)
##' @format A \code{data.frame} containing 6 columns with year (Year),
##' @format **Year** : Year,
##' @format **Month** : Month,
##' @format **Day** : Day,
##' @format **TEM** : Temperature,
##' @format **soilM** : Soil moisture,
##' @format **VPD** : vapor Pressure Deficit.
##'
'Climdata'

##' Meteorological data 2 , easily accessible from meteorological stations , input data for the BTR model from Liangshui,China
##'
##' Model primary inputs include date, temperature, VPD soil moisture at Liangshui, Yichun, Heilongjiang, China.
##' @source A process-based model of climate-driven xylogenesis and tree-ring formation in broad-leaved trees (BTR).[DOI:10.1093/treephys/tpae127]
##'
##' @docType data
##' @keywords datasets
##' @name Climdata.2
##' @usage data(Climdata.2)
##' @format A \code{data.frame} containing 6 columns with year (Year),
##' @format **Year** : Year,
##' @format **Month** : Month,
##' @format **Day** : Day,
##' @format **TEM** : Temperature,
##' @format **RH** : Relative humidity,
##' @format **PRE** : Precipitation.
##'
'Climdata.2'
