#### try try
library(devtools)
library(data.table)
load('data/Climdata.rda')
load('data/BPparam.rda')
load('data/BPage.rda')

source('R/1_MicroClimNeo.R')
source('R/5_SummaryRes.R')
source('R/compute_gR.R')

clim <- Compute_clim( climdata = Climdata,lat = 47.13)

climdata <- Climdata

parameters <- BPparam ; age <- BPage; start.year = 1998; end.year = 2000
write.result = F; save.daily = F; gR.method = "Jonhson"
# ; mode = 'mix'Dcase = "min";
division = 'dynamic'; Named = NULL; set.layers = 300

##
document()


usethis::use_r('btr')
usethis::use_r('microClimata')
