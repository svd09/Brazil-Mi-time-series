#---------------------------------------#
# ETE AND ECE MODELS FOR THE WHOLE DATA #
#---------------------------------------#
# THIS IS SCRIPT FOR THE PAPER.
# FIT ETE AND ECE MODELS FOR THE WHOLE DATA.
# GET THE OR AND REPORT THEM IN TABLE.
# load libraries.
rm(list = ls()) # clean env
options(scipen = 999) # no scientific notation

library(easypackages)
libraries(c(
  "tidyverse", "tidylog", "vroom", "lubridate",
  "stringr", "purrr", "lubridate", "tidyr","sf","geobr",
  "readxl","ggrepel","weathermetrics","zoo","survival",
  "dlnm","splines","mvmeta","metafor","gnm","broom","ggh4x",
  "rlist","mixmeta",'ckbplotr',"janitor",
  "bannerCommenter","qs","geobr","data.table"))


# get the data
data <- qread("data/whole_data_mi_mortality_cpr_model.qs")
glimpse(data)
length(unique(data$loc6digit))
# 2443 locations.
# all data formatted for models.
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")

brazil_res = data.frame(
  event = heat_vars,
  est = rep(NA,length(heat_vars)),
  high = rep(NA,length(heat_vars)),
  low = rep(NA,length(heat_vars))
)

for(j in 1:length(heat_vars)){
  cat(j,"")
  var = heat_vars[j]
  cb <- crossbasis(data[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data$stratum),
               data = data)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  brazil_res[j,2] <- or
  brazil_res[j,3] <- high
  brazil_res[j,4] <- low
  # remove big stuff
  rm(model,cb,cp)
}

brazil_res
write_csv(brazil_res,"results/brazil_heat_res.csv")
# ==============================================================================
#------#
# COLD #
#------#

cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
              "EHE_01","EHE_01_2","EHE_01_3")
brazil_cold = data.frame(
  event = cold_vars,
  est = rep(NA,length(cold_vars)),
  high = rep(NA,length(cold_vars)),
  low = rep(NA,length(cold_vars))
)

for(j in 1:length(cold_vars)){
  cat(j,"")
  
  var = cold_vars[j]
  cb <- crossbasis(data[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data$stratum),
               data = data)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  brazil_cold[j,2] <- or
  brazil_cold[j,3] <- high
  brazil_cold[j,4] <- low
}
brazil_cold
write_csv(brazil_cold,"results/brazil_cold_res.csv")

#- END #========================================================================