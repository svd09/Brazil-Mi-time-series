#-------------------------------#
# MI MORTALITY AND TEMPERATURE  #
#-------------------------------#
# PAPER SCRIPT TO MODEL MI MORTALITY AND TEMPERATURE IN BRAZIL AND REGIONS.
# SUBCOHORT ANALYSIS BY SEX AND AGE GROUPS.

# load libraries.
rm(list = ls())
options(scipen = 999)

library(easypackages)
libraries(c(
  "tidyverse", "tidylog", "vroom", "lubridate",
  "stringr", "purrr", "lubridate", "tidyr","sf","geobr",
  "readxl","ggrepel","zoo","survival",
  "dlnm","splines","metafor","gnm","broom",
  "rlist","mixmeta","janitor","qs","geobr","data.table","pbapply",
  "MASS","Matrix","weathermetrics","dtplyr"))
# get the data.
path = "/mnt/vstor/SOM_CVRI_JXD101/DATA_SCIENCE/brazil mortality data/"
dlist <-  qread(paste0(path,"comb_data_subcohorts.qs"))
# see the col for elements.
names(dlist) # no names.
# so from the prior script.
names(dlist) <- c("males","females","young","old")
colnames(dlist[[1]])
dlist$males$data_name = "males"
dlist$females$data_name = "females"
dlist$young$data_name = "young"
dlist$old$data_name = "old"

# convert NA in counts to 0.
for(k in 1:4){
  df = dlist[[k]]
  df$daily_counts[is.na(df$daily_counts)] = 0
  dlist[[k]] = df
}
# look at 1 dataset.
str(dlist$males)
dlist[[1]]
# make the app_temp variable.
# - HEAT.
# make heatwave variables & stratum variable.
# best to all join them together and then run across the whole dataset.
data = bind_rows(dlist)
data$umidade[data$umidade > 100] = 100
data <-  data %>% dtplyr::lazy_dt() %>%
  mutate(
    year = year(date),
    month = month(date),
    dow = weekdays(date),
    stratum = factor(year):factor(month):factor(dow):factor(loc6digit),
    tmean = weathermetrics::heat.index(
      t = temperatura, 
      rh = umidade,
      temperature.metric = "celcius")) %>%
  as_tibble()

# make categorical variables.
data <- data  %>% mutate(
  EHE_90 = ifelse(tmean >= quantile(tmean,0.90) ,1,0),
  EHE_95 = ifelse(tmean >= quantile(tmean,0.95) ,1,0),
  EHE_99 = ifelse(tmean >= quantile(tmean,0.99,na.rm = TRUE),1,0),
  EHE_10 = ifelse(tmean <= quantile(tmean,0.10) ,1,0),
  EHE_5 = ifelse(tmean <= quantile(tmean,0.05) ,1,0),
  EHE_01 = ifelse(tmean <= quantile(tmean,0.01,na.rm = TRUE),1,0),
  .by=loc6digit
) 

data <- data %>%  mutate(
  EHE_90_2 = case_when(
    EHE_90 == 1 & lead(EHE_90, 1, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  EHE_90_3 = case_when(
    EHE_90 == 1 & 
      lead(EHE_90, 1, default = 0) == 1 & 
      lead(EHE_90, 2, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  EHE_99_2 = case_when(
    EHE_99 == 1 & lead(EHE_99, 1, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  EHE_99_3 = case_when(
    EHE_99 == 1 & 
      lead(EHE_99, 1, default = 0) == 1 & 
      lead(EHE_99, 2, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  EHE_10_2 = case_when(
    EHE_10 == 1 & lead(EHE_10, 1, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  EHE_10_3 = case_when(
    EHE_10 == 1 & 
      lead(EHE_10, 1, default = 0) == 1 & 
      lead(EHE_10, 2, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  EHE_01_2 = case_when(
    EHE_01 == 1 & lead(EHE_01, 1, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  EHE_01_3 = case_when(
    EHE_01 == 1 & 
      lead(EHE_01, 1, default = 0) == 1 & 
      lead(EHE_01, 2, default = 0) == 1 ~ 1,
    TRUE ~ 0
  ),
  .by=loc6digit
) 

# fit categorical models for each group.
# keep only data in the env to clear memory.
rm(list = setdiff(ls(), "data"))
saveRDS(data, 
        paste0(path,"subcohort_data.rds"))
# this data saved for further models.

data_m = data %>% filter(data_name == 'males')
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")
cb_list_m <- vector("list",length(heat_vars))

males_res = data.frame(
  event = heat_vars,
  est = rep(NA,length(heat_vars)),
  high = rep(NA,length(heat_vars)),
  low = rep(NA,length(heat_vars))
)


for(j in 1:length(heat_vars)){
  cat(j,"")
  pb$tick()
  var = heat_vars[j]
  cb <- crossbasis(data_m[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_m$stratum),
               data = data_m)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  males_res[j,2] <- or
  males_res[j,3] <- high
  males_res[j,4] <- low
}
males_res
write_csv(males_res,
          "results/males_or.csv")

#---------#
# FEMALES #
#---------#

data_f = data %>% filter(data_name == 'females')
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")

females_res = data.frame(
  event = heat_vars,
  est = rep(NA,length(heat_vars)),
  high = rep(NA,length(heat_vars)),
  low = rep(NA,length(heat_vars))
)

for(j in 1:length(heat_vars)){
  cat(j,"")
  
  var = heat_vars[j]
  cb <- crossbasis(data_f[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_f$stratum),
               data = data_f)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  females_res[j,2] <- or
  females_res[j,3] <- high
  females_res[j,4] <- low
}

females_res

write_csv(females_res,
          "/home/svd14/brazil mi mortality/results/females_or.csv")


#-------#
# YOUNG #
#-------#

data_y = data %>% filter(data_name == 'young')
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")

young_res = data.frame(
  event = heat_vars,
  est = rep(NA,length(heat_vars)),
  high = rep(NA,length(heat_vars)),
  low = rep(NA,length(heat_vars))
)

for(j in 1:length(heat_vars)){
  cat(j,"")
  
  var = heat_vars[j]
  cb <- crossbasis(data_y[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_y$stratum),
               data = data_y)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  young_res[j,2] <- or
  young_res[j,3] <- high
  young_res[j,4] <- low
}
young_res
write_csv(young_res,
          "/home/svd14/brazil mi mortality/results/young_or.csv")


#-----#
# OLD #
#-----#

data_o = data %>% filter(data_name == 'old')
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")

old_res = data.frame(
  event = heat_vars,
  est = rep(NA,length(heat_vars)),
  high = rep(NA,length(heat_vars)),
  low = rep(NA,length(heat_vars))
)

for(j in 1:length(heat_vars)){
  cat(j,"")
  
  var = heat_vars[j]
  cb <- crossbasis(data_o[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_o$stratum),
               data = data_o)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  old_res[j,2] <- or
  old_res[j,3] <- high
  old_res[j,4] <- low
}
old_res
write_csv(old_res,
          "/home/svd14/brazil mi mortality/results/old_or.csv")


#- END =============================================================
# THE SCRIPT FOR HEAT VARS IS OVER.
# NOW TO DO THE SAME FOR THE COLD VARS.

data_m = data %>% filter(data_name == "males")

cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
"EHE_01","EHE_01_2","EHE_01_3")

males_res_cold = data.frame(
  event = cold_vars,
  est = rep(NA,length(cold_vars)),
  high = rep(NA,length(cold_vars)),
  low = rep(NA,length(cold_vars))
)

for(j in 1:length(cold_vars)){
  cat(j,"")
 
  var = cold_vars[j]
  cb <- crossbasis(data_m[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_m$stratum),
               data = data_m)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  males_res_cold[j,2] <- or
  males_res_cold[j,3] <- high
  males_res_cold[j,4] <- low
}
males_res_cold
write_csv(males_res_cold,
          "/home/svd14/brazil mi mortality/results/males_res_cold.csv")

#---------#
# FEMALES #
#---------#

data_f = data %>% filter(data_name == "females")

cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
              "EHE_01","EHE_01_2","EHE_01_3")

females_res_cold = data.frame(
  event = cold_vars,
  est = rep(NA,length(cold_vars)),
  high = rep(NA,length(cold_vars)),
  low = rep(NA,length(cold_vars))
)

for(j in 1:length(cold_vars)){
  cat(j,"")
  
  var = cold_vars[j]
  cb <- crossbasis(data_f[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_f$stratum),
               data = data_f)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  females_res_cold[j,2] <- or
  females_res_cold[j,3] <- high
  females_res_cold[j,4] <- low
}

females_res_cold

write_csv(females_res_cold,
  "/home/svd14/brazil mi mortality/results/females_res_cold.csv")
#-------#
# YOUNG #
#-------#
data_y = data %>% filter(data_name == "young")
cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
              "EHE_01","EHE_01_2","EHE_01_3")
young_res_cold = data.frame(
  event = cold_vars,
  est = rep(NA,length(cold_vars)),
  high = rep(NA,length(cold_vars)),
  low = rep(NA,length(cold_vars))
)

for(j in 1:length(cold_vars)){
  cat(j,"")
  
  var = cold_vars[j]
  cb <- crossbasis(data_y[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_y$stratum),
               data = data_y)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  young_res_cold[j,2] <- or
  young_res_cold[j,3] <- high
  young_res_cold[j,4] <- low
}

young_res_cold
write_csv(young_res_cold,
  "/home/svd14/brazil mi mortality/results/young_res_cold.csv")

#-----#
# OLD #
#-----#
data_o = data %>% filter(data_name == "old")
cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
              "EHE_01","EHE_01_2","EHE_01_3")
old_res_cold = data.frame(
  event = cold_vars,
  est = rep(NA,length(cold_vars)),
  high = rep(NA,length(cold_vars)),
  low = rep(NA,length(cold_vars))
)
for(j in 1:length(cold_vars)){
  cat(j,"")
  
  var = cold_vars[j]
  cb <- crossbasis(data_o[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_o$stratum),
               data = data_o)
  cp <- crosspred(cb,model,at=1,cen=0)
  or <- cp$allRRfit
  low <- cp$allRRlow
  high <- cp$allRRhigh
  est = paste0(or,"(",low,"-",high,")")
  old_res_cold[j,2] <- or
  old_res_cold[j,3] <- high
  old_res_cold[j,4] <- low
}
old_res_cold
write_csv(old_res_cold,
  "/home/svd14/brazil mi mortality/results/old_res_cold.csv")
# END ========================================================
# ALL MODELS FOR SUBCOHORTS DONE.
#########################===##################################.







#- junk code below this.


cp = crosspred(cb,model,at=1,cen=0)



names(model)[1:4] = colnames(cb)
coef = coef(model)
vcov = vcov(model)
length(grep("^cbtemp_", names(coef(model))))


















names(cb_list_m) = heat_vars
for(j in 1:length(heat_vars)){
  cat(j,"")
  var <- heat_vars[j]
  cb <- crossbasis(data_m[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  colnames(cb) = paste0("cbtemp_", seq_len(ncol(cb_temp)))
  cb_list_m[[j]] <- cb
  rm(cb)
}

model_list_m <- vector("list", length(heat_vars))
names(model_list_m) <- heat_vars

for(n in 1:length(model_list_m)){
  cat(n,"")
  cb <- cb_list_m[[n]]
  formula <- as.formula(paste("daily_counts ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(data_m$stratum),
               data = data_m)
  model_list_m[[n]] <- model
  rm(model)
}
model_list_m[[1]]

crosspred(cb_list_m[[1]],model_list_m[[1]],at=1,cen=0)

gc() # clean memory that is not needed now.
# do crosspred to get the cumulative RR.
crosspred_list_m = vector("list",length(heat_vars))
names(crosspred_list_m) = heat_vars
for(j in 1:length(heat_vars)){
  cat(j,"")
  cb = cb_list_m[[j]]
  model = model_list_m[[j]]
  cp = crosspred(cb,model,at=1,cen=0,cumul=T)
  crosspred_list_m[[j]] <- cp
}
crosspred_list_m
#- snippet for progress bar.
library(progress)

n <- 100
pb <- progress_bar$new(
  format = "  doing [:bar] :percent in :elapsed",
  total = n, clear = FALSE, width = 60
)

for (i in 1:n) {
  Sys.sleep(0.05)
  pb$tick()
}

#- end snippet.

# extract the OR from the crosspred for the cumulative effect.
males_or 