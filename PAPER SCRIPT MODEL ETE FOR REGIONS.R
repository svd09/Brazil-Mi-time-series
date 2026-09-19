#-------------------------#
# regional models for ETE #
#-------------------------#
# PAPER SCRIPT
# MODEL ETE FOR EACH REGION AND SAVE OBJ.
# USE THEM TO GET AF/AN LATER.

#------#
# COLD #
#------#


rm(list = ls()) # clean env
options(scipen = 999) # no scientific notation

library(easypackages)
libraries(c(
  "tidyverse", "tidylog", "vroom", "lubridate",
  "stringr", "purrr", "lubridate", "tidyr","sf","geobr",
  "readxl","ggrepel","weathermetrics","zoo","survival",
  "dlnm","splines","mvmeta","metafor","gnm","broom","ggh4x",
  "rlist","mixmeta",'ckbplotr',"janitor",
  "bannerCommenter","qs","geobr","data.table","pbapply",
  "MASS","Matrix"))

# North
north = qread("data/north_mi_mortality_with_exp_overall.qs")
# we need to get the ETE variables separately for each region.
# ETE will change for each region based on their tmean.
glimpse(north)
north = north %>% dplyr::select(-starts_with('EHE'))
glimpse(north)

# make a list of the data.
dlist <- split(north, north$loc6digit)
names(dlist) <- sapply(dlist, function(x) unique(x$state_loc.x))

# - HEAT.
# make heatwave variables.
for(k in 1:length(dlist)){
  df = dlist[[k]]
  df2 = df %>% mutate(
    EHE_90 = ifelse(tmean >= quantile(tmean,0.90) ,1,0),
    EHE_95 = ifelse(tmean >= quantile(tmean,0.95) ,1,0),
    EHE_99 = ifelse(tmean >= quantile(tmean,0.99,na.rm = TRUE),1,0))
  
  dlist[[k]] <- df2
}

for(k in 1:length(dlist)){
  df = dlist[[k]]
  df2 = df %>% mutate(
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
    )
  )
  
  dlist[[k]] <- df2
}

#- COLD.
# for cold variables choose 1%, 5%, and 10%.

for(k in 1:length(dlist)){
  df = dlist[[k]]
  df2 = df %>% mutate(
    EHE_10 = ifelse(tmean <= quantile(tmean,0.10) ,1,0),
    EHE_5 = ifelse(tmean <= quantile(tmean,0.05) ,1,0),
    EHE_01 = ifelse(tmean <= quantile(tmean,0.01,na.rm = TRUE),1,0))
  
  dlist[[k]] <- df2
}

for(k in 1:length(dlist)){
  df = dlist[[k]]
  df2 = df %>% mutate(
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
    )
  )
  
  dlist[[k]] <- df2
}

north <- bind_rows(dlist)
glimpse(north)
# fit models.
# MODELS for COLD.
cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
               "EHE_01","EHE_01_2","EHE_01_3")


cb_list <- vector("list", length(cold_vars))
names(cb_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  var <- cold_vars[n]
  cb <- crossbasis(north[[var]], lag = 7,
        argvar = list(fun = "lin"),
        arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(cold_vars))
names(model_list) <- cold_vars

for(n in 1:length(model_list)){
  cat(n,"")
  cb <- cb_list[[n]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(north$stratum),
               data = north)
  model_list[[n]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_list = vector("list",length(cold_vars))
names(cp_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  cb <- cb_list[[n]]
  model <- model_list[[n]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_list[[n]] <- cp
  rm(cp)
}
north_cold_res = data.frame(
  var = cold_vars,
  RR = rep(NA, length(cold_vars)),
  L  = rep(NA, length(cold_vars)),
  U  = rep(NA, length(cold_vars))
)
# get the cumulative RR from each model.

for(j in 1:6){
  cp = cp_list[[j]]
  north_cold_res[j,2] = cp$allRRfit
  north_cold_res[j,3] = cp$allRRlow
  north_cold_res[j,4] = cp$allRRhigh
}
# save this 
saveRDS(north_cold_res,
        "model objects/north_cold_results.rds")

write_csv(north_cold_res,
          "results/north_cold_res.csv")
# remove ENV.
rm(list = ls())
#=========================================================
# get NE data.
ne = qread("data/ne_mi_mortality_with_exp_overall.qs")
ne = ne %>% dplyr::select(-starts_with('EHE'))

ne <- ne  %>% mutate(
  EHE_90 = ifelse(tmean >= quantile(tmean,0.90) ,1,0),
  EHE_95 = ifelse(tmean >= quantile(tmean,0.95) ,1,0),
  EHE_99 = ifelse(tmean >= quantile(tmean,0.99,na.rm = TRUE),1,0),
  EHE_10 = ifelse(tmean <= quantile(tmean,0.10) ,1,0),
  EHE_5 = ifelse(tmean <= quantile(tmean,0.05) ,1,0),
  EHE_01 = ifelse(tmean <= quantile(tmean,0.01,na.rm = TRUE),1,0),
  .by=loc6digit
) 

ne <- ne %>%  mutate(
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



glimpse(ne)
# fit models.
# MODELS for COLD.
cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
              "EHE_5","EHE_01","EHE_01_2","EHE_01_3")


cb_list <- vector("list", length(cold_vars))
names(cb_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  var <- cold_vars[n]
  cb <- crossbasis(ne[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(cold_vars))
names(model_list) <- cold_vars

for(n in 1:length(model_list)){
  cat(n,"")
  cb <- cb_list[[n]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(ne$stratum),
               data = ne)
  model_list[[n]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_list = vector("list",length(cold_vars))
names(cp_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  cb <- cb_list[[n]]
  model <- model_list[[n]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_list[[n]] <- cp
  rm(cp)
}
# save this for north.
ne_res = list(
  cb = cb_list,
  model = model_list,
  cp = cp_list
)
# get crossreduce and then get the cumulative OR.
ne_cold_res = data.frame(
  var = cold_vars,
  RR = rep(NA, length(cold_vars)),
  L  = rep(NA, length(cold_vars)),
  U  = rep(NA, length(cold_vars))
)

for(j in 1:7){
  cb = cb_list[[j]]
  model = model_list[[j]]
  cr = crossreduce(cb,model,at=1,cen=0)
  ne_cold_res[j,2] = cr$RRfit
  ne_cold_res[j,3] = cr$RRlow
  ne_cold_res[j,4] = cr$RRhigh
}

ne_cold_res
write_csv(ne_cold_res,"results/ne_cold_res.csv")

# save this 
qs::qsave(ne_res,
        "model objects/ne_cold_results.qs")
rm(ne,ne_res,model_list,cp_list,cb_list)
rm(dlist)
rm(cb,df,df2,model)
rm(list = ls())
#==================================================================
# SOUTHEAST
se = qread("data/se_mi_mortality_with_exp_overall.qs")
se = se %>% dplyr::select(-starts_with('EHE'))
# get the SE data.
se <- se  %>% mutate(
  EHE_90 = ifelse(tmean >= quantile(tmean,0.90) ,1,0),
  EHE_95 = ifelse(tmean >= quantile(tmean,0.95) ,1,0),
  EHE_99 = ifelse(tmean >= quantile(tmean,0.99,na.rm = TRUE),1,0),
  EHE_10 = ifelse(tmean <= quantile(tmean,0.10) ,1,0),
  EHE_5 = ifelse(tmean <= quantile(tmean,0.05) ,1,0),
  EHE_01 = ifelse(tmean <= quantile(tmean,0.01,na.rm = TRUE),1,0),
  .by=loc6digit
) 

se <- se %>%  mutate(
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
glimpse(se)
# fit models.
# MODELS for COLD.
cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
              "EHE_01","EHE_01_2","EHE_01_3")


cb_list <- vector("list", length(cold_vars))
names(cb_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  var <- cold_vars[n]
  cb <- crossbasis(se[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(cold_vars))
names(model_list) <- cold_vars

for(n in 1:length(model_list)){
  cat(n,"")
  cb <- cb_list[[n]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(se$stratum),
               data = se)
  model_list[[n]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_list = vector("list",length(cold_vars))
names(cp_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  cb <- cb_list[[n]]
  model <- model_list[[n]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_list[[n]] <- cp
  rm(cp)
}
# save this for north.
se_res = list(
  cb = cb_list,
  model = model_list,
  cp = cp_list
)
# make container for results.
se_cold_res = data.frame(
  var = cold_vars,
  RR = rep(NA, length(cold_vars)),
  L  = rep(NA, length(cold_vars)),
  U  = rep(NA, length(cold_vars))
)
for(j in 1:6){
  cb = cb_list[[j]]
  model = model_list[[j]]
  cr = crossreduce(cb,model,at=1,cen=0)
  se_cold_res[j,2] = cr$RRfit
  se_cold_res[j,3] = cr$RRlow
  se_cold_res[j,4] = cr$RRhigh
}
se_cold_res
write_csv(se_cold_res,"results/se_cold_res.csv")

# save this 
qsave(se_res,
        "model objects/se_cold_results.qs")
rm(se,df,df2,cb,dlist,
   se_res,model_list,cp_list,cb_list)
rm(model)
rm(list = ls())
#==================================================================
# SOUTH
sul = qread("data/sul_mi_mortality_with_exp_overall.qs")
sul = sul %>% dplyr::select(-starts_with('EHE'))
sul <- sul  %>% mutate(
  EHE_90 = ifelse(tmean >= quantile(tmean,0.90) ,1,0),
  EHE_95 = ifelse(tmean >= quantile(tmean,0.95) ,1,0),
  EHE_99 = ifelse(tmean >= quantile(tmean,0.99,na.rm = TRUE),1,0),
  EHE_10 = ifelse(tmean <= quantile(tmean,0.10) ,1,0),
  EHE_5 = ifelse(tmean <= quantile(tmean,0.05) ,1,0),
  EHE_01 = ifelse(tmean <= quantile(tmean,0.01,na.rm = TRUE),1,0),
  .by=loc6digit
) 

sul <- sul %>%  mutate(
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




glimpse(sul)
# fit models.
# MODELS for COLD.
cold_vars = c("EHE_10","EHE_10_2","EHE_10_3",
              "EHE_01","EHE_01_2","EHE_01_3")


cb_list <- vector("list", length(cold_vars))
names(cb_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  var <- cold_vars[n]
  cb <- crossbasis(sul[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(cold_vars))
names(model_list) <- cold_vars

for(n in 1:length(model_list)){
  cat(n,"")
  cb <- cb_list[[n]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(sul$stratum),
               data = sul)
  model_list[[n]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_list = vector("list",length(cold_vars))
names(cp_list) <- cold_vars
for(n in 1:length(cold_vars)){
  cat(n,"")
  cb <- cb_list[[n]]
  model <- model_list[[n]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_list[[n]] <- cp
  rm(cp)
}
# save this for north.
sul_res = list(
  cb = cb_list,
  model = model_list,
  cp = cp_list
)
# get cumulative RR results.
sul_cold_res = data.frame(
  var = cold_vars,
  RR = rep(NA, length(cold_vars)),
  L  = rep(NA, length(cold_vars)),
  U  = rep(NA, length(cold_vars))
)
for(j in 1:6){
  cb = cb_list[[j]]
  model = model_list[[j]]
  cr = crossreduce(cb,model,at=1,cen=0)
  sul_cold_res[j,2] = cr$RRfit
  sul_cold_res[j,3] = cr$RRlow
  sul_cold_res[j,4] = cr$RRhigh
}
sul_cold_res
write_csv(sul_cold_res,"results/sul_cold_res.csv")

# save this 
qsave(sul_res,
      "model objects/sul_cold_results.qs")
rm(sul,df,df2,cb,dlist,
   sul_res,model_list,cp_list,cb_list)
rm(model)
#================================================================
