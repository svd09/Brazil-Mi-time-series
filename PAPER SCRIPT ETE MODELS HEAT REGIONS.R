#-------------------------#
# regional models for ETE #
#-------------------------#
# PAPER SCRIPT
# MODEL ETE FOR EACH REGION AND SAVE OBJ.
# USE THEM TO GET AF/AN LATER.

#------#
# HEAT #
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
glimpse(dlist[[1]])
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

north <- bind_rows(dlist)
glimpse(north)

# MODELS for HEAT
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")


cb_list <- vector("list", length(heat_vars))
names(cb_list) <- heat_vars
for(n in 1:length(heat_vars)){
  cat(n,"")
  var <- heat_vars[n]
  cb <- crossbasis(north[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(heat_vars))
names(model_list) <- heat_vars

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
cp_list = vector("list",length(heat_vars))
names(cp_list) <- heat_vars
for(n in 1:length(heat_vars)){
  cat(n,"")
  cb <- cb_list[[n]]
  model <- model_list[[n]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_list[[n]] <- cp
  rm(cp)
}
# save this for north.
north_res = list(
  cb = cb_list,
  model = model_list,
  cp = cp_list
)
# get the cumulative RR for heat.
north_res_heat = data.frame(
  var = heat_vars,
  rr = rep(NA,length(heat_vars)),
  rr_l = rep(NA,length(heat_vars)),
  rr_u = rep(NA,length(heat_vars))
)
for(n in 1:length(heat_vars)){
  cat(n,"")
  cp <- cb_list[[n]]
  model = model_list[[n]]
  cr = crossreduce(cb,model,at=1,cen=0)
  north_res_heat[n,2] = cr$RRfit
  north_res_heat[n,3] = cr$RRlow
  north_res_heat[n,4] = cr$RRhigh
}
north_res_heat
write_csv(north_res_heat,
          "results/north_heat_results.csv")
# save this 
qs::qsave(north_res,
          "model objects/north_heat_results.qs")
rm(north,north_res,model_list,cp_list,cb_list)
rm(dlist)
rm(cb,df,df2,model)
rm(list = ls())
#===============================================================================
# NORTHEAST.
ne = qread("data/ne_mi_mortality_with_exp_overall.qs")
# we need to get the ETE variables separately for each region.
# ETE will change for each region based on their tmean.
glimpse(ne)
ne = ne %>% dplyr::select(-starts_with('EHE'))
glimpse(ne)

# make a list of the data.
dlist <- split(ne, ne$loc6digit)
names(dlist) <- sapply(dlist, function(x) unique(x$state_loc.x))
glimpse(dlist[[1]])
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

ne <- bind_rows(dlist)
glimpse(ne)

# MODELS for HEAT
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")


cb_list <- vector("list", length(heat_vars))
names(cb_list) <- heat_vars
for(n in 1:length(heat_vars)){
  cat(n,"")
  var <- heat_vars[n]
  cb <- crossbasis(ne[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(heat_vars))
names(model_list) <- heat_vars

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
cp_list = vector("list",length(heat_vars))
names(cp_list) <- heat_vars
for(n in 1:length(heat_vars)){
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
ne_res_heat = data.frame(
  var = heat_vars,
  rr = rep(NA,length(heat_vars)),
  rr_l = rep(NA,length(heat_vars)),
  rr_u = rep(NA,length(heat_vars))
)
for(n in 1:length(heat_vars)){
  cat(n,"")
  cp <- cb_list[[n]]
  model = model_list[[n]]
  cr = crossreduce(cb,model,at=1,cen=0)
  ne_res_heat[n,2] = cr$RRfit
  ne_res_heat[n,3] = cr$RRlow
  ne_res_heat[n,4] = cr$RRhigh
}
ne_res_heat
write_csv(ne_res_heat,
          "results/ne_heat_results.csv")

# save this 
qs::qsave(ne_res,
          "model objects/ne_heat_results.qs")
rm(ne,ne_res,model_list,cp_list,cb_list)
rm(dlist)
rm(cb,df,df2,model)
gc()
rm(list = ls())
#===============================================================================
# SOUTHEAST
se = qread("data/se_mi_mortality_with_exp_overall.qs")
# we need to get the ETE variables separately for each region.
# ETE will change for each region based on their tmean.
glimpse(se)
se = se %>% dplyr::select(-starts_with('EHE'))
glimpse(se)

# make a list of the data.
dlist <- split(se, se$loc6digit)
names(dlist) <- sapply(dlist, function(x) unique(x$state_loc.x))
glimpse(dlist[[1]])
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

se <- bind_rows(dlist)
glimpse(se)

# MODELS for HEAT
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")



cb_list <- vector("list", length(heat_vars))
names(cb_list) <- heat_vars
for(n in 1:length(heat_vars)){
  cat(n,"")
  var <- heat_vars[n]
  cb <- crossbasis(se[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(heat_vars))
names(model_list) <- heat_vars

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
cp_list = vector("list",length(heat_vars))
names(cp_list) <- heat_vars
for(n in 1:length(heat_vars)){
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
se_res_heat = data.frame(
  var = heat_vars,
  rr = rep(NA,length(heat_vars)),
  rr_l = rep(NA,length(heat_vars)),
  rr_u = rep(NA,length(heat_vars))
)
for(n in 1:length(heat_vars)){
  cat(n,"")
  cp <- cb_list[[n]]
  model = model_list[[n]]
  cr = crossreduce(cb,model,at=1,cen=0)
  se_res_heat[n,2] = cr$RRfit
  se_res_heat[n,3] = cr$RRlow
  se_res_heat[n,4] = cr$RRhigh
}
se_res_heat
write_csv(se_res_heat,
          "results/se_heat_results.csv")

# save this 
qs::qsave(se_res,
          "model objects/se_heat_results.qs")
rm(se,se_res,model_list,cp_list,cb_list)
rm(dlist)
rm(cb,df,df2,model)
gc()
rm(list = ls())
#======================================================================
# SOUTH
sul = qread("data/sul_mi_mortality_with_exp_overall.qs")
# we need to get the ETE variables separately for each region.
# ETE will change for each region based on their tmean.
glimpse(sul)
sul = sul %>% dplyr::select(-starts_with('EHE'))
glimpse(sul)

# make a list of the data.
dlist <- split(sul, sul$loc6digit)
names(dlist) <- sapply(dlist, function(x) unique(x$state_loc.x))
glimpse(dlist[[1]])
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

sul <- bind_rows(dlist)
glimpse(sul)

# MODELS for HEAT
heat_vars = c("EHE_90","EHE_90_2","EHE_90_3",
              "EHE_99","EHE_99_2","EHE_99_3")



cb_list <- vector("list", length(heat_vars))
names(cb_list) <- heat_vars
for(n in 1:length(heat_vars)){
  cat(n,"")
  var <- heat_vars[n]
  cb <- crossbasis(sul[,var], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_list[[n]] <- cb
  rm(cb)
}

model_list <- vector("list", length(heat_vars))
names(model_list) <- heat_vars

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
cp_list = vector("list",length(heat_vars))
names(cp_list) <- heat_vars
for(n in 1:length(heat_vars)){
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
sul_res_heat = data.frame(
  var = heat_vars,
  rr = rep(NA,length(heat_vars)),
  rr_l = rep(NA,length(heat_vars)),
  rr_u = rep(NA,length(heat_vars))
)
for(n in 1:length(heat_vars)){
  cat(n,"")
  cp <- cb_list[[n]]
  model = model_list[[n]]
  cr = crossreduce(cb,model,at=1,cen=0)
  sul_res_heat[n,2] = cr$RRfit
  sul_res_heat[n,3] = cr$RRlow
  sul_res_heat[n,4] = cr$RRhigh
}
sul_res_heat
write_csv(sul_res_heat,
          "results/sul_heat_results.csv")

# save this 
qs::qsave(sul_res,
          "model objects/sul_heat_results.qs")
rm(sul,sul_res,
   model_list,cp_list,cb_list)
rm(dlist)
rm(cb,df,df2,model)
gc()
# END.
# SCRIPT DONE. ALL REGIONS MODELS FOR HEAT DONE AND SAVED.
#=================================================================