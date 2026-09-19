##################################################################
##                    STEMI mortality models                    ##
##################################################################
# Salil V Deo 
# SCRIPT FOR THE PAPER
# MODEL TEMP AS ETE CATEGORICAL AND REPORT OR.
# ALSO OBTAIN AF AND AN FOR THE ETE.
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

# get function for AF/AN calculations.
source("codes/AF and AN for conditional poisson models function.R")

# get the data.
# WHOLE DATA
whole_data <- readRDS("data/mi_mortality_wf.rds")
# MODELS for COLD.
# ETE 01.
cb <- crossbasis(whole_data$EHE_01, lag = 7,
                 argvar = list(fun = "lin"),
                 arglag = list(fun = "ns", knots = logknots(7, 2)))
cb
formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
# removed the other variables as mem does not allow full model.
# this should be fine for the paper.
model <- gnm(formula = formula,
             family = quasipoisson(),
             eliminate = factor(whole_data$stratum),
             data = whole_data)
cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
# see the trends in OR across lag terms.
# get the final value.

fit <- as.numeric(cp$cumRRfit[1,])
low <- as.numeric(cp$cumRRlow[1,])
high <- as.numeric(cp$cumRRhigh[1,])
lag <- as.factor(seq(0, 7, 1))

df_ehe01 <- data.frame(lag, fit, low, high)
write_csv(df_ehe01,
          "results/ehe_01_overall.csv")
# calculating AF and AN
summary(model)

L = 7
lags <- 0:L
nsim <- 100
set.seed(1974)

#- get the exposure and outcome data.
x <- whole_data$EHE_01
y <- whole_data$events

# get coef.
coef_all <- coef(model)
vcov_all <- vcov(model)

#- get only the EHE coefficients.
ix = grep("^cb", names(coef_all))
stopifnot(length(ix)> 0)

b <- coef_all[ix]
V <- vcov_all[ix,ix,drop = F]

al = attr(cb,"arglag")
L = attr(cb,"lag")

B <- do.call(dlnm::onebasis, c(list(x = lags), al))

b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V)

AF_draws <- pblapply(seq_len(nrow(b_draws)),function(i){
  bi <- b_draws[i,]
  beta_l <- as.numeric(B %*% bi)
  af_an_from_beta(beta_l, x, y, trunc_excess = TRUE)
})

results = do.call(rbind, AF_draws)
summary(results)
# AF               AN_tot    
# Min.   :0.004369   Min.   :6190  
# 1st Qu.:0.004547   1st Qu.:6442  
# Median :0.004654   Median :6593  
# Mean   :0.004652   Mean   :6591  
# 3rd Qu.:0.004753   3rd Qu.:6734  
# Max.   :0.004972   Max.   :7044   


# example plot.
ggplot(data = df,aes(x = lag, y = fit, color=lag)) +
  geom_point() +
  geom_errorbar(aes(ymin = low, ymax = high), alpha = 0.3,width=0,linewidth=1.5) +
  labs(
    title = "Cumulative OR for EHE_01",
    x = "Lag (days)",
    y = "Cumulative OR"
  ) +
  theme_minimal()

rm(cb,model,cp)

# EHE 10
cb <- crossbasis(whole_data$EHE_10, lag = 7,
                 argvar = list(fun = "lin"),
                 arglag = list(fun = "ns", knots = logknots(7, 2)))
formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
model <- gnm(formula = formula,
             family = quasipoisson(),
             eliminate = factor(whole_data$stratum),
             data = whole_data)
cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)

fit <- as.numeric(cp$cumRRfit[1,])
low <- as.numeric(cp$cumRRlow[1,])
high <- as.numeric(cp$cumRRhigh[1,])
lag <- as.factor(seq(0, 7, 1))

df_ehe10 <- data.frame(lag, fit, low, high)
write_csv(df_ehe10,
          "results/ehe_10_overall.csv")
# Calculating AF/AN.

L = 7
lags <- 0:L
nsim <- 100
set.seed(1974)

#- get the exposure and outcome data.
x <- whole_data$EHE_10
y <- whole_data$events

# get coef.
coef_all <- coef(model)
vcov_all <- vcov(model)

#- get only the EHE coefficients.
ix = grep("^cb", names(coef_all))
stopifnot(length(ix)> 0)

b <- coef_all[ix]
V <- vcov_all[ix,ix,drop = F]

al = attr(cb,"arglag")
L = attr(cb,"lag")

B <- do.call(dlnm::onebasis, c(list(x = lags), al))
nsim = 100
b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V)

AF_draws <- pblapply(seq_len(nrow(b_draws)),function(i){
  bi <- b_draws[i,]
  beta_l <- as.numeric(B %*% bi)
  af_an_from_beta(beta_l, x, y, trunc_excess = TRUE)
})
results = do.call(rbind, AF_draws)
summary(results)
# AF              AN_tot     
# Min.   :0.01653   Min.   :23419  
# 1st Qu.:0.01703   1st Qu.:24132  
# Median :0.01733   Median :24559  
# Mean   :0.01733   Mean   :24553  
# 3rd Qu.:0.01761   3rd Qu.:24951  
# Max.   :0.01824   Max.   :25847  

rm(cp,cb,model) # clear the env.
#------------------------------------------------------------------

# consecutive days of cold temp.
# EHE 10
model_2 = c("EHE_10_2","EHE_10_3")
cb_2 <- vector("list", length(model_2))
names(cb_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- crossbasis(whole_data[[var]], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_2[[n]] <- cb
  rm(cb)
}
model_list_2 <- vector("list", length(model_2))
names(model_list_2) <- model_2

for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(whole_data$stratum),
               data = whole_data)
  model_list_2[[var]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_2 = vector("list",length(model_2))
names(cp_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  model <- model_list_2[[var]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_2[[var]] <- cp
  rm(cp)
}

cp_2[[1]]
cp_2[[2]]
# save this list.
results_ehe_10_2_3 <- list(cb = cb_2,
                           model = model_list_2,
                           cp = cp_2)
saveRDS(results_ehe_10_2_3,
        "model objects/results_ehe_10_2_3.rds")
rm(model_list_2,cb_2,cp_2) # remove these so that we can run this again.
# am going to calculate the AF and AN later for all these objects together.
#--------------------------------------------------
# EHE 01 2 and 3

model_2 = c("EHE_01_2","EHE_01_3")
cb_2 <- vector("list", length(model_2))
names(cb_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- crossbasis(whole_data[[var]], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_2[[n]] <- cb
  rm(cb)
}
model_list_2 <- vector("list", length(model_2))
names(model_list_2) <- model_2

for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(whole_data$stratum),
               data = whole_data)
  model_list_2[[var]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_2 = vector("list",length(model_2))
names(cp_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  model <- model_list_2[[var]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_2[[var]] <- cp
  rm(cp)
}

cp_2[[1]]
cp_2[[2]]
# save this list.
results_ehe_01_2_3 <- list(cb = cb_2,
                           model = model_list_2,
                           cp = cp_2)
saveRDS(results_ehe_01_2_3,
        "model objects/results_ehe_01_2_3.rds")
rm(cb_2,cp_2,model,model_list_2,cb,results_ehe_01_2_3)
#-- ALL COLD MODELS DONE.
#-------------------------------------------------------
#--------------#
# HEAT MODELS  #
#--------------#
# EHE_90, EHE_99
model_h1 = c("EHE_90","EHE_99")
cb_h1 <- vector("list", length(model_2))
names(cb_h1) <- model_h1
for(n in 1:length(model_h1)){
  cat(n,"")
  var <- model_h1[n]
  cb <- crossbasis(whole_data[[var]], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_h1[[n]] <- cb
  rm(cb)
}
model_list_h1 <- vector("list", length(model_h1))
names(model_list_h1) <- model_h1

for(n in 1:length(model_h1)){
  cat(n,"")
  var <- model_h1[n]
  cb <- cb_h1[[var]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(whole_data$stratum),
               data = whole_data)
  model_list_h1[[var]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_h1 = vector("list",length(model_h1))
names(cp_h1) <- model_h1
for(n in 1:length(model_h1)){
  cat(n,"")
  var <- model_h1[n]
  cb <- cb_h1[[var]]
  model <- model_list_h1[[var]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_h1[[var]] <- cp
  rm(cp)
}
cp_h1[[1]]
cp_h1[[2]]
# save this list.
results_ehe_90_99 <- list(cb = cb_h1,
                           model = model_list_h1,
                           cp = cp_h1)
saveRDS(results_ehe_90_99,
        "model objects/results_ehe_90_99.rds")
rm(cb_h1,cp_h1,model,model_list_h1,cb)
#---------------------------------------------------------------
# EHE 90 2 and 3
model_2 = c("EHE_90_2","EHE_90_3")
cb_2 <- vector("list", length(model_2))
names(cb_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- crossbasis(whole_data[[var]], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_2[[n]] <- cb
  rm(cb)
}
model_list_2 <- vector("list", length(model_2))
names(model_list_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(whole_data$stratum),
               data = whole_data)
  model_list_2[[var]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_2 = vector("list",length(model_2))
names(cp_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  model <- model_list_2[[var]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_2[[var]] <- cp
  rm(cp)
}
cp_2[[1]]
cp_2[[2]]
# save this list.
results_ehe_90_2_3 <- list(cb = cb_2,
                           model = model_list_2,
                           cp = cp_2)
saveRDS(results_ehe_90_2_3,
        "model objects/results_ehe_90_2_3.rds")
rm(cb_2,cp_2,model,model_list_2,cb,results_ehe_90_2_3)
#-----------------------------------------------------------
# EHE 99 2 and 3
model_2 = c("EHE_99_2","EHE_99_3")
cb_2 <- vector("list", length(model_2))
names(cb_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- crossbasis(whole_data[[var]], lag = 7,
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "ns", knots = logknots(7, 2)))
  cb_2[[n]] <- cb
  rm(cb)
}
model_list_2 <- vector("list", length(model_2))
names(model_list_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
  # removed the other variables as mem does not allow full model.
  # this should be fine for the paper.
  model <- gnm(formula = formula,
               family = quasipoisson(),
               eliminate = factor(whole_data$stratum),
               data = whole_data)
  model_list_2[[var]] <- model
  rm(model)
}
gc()
# get only the cumulative RR from each model.
cp_2 = vector("list",length(model_2))
names(cp_2) <- model_2
for(n in 1:length(model_2)){
  cat(n,"")
  var <- model_2[n]
  cb <- cb_2[[var]]
  model <- model_list_2[[var]]
  cp <- crosspred(cb, model, cen = 0, at = 1, cumul = TRUE)
  cp_2[[var]] <- cp
  rm(cp)
}

results_ehe_99_2_3 = list(
  cb_2,model_list_2,cp_2
)

saveRDS(results_ehe_99_2_3,
        "model objects/results_ehe_99_2_3.rds")
rm(cp_2,model_list_2,cb_2,model,cb,results_ehe_99_2_3)
#-----------------------------------------------------------#
# END                                                       #
# next script get region specific models for each location. #
# run them and save all the results together.               #
# END #-----------------------------------------------------#
# ne = qread("data/ne_mi_mortality_with_exp_overall.qs")
# # 
# 
# 
# north = qread("data/north_mi_mortality_with_exp_overall.qs")
# se = qread("data/se_mi_mortality_with_exp_overall.qs")
# sul = qread("data/sul_mi_mortality_with_exp_overall.qs")