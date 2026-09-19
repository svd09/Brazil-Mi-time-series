##################################################################
##                    STEMI mortality models                    ##
##################################################################
# Salil V Deo 
# SCRIPT FOR THE PAPER
# REVISED SCRIPT - GETTING AF & AN FOR EACH EVENT DAY OF THE EXTREME EVENT
# THIS IS BECAUSE CUMULATIVE EXTREME DAYS ARE FEWER HENCE THE RATE & NUMBER MAY BE LOWER FOR THESE 
# EVENTS BECAUSE OF THAT - ADJUSTING FOR THE NUMBER OF DAYS CORRECTS THIS.
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
# load function
source("codes/AN ON EVENT DAYS FUNCTION.R")
# run the code from before using this function instead.

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
  af_an_from_beta_event_days(beta_l, x, y, trunc_excess = TRUE)
})

results = do.call(rbind, AF_draws)
summary(results)
# > summary(results)
# AF_all           AN_tot_all     n_days_all        n_event_days    AN_sum_event_days
# Min.   :0.004369   Min.   :6190   Min.   :18737810   Min.   :500089   Min.   :2197     
# 1st Qu.:0.004547   1st Qu.:6442   1st Qu.:18737810   1st Qu.:500089   1st Qu.:2295     
# Median :0.004654   Median :6593   Median :18737810   Median :500089   Median :2351     
# Mean   :0.004652   Mean   :6591   Mean   :18737810   Mean   :500089   Mean   :2364     
# 3rd Qu.:0.004753   3rd Qu.:6734   3rd Qu.:18737810   3rd Qu.:500089   3rd Qu.:2428     
# Max.   :0.004972   Max.   :7044   Max.   :18737810   Max.   :500089   Max.   :2627     
# AN_per_event_day   AF_on_event_days 
# Min.   :0.004394   Min.   :0.06017  
# 1st Qu.:0.004590   1st Qu.:0.06284  
# Median :0.004702   Median :0.06438  
# Mean   :0.004728   Mean   :0.06474  
# 3rd Qu.:0.004854   3rd Qu.:0.06647  
# Max.   :0.005253   Max.   :0.07193 

write_csv(results,
          "results/ehe_01_af_an_event_days_per_day.csv")

#-------------------------------------------------------------------------------

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
  af_an_from_beta_event_days(beta_l, x, y, trunc_excess = TRUE)
})
results = do.call(rbind, AF_draws)
summary(results)
# > summary(results)
# AF_all          AN_tot_all      n_days_all        n_event_days     AN_sum_event_days
# Min.   :0.01653   Min.   :23419   Min.   :18737810   Min.   :3084152   Min.   :14994    
# 1st Qu.:0.01703   1st Qu.:24132   1st Qu.:18737810   1st Qu.:3084152   1st Qu.:15385    
# Median :0.01733   Median :24559   Median :18737810   Median :3084152   Median :15635    
# Mean   :0.01733   Mean   :24553   Mean   :18737810   Mean   :3084152   Mean   :15653    
# 3rd Qu.:0.01761   3rd Qu.:24951   3rd Qu.:18737810   3rd Qu.:3084152   3rd Qu.:15912    
# Max.   :0.01824   Max.   :25847   Max.   :18737810   Max.   :3084152   Max.   :16530    
# AN_per_event_day   AF_on_event_days 
# Min.   :0.004862   Min.   :0.06302  
# 1st Qu.:0.004988   1st Qu.:0.06467  
# Median :0.005069   Median :0.06572  
# Mean   :0.005075   Mean   :0.06579  
# 3rd Qu.:0.005159   3rd Qu.:0.06688  
# Max.   :0.005360   Max.   :0.06948  

write_csv(results,
          "results/ehe_10_af_an_event_days_per_day.csv")
#-------------------------------------------------------------------------------
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
# CUMULATIVE EHE_10_2 & EHE_10_3
source("codes/AN ON EVENT DAYS FUNCTION.R")
whole_data <- readRDS("data/mi_mortality_wf.rds")

colnames(whole_data)

# MODELS for COLD.
# EHE_01_3
cb_ehe_01_2 <- crossbasis(whole_data$EHE_01_2, lag = 7,
                            argvar = list(fun = "lin"),
                            arglag = list(fun = "ns", knots = logknots(7, 2)))
formula <- as.formula(paste("events ~ cb_ehe_01_2 + so2 + co + no2 + go3 + pm25"))
model_ehe_01_2 <- gnm(formula = formula,
                         family = quasipoisson(),
                         eliminate = factor(whole_data$stratum),
                         data = whole_data)

# EHE_01_3
cb_ehe_01_3 <- crossbasis(whole_data$EHE_01_3, lag = 7,
                          argvar = list(fun = "lin"),
                          arglag = list(fun = "ns", knots = logknots(7, 2)))
formula <- as.formula(paste("events ~ cb_ehe_01_3 + so2 + co + no2 + go3 + pm25"))
model_ehe_01_3 <- gnm(formula = formula,
                      family = quasipoisson(),
                      eliminate = factor(whole_data$stratum),
                      data = whole_data)


#- get the AF/AN so that can remove these large obj.
# EHE_01_2
L = 7
lags <- 0:L
nsim <- 100
set.seed(1974)

#- get the exposure and outcome data.
x <- whole_data$EHE_01_2
y <- whole_data$events

# get coef.
coef_all <- coef(model_ehe_01_2)
vcov_all <- vcov(model_ehe_01_2)

#- get only the EHE coefficients.
ix = grep("^cb", names(coef_all))
stopifnot(length(ix)> 0)

b <- coef_all[ix]
V <- vcov_all[ix,ix,drop = F]

al = attr(cb_ehe_01_2,"arglag")
L = attr(cb_ehe_01_2,"lag")

B <- do.call(dlnm::onebasis, c(list(x = lags), al))
nsim = 100
b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V)

AF_draws_ehe_01_2 <- pblapply(seq_len(nrow(b_draws)),function(i){
  bi <- b_draws[i,]
  beta_l <- as.numeric(B %*% bi)
  af_an_from_beta_event_days(beta_l, x, y, trunc_excess = TRUE)
})
results_01_2 = do.call(rbind, AF_draws_ehe_01_2)
summary(results_01_2)
rm(list = setdiff(ls(),c("results_01_2","AF_draws_ehe_01_2",
                              "model_ehe_01_2","cb_ehe_01_2",
                              "models_ehe","cb_ehe")))
gc()







# EHE_10_3
cb_ehe_10_3 <- crossbasis(whole_data$EHE_10_3, lag = 7,
                          argvar = list(fun = "lin"),
                          arglag = list(fun = "ns", knots = logknots(7, 2)))
formula <- as.formula(paste("events ~ cb_ehe_10_3 + so2 + co + no2 + go3 + pm25"))
model_ehe_10_3 <- gnm(formula = formula,
                      family = quasipoisson(),
                      eliminate = factor(whole_data$stratum),
                      data = whole_data)

# EHE_10_2
cb_ehe_10_2 <- crossbasis(whole_data$EHE_10_2, lag = 7,
                          argvar = list(fun = "lin"),
                          arglag = list(fun = "ns", knots = logknots(7, 2)))
formula <- as.formula(paste("events ~ cb_ehe_10_2 + so2 + co + no2 + go3 + pm25"))
model_ehe_10_2 <- gnm(formula = formula,
                      family = quasipoisson(),
                      eliminate = factor(whole_data$stratum),
                      data = whole_data)

#- now we only need these two obj.
# so remove every thing else from the global env.
keep = c("models_ehe","cb_ehe")
rm(list = setdiff(ls(envir = .GlobalEnv),keep),
   envir = .GlobalEnv)
gc()
#load the function
source("codes/AN ON EVENT DAYS FUNCTION.R")


rm(list = ls())

cb <- crossbasis(whole_data$EHE_01_2, lag = 7,
                 argvar = list(fun = "lin"),
                 arglag = list(fun = "ns", knots = logknots(7, 2)))
cb
formula <- as.formula(paste("events ~ cb + so2 + co + no2 + go3 + pm25"))
# removed the other variables as mem does not allow full model.
# this should be fine for the paper.
model_01_2 <- gnm(formula = formula,
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
  af_an_from_beta_event_days(beta_l, x, y, trunc_excess = TRUE)
})

results = do.call(rbind, AF_draws)
summary(results)