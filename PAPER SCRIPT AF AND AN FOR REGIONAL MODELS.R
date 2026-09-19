#---------------------------#
# AF/AN for regional models #
#---------------------------#
# PAPER SCRIPT
# GET AF/AN FROM THE REGIONAL MODELS SAVED FROM THE EARLIER SCRIPT.
# load libaries.
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
# get the AF and AN function code.
source("codes/AF and AN for conditional poisson models function.R")
#-------#
# NORTH.#
#-------#

north <- readRDS("model objects/north_cold_results.rds")
n_data <- qread("data/north_mi_mortality_with_exp_overall.qs")

cold_vars <- names(north$cb)
# make an empty list to hold results.


# # EHE_10
# L <- 7
# lags <- 0:L
# nsim <- 1000
# set.seed(1974)
# # get the exposure and outcome data.
# x <- n_data$EHE_10
# y <- n_data$events
# # get coef.
# model <- north$model$EHE_10
# coef_all <- coef(model)
# vcov_all <- vcov(model)
# # get only the EHE coefficients.
# ix = grep("^cb", names(coef_all))
# stopifnot(length(ix)> 0)
# b <- coef_all[ix]
# V <- vcov_all[ix,ix,drop = F]
# # define cb
# cb <- north$cb$EHE_10
# al = attr(cb,"arglag")
# L = attr(cb,"lag")
# B <- do.call(dlnm::onebasis, c(list(x = lags), al))
# b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V)
# AF_draws <- pblapply(seq_len(nrow(b_draws)),function(i){
#   bi <- b_draws[i,]
#   beta_l <- as.numeric(B %*% bi)
#   af_an_from_beta(beta_l, x, y, trunc_excess = TRUE)
# })
# bind_rows(Af_draws)
# 
# Af_draws <- pblapply(seq_len(nrow(b_draws)),function(i){
#   bi <- b_draws[i,]
#   beta_l <- as.numeric(B %*% bi)
#   af_an_from_beta(beta_l, x, y, trunc_excess = TRUE)
# })

results_df = NULL
gc()
for(k in 1:length(cold_vars)){
  cat(k,"working")
  var <- cold_vars[k]
  L <- 7
  lags <- 0:L
  nsim <- 1000
  set.seed(1974)
  # get the exposure and outcome data.
  x <- n_data[[var]]
  x <- as.numeric(x) 
  y <- n_data$events
  # get coef.
  model <- north$model[[var]]
  coef_all <- coef(model)
  vcov_all <- vcov(model)
  # get only the EHE coefficients.
  ix = grep("^cb", names(coef_all))
  stopifnot(length(ix)> 0)
  b <- coef_all[ix]
  V <- vcov_all[ix,ix,drop = F]
  # define cb
  cb <- north$cb[[var]]
  al = attr(cb,"arglag")
  L = attr(cb,"lag")
  B <- do.call(dlnm::onebasis, c(list(x = lags), al))
  b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V)
  AF_draws <- pblapply(seq_len(nrow(b_draws)),function(i){
    bi <- b_draws[i,]
    beta_l <- as.numeric(B %*% bi)
    af_an_from_beta(beta_l, x, y, trunc_excess = TRUE)
  })
  res_df <- bind_rows(AF_draws)
  colnames(res_df) = c(paste0("AF","_",var),
                       paste0("AN","_",var))
  if(is.null(results_df)){
    results_df = res_df
  } else {
    results_df = cbind(results_df,res_df)
  }
  rm(model,cb) # remove the large objects to save memory.
}

glimpse(results_df)
summary(results_df)
# check the number per day of event.
sum(n_data$EHE_10)

(435.7/sum(n_data$EHE_10))*1000
(374.2/sum(n_data$EHE_10_2))*1000
(128.41/sum(n_data$EHE_10_3))*1000


write_csv(results_df,
          "results/north_cold_af_an_results.csv")
gc()
rm(north,n_data,results_df) # clear the env.
#-----#
# SE  #
#-----#

# get SE region.
rm(list = ls())
se_data = qread("data/se_mi_mortality_with_exp_overall.qs")
se <- qread("model objects/se_cold_results.qs")
source("codes/AF and AN for conditional poisson models function.R")

cold_vars <- names(se$cb)
gc()

results_df = NULL
gc()
for(k in 1:length(cold_vars)){
  cat(k,"working")
  var <- cold_vars[k]
  L <- 7
  lags <- 0:L
  nsim <- 1000
  set.seed(1974)
  # get the exposure and outcome data.
  x <- se_data[[var]]
  x <- as.numeric(x) 
  y <- se_data$events
  # get coef.
  model <- se$model[[var]]
  coef_all <- coef(model)
  vcov_all <- vcov(model)
  # get only the EHE coefficients.
  ix = grep("^cb", names(coef_all))
  stopifnot(length(ix)> 0)
  b <- coef_all[ix]
  V <- vcov_all[ix,ix,drop = F]
  # define cb
  cb <- se$cb[[var]]
  al = attr(cb,"arglag")
  L = attr(cb,"lag")
  B <- do.call(dlnm::onebasis, c(list(x = lags), al))
  b_draws <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V)
  AF_draws <- pblapply(seq_len(nrow(b_draws)),function(i){
    bi <- b_draws[i,]
    beta_l <- as.numeric(B %*% bi)
    af_an_from_beta(beta_l, x, y, trunc_excess = TRUE)
  })
  res_df <- bind_rows(AF_draws)
  colnames(res_df) = c(paste0("AF","_",var),
                       paste0("AN","_",var))
  if(is.null(results_df)){
    results_df = res_df
  } else {
    results_df = cbind(results_df,res_df)
  }
  rm(model,cb) # remove the large objects to save memory.
}
results_df
write_csv(results_df,
          "results/se_cold_af_an_results.csv")

rm(se_data,se,results_df) # clear the env.
rm(list = ls())
