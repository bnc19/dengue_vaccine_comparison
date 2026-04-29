# ------------------------------------------------------------------------------
# run_model_TR
#
# End-to-end wrapper for fitting the multi-trial Stan model. This function
# preprocesses Qdenga, Dengvaxia, and Butantan-DV trial datasets, constructs
# the combined Stan data list via functions located in "R/format_stan_data.R", 
# compiles and runs MCMC sampling  with cmdstanr given the specified Stan model, 
# "stan_model" and saves posterior summaries, derived quantities, and diagnostics. 
# 
# Saves posterior estimates and optionally evaluation metrics to file path given
# in "folder" 
#
# Intended as a high-level execution pipeline for reproducible model fitting 
# by switching parameters on and off using the model flags. 
# ------------------------------------------------------------------------------

run_model_TR = function(n_it = 10000,
                          adapt_delta = 0.8,
                          n_chains = 4,
                          stan_model,
                          folder,
                          init_values = function() {
                            list(
                              hs = replicate(2, runif(3, 2, 20)),
                              hl = runif(3, 60, 132),
                              ts = replicate(2, runif(3, 0, 12)),
                              lambda_D_Q = replicate(6, runif(4, 0, 0.02)),
                              lambda_D_De = runif(4, 0, 0.02),
                              lambda_KD_Bu = replicate(2, runif(2, 0, 0.02)), 
                              p = replicate(3, runif(3, 0.1, 0.9)),
                              pK3 = replicate(4, runif(3, 0.1, 0.9)),
                              gamma = runif(1, 0.1, 0.9),
                              rho = runif(4, 0, 1),
                              phi = runif(1, 0.1, 0.5),
                              epsilon = runif(1, 1, 4),
                              delta = replicate(4, runif(3, 0.1, 0.8)),
                              L = replicate(4, runif(3, 0, 4)),
                              w = replicate(4, runif(3, 1, 4)),
                              alpha = replicate(4, runif(3, 0, 2)),
                              beta = replicate(2, runif(3, 0, 2)),
                              kappa =  runif(3, 0, 2),
                              omega =  runif(3, 0, 2),
                              tau = replicate(4, runif(3, 1, 5)),
                              lc = replicate(4, replicate(3, runif(3, 2, 8))),
                              sens = c(runif(1, .8, 1), runif(1, .6, 1),runif(1, .8, 1)),
                              spec = c(runif(1, 0.95, 1), runif(1, .7, 1), runif(1, 0.95, 1))
                            )
                          },
                          
                          VCD_years_Q = c(12, 18, 24, 36, 48, 54) / 12,
                          time_Q = 1:54,
                          mu_Q,
                          baseline_SP_Q,
                          VCD_Q,
                          hosp_Q,
                          baseline_SP_De,
                          cases_De,
                          mu_De,
                          cases_Bu,
                          baseline_SP_Bu,
                          mu_Bu,
                          VCD_years_Bu = c(24, 60),
                          VCD_years_De = c(13, 24, 36, 60),
                          start_time_Bu = 1,
                          end_time_Bu = 60,
                          start_time_De = 13,
                          end_time_De = 72,
                          B = 2,
                          K = 4,
                          V = 2,
                          R = 2,
                          TR = 3,
                          J_Q = 3,
                          J_De = 2,
                          J_Bu = 3,
                          C = 3,
                          HI = 12,
                          metric = c("VCD_BVKD",
                                   "VCD_BVJA",
                                   "VCD_KJ2",
                                   "HOSP_BVK4" ,
                                   "HOSP_BVJA" ,
                                   "HOSP_KJ2"),
                          BF = F,
                          diagnostics = F,
                        
                          # Model Flags 
                          share_alpha = 0, 
                          share_omega_kappa = 0, 
                          include_pK3 = 0,
                          include_eps = 0,
                          include_beta = c(0, 0, 0), # trial specific
                          mono_lc = 0, 
                          mono_lc_SN = 0,
                          mono_lc_MU = 0,
                          mono_lc_MO = 0,
                          rho_K = 0,
                          L_K = c(0,0,0),  # trial specific
                          delta_KJ = c(0,0,0),  # trial specific
                          inc_FOIJ = c(0,0,0),  # trial specific
                          w_CK = 0,
                          alpha_CK = 0,
                          uniform = 0,
                          tau_K = c(0,0,0),
                          L_sd = 1,
                          L_mean = 0,
                          psi_J = 0,
                          enhancement = c(1,1,1),
                          share_n_param = 0,
                          average_mu = 0) {
  
# check


if(mono_lc_MO == 1 & mono_lc_SN != 1) stop("mono lcMO means mono lcSN")
if(mono_lc_MO == 1 & mono_lc_MU != 1) stop("mono lcMO means mono lcMU")
  
# set up -----------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)  
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
options(mc.cores = parallel::detectCores())

# source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)

comp_model = cmdstan_model(paste0("models/", stan_model),stanc_options = list("O1"))
file_path = (paste0("output/", folder))
dir.create(file_path, recursive = T, showWarnings = F)

# data -------------------------------------------------------------------------
hosp_Q = factor_VCD_Q(hosp_Q)
VCD_Q = factor_VCD_Q(VCD_Q)

time_De = start_time_De : end_time_De
time_De = time_De - time_De[1] + 1 # start at 1

time_Bu = start_time_Bu : end_time_Bu

if(average_mu == 1){ # titres are not serotype-specific 
  mean_Q = exp(apply(log(mu_Q), 1, mean))
  mean_De = exp(apply(log(mu_De[,2:5]), 1, mean))
  mean_Bu = exp(apply(log(mu_Bu), 1, mean))

  mu_Q[1,] = mean_Q[1]
  mu_Q[2,] = mean_Q[2]
  mu_De[1,2:5] = mean_De[1]
  mu_De[2,2:5] = mean_De[2]
  mu_Bu[1,] = mean_Bu[1]
  mu_Bu[2,] = mean_Bu[2]
}

# fit Stan model ---------------------------------------------------------------
list_data = format_stan_data_Q_De(
  baseline_SP_Q = baseline_SP_Q,
  VCD_Q = VCD_Q,
  hosp_Q = hosp_Q,
  mu_Q = mu_Q,
  VCD_years_Q = VCD_years_Q,
  time_Q = time_Q,
  baseline_SP_De = baseline_SP_De,
  cases_De = cases_De,
  mu_De = mu_De,
  VCD_years_Bu = VCD_years_Bu, 
  VCD_years_De = VCD_years_De,
  cases_Bu = cases_Bu,
  baseline_SP_Bu = baseline_SP_Bu,
  mu_Bu = mu_Bu,
  time_De = time_De,
  time_Bu = time_Bu, 
  C = C,
  B = B,
  K = K,
  V = V,
  J_Q = J_Q,
  J_De = J_De,
  J_Bu = J_Bu, 
  R = R,
  TR = TR,
  HI = HI,
  uniform = uniform,
  include_pK3 = include_pK3,
  include_eps = include_eps,
  L_mean = L_mean,
  L_sd = L_sd,
  enhancement = enhancement,
  mono_lc = mono_lc, 
  mono_lc_MO = mono_lc_MO,
  mono_lc_SN = mono_lc_SN,
  mono_lc_MU = mono_lc_MU,
  rho_K = rho_K,
  w_CK = w_CK,
  alpha_CK = alpha_CK,
  tau_K = tau_K,
  include_beta = include_beta,
  L_K = L_K,
  delta_KJ = delta_KJ,
  psi_J = psi_J,
  inc_FOIJ = inc_FOIJ,
  share_n_param = share_n_param,
  share_omega_kappa = share_omega_kappa,
  share_alpha=share_alpha
)

start = Sys.time()
stan_fit = comp_model$sample(
  data = list_data,
  chains = n_chains,
  parallel_chains = n_chains,
  iter_warmup = floor(n_it/2),
  iter_sampling = floor(n_it/2),
  thin = 1,
  init = init_values,
  seed = 140,
  refresh = 500,
  adapt_delta = adapt_delta)
end = Sys.time()
print(end - start)

# save posts -------------------------------------------------------------------

fit_ext = stan_fit$draws(format = "df")

# save up to log likelihood 
i1 = which(names(fit_ext) == "ll_Q")
names_select = names(fit_ext)[1:i1]
posterior_chains = fit_ext[names_select]
write.csv(posterior_chains, paste0(file_path, "/posterior_chains.csv"))

posterior = summarise_draws(posterior_chains)
write.csv(posterior, paste0(file_path, "/posterior.csv"))

if (BF == T & L_K[1] == 0 & L_K[2] == 0 & L_K[3] == 0){ # only set up to compare single enhancement param
  calculate_bayes(
    stan_fit = stan_fit,
    L_mean = L_mean,
    L_sd = L_sd,
    file_path = file_path,
    lower = 0
  )
}

# save posterior ---------------------------------------------------------------
AR_Q = which(grepl("AR", names(fit_ext)) & grepl("Q", names(fit_ext)))
AR_out_Q = fit_ext[AR_Q] %>%  as.data.frame()
saveRDS(AR_out_Q, paste0(file_path, "/AR_Q.RDS"))

AR_De = which(grepl("AR", names(fit_ext)) & grepl("De", names(fit_ext)))
AR_out_De = fit_ext[AR_De] %>%  as.data.frame()
saveRDS(AR_out_De, paste0(file_path, "/AR_De.RDS"))

AR_Bu = which(grepl("AR", names(fit_ext)) & grepl("Bu", names(fit_ext)))
AR_out_Bu = fit_ext[AR_Bu] %>%  as.data.frame()
saveRDS(AR_out_Bu, paste0(file_path, "/AR_Bu.RDS"))

VE_Q = which(grepl("VE", names(fit_ext)) & grepl("Q", names(fit_ext)))
VE_out_Q = fit_ext[VE_Q] %>%  as.data.frame()
saveRDS(VE_out_Q, paste0(file_path, "/VE_Q.RDS"))

VE_De = which(grepl("VE", names(fit_ext)) & grepl("De", names(fit_ext)))
VE_out_De = fit_ext[VE_De] %>%  as.data.frame()
saveRDS(VE_out_De, paste0(file_path, "/VE_De.RDS"))

VE_Bu = which(grepl("VE", names(fit_ext)) & grepl("Bu", names(fit_ext)))
VE_out_Bu = fit_ext[VE_Bu] %>%  as.data.frame()
saveRDS(VE_out_Bu, paste0(file_path, "/VE_Bu.RDS"))


n_Q = which(grepl("n_Q", names(fit_ext)))
n_out_Q = fit_ext[n_Q] %>%  as.data.frame()
saveRDS(n_out_Q, paste0(file_path, "/n_Q.RDS"))

n_De = which(grepl("n_De", names(fit_ext)))
n_out_De = fit_ext[n_De] %>%  as.data.frame()
saveRDS(n_out_De, paste0(file_path, "/n_De.RDS"))

n_Bu = which(grepl("n_Bu", names(fit_ext)))
n_out_Bu = fit_ext[n_Bu] %>%  as.data.frame()
saveRDS(n_out_Bu, paste0(file_path, "/n_Bu.RDS"))


# run diagnostics --------------------------------------------------------------
if(diagnostics == T) diagnose_stan_fit(stan_fit, file_path, pars)

# WAIC 
loglik = stan_fit$draws("log_lik")
WAIC = waic(loglik)
saveRDS(WAIC, file = paste0(file_path, "/WAIC.RDS"))
saveRDS(loglik, file = paste0(file_path, "/loglik.RDS"))

# # PPC -------------------------------------------------------------------------- 
# y = y_rep = out = out2 = list()
# for (i in 1:length(metric)) {
#   y[[i]] = c(list_data[[metric[i]]])
#   y_rep[[i]] = as.matrix(fit_ext[, grepl(paste0("pred_" , metric[i]), names(fit_ext))])
#   out[[i]] = ppc_dens_overlay(y[[i]], y_rep[[i]][1:50,])
#   ggsave(out[[i]], file = paste0(file_path, "/PPC_", metric[i], ".jpg"))
#   out2[[i]] = ppc_hist(y[[i]], y_rep[[i]][1:10,], binwidth = 5)
#   ggsave(out2[[i]], file = paste0(file_path, "/PPC2_", metric[i], ".jpg"))
# }

}
