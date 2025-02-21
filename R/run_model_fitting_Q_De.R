run_model_Q_DE = function(n_it = 10000,
                          adapt_delta = 0.8,
                          n_chains = 4,
                          stan_model,
                          folder,
                          init_values = function() {
                            list(
                              hs = runif(2, 2, 20),
                              hl = runif(1, 60, 132),
                              ts = runif(2, 0, 12),
                              lambda_D_Q = replicate(6, runif(4, 0, 0.02)),
                              lambda_D_De = runif(4, 0, 0.02),
                              p = replicate(3, runif(2, 0.1, 0.9)),
                              pK3 = replicate(4, runif(2, 0.1, 0.9)),
                              gamma = runif(1, 0.1, 0.9),
                              rho = runif(4, 1, 5),
                              phi = runif(1, 0.1, 0.5),
                              epsilon = runif(1, 1, 4),
                              delta = replicate(4, runif(2, 0.1, 0.8)),
                              L = replicate(4, runif(2, 0, 4)),
                              w = replicate(4, runif(2, 1, 4)),
                              alpha = replicate(4, runif(2, 0, 2)),
                              beta = replicate(2, runif(2, 0, 2)),
                              tau = replicate(4, runif(2, 1, 5)),
                              lc = replicate(4, replicate(3, runif(2, 2, 8))),
                              sens = c(runif(1, .8, 1), runif(1, .6, 1)),
                              spec = c(runif(1, 0.95, 1), runif(1, .7, 1))
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
                          VCD_years_De = c(13, 24, 36, 60),
                          start_time_De = 13,
                          end_time_De = 72,
                          B = 2,
                          K = 4,
                          V = 2,
                          R = 2,
                          TR = 2,
                          J_Q = 3,
                          J_De = 2,
                          C = 3,
                          HI = 12,
                          include_pK3 = 0,
                          include_eps = 0,
                          include_beta = c(0, 0),
                          # trial specific
                          mono_lc_SN = 0,
                          mono_lc_MU = 0,
                          rho_K = 0,
                          L_K = c(0, 0),
                          # trial specific
                          delta_KJ = c(0, 0),
                          # trial specific
                          w_CK = 0,
                          alpha_CK = 0,
                          uniform = 0,
                          tau_K = 0,
                          L_sd = 1,
                          L_mean = 0,
                          psi_J = 0,
                          enhancement = 1,
                          diagnostics = F,
                          metric = c("VCD_BVKD",
                                     "VCD_BVJA",
                                     "VCD_KJ2",
                                     "HOSP_BVK4" ,
                                     "HOSP_BVJA" ,
                                     "HOSP_KJ2"),
                          BF = F) {
  
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
  VCD_years_De = VCD_years_De,
  time_De = time_De,
  C = C,
  B = B,
  K = K,
  V = V,
  J_Q = J_Q,
  J_De = J_De,
  R = R,
  TR = TR,
  HI = HI,
  uniform = uniform,
  include_pK3 = include_pK3,
  include_eps = include_eps,
  L_mean = L_mean,
  L_sd = L_sd,
  enhancement = enhancement,
  mono_lc_SN = mono_lc_SN,
  mono_lc_MU = mono_lc_MU,
  rho_K = rho_K,
  w_CK = w_CK,
  alpha_CK = alpha_CK,
  tau_K = tau_K,
  include_beta = include_beta,
  L_K = L_K,
  delta_KJ = delta_KJ,
  psi_J = psi_J
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
  seed = 14,
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

if (BF == T & L_K[1] == 0 & L_K[2] == 0){
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

VE = which(grepl("VE" , names(fit_ext)))
VE_out = fit_ext[VE] %>%  as.data.frame()
saveRDS(VE_out, paste0(file_path, "/VE.RDS"))

VE_Q = which(grepl("VE", names(fit_ext)) & grepl("Q", names(fit_ext)))
VE_out_Q = fit_ext[VE_Q] %>%  as.data.frame()
saveRDS(VE_out_Q, paste0(file_path, "/VE_Q.RDS"))

VE_De = which(grepl("VE", names(fit_ext)) & grepl("De", names(fit_ext)))
VE_out_De = fit_ext[VE_De] %>%  as.data.frame()
saveRDS(VE_out_De, paste0(file_path, "/VE_De.RDS"))

# n = which(grepl("n" , names(fit_ext)))
# n_out = fit_ext[n] %>%  as.data.frame()
# saveRDS(n_out, paste0(file_path, "/n.RDS"))


# run diagnostics --------------------------------------------------------------
if(diagnostics == T) diagnose_stan_fit(stan_fit, file_path, pars)

# WAIC 
WAIC = waic(stan_fit$draws("log_lik"))
saveRDS(WAIC, file = paste0(file_path, "/WAIC.RDS"))

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
