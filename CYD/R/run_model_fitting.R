run_CYD_model = function(n_it = 10000,
                         adapt_delta = 0.8,
                         init_values = function() {
                           list(
                             hs = runif(2, 2, 30),
                             hl = runif(1, 60, 132),
                             ts = runif(1, -12, 12),
                             lambda_D = runif(4, 0, 0.02),
                             p = runif(2, 0.1, 0.9),
                             pK3 = runif(4, 0.1, 0.9),
                             gamma = runif(1, 0.1, 0.9),
                             rho = runif(4, 1, 5),
                             phi = runif(1, 0.1, 0.5),
                             delta = runif(4, 0.1, 0.8),
                             L = runif(4, 0, 4),
                             w = runif(4, 1, 4),
                             alpha = runif(4, 0, 2),
                             beta = runif(2, 0, 2),
                             tau = runif(4, 1, 5),
                             lc = replicate(4, runif(3, 2, 8)),
                             sens = runif(1, .6, 1),
                             spec = runif(1, .7, 1),
                             epsilon = runif(1, 1, 4)
                           )
                         },
                         n_chains = 4,
                         stan_model,
                         folder,
                         pars = NULL,
                         VCD_years = c(13, 24, 36, 60),
                         start_time = 13,
                         end_time = 72,
                         BF = F,
                         B = 2,
                         K = 4,
                         V = 2,
                         R = 2,
                         J = 2,
                         C = 3,
                         HI = 12,
                         include_pK2 = 0,
                         include_eps = 0,
                         include_beta = 0,
                         mono_lc_SN = 0,
                         mono_lc_MU = 0,
                         rho_K = 0,
                         L_K = 0,
                         w_CK = 0,
                         alpha_CK = 0,
                         tau_K = 0,
                         delta_KJ = 0,
                         psi_J = 0,
                         FOI_J = 0,
                         L_mean = 0,
                         L_sd = 1,
                         lower_bound_L = 0,
                         MU_test_SN = 0,
                         MU_symp = 1,
                         enhancement = 1,
                         mu,
                         baseline_SP,
                         severe = F,
                         cases
) {

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
file.sources = paste0("CYD/R/", list.files(path = "CYD/R/"))
sapply(file.sources, source)
  
comp_model = cmdstan_model(paste0("CYD/models/", stan_model),stanc_options = list("O1"))

file_path = (paste0("CYD/output/", folder))
dir.create(file_path)

# fit Stan model -------------------------------------------------------------

time = start_time : end_time
time = time - time[1] + 1 # start at 1

list_data = format_CYD_stan_data(
  baseline_SP = baseline_SP,
  cases = cases,
  B = B,
  R = R,
  K = K,
  V = V,
  J = J,
  C = C,
  HI = HI,
  VCD_years = VCD_years,
  time = time,
  mu = mu,
  include_pK2 = include_pK2,
  include_eps = include_eps,
  include_beta = include_beta,
  mono_lc_SN = mono_lc_SN,
  mono_lc_MU = mono_lc_MU,
  rho_K = rho_K,
  L_K = L_K,
  w_CK = w_CK,
  alpha_CK = alpha_CK,
  tau_K = tau_K,
  delta_KJ = delta_KJ,
  psi_J = psi_J,
  FOI_J = FOI_J,
  L_mean = L_mean,
  L_sd = L_sd,
  lower_bound_L = lower_bound_L,
  MU_test_SN = MU_test_SN,
  MU_symp = MU_symp,
  enhancement = enhancement
) 

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

# WAIC 
WAIC = waic(stan_fit$draws("log_lik"))
saveRDS(WAIC, file = paste0(file_path, "/WAIC.RDS"))

# save posts -------------------------------------------------------------------
  
  fit_ext = stan_fit$draws(format = "df")
   
  # save up to log likelihood 
  i1 = which(names(fit_ext) == "ll")
  names_select = names(fit_ext)[1:i1]
  posterior_chains = fit_ext[names_select]
  write.csv(posterior_chains,paste0(file_path, "/posterior_chains.csv") )
  
  posterior = summarise_draws(posterior_chains)
  write.csv(posterior, paste0(file_path, "/posterior.csv"))
  
  if (BF == T)
    calculate_CYD_bayes(
      stan_fit = stan_fit,
      L_mean = L_mean,
      L_sd = L_sd,
      file_path = file_path,
      lower = lower_bound_L
    )
  
# plot fit ---------------------------------------------------------------------  

  AR = which(grepl( "AR" , names( fit_ext ) ))
  AR_out = fit_ext[AR] %>%  as.data.frame()
  saveRDS(AR_out, paste0(file_path, "/AR.RDS"))
  
  VE = which(grepl( "VE" , names( fit_ext ) ))
  VE_out = fit_ext[VE]%>%  as.data.frame()
  saveRDS(VE_out, paste0(file_path, "/VE.RDS"))
  
  n = which(grepl( "n" , names( fit_ext ) ))
  n_out = fit_ext[n]%>%  as.data.frame()
  saveRDS(n_out, paste0(file_path, "/n.RDS"))

  plot_CYD_output(file_path = file_path,
                  severe = severe)
  
}
