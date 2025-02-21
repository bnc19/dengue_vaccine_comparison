  
# source functions 
rm(list = ls())
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
n_it = 4000

# M1 - start serotype hosp, no tau, no age, single n50 for MU and SN  ----------
run_model_Q_DE (
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 0,
  delta_KJ = c(1,1), 
  L_K = c(0,0), 
  include_beta = c(0,0),
  folder = "M1",
  n_it = n_it,
  adapt_delta = 0.75,
  stan_model = "D_Q_model_current.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
)
