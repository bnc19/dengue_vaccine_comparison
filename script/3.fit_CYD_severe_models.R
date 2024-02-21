
# source functions 
rm(list = ls())
file.sources = paste0("CYD/R/", list.files(path = "CYD/R/"))
sapply(file.sources, source)
n_it = 10000

# M0 - start with best fitting CYD-TDV model without severe data ---------------
run_CYD_model (
  mono_lc_SN = 1,
  L_K = 1, 
  tau_K = 0,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M0",
  n_it = n_it,
  adapt_delta = 0.8,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)

# M1 - M1 but single L ---------------------------------------------------------
run_CYD_model (
  L_K = 0, 
  mono_lc_SN = 1,
  tau_K = 0,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M1",
  n_it = n_it,
  adapt_delta = 0.99,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)

# M2 - M1 but age-group offset for RR_severe -----------------------------------
run_CYD_model (
  include_beta = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M2",
  n_it = n_it,
   adapt_delta = 0.99,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)

# M3 - M1 but age-specific psi -------------------------------------------------
run_CYD_model (
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M3",
  n_it = n_it,
   adapt_delta = 0.99,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)

# M4 - M3 but include epsilon --------------------------------------------------
run_CYD_model (
  include_eps = 1,
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M4",
  n_it = n_it,
   adapt_delta = 0.99,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)


# M5 - M4 but include single beta ----------------------------------------------
run_CYD_model (
  include_beta = 1,
  include_eps = 1,
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M5",
  n_it = n_it,
   adapt_delta = 0.99,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)


# M6 - M5 but include outcome specific beta ------------------------------------
run_CYD_model (
  include_beta = 2,
  include_eps = 1,
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M6",
  n_it = n_it,
  adapt_delta = 0.99,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)


# M7 - M5 but no epsilon -------------------------------------------------------
run_CYD_model (
  include_eps = 0,
  include_beta = 1,
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M7",
  n_it = n_it,
   adapt_delta = 0.99,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)


# M8 - M6 but no epsilon -------------------------------------------------------
run_CYD_model (
  include_beta = 2,
  include_eps = 0,
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M8",
  n_it = n_it,
  adapt_delta = 0.9,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)


# M9 - M6 but serotype specific L ----------------------------------------------
run_CYD_model (
  L_K = 1, 
  include_beta = 2,
  include_eps = 1,
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M9",
  n_it = n_it,
  adapt_delta = 0.9,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)




# M10 - M6 but serotype specific tau ----------------------------------------------
run_CYD_model (
  tau_K = 1, 
  include_beta = 2,
  include_eps = 1,
  psi_J = 1,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M10",
  n_it = n_it,
  adapt_delta = 0.9,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)

