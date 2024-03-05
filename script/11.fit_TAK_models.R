# Script to run TAK models starting with best fitting from impact modelling paper
rm(list=ls())

# source functions
file.sources = paste0("TAK/R/", list.files(path = "TAK/R/"))
sapply(file.sources, source)
n_it = 10000

# M1 - best fitting TAK model --------------------------------------------------
run_TAK_model (
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M1",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("TAK/data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("TAK/data/processed/vcd_data.csv"),
  hosp = read.csv("TAK/data/processed/hosp_data.csv") ,
  mu =  array(read.csv("TAK/data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M2 - as M1 MU offset from MO -------------------------------------------------
run_TAK_model (
  BF = T,
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M2",
  n_it = n_it,
  adapt_delta = 0.75,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("TAK/data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("TAK/data/processed/vcd_data.csv"),
  hosp = read.csv("TAK/data/processed/hosp_data.csv") ,
  mu =  array(read.csv("TAK/data/processed/n0_new.csv")$mean, dim = c(2, 4))
)


# M3 - as M1 but SN offset from MO ---------------------------------------------
run_TAK_model (
  mono_lc_SN = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M3",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("TAK/data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("TAK/data/processed/vcd_data.csv"),
  hosp = read.csv("TAK/data/processed/hosp_data.csv") ,
  mu =  array(read.csv("TAK/data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# M4 - as M3 but MU offset from MO ---------------------------------------------
run_TAK_model (
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M4",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("TAK/data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("TAK/data/processed/vcd_data.csv"),
  hosp = read.csv("TAK/data/processed/hosp_data.csv") ,
  mu =  array(read.csv("TAK/data/processed/n0_new.csv")$mean, dim = c(2, 4))
)


