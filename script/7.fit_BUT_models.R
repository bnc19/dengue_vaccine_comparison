# Script to fit the model to data from the BUTANTAN-DV phase III clinical trial 

rm(list=ls())

# source functions 
file.sources = paste0("BUT/R/", list.files(path = "BUT/R/"))
sapply(file.sources, source)
n_it = 10000

# M1 - start with best fitting TAK-003 model (without delta and tau) -----------
# Fixed antibody parameters, sens and spec 
run_BUT_model (
  include_beta = 2,  # youngest only 
  rho_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M1",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

# M2 - M1 but no enhancement parameter L  --------------------------------------
run_BUT_model (
  enhancement = 0, 
  include_beta = 2,  # youngest only 
  rho_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M2",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

# M3 - M2 but no beta  ---------------------------------------------------------
run_BUT_model (
  enhancement = 0, 
  include_beta = 0,  
  rho_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 1,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M3",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

# M4 - M3 but MU lc50 offset from MO -------------------------------------------
run_BUT_model (
  mono_lc_MU = 2,
  mono_lc_SN = 1,
  enhancement = 0, 
  include_beta = 0,  
  rho_K = 0,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M4",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

# M5 - M4 but SN lc50 offset from MO -------------------------------------------
run_BUT_model (
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  enhancement = 0, 
  include_beta = 0,  
  rho_K = 0,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M5",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

# M6 - M4 but include beta youngest only  ----------------------------------------
run_BUT_model (
  include_beta = 2, 
  mono_lc_MU = 2,
  mono_lc_SN = 1,
  enhancement = 0, 
  rho_K = 0,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M6",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

# M7 - M4 but age specific FOI -------------------------------------------------
run_BUT_model (
  inc_FOIJ = 1, 
  mono_lc_MU = 2,
  mono_lc_SN = 1,
  enhancement = 0, 
  include_beta = 0,  
  rho_K = 0,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M7",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

# M8 - M7 but single lc50 value ------------------------------------------------
run_BUT_model (
  single_lc = 1, 
  inc_FOIJ = 1, 
  enhancement = 0, 
  include_beta = 0,  
  rho_K = 0,
  MU_test_SN = 1,
  stan_model = "BUT_model.stan",
  folder = "M8",
  n_it = n_it,
  adapt_delta =0.75,
  sens = 0.91,
  spec = 0.99,
  ts = c(-2.1,0.31),
  hs = c(1.93,4.34),
  hl = 72.14,
  baseline_seropos = readRDS("BUT/data/processed/baseline_SP.RDS"),
  cases = readRDS("BUT/data/processed/cases_stan_format.RDS"),
  titres = readRDS("BUT/data/processed/mu.RDS")
)

