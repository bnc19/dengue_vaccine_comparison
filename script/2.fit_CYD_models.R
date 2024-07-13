  
  # source functions 
  rm(list = ls())
  file.sources = paste0("CYD/R/", list.files(path = "CYD/R/"))
  sapply(file.sources, source)
  n_it = 10000

# M1 - start with best fitting TAK-003 model -----------------------------------
run_CYD_model (
  include_pK2 = 0,
  include_beta = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  delta_KJ = 1, 
  folder = "M1",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M2 - M1 but serotype-specific p ----------------------------------------------
run_CYD_model (
  include_pK2 = 1,
  include_beta = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  delta_KJ = 1, 
  folder = "M2",
  n_it = n_it,
  adapt_delta =0.77,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M3 - M1 but single tau -------------------------------------------------------
run_CYD_model (
  tau_K = 0,
  include_pK2 = 0,
  include_beta = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  MU_test_SN = 1,
  delta_KJ = 1, 
  folder = "M3",
  n_it = n_it,
  adapt_delta =0.77,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M4 - M1 but single delta -----------------------------------------------------
run_CYD_model (
  delta_KJ = 0, 
  include_pK2 = 0,
  include_beta = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M4",
  n_it = n_it,
  adapt_delta =0.99,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M5 - M1 but age delta -----------------------------------------------------
run_CYD_model (
  delta_KJ = 2, 
  include_pK2 = 0,
  include_beta = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M5",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M6 - M5 but monotypic tau  ---------------------------------------------------
run_CYD_model (
  tau_K = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  include_beta = 1,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  MU_test_SN = 1,
  folder = "M6",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)


# M7 - M5 but no beta ----------------------------------------------------------
run_CYD_model (
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M7",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M8 - M7 but lc MU is offset from lc MO  --------------------------------------
run_CYD_model (
  mono_lc_MU = 2,
  mono_lc_SN = 1,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M8",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M9 - M8 but lc SN is offset from lc MO  --------------------------------------
run_CYD_model (
  mono_lc_SN = 2,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "M9",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M10 - M9 but mono tau   ------------------------------------------------------
run_CYD_model (
  tau_K = 0,
  mono_lc_SN = 2,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "M10",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)


# M11 - M10 but serotype-specific L  -------------------------------------------
run_CYD_model (
  L_K = 1, 
  tau_K = 0,
  mono_lc_SN = 2,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "M11",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)

# M12 - M11 but mono lc50 SN   -------------------------------------------------
run_CYD_model (
  mono_lc_SN = 1,
  L_K = 1, 
  tau_K = 0,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "M12",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  BF = T
)

# M13 - M12 but single L -------------------------------------------------
run_CYD_model (
  L_K = 0, 
  mono_lc_SN = 1,
  tau_K = 0,
  mono_lc_MU = 2,
  include_beta = 0,
  delta_KJ = 2, 
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "M13",
  n_it = n_it,
  adapt_delta =0.8,
  stan_model = "CYD_model.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv")
)


# SA 
  
#  SA1 - M12 but end_time = 60 (minimum)
  run_CYD_model (
    end_time= 60,
    mono_lc_SN = 1,
    L_K = 1,
    tau_K = 0,
    mono_lc_MU = 2,
    include_beta = 0,
    delta_KJ = 2,
    include_pK2 = 0,
    MU_test_SN = 1,
    folder = "SA1",
    n_it = n_it,
    adapt_delta = 0.8,
    stan_model = "CYD_model.stan",
    baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
    cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
    mu =  read.csv("CYD/data/processed/mu.csv"),
    BF = T
  )
  