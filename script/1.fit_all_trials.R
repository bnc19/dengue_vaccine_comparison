  
# source functions 
rm(list = ls())
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
n_it = 6000

# M1 - start serotype hosp, no tau, no age, single n50 for MU and SN  ----------
run_model_TR (
  share_n_param = 2,
  uniform = 1, 
  mono_lc_MU = 1,
  mono_lc_SN = 1,
  tau_K = c(0,0,0),
  delta_KJ = c(1,1,1), 
  L_K = c(0,0,0), 
  include_beta = c(0,0,0),
  folder = "M1",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M2 - M1 but scale n50 for MU and SN  -----------------------------------------
run_model_TR (
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1, 
  tau_K = c(0,0,0),
  delta_KJ = c(1,1,1), 
  L_K = c(0,0,0), 
  include_beta = c(0,0,0),
  folder = "M2",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M3 - M2 but delta J for Dengvaxia  -----------------------------------------
run_model_TR (
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1, 
  tau_K = c(0,0,0),
  L_K = c(0,0,0), 
  include_beta = c(0,0,0),
  folder = "M3",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M4 - M3 but age VE for Qdenga  -----------------------------------------------
run_model_TR (
  include_beta = c(1,0,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  L_K = c(0,0,0), 
  folder = "M4",
  n_it = n_it,
  adapt_delta = 0.8,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)



# M5 - M4 but age VE for Qdenga only youngest  ---------------------------------
run_model_TR (
  include_beta = c(3,0,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  L_K = c(0,0,0), 
  folder = "M5",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M6 - M5 but age VE for Dengvaxia  --------------------------------------------
run_model_TR (
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  L_K = c(0,0,0), 
  folder = "M6",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M7 - M6 but age VE for Butantan-DV -------------------------------------------
run_model_TR (
  include_beta = c(3,1,1),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  L_K = c(0,0,0), 
  folder = "M7",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M8 - M6 but serotype specific L for Qdenga   ---------------------------------
run_model_TR (
  L_K = c(1,0,0), 
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M8",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M9 - M6 but serotype specific L for Dengvaxia   ---------------------------------
run_model_TR (
  L_K = c(0,1,0), 
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M9",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M10 - M9 but serotype specific L for But   ---------------------------------
run_model_TR (
  L_K = c(0,1,1), 
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M10",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M11 - M9 but serotype specific tau for Q -------------------------------------
run_model_TR (
  tau_K = c(1,0,0),
  L_K = c(0,1,0), 
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  folder = "M11",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M12 - M9 but serotype specific tau for D -------------------------------------
run_model_TR (
  tau_K = c(0,1,0),
  L_K = c(0,1,0), 
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  folder = "M12",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M13 - M9 but inc_FOIJ for B -------------------------------------------------

run_model_TR (  
  BF = T, 
  inc_FOIJ = c(0,0,1),
  L_K = c(0,1,0), 
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M13",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M14 - M13 but rho K ----------------------------------------------------------
run_model_TR (
  rho_K = 1, 
  inc_FOIJ = c(0,0,1),
  L_K = c(0,1,0), 
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M14",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M15 - M13 but single Butantan age VE ------------------------------------
run_model_TR ( 
  include_beta = c(3,1,2),
  inc_FOIJ = c(0,0,1),
  L_K = c(0,1,0), 
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M15",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M16 - M13 but serotype-specific L for Qdenga ------------------------------------

run_model_TR (  
  L_K = c(1,1,0), 
  BF = T, 
  inc_FOIJ = c(0,0,1),
  include_beta = c(3,1,0),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M16",
  n_it = n_it,
  adapt_delta = 0.78,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M17 - M13 but no Qdenga age effects -------------------------------------------------

run_model_TR (  
  include_beta = c(0,1,0),
  BF = T, 
  inc_FOIJ = c(0,0,1),
  L_K = c(0,1,0), 
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M17",
  n_it = n_it,
  adapt_delta = 0.78,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)

# M18 - M13 but no Dengvaxia age effects -------------------------------------------------

run_model_TR (  
  include_beta = c(3,0,0),
  BF = T, 
  inc_FOIJ = c(0,0,1),
  L_K = c(0,1,0), 
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  tau_K = c(0,0,0),
  folder = "M18",
  n_it = n_it,
  adapt_delta = 0.78,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M19 - M18 but single Dengvaxia L and serotype tau --------------------------------------

run_model_TR ( 
  L_K = c(0,0,0), 
  tau_K = c(0,1,0),
  include_beta = c(3,0,0),
  inc_FOIJ = c(0,0,1),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  folder = "M19",
  n_it = n_it,
  adapt_delta = 0.78,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M20 - M19 but Qdenga serotype tau --------------------------------------

run_model_TR ( 
  tau_K = c(1,1,0),
  L_K = c(0,0,0), 
  include_beta = c(3,0,0),
  inc_FOIJ = c(0,0,1),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  folder = "M20",
  n_it = n_it,
  adapt_delta = 0.78,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# M21 - M19 but single Dengvaxia tau --------------------------------------

run_model_TR ( 
  BF = T, 
  L_K = c(0,0,0), 
  tau_K = c(0,0,0),
  include_beta = c(3,0,0),
  inc_FOIJ = c(0,0,1),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  folder = "M21",
  n_it = n_it,
  adapt_delta = 0.78,
  stan_model = "final_model.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


# MS1 - M21 but severe data and age-specific psi  ------------------------------
run_model_TR (
  psi_J = 1, 
  L_K = c(0,0,0), 
  tau_K = c(0,0,0),
  include_beta = c(3,0,0),
  inc_FOIJ = c(0,0,1),
  delta_KJ = c(1,2,1), 
  mono_lc_MU = 2,
  mono_lc_SN = 2,
  share_n_param = 2,
  uniform = 1,
  folder = "MS1",
  n_it = n_it,
  adapt_delta = 0.81,
  stan_model = "final_model_severe.stan",
  mu_Q  = array(read.csv("data/processed_Q/n0_new.csv")$mean, dim = c(2, 4)),
  baseline_SP_Q = read.csv("data/processed_Q/seropositive_by_age_baseline.csv"),
  VCD_Q = read.csv("data/processed_Q/vcd_data.csv"),
  hosp_Q = read.csv("data/processed_Q/hosp_data.csv") ,
  baseline_SP_De = read.csv("data/processed_De/baseline_SP.csv"),
  cases_De = readRDS("data/processed_De/cases_stan_format.RDS"),
  mu_De = read.csv("data/processed_De/mu.csv"),
  cases_Bu  = readRDS("data/processed_Bu/cases_stan_format.RDS"),
  baseline_SP_Bu = readRDS("data/processed_Bu/baseline_SP.RDS"),
  mu_Bu = readRDS("data/processed_Bu/mu.RDS")
)


