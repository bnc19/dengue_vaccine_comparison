# M2 year 2  - M2 but only fit to the first 2 years of symptomatic data to match Butantan-DV 

VCD_years = c(12, 18, 24) / 12
time =1:24
mono_lc_MU = 2
include_pK3 = 0
rho_K = 0
include_beta = 3
mono_lc_SN = 1
tau_K = 1
MU_test_SN = 1
folder = "M2_year2"
n_it = 10000
n_chains = 4
adapt_delta = 0.77
stan_model = "TAK_model_year2.stan"
baseline_SP = read.csv("TAK/data/raw/seropositive_by_age_baseline.csv")
VCD =  read.csv("TAK/data/processed/vcd_data.csv")
hosp = read.csv("TAK/data/processed/hosp_data.csv") 
mu =  array(read.csv("TAK/data/processed/n0_new.csv")$mean, dim = c(2, 4))

init_values = function(){
  list(
    hs = runif(2, 2,4),
    hl = runif(1, 108,132),
    ts = runif(2, -2,2),
    lambda_D = replicate(3, runif(4,0,0.02)),
    p = runif(3,0.1,0.9), 
    pK3 = runif(4,0.1,0.9),
    gamma = runif(1,0.1,0.9),
    rho = runif(4,1,5), 
    phi = runif(1,0.1,0.5),
    delta = runif(4,0.1,0.8),
    L = runif(4,0,4),
    w = runif(4,1,4),
    alpha = runif(4,0,2),
    beta = runif(2,0,2),
    tau = runif(4,1,5),
    lc = replicate(4,runif(3,2,8)),
    sens = runif(1,.8,1),
    spec = runif(1, 0.95,1),
    epsilon = runif(1,1,4)
  )}

B = 2
K = 4
V = 2
R = 2
J = 3
C = 3
HI = 12
include_eps = 0
rho_K = 0
L_K = 0
w_CK = 0
alpha_CK = 0
L_sd = 1
L_mean =0
lower_bound_L = 0
MU_test_SN = 1
MU_symp = 1
enhancement = 1

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
file.sources = paste0("TAK/R/", list.files(path = "TAK/R/"))
sapply(file.sources, source)

comp_model = cmdstan_model(paste0("TAK/models/", stan_model),stanc_options = list("O1"))
file_path = (paste0("TAK/output/", folder))
dir.create(file_path)

# data -----------------------------------------------------------------------
hosp_f = factor_TAK_VCD(hosp)
VCD_f = factor_TAK_VCD(VCD)  

VCD2 = filter(VCD_f, month <= 24)
HOSP2 = filter(hosp_f, month <= 24)

#  Stan data 
T = length(time)
D = length(VCD_years)

  
# calculate # VCD over time
VCD_D = VCD2 %>%
    filter(serostatus == "both",
           age == "all",
           serotype == "all",
           year %in% VCD_years) %>%
    group_by(year) %>%
    summarise(Y = sum(Y)) %>%
    select(Y)
  
# calculate # hosp over time

hosp_D = HOSP2 %>%
  filter(age != "all", trial != "both") %>%
  group_by(year) %>%
  summarise(Y = sum(Y)) %>%
  select(Y)
  
# calculate # pop by trial arm, serostatus, age group, over time
  N_pop_BVJD = VCD2 %>%
    filter(serostatus != "both",
           age != "all",
           serotype == "all",
           trial != "both") %>%
    arrange(year, age, trial , serostatus)
  
pop_BVJD = array(N_pop_BVJD$N, dim = c(B, V, J, D))
  
# calculate # hosp by serostatus + trial arm + serotype, over time (36-54 months)
N_hosp_BVK = HOSP2 %>%
    filter(serotype != "all", trial != "both") %>%
    arrange(year, serostatus, trial , serotype)
  
hosp_BVK_m = array(N_hosp_BVK$Y, dim = c(B * K * V))  # (D=4)
  
# calculate # VCD by serostatus + trial arm + serotype, over time
# for multinomial likelihood

N_VCD_BVKD_m = VCD2 %>%
    group_by(serostatus, trial, age, serotype, year) %>%
    summarise(Y = sum(Y)) %>%
    filter(serostatus != "both",
           age == "all",
           serotype != "all",
           trial != "both",
           year %in% VCD_years) %>%
    arrange(year, serostatus, trial , serotype)
  
VCD_BVKD_m = array(N_VCD_BVKD_m$Y, dim = c(B * K * V, D))
  
# calculate # hosp by serostatus + trial arm + age-group, over time
 N_hosp_BVJA_m =  HOSP2 %>%
    filter(serotype == "all", trial != "both", age != "all") %>%
    arrange(year, serostatus, trial , age)
  
hosp_BVJA_m = array(N_hosp_BVJA_m$Y, dim = c(B * V * J, D))
  
# calculate # VCD by serostatus + trial arm + age-group, over time
# for multinomial likelihood

N_VCD_BVJA_m = VCD2 %>%
    group_by(serostatus, trial, age, serotype, year) %>%
    summarise(Y = sum(Y)) %>%
    filter(serostatus != "both",
           age != "all",
           serotype == "all",
           trial != "both") %>%
    arrange(year, serostatus, trial , age)
  
  VCD_BVJA_m = array(N_VCD_BVJA_m$Y, dim = c(B * V * J, D))
  
# VCD age and serotype year and 2
  N_VCD_KJ2_m = VCD2 %>%
    filter(age != "all", serotype != "all") %>%
    arrange(year, serotype, age)
  
VCD_KJ2_m = array(N_VCD_KJ2_m$Y, dim = c(K * J, 2))
  
# hosp age and serotype year and 2
  N_hosp_KJ2_m = HOSP2 %>%
    filter(age != "all", serotype != "all") %>%
    arrange(year, serotype, age)
  
hosp_KJ2_m = array(N_hosp_KJ2_m$Y, dim = c(K * J, 2))

# data
  
  stan_data = list(
    time = time,
    T = T,
    J = J,
    K = K,
    B = B,
    V = V,
    D = D,
    C = C,
    HI = HI,
    R = R,
    include_pK3 = include_pK3,
    include_eps = include_eps,
    include_beta = include_beta,
    mono_lc_SN = mono_lc_SN,
    mono_lc_MU = mono_lc_MU,
    rho_K = rho_K,
    L_K = L_K,
    w_CK = w_CK,
    alpha_CK = alpha_CK,
    tau_K = tau_K,
    L_sd = L_sd,
    L_mean = L_mean,
    lower_bound_L = lower_bound_L,
    MU_test_SN = MU_test_SN,
    MU_symp = MU_symp,
    enhancement = enhancement, 
    SP_J = baseline_SP$Seropositive,
    pop_J = baseline_SP$Total,
    VCD_D = as.numeric(VCD_D$Y),
    pop = pop_BVJD,
    VCD_BVKD = VCD_BVKD_m,
    VCD_BVJD = VCD_BVJA_m,
    mu = mu,
    HOSP_BVJD = hosp_BVJA_m,
    HOSP_BVK = hosp_BVK_m,
    HOSP_D = as.numeric(hosp_D$Y),
    HOSP_KJ2 = hosp_KJ2_m,
    VCD_KJ2 = VCD_KJ2_m
  )

  # fit model ------------------------------------------------------------------
start = Sys.time()
stan_fit = comp_model$sample(
  data = stan_data,
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
i1 = which(names(fit_ext) == "ll")
names_select = names(fit_ext)[1:i1]
posterior_chains = fit_ext[names_select]
write.csv(posterior_chains, paste0(file_path, "/posterior_chains.csv"))

posterior = summarise_draws(posterior_chains)
write.csv(posterior, paste0(file_path, "/posterior.csv"))


# plot fit ---------------------------------------------------------------------  
AR = which(grepl("AR" , names(fit_ext)))
AR_out = fit_ext[AR] %>%  as.data.frame()
saveRDS(AR_out, paste0(file_path, "/AR.RDS"))

VE = which(grepl("VE" , names(fit_ext)))
VE_out = fit_ext[VE] %>%  as.data.frame()
saveRDS(VE_out, paste0(file_path, "/VE.RDS"))

n = which(grepl("n" , names(fit_ext)))
n_out = fit_ext[n] %>%  as.data.frame()
saveRDS(n_out, paste0(file_path, "/n.RDS"))

plot_TAK_output(file_path = file_path)


}
