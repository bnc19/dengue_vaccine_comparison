run_BUT_model = function(n_it = 4000,
                     adapt_delta = 0.8,
                     init_values = function(){
                       list(
                         hs = runif(2, 2,4),
                         hl = runif(1, 108,132),
                         ts = runif(2, -2,2),
                         lambda_K = runif(2,0,0.02),
                         p = runif(3,0.1,0.9), 
                         pK3 = runif(4,0.1,0.9),
                         gamma = runif(1,0.1,0.9),
                         rho = runif(2,1,5), 
                         phi = runif(1,0.1,0.5),
                         L = runif(2,0,4),
                         w = runif(3,1,4),
                         beta = runif(2,0,2),
                         lc = replicate(2,runif(3,2,8)),
                         sens = runif(1,.8,1),
                         spec = runif(1, 0.95,1))},
                     start_time = 1,
                     end_time = 24, 
                     n_chains = 4,
                     stan_model,
                     folder,
                     pars = NULL,
                     BF=F,
                     B = 2,
                     K = 2,
                     V = 2,
                     J = 3,
                     C = 3,
                     M = 4,
                     HI = 12,
                     include_pK3=0,
                     include_beta=0,
                     mono_lc_SN=0,
                     mono_lc_MU=0,
                     rho_K=0,
                     L_K=0,
                     w_CK=0,
                     L_mean=0,
                     L_sd=1,
                     lower_bound_L=0,
                     MU_test_SN=0,
                     MU_symp=1,
                     enhancement=1,
                     sens = 0.91,
                     spec = 0.99,
                     ts = c(-2.12,0.31),
                     hs = c(1.9,4.3),
                     hl = 73.5,
                     diagnostics = F,
                     serotype_serostatus_cases,
                     age_cases,
                     baseline_seropos,
                     titres,
                     pop, 
                     serostatus = c("seronegative", "monotypic",
                                    "multitypic")
) {
# set up -----------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(loo)
library(readxl)
options(mc.cores = parallel::detectCores())
  
# source functions
file.sources = paste0("BUT/R/", list.files(path = "BUT/R/"))
sapply(file.sources, source)
  
comp_model = cmdstan_model(paste0("BUT/models/", stan_model),stanc_options = list("O1"))
  
file_path = (paste0("BUT/output/", folder))
dir.create(file_path)
# data -----------------------------------------------------------------------
serotype_serostatus_cases = factor_BUT_cases(serotype_serostatus_cases)
age_cases = factor_BUT_cases(age_cases)
titres = factor_BUT_cases(titres)
pop = factor_BUT_cases(pop)

mu = titres %>% 
  filter(Serotype == "DENV1" | Serotype == "DENV2") %>% 
  arrange(Serostatus, Serotype) %>% 
  pivot_wider(names_from = Serotype, values_from = "Titre") %>%  
  select(DENV1, DENV2) %>% 
  as.matrix()

VCD_BVK = serotype_serostatus_cases %>% 
  arrange(Serostatus, Arm, Serotype) %>%  
  pull(Y)

VCD_VJ = age_cases %>% 
  arrange(Arm, Age) %>%  
  pull(Y)

pop_BV = serotype_serostatus_cases %>% 
  filter(Serotype == "DENV2") %>%  
  select(Serostatus, Arm, N) %>% 
  arrange(Serostatus, Arm) %>% 
  pivot_wider(names_from = Arm, values_from = N) %>%  
  select(placebo,vaccine) %>%  
  as.matrix()
  
pop_VJ = age_cases %>% 
  select(Age, Arm, N) %>% 
  arrange(Arm, Age) %>% 
  pivot_wider(names_from = Age, values_from = N) %>% 
  select(-Arm) %>%  
  as.matrix()

# fit Stan model ---------------------------------------------------------------
time = start_time : end_time
T = length(time)
  list_data = list(
    M = M,
    B = B,
    K = K,
    V = V,
    J = J,
    C = C,
    T=T,
    HI = HI,
    time = time,
    mu = mu,
    sens = sens,
    spec = spec,
    ts = ts,
    hs = hs,
    hl = hl,
    include_pK3=include_pK3,
    include_beta=include_beta,
    mono_lc_SN=mono_lc_SN,
    mono_lc_MU=mono_lc_MU,
    rho_K=rho_K,
    L_K=L_K,
    w_CK=w_CK,
    L_mean=L_mean,
    L_sd=L_sd,
    lower_bound_L=lower_bound_L,
    MU_test_SN=MU_test_SN,
    MU_symp=MU_symp,
    enhancement=enhancement,
    SP = sum(baseline_seropos$SP), 
    pop_J = pull(pop, N), 
    VCD = sum(age_cases$Y), 
    VCD_BVK = VCD_BVK, 
    VCD_VJ = VCD_VJ,
    pop_BV = pop_BV,
    pop_VJ = pop_VJ, 
    pop = sum(pop_VJ)) 

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
  
i1 = which(names(fit_ext) == "ll")
names_select = names(fit_ext)[1:i1]
posterior_chains = fit_ext[names_select]
write.csv(posterior_chains,paste0(file_path, "/posterior_chains.csv") )
  
posterior = summarise_draws(posterior_chains)
write.csv(posterior, paste0(file_path, "/posterior.csv"))
  
  
if(BF== T) calculate_BUT_bayes(stan_fit =stan_fit, L_mean =L_mean, L_sd=L_sd,
                           file_path=file_path, lower= (lower_bound_L-1))

# WAIC 
WAIC = waic(stan_fit$draws("log_lik"))
saveRDS(WAIC, file = paste0(file_path, "/WAIC.RDS"))

  
 # plot fit --------------------------------------------------------------------
AR = which(grepl( "AR" , names( fit_ext ) ))
AR_out = fit_ext[AR] %>%  as.data.frame()
saveRDS(AR_out, paste0(file_path, "/AR.RDS"))
  
VE = which(grepl( "VE" , names( fit_ext ) ))
VE_out = fit_ext[VE]%>%  as.data.frame()
saveRDS(VE_out, paste0(file_path, "/VE.RDS"))
  
n = which(grepl( "n" , names( fit_ext ) ))
n_out = fit_ext[n]%>%  as.data.frame()
saveRDS(n_out, paste0(file_path, "/n.RDS"))

plot_BUT_output(file_path =file_path,
            serotype_serostatus_cases=serotype_serostatus_cases,
            age_cases=age_cases,
            include_beta = include_beta)
  

}