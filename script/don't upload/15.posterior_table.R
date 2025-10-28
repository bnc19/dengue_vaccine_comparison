# Format CYD distribution for table in manuscript 
library(tidyverse)

# source data 
model = "M21"
posterior = read.csv(paste0("output/", model, "/posterior.csv"))
posterior_chains = read.csv(paste0("output/", model, "/posterior_chains.csv"))

# format main param exlusing lambda 
parameters = posterior %>% 
  select(variable, mean, q5, q95) %>%  
  filter(grepl("p\\[", variable) |
         grepl("delta\\[1", variable) |
         grepl("delta\\[2", variable) |
         grepl("sens", variable)|
         grepl("spec", variable)|
         grepl("omega", variable)|
         grepl("kappa", variable)|
         variable %in% c("L[1,1]", "L[2,1]", "L[3,1]", 
                         "w[1,1]", "w[2,1]", "w[3,1]", 
                         "lc[1,2,1]", "lc[2,2,1]", "lc[3,2,1]",
                         "lc[1,2,2]", "lc[2,2,2]", "lc[3,2,2]",
                         "lc[1,2,3]", "lc[2,2,3]", 
                         "lc[1,2,4]", "lc[2,2,4]", 
                         "tau[1,1]",  "tau[2,1]",
                         "alpha[1,1]", "alpha[2,1]",
                         "beta[1,1]",
                         "gamma", "rho[1]", "phi",
                         "hs[1,1]", "hs[1,2]", "hl[1]", "ts[1,1]",  "ts[1,2]"  ,  
                         "hs[2,1]", "hs[2,2]", "hl[2]", "ts[2,1]",  "ts[2,2]"  
                         )
         ) %>% 
  filter(variable != "delta[2,4]", variable != "delta[2,3]")  %>% 
  mutate_if(is.numeric, ~sprintf("%.2f", .)) %>%
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2", "dep3")) %>% 
  arrange(variable, dep) %>% 
  pivot_wider(names_from = dep, values_from = X)
         
write.csv(parameters, paste0("output/", model, "/posterior_formatted.csv"))  

# 61 parameters 

# foi 

lamba1 = posterior %>% select(variable, mean, q5, q95) %>% 
  filter(grepl("theta", variable)| variable == "FOI_J1[3]") %>%  
  mutate_if(is.numeric, round, 2) %>% 
  mutate_if(is.numeric, as.character)
  
  
  
lambda = posterior %>%
  select(variable, mean, q5, q95) %>%  
  filter(grepl("lambda", variable))%>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2) %>% 
  bind_rows(lamba1) %>% 
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2", "dep3")) 

# 35 parameters 

write.csv(lambda, paste0("output/", model, "/lambda_formatted.csv"))


# get rho etc. 
posterior_chains %>% 
  as.data.frame() %>% 
  select("rho.1.", "gamma", "phi") %>% 
  mutate(gamma = rho.1. * gamma) %>% 
  mutate(phi = gamma * phi) %>% 
  pivot_longer(everything()) %>% 
  group_by(name) %>% 
  mutate(value = value * 100) %>% 
  summarise(mean = mean(value),
            lower = quantile(value, p = 0.025),
            upper = quantile(value, p = 0.975)) %>% 
  mutate_if(is.numeric, ~sprintf("%.0f", .))
  

# compare FOI

lambda_chains = posterior_chains %>% 
  as.data.frame() %>% 
  select(contains("lambda"))

lambda_chains %>% 
  select(contains("lambda_K_Bu")) %>% 
  mutate(foi =( lambda_K_Bu.1. + lambda_K_Bu.2.) * 12) %>% # monthly foi 
  summarise(mean(foi))

Q_weights = data.frame(time = as.character(1:6), 
                       W = c(12,6,6,12,12,6))

lambda_chains %>% 
  select(contains("lambda_D_Q")) %>% 
  rownames_to_column(var = "ni") %>% 
  pivot_longer(-ni) %>%  
  separate(name, into = c("var", "serotype", "time"), sep = "\\.") %>% 
  group_by(time, ni) %>%  
  summarise(cum_FOI = sum(value)) %>%  # cumulative FOI 
  left_join(Q_weights) %>% 
  mutate(w_mean = cum_FOI * W) %>%  # monthly FOI multiplied by number of months 
  group_by(ni) %>% 
  summarise(FOI = sum(w_mean)) %>% # sum over t_d,
  mutate(FOI = round(FOI / 4.5, 2 )) %>% 
  ungroup() %>% 
  summarise(mean(FOI))
  
 
D_weights = data.frame(time = as.character(1:4), 
                       W = c(length(1:13), length(14:24) ,length(25:36), length(37:60)))

lambda_chains %>% 
  select(contains("lambda_D_De")) %>% 
  rownames_to_column(var = "ni") %>% 
  pivot_longer(-ni) %>%  
  separate(name, into = c("var", "time"), sep = "\\.") %>% 
  left_join(D_weights) %>% 
  mutate(w_mean = value * W) %>%  # monthly FOI multiplied by number of months 
  group_by(ni) %>% 
  summarise(FOI = sum(w_mean)) %>% # sum over t_d,
  mutate(FOI = round(FOI / 5, 2 )) %>% 
  ungroup() %>% 
  summarise(mean(FOI))

 