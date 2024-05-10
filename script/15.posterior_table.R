# Format CYD distribution for table in manuscript 
library(tidyverse)

# source data 

posterior_CYD = read.csv("CYD/output/M12/posterior.csv")

tidy_post_CYD = posterior_CYD %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__','pK2[1]','pK2[2]','pK2[3]','pK2[4]',
                          "delta[3]", "delta[4]", 
                          "epsilon", 'rho[2]','rho[3]','rho[4]',
                          'w[2]','w[3]','w[4]', "tau[2]", "tau[3]", "tau[4]",
                          "lc[1,2]", "lc[1,3]", "lc[1,4]",
                          "lc[3,1]", "lc[3,2]", "lc[3,3]","lc[3,4]",
                          "kappa",  "beta[1]", "beta[2]","beta[3]", "beta[4]",
                          "scale_FOI",
                          'alpha[2]','alpha[3]','alpha[4]',
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 

# only present min and max FOI and present in scientific notation 
lambda_post_CYD = tidy_post_CYD %>%
  filter(grepl("lambda", variable)) %>%
  filter(mean == max(mean) |
           mean == min(mean)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)


out_CYD = tidy_post_CYD %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, as.character) %>% 
  bind_rows(lambda_post_CYD) %>% 
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 

write.csv(out_CYD, "CYD/output/M12/posterior_formatted.csv")  


# Format Dengvaxia severe fit for table in manuscript

posterior_CYD_sev = read.csv("CYD/output/severe/M9/posterior.csv")

tidy_post_CYD_sev = posterior_CYD_sev %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__','pK2[1]','pK2[2]','pK2[3]','pK2[4]',
                          "delta[3]", "delta[4]", 
                          'rho[2]','rho[3]','rho[4]',
                          'w[2]','w[3]','w[4]', "tau[2]", "tau[3]", "tau[4]",
                          "lc[1,2]", "lc[1,3]", "lc[1,4]",
                          "lc[3,1]", "lc[3,2]", "lc[3,3]","lc[3,4]",
                          "kappa", "beta[3]", "beta[4]",
                          "scale_FOI",
                          'alpha[2]','alpha[3]','alpha[4]',
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 


# only present min and max FOI and present in scientific notation 
lambda_post_CYD_sev = tidy_post_CYD_sev %>%
  filter(grepl("lambda", variable)) %>%
  filter(mean == max(mean) |
           mean == min(mean)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)


out_CYD_sev = tidy_post_CYD_sev %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, as.character) %>% 
  bind_rows(lambda_post_CYD_sev) %>% 
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 

write.csv(out_CYD_sev, "CYD/output/severe/M9/posterior_formatted.csv")  


# Format Butantan-DV posterior distribution for table in manuscript 

# source data 

posterior_BUT = read.csv("BUT/output/M7/posterior.csv")

tidy_post_BUT = posterior_BUT %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__', "epsilon", 'rho[2]',  'L[1]',
                          'L[2]', 'w[2]','w[3]', "lc[1,2]",  
                          "lc[3,1]", "lc[3,2]", "beta[1]","beta[2]", "kappa",
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 

# FOI  in scientific notation 
lambda_post_BUT = tidy_post_BUT %>%
  filter(grepl("lambda", variable)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)

out_BUT = tidy_post_BUT %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, as.character) %>% 
  bind_rows(lambda_post_BUT) %>% 
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 


write.csv(out_BUT, "BUT/output/M7/posterior_formatted.csv")  

# Format Qdenga posterior distribution for table in manuscript 

posterior_TAK = read.csv("TAK/output/M2/posterior.csv")

tidy_post_TAK = posterior_TAK %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__','pK3[1]','pK3[2]','pK3[3]','pK3[4]',
                          "epsilon", 'rho[2]','rho[3]','rho[4]',
                          'L[2]','L[3]','L[4]', 'w[2]','w[3]','w[4]',
                          "lc[1,2]", "lc[3,1]", "lc[3,2]", "lc[1,3]", "lc[3,3]", 
                          "lc[1,4]", "lc[3,4]", 'w[2]','w[3]','w[4]',
                          'alpha[2]','alpha[3]','alpha[4]', "beta[2]",
                          "kappa",
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 

# only present min and max FOI and present in scientific notation 
lambda_post_TAK = tidy_post_TAK %>%
  filter(grepl("lambda", variable)) %>%
  filter(mean == max(mean) |
           mean == min(mean)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)

tidy_post_TAK %>%
  filter(grepl("lambda", variable)) %>% 
  separate(variable, into  = c("x", "y"), sep = "\\[") %>% 
  separate(y, into = c("y", "z"), sep = ",") %>% 
  group_by(z) %>%  
  summarise(cum_FOI = sum(mean)) %>% 
  filter(cum_FOI == max(cum_FOI) |
           cum_FOI == min(cum_FOI))

out_TAK = tidy_post_TAK %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, formatC, 2, format="f") %>%  
  mutate_if(is.numeric, as.factor) %>% 
  bind_rows(lambda_post_TAK) %>% 
  unite("X", mean:q5, sep = " (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 


write.csv(out_TAK, "TAK/output/M2/posterior_formatted.csv")  



# Calculate average yearly FOI over each trial 
Deng_weights = c(length(1:13), length(14:24) ,length(25:36), length(37:60))
Tak_weights = c(12,6,6,12,12,6) 
But_weights = 24

tidy_post_CYD %>%
  filter(grepl("lambda", variable)) %>%  
  mutate(weight = Deng_weights) %>% 
  mutate(FOI = mean * Deng_weights) %>% # monthly FOI multiplied by number of months 
  summarise(FOI = sum(FOI)) %>% # sum over t_d,
  mutate(FOI = round(FOI / 5, 2 ))  # divide by number of years for yearly FOI 


tidy_post_BUT %>%
  filter(grepl("lambda", variable))  %>% 
  summarise(cum_FOI = sum(mean)) %>%  # cumulative FOI 
  mutate(cum_FOI = cum_FOI * But_weights) %>%  # monthly FOI multiplied by number of months 
  mutate(cum_FOI = round(cum_FOI / 2,2 ))  # divide by number of years for yearly FOI 

tidy_post_TAK %>%
  filter(grepl("lambda", variable)) %>% 
  separate(variable, into  = c("x", "y"), sep = "\\[") %>% 
  separate(y, into = c("y", "z"), sep = ",") %>% 
  group_by(z) %>%  
  summarise(cum_FOI = sum(mean)) %>%  # cumulative FOI 
  mutate(weight = Tak_weights) %>% 
  mutate(w_mean = cum_FOI * weight) %>%  # monthly FOI multiplied by number of months 
  summarise(FOI = sum(w_mean)) %>% # sum over t_d,
  mutate(FOI = round(FOI / 4.5, 2 ))  # divide by number of years for yearly FOI 

  
