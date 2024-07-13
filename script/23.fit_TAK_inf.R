# Script to run TAK models starting with best fitting from impact modelling paper
rm(list=ls())


# load functions 
library(tidyverse)
library(Hmisc)
library(wesanderson)

# source functions
file.sources = paste0("TAK/R/", list.files(path = "TAK/R/"))
sapply(file.sources, source)
n_it = 10000

# colours 
age_fill_VE = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
age_fill_AR = scales::brewer_pal(palette = "Blues")(4)[2:4]
serotype_fill = c( "#FA9FB5", "#A6BDDB", "#C51B8A", "#1C9099") 
trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]


# best fitting model plus inf 
run_TAK_model (
  chi_C = 0, 
  BF = T,
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "INF_1",
  n_it = n_it,
  adapt_delta = 0.75,
  stan_model = "TAK_model_infection.stan",
  baseline_SP = read.csv("TAK/data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("TAK/data/processed/vcd_data.csv"),
  hosp = read.csv("TAK/data/processed/hosp_data.csv") ,
  infections = read.csv("TAK/data/processed/asymptomatic_cases.csv"),
  mu =  array(read.csv("TAK/data/processed/n0_new.csv")$mean, dim = c(2, 4))
)


run_TAK_model (
  chi_C = 1, 
  BF = T,
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "INF_2",
  n_it = n_it,
  adapt_delta = 0.75,
  stan_model = "TAK_model_infection.stan",
  baseline_SP = read.csv("TAK/data/raw/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("TAK/data/processed/vcd_data.csv"),
  hosp = read.csv("TAK/data/processed/hosp_data.csv") ,
  infections = read.csv("TAK/data/processed/asymptomatic_cases.csv"),
  mu =  array(read.csv("TAK/data/processed/n0_new.csv")$mean, dim = c(2, 4))
)


# plot TAK inf 

path = "TAK/output/INF_2/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
VCD =  read.csv("TAK/data/processed/vcd_data.csv")
hosp = read.csv("TAK/data/processed/hosp_data.csv")
asymp = read.csv("TAK/data/processed/asymp_data.csv")

# Tidy data 
VCD = factor_TAK_VCD(VCD)
hosp = factor_TAK_VCD(hosp)

asymp = asymp %>% 
  mutate(serostatus = factor(Serostatus,
                             levels = c("SN", "SP", "both"),
                             labels = c("seronegative", "seropositive", "both")),
         trial = factor(Arm, levels = c("Placebo", "TAK","both"), 
                        labels = c("placebo", "vaccine", "both")),
         year = Month,
         outcome = "asymp") %>% 
  select(-Serostatus, -Arm, -Outcome, -Month) 

# summarise over iterations 
VE_model = extract_TAK_model_results(VE)  
AR_model = extract_TAK_model_results(AR)  

# plot VE inf fit 

AR_model %>%
  filter(group == "AR_As_BV3") %>%
  separate(name, into = c("serostatus", "trial", "year")) %>%
  mutate(
    trial = factor(trial, levels = c(1,2), labels = c("placebo", "vaccine")),
    year = ifelse(year == 1, 6,
                  ifelse(year == 2, 12, 24)),
    serostatus = factor(serostatus, levels = c(1,2), labels = c("seronegative", "seropositive"))) %>%
  bind_rows(asymp_AR) %>% 
  ggplot(aes(x = year, y = mean )) +
  geom_point(aes(shape = type, color = serostatus, group = interaction(type, serostatus)),
             position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serostatus),
      linetype = type,
      color = serostatus
    ),
    position = position_dodge(width =  0.7),
    width =  0.4,
    linewidth = 1
  )  +
  facet_wrap(~trial)


# calculate attack rates from data 
AR_data = calc_TAK_attack_rates(VCD)
HR_data = calc_TAK_hosp_rates(VCD=VCD, hosp = hosp)

asymp_AR = asymp %>% 
  mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
         lower = binconf(Y,N, method = "exact")[,2] * 100, 
         upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
  select(- c(Y,N)) %>% 
  mutate(type = "data") 

data = bind_rows(AR_data, HR_data)

theme_set(
  theme_light() +
    theme(
      text = element_text(size = 14),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      # legend.position = c(0.86,0.25),  CHANGE TO THIS IS PLOTTING VE AND AR TOGETHER 
      legend.position = c(0.9,0.75), 
      legend.title = element_blank()
    ))

# plot VE 
VE_plot =  VE_model %>%
  filter(group == "VE") %>%
  separate(name, into = c("Serostatus", "Serotype", "Age","Outcome", "Month")) %>%
  mutate(
    Age = ifelse(Age ==1, 1, 2), # age groups 2 and 3 have same VE 
    Age = factor(Age, labels = c("4-5yrs","6-16yrs")),
    Outcome = factor(Outcome, labels = c("symptomatic", "hospitalised")),
    Serotype = factor(Serotype, labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    Month = as.numeric(Month),
    Serostatus = factor(Serostatus, labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = Month , y = mean)) +
  geom_line(aes(color = Age)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Age), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 54,12)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  facet_grid(Serostatus+Outcome ~ Serotype, scale= "free" ) +
  theme(legend.position = c(0.07,0.07)) +
  scale_color_manual(values = age_fill_VE)+
  scale_fill_manual(values = age_fill_VE)


# plot VE inf 
VE_plot_inf =  VE_model %>%
  filter(group == "VE_asymp") %>%
  separate(name, into = c("Serostatus", "Serotype", "Age", "Month")) %>%
  mutate(
    Age = ifelse(Age ==1, 1, 2), # age groups 2 and 3 have same VE 
    Age = factor(Age, labels = c("4-5yrs","6-16yrs")),
    Serotype = factor(Serotype, labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    Month = as.numeric(Month),
    Serostatus = factor(Serostatus, labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = Month , y = mean)) +
  geom_line(aes(color = Serotype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Serotype), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 54,12)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  facet_grid(Serotype ~Serostatus+Age, scale= "free" ) +
  theme(legend.position = "none") +
  scale_color_manual(values = serotype_fill)+
  scale_fill_manual(values = serotype_fill)


# Plot AR  across time 
AR_VRD = data %>%
  filter(serostatus == "both",
         age == "all",
         trial != "both",
         serotype == "all") %>%
  mutate(trial = ifelse(trial == "TAK-003", "vaccine", trial)) %>%
  select(-serotype,-age,-serostatus,-X, -year) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  ))

time_plot =  AR_model %>%
  filter(group == "AR_VRD") %>%
  separate(name, into = c("trial", "outcome", "month")) %>%
  mutate(
    trial = factor(trial, labels = c("placebo", "vaccine")),
    month = ifelse(month == 1, 12,
                   ifelse(
                     month == 2, 18,
                     ifelse(month == 3, 24,
                            ifelse(month == 4, 36,
                                   ifelse(month == 5, 48, 54))))),
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))) %>%
  bind_rows(AR_VRD) %>%
  mutate(month = as.numeric(month)) %>%
  ggplot(aes(x = month, y = mean)) +
  geom_point(aes(color = trial, shape = type), position = position_dodge(width = 5), size = 3) +
  geom_errorbar( aes(ymin = lower ,ymax = upper ,color = trial, linetype = type),
                 position = position_dodge(width = 5),width =  0.4,linewidth =1) +
  labs(x = "Month", y = "Attack rate (%)") +
  scale_x_continuous(breaks = c(12,18,24,36,48,54)) +
  scale_color_manual(values = trial_fill) +
  facet_wrap(~ outcome, ncol = 2) # CHANGE TO 1 IF PLOTTING AR AND VE TOGETHER 


# Plot AR by age and trial and serostatus 

AR_BVJR = data %>%
  filter(serotype == "all" |
           age != "all", serostatus != "both", month == 1000) %>%
  mutate(trial = ifelse(trial == "TAK-003", "vaccine", trial)) %>%
  select(-month,-X, -year,-serotype) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  ))

age_plot =  AR_model %>%
  filter(group == "AR_BVJR") %>%
  separate(name, into = c("serostatus", "trial", "age", "outcome")) %>% 
  mutate(outcome = factor(outcome,labels = c("symptomatic", "hospitalised"))) %>% 
  mutate(serostatus = factor(serostatus, labels=c("seronegative", "seropositive")),
         trial = factor(trial, labels=c("placebo", "vaccine")),
         age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs"),
                      levels = c(1,2,3))) %>% 
  bind_rows(AR_BVJR) %>%
  unite(c(trial, serostatus), col = "x") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(
    aes( shape = type,  color = age, group = interaction(type, age)),
    position = position_dodge(width = 0.5), size = 3 ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, group = interaction(type, age),
        linetype = type,color = age),
    position = position_dodge(width =  0.5),
    width =  0.4,
    linewidth = 1
  ) +
  # labs(x = " ", y ="" ) +  # CHANGE IF PLOTTING AR AND VE TOGETHER 
  labs(x = " ", y = "Attack rate (%)") +
  scale_color_manual(values = age_fill_AR) +
  scale_x_discrete(
    labels  = c(
      "placebo \nseronegative",
      "placebo \nseropositive",
      "vaccine \nseronegative",
      "vaccine \nseropositive"
    )) +
  guides(shape = "none",
         linetype = "none")+
  facet_wrap(~ outcome, ncol=2) # CHANGE TO 1 IF PLOTTING AR AND VE TOGETHER 


# Plot AR by serotype and trial and serostatus 

AR_BVKR = data %>%
  filter(serotype != "all" |
           age == "all", serostatus != "both", month == 1000) %>%
  mutate(trial = ifelse(trial == "TAK-003", "vaccine", trial)) %>%
  select(-month,-X, -year,-age) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  )) %>%
  mutate(serotype = factor(serotype, labels = c("DENV1", "DENV2", "DENV3", "DENV4")))


serotype_plot = AR_model %>%
  filter(group == "AR_BVKR") %>%
  separate(name, into = c("serostatus", "trial", "serotype", "outcome")) %>%
  mutate(outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))) %>%
  mutate(
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
    trial = factor(trial, labels = c("placebo", "vaccine")),
    serotype = factor(serotype, labels = c("DENV1", "DENV2", "DENV3", "DENV4"))) %>%
  bind_rows(AR_BVKR) %>% 
  unite(c(trial, serostatus), col = "x") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serotype),
      linetype = type,
      color = serotype
    ),
    position = position_dodge(width =  0.7),
    width =  0.4,
    linewidth = 1
  ) +
  # labs(x = "Month", y = "Attack rate (%)") + # CHANGE IF PLOTTING AR AND VE TOGETHER 
  labs(x = " ", y = "Attack rate (%)") +
  scale_color_manual(values = serotype_fill) +
  scale_x_discrete(
    labels  = c(
      "placebo \nseronegative",
      "placebo \nseropositive",
      "vaccine \nseronegative",
      "vaccine \nseropositive")) +
  guides(shape = "none",
         linetype = "none")+
  facet_wrap(~ outcome, ncol=2) # CHANGE TO 1 IF PLOTTING AR AND VE TOGETHER 

# combine all plots 
g1 = cowplot::plot_grid(time_plot, age_plot, 
                        serotype_plot, ncol=3, 
                        axis = "tblr", align = "h",
                        labels = c("a", "b", "c"))

g2 = cowplot::plot_grid(g1, VE_plot, rel_heights = c(1,1.8),
                        ncol =1, labels = c("", "d"))  

