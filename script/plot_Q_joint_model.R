# Script to plot the main TAK figure including model fit and the vaccine  
# efficacy estimate 
rm(list=ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(wesanderson)

# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
age_fill_AR = scales::brewer_pal(palette = "Blues")(4)[c(2:4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59", "#111111")
trial_fill = c("#C51B8A", "#99CC99")
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")


# source files 
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
path = "output/M1/"

# Data 
VE = readRDS(paste0(path, "VE_Q.RDS"))
AR = readRDS(paste0(path, "AR_Q.RDS"))
VCD =  read.csv("data/processed_Q/vcd_data.csv")
hosp = read.csv("data/processed_Q/hosp_data.csv")

# Tidy data 
VCD = factor_VCD_Q(VCD)
hosp = factor_VCD_Q(hosp)

# summarise over iterations 
VE_model = extract_model_results_Q(VE)  
AR_model = extract_model_results_Q(AR)  

# calculate attack rates from data 
AR_data = calc_attack_rates_Q(VCD)
HR_data = calc_hosp_rates_Q(VCD=VCD, hosp = hosp)
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
  filter(group == "VE_Q") %>%
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
  scale_color_manual(values = age_fill)+
  scale_fill_manual(values = age_fill)
  
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
  filter(group == "AR_VRD_Q") %>%
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
    filter(group == "AR_BVJR_Q") %>%
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
  filter(group == "AR_BVKR_Q") %>%
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


# save sep plots 
g = cowplot::plot_grid(time_plot, age_plot, 
                        serotype_plot, ncol=1, 
                        axis = "tblr", align = "h",
                        labels = c("a", "b", "c"))

ggsave(
  plot = g,
  filename =  paste0(path, "AR_plot_Q.png"),
  height = 35,
  width = 32,
  units = "cm",
  dpi = 600,
  scale = 0.75
)

ggsave(
  plot = VE_plot,
  filename =  paste0(path, "VE_plot_Q.png"),
  height = 40,
  width = 45,
  units = "cm",
  dpi = 600,
  scale = 0.75
)

