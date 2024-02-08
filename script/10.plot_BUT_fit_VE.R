# Script to plot the main figure including model fit and the vaccine efficacy 
# estimate 
rm(list=ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(wesanderson)
library(readxl)
library(patchwork)

theme_set(
  theme_light() +
    theme(
      text = element_text(size = 16),
      legend.position ="top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))

# source files 
file.sources = paste0("BUT/R/", list.files(path = "BUT/R/"))
sapply(file.sources, source)
path = "BUT/output//M4/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
cases = readRDS("BUT/data/processed/cases_stan_format.RDS")
# plot attack rates 
# add aggregated populations to data and calculate attack rates

AR_age_data = calc_BUT_attack_rates(cases$Sy_BVJ)
AR_serotype_data = calc_BUT_attack_rates(cases$Sy_BVK)
AR_model = extract_BUT_model_results(AR)

# plot serotype, serostatus attack rate 
AR_plot_BVK = AR_model %>%
  filter(group == "AR_BVK") %>%
  separate(name, into = c("serostatus", "arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         serotype = factor(serotype, labels = c(paste0("DENV", 1:2)))) %>% 
  bind_rows(AR_serotype_data) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type,color = serotype,group = interaction(type, serotype)),
    position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
      linetype = type, color = serotype),
    position = position_dodge(width =  0.5),width =  0.4,linewidth = 1) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  facet_wrap(~ serostatus) + theme_light() +
  theme(legend.position =c(0.87,0.8),
        text = element_text(size = 18),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0)) +
  scale_color_brewer(palette = "Set2") 
  

# plot symp attack rate by age and trial arm

AR_plot_BVJ = AR_model %>%
  filter(group == "AR_BVJ") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"))) %>% 
  bind_rows(AR_age_data) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
    position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
      linetype = type, color = age),
    position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "") +
  scale_color_brewer(palette = "Accent") +  theme_light() +
  guides(shape = "none",linetype = "none") +
  theme(legend.position =c(0.87,0.8),
        text = element_text(size = 18),
        legend.title = element_blank()) +
  facet_wrap(~ serostatus) 

# plot VE 

VE_model = extract_BUT_model_results(VE)

VE_plot =  VE_model %>%
  separate(name, into = c("serostatus", "serotype", "age", "month")) %>%
  filter(age == 1) %>% # all the same 
  mutate(
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = serotype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = serotype), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 24,6)) +
  facet_wrap(~serostatus) + theme_light() + 
  theme(legend.position = "none",
        text = element_text(size = 18),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") 

# combine all plots 
g1 = (AR_plot_BVK + AR_plot_BVJ )/  VE_plot + plot_annotation(tag_levels = 'a')
      

ggsave(
  plot = g1,
  filename =  "BUT/output/figures/main_fit_ve_fig.png",
  height = 30,
  width = 40,
  units = "cm",
  dpi = 600,
  scale = 0.8
)

