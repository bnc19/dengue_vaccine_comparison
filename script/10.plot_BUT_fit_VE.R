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
serotype_serostatus_cases = read_excel("BUT/data/cases.xlsx")
age_cases = read_excel("BUT/data/cases.xlsx", sheet = 2)

# Tidy data 
serotype_serostatus_cases = factor_BUT_cases(serotype_serostatus_cases)
age_cases = factor_BUT_cases(age_cases)

# plot attack rates 
# add aggregated populations to data and calculate attack rates

AR_age_data = calc_BUT_attack_rates(age_cases)
AR_serotype_data = calc_BUT_attack_rates(serotype_serostatus_cases)
AR_model = extract_BUT_model_results(AR)

# plot serotype serostatus attack rate 
AR_plot_BVK = AR_model %>%
  filter(group == "AR_BVK") %>%
  separate(name, into = c("Serostatus", "Arm", "Serotype")) %>% 
  mutate(Arm = factor(Arm, labels = c("placebo", "vaccine")),
         Serostatus = factor(Serostatus, labels = c("seronegative", "seropositive")),
         Serotype = factor(Serotype, labels = c(paste0("DENV", 1:2)))) %>% 
  bind_rows(AR_serotype_data) %>%
  ggplot(aes(x = Arm, y = mean)) +
  geom_point(aes(shape = type,color = Serotype,group = interaction(type, Serotype)),
    position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, Serotype),
      linetype = type, color = Serotype),
    position = position_dodge(width =  0.5),width =  0.4,linewidth = 1) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  facet_wrap(~ Serostatus) + theme_light() +
  theme(legend.position =c(0.87,0.8),
        text = element_text(size = 18),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0)) +
  scale_color_brewer(palette = "Set2") 
  

# plot symp attack rate by age and trial arm

AR_plot_VJ = AR_model %>%
  filter(group == "AR_VJ") %>%
  separate(name, into = c("Arm", "Age")) %>% 
  mutate(Arm = factor(Arm, labels = c("placebo", "vaccine")),
         Age = factor(Age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"))) %>% 
  bind_rows(AR_age_data) %>%
  ggplot(aes(x = Arm, y = mean)) +
  geom_point(aes(shape = type, color = Age, group = interaction(type, Age)),
    position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, Age),
      linetype = type, color = Age),
    position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "") +
  scale_color_brewer(palette = "Accent") +  theme_light() +
  guides(shape = "none",linetype = "none") +
  theme(legend.position =c(0.87,0.8),
        text = element_text(size = 18),
        legend.title = element_blank()) 

# plot VE 

VE_model = extract_BUT_model_results(VE)

VE_plot =  VE_model %>%
  separate(name, into = c("Serostatus", "Serotype", "Age", "Month")) %>%
  filter(Age == 1) %>% # all the same 
  mutate(
    Serotype = factor(Serotype, 
                      labels = c("DENV1", "DENV2")),
    Month = as.numeric(Month),
    Serostatus = factor(Serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = Month , y = mean)) +
  geom_line(aes(color = Serotype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Serotype), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 24,6)) +
  facet_wrap(~Serostatus) + theme_light() + 
  theme(legend.position = "none",
        text = element_text(size = 18),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") 

# combine all plots 
g1 = (AR_plot_BVK + AR_plot_VJ )/  VE_plot + plot_annotation(tag_levels = 'a')
      

ggsave(
  plot = g1,
  filename =  "BUT/output/figures/main_fit_ve_fig.png",
  height = 30,
  width = 40,
  units = "cm",
  dpi = 600,
  scale = 0.8
)

