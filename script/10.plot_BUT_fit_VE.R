# Script to plot the main figure including model fit and the vaccine efficacy 
# estimate 
rm(list=ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(wesanderson)
library(readxl)
library(cowplot)

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

# colours
age_fill = scales::brewer_pal(palette = "Blues")(4)[2:4]
serotype_fill = c(scales::brewer_pal(palette = "RdPu")(6)[2:5], "#CCCCCC") 


# source files 
file.sources = paste0("BUT/R/", list.files(path = "BUT/R/"))
sapply(file.sources, source)
path = "BUT/output/M7/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
serotype_serostatus_cases = read_excel("BUT/data/raw/cases.xlsx")
age_cases = read_excel("BUT/data/raw/cases.xlsx", sheet = 2)

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
  separate(name, into = c("serostatus", "arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         serotype = factor(serotype, labels = c(paste0("DENV", 1:2)))) %>% 
  bind_rows(AR_serotype_data) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type,color = serotype,
                 group = interaction(type, serotype)),
    position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    group = interaction(type, serotype),
      linetype = type, color = serotype),
    position = position_dodge(width =  0.5),width =  0.4,linewidth = 1) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  facet_wrap(~ serostatus) + theme_light() +
  theme(legend.position =c(0.87,0.73),
        text = element_text(size = 18),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0)) +
  scale_color_manual(values = serotype_fill)
  

# plot symp attack rate by age and trial arm
AR_plot_BVJ = AR_model %>%
  filter(group == "AR_BVJ") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"))) %>% 
  bind_rows(AR_age_data) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = age,
                 group = interaction(type, age)),
    position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    group = interaction(type, age),
      linetype = type, color = age),
    position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = " ") + # DELETE Y LAB IF COMBINED AR AND VE FIGURE 
  guides(shape = "none",linetype = "none") +
  theme(legend.position =c(0.85,0.8),
        text = element_text(size = 18),
        legend.title = element_blank()) +
  facet_wrap(~ serostatus)  + 
  scale_color_manual(values = age_fill)


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
  facet_wrap(~serostatus) + theme_light() + 
  theme(
        # legend.position = "none", # INCLUDE IS PLOTTING AR AND VE TOGETHER 
        legend.position = c(0.87,0.11), 
        text = element_text(size = 18),
        legend.title = element_blank()) +
  scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = serotype_fill) +
  scale_fill_manual(values = serotype_fill)


# plot sep 

g2 = plot_grid(AR_plot_BVK,  AR_plot_BVJ,
               labels = c("a", "b"),
               rel_widths = c(1.1,1))

ggsave(
  plot = g2,
  filename =  "BUT/output/figures/main_fit_B.png",
  height = 9,
  width = 29,
  units = "cm",
  dpi = 600,
  scale = 0.9
)

ggsave(
  plot = VE_plot,
  filename =  "BUT/output/figures/ve_B.png",
  height = 15,
  width = 30,
  units = "cm",
  dpi = 600,
  scale = 0.8
)
