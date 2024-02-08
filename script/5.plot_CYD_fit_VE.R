# Script to plot the main CYD figure including model fit and the vaccine  
# efficacy estimate 
rm(list = ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(wesanderson)
library(readxl)
library(cowplot)

# source files 
file.sources = paste0("CYD/R/", list.files(path = "CYD/R/"))
sapply(file.sources, source)
path = "CYD/output/M12/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
VCD =  readRDS("CYD/data/processed/cases_stan_format.RDS")

theme_set(
  theme_light() +
    theme(
      text = element_text(size = 16),
      legend.position = c(0.85,0.74),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))

# plot attack rates 
AR_plot = plot_CYD_attack_rate(
  VCD =  VCD,
  file_path = path,
  AR = AR
)

# plot VE 
VE_model = extract_CYD_model_results(VE)

VE_plot =  VE_model %>%
  filter(group == "VE_BKRT") %>% 
  separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>%
  mutate(
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = serostatus)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = serostatus), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(serotype~outcome) + theme_light() + 
  theme(legend.position = "none",
        text = element_text(size = 18),
        legend.title = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1)
  
# combine all plots 
g1 = plot_grid(AR_plot, VE_plot, labels = c("", "e"), rel_widths = c(1.2,1))

ggsave(
  plot = g1,
  filename =  "CYD/output/figures/main_fit_ve_fig.png",
  height = 30,
  width = 55,
  units = "cm",
  dpi = 600,
  scale = 0.9
)

