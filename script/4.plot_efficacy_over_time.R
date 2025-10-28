# Script to plot VE of Dengvaxia and Qdenga up to 54 months. Script outputs Fig.3 and Fig.4. 

rm(list = ls())
library(tidyverse)
source("R/plotting_functions.R")

# colours 
cols = c("#E6D192", "#A0C4D6","#8495B0")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.subtitle = element_text(size = 8),
      plot.caption = element_text(size = 8),
      legend.text = element_text(size = 6.5),
      legend.title = element_blank()
    ))

# files 
input = "output/M21/"

VE_Q = readRDS(paste0(input, "VE_Q.RDS"))  # without age 
VE_D = readRDS(paste0(input, "VE_De.RDS"))


# VE 
VE_D_model = extract_model_results(VE_D)  
VE_Q_model = extract_model_results(VE_Q)  
VE_D_model$vaccine = "Dengvaxia"
VE_Q_model$vaccine = "Qdenga"

# combine both vaccines
comb_VE = VE_Q_model %>%
  bind_rows(VE_D_model) %>% 
  filter(group == "VE_Q" | group == "VE_De") %>%  
  separate(name, into = c("serostatus", "serotype","age", "outcome","time")) %>%  
  filter(age == 1 | age == 2 & group == "VE_Q") %>% # only Qdenga has age
  mutate(
    month = as.numeric(time),
    serostatus = factor(serostatus, levels = 1:3, labels = c("seronegative", "monotypic", "multitypic")),
    serotype = factor(serotype, levels = 1:4, labels = c(paste0("DENV", 1:4))),
    outcome = factor(outcome, levels = 1:2,  labels = c("symptomatic", "hospitalised"))
    ) %>% 
  filter(month <=54) %>% 
  mutate(age = ifelse(vaccine == "Qdenga" & age == 1, "4-5yrs", 
                      ifelse(vaccine == "Qdenga" & age == 2, "6-16yrs",
                             ifelse(vaccine == "Dengvaxia", "2-16yrs", NA))))

plot_ve_symp = comb_VE %>%   
  filter(outcome == "symptomatic") %>% 
  ggplot(aes(x = month, y = mean, group = interaction(vaccine, age, sep = " "))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(vaccine, age, sep = " ")), alpha = 0.3) +
  geom_line(aes(color = interaction(vaccine, age, sep = " ")), size = 1.1) +
  labs(x = "Months post final dose", 
       y = "Vaccine efficacy\nagainst symptomatic disease (%)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=1) +
  facet_grid(serostatus ~ serotype) + 
  theme(legend.position = "top") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(0,12,24,36,48))

ggsave(
  plot =  plot_ve_symp,
  filename = "output/figures/Fig3.pdf",
  height = 14,
  width = 18,
  units = "cm",
  dpi = 300
)



plot_ve_hosp = comb_VE %>%   
  filter(outcome != "symptomatic") %>% 
  ggplot(aes(x = month, y = mean, group = interaction(vaccine, age, sep = " "))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = interaction(vaccine, age, sep = " ")), alpha = 0.3) +
  geom_line(aes(color = interaction(vaccine, age, sep = " ")), size = 1.1) +
  labs(x = "Months post final dose", 
       y = "Vaccine efficacy\nagainst hospitalisation (%)") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=1) +
  facet_grid(serostatus ~ serotype) + 
  theme(legend.position = "top") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(breaks = c(0,12,24,36,48))


ggsave(
  plot =  plot_ve_hosp,
  filename = "output/figures/Fig4.pdf",
  height = 14,
  width = 18,
  units = "cm",
  dpi = 300
)

