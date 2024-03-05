# Plot VE of Dengvaxia and Qdenga side by side
rm(list = ls())
library(tidyverse)
source("compare_vaccines/R/plotting_functions.R")

# colours 
mycols = c( "#FFCC99", "#9999FF")

# files 
TAK_file = "TAK/output/M2/VE.RDS"
CYD_file = "CYD/output/M12/VE.RDS"

set.seed(123)

CYD_VE = readRDS(CYD_file)
TAK_VE = readRDS(TAK_file)

# VE 
C_VE_model = extract_model_results(CYD_VE)  
T_VE_model = extract_model_results(TAK_VE)  
C_VE_model$vaccine = "Dengvaxia"
T_VE_model$vaccine = "Qdenga"

# combine both vaccines
comb_VE = C_VE_model %>%
  bind_rows(T_VE_model) %>% 
  filter(group == "VE_BKRT") %>%  
  separate(name, into = c("serostatus", "serotype","outcome","time")) %>%  
  mutate(
    month = as.numeric(time),
    serostatus = factor(serostatus, levels = 1:3, labels = c("seronegative", "monotypic", "multitypic")),
    serotype = factor(serotype, levels = 1:4, labels = c(paste0("DENV", 1:4))),
    outcome = factor(outcome, levels = 1:2,  labels = c("symptomatic", "hospitalised"))
    )

plot_ve_over_time = comb_VE %>%   
  filter(month < 55) %>% 
  ggplot(aes(x = month, y = mean)) +
  geom_line(aes(color = vaccine)) +
  geom_ribbon(aes( ymin = lower, ymax = upper, fill = vaccine), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  geom_hline(yintercept = 0, linetype = "dashed",color = "black", linewidth=1) +
  facet_grid(serostatus + outcome ~ serotype, scales = "free") + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5), legend.position = "top")+
  scale_color_manual(values = mycols) +
  scale_fill_manual(values = mycols)


ggsave(
  plot = plot_ve_over_time,
  filename = "compare_vaccines/output/VE_overtime.jpg",
  height = 30,
  width = 30,
  units = "cm",
  dpi = 300,
  scale = 0.6
)

