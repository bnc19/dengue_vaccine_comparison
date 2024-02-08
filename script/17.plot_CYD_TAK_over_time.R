# Plot VE of Dengvaxia and Qdenga side by side
library(tidyverse)
source("R/plotting_functions.R")
n= 100

TAK_file = "/Users/bethancracknelldaniels/Desktop/R projects/Qdenga_efficacy/output/final/M32/VE.RDS"
CYD_file = "/Users/bethancracknelldaniels/Desktop/R projects/CYD_Efficacy/output/M10/VE.RDS"
set.seed(123)
CYD_VE = as.data.frame(sapply(readRDS(CYD_file), sample, n))
TAK_VE = as.data.frame(sapply(readRDS(TAK_file), sample, n))

# VE 
C_VE_model = extract_model_results(CYD_VE)  
T_VE_model = extract_model_results(TAK_VE)  
C_VE_model$vaccine = "CYD"
T_VE_model$vaccine = "TAK"

# combine both vaccines

comb_VE = C_VE_model %>%
  bind_rows(T_VE_model) %>% 
  filter(group == "VE_BKRT") %>%  
  separate(name, into = c("serostatus", "serotype","outcome","time")) %>%  
  mutate(month= as.numeric(time),
         serostatus = factor(serostatus, levels = 1:3, 
                             labels = c("seronegative", "monotypic", "multitypic")),
         serotype = factor(serotype, levels = 1:4, 
                           labels = c(paste0("DENV", 1:4))),
         outcome = factor(outcome, levels = 1:2, 
                          labels = c("symptomatic", "hospitalised")), 
         vaccine = factor(vaccine, levels = c("CYD", "TAK"), 
                          labels = c("Dengvaxia", "Qdenga")))


mycols = c("#9999FF", "#FF9900")

plot_ve_over_time = comb_VE %>%   
  filter(month < 55) %>% 
  filter(serostatus != "multitypic") %>% 
  ggplot(aes(x = month, y = mean)) +
  geom_line(aes(color = vaccine)) +
  geom_ribbon(
    aes(
      ymin = lower ,
      ymax = upper,
      fill = vaccine),
    alpha = 0.5
  ) +
  labs(x = "Month", y = "Vaccine Efficacy (%)") +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  facet_grid(serotype~outcome+serostatus) + 
  theme_light() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")+
  scale_color_manual(values= mycols) +
  scale_fill_manual(values= mycols)


ggsave(
  plot = plot_ve_over_time,
  filename = "output/VE_overtime.jpg",
  height = 20,
  width = 30,
  units = "cm",
  dpi = 600,
  scale = 0.6
)

