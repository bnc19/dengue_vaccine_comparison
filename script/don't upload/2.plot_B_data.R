rm(list = ls())

library(patchwork)
library(tidyverse)

cases = readRDS("data/processed_bu/cases_stan_format.RDS")

# colours 
serotype_fill = c("#8EAFBF", "#D55E00", "#CC79A7", "#016C59", "darkgrey")
age_fill = c("#5F9E9C", "#9A92B2", "#B07B4C")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1,0,0,0), "cm"),
      legend.margin = margin(0, 0, 0, 0)
    ) 
  
) 

# plot case data ---------------------------------------------------------------

Sy_BVJ_plot  = cases$Sy_BVJ %>%   
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab("Symptomatic cases") + xlab("") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_fill_manual(values = age_fill) +
  theme(legend.position = c(0.92,0.7)) +
  facet_grid(~serostatus) 

Sy_BVK_plot = cases$Sy_BVK %>%   
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = serotype), position = "stack", stat = "identity") +
  ylab("Symptomatic cases") + xlab(" ") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_fill_manual(values = serotype_fill) +
  theme(legend.position = c(0.93,0.8)) +
  facet_grid(~serostatus) 

final_plot = Sy_BVK_plot / Sy_BVJ_plot +
  plot_annotation(tag_levels = 'a') 

ggsave(final_plot, file = "output/figures/SupFig3.jpg",
       height = 12, width = 18, unit = "cm" )

