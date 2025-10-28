# Script to tidy the CYD case data 
rm(list = ls())
library(tidyverse)
library(patchwork)

cases = readRDS("data/processed_De/cases_stan_format.RDS")

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

Ho_BVJD_plot = cases$Ho_BVJD %>% 
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab("Hospitalised cases") + xlab("Month") +
  scale_fill_manual(values = age_fill) +
  facet_grid(arm~serostatus) +
  theme(legend.position = "none") +
  scale_x_discrete(labels  = c("13", 
                               "24",
                               "36", 
                               "60"))
  
Ho_BVKJ_plot = cases$Ho_BVKJ %>% 
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = serotype), position = "stack", stat = "identity") +
  ylab(" ") + 
  scale_fill_manual(values = serotype_fill) +
  facet_grid(age~serostatus) +
  theme(legend.position = "none", 
        axis.title.x = element_blank())  +
  guides(fill = "none")

Sy_BVJ_plot  = cases$Sy_BVJ %>%   
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab("Symptomatic cases") +
  scale_fill_manual(values = age_fill) +
  facet_grid(~serostatus) +
  theme(legend.position = c(0.84,0.84), 
        axis.title.x = element_blank())

Sy_VK_plot = cases$Sy_VK %>%   
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = serotype), position = "stack", stat = "identity") +
  ylab(" ") +
  scale_fill_manual(values = serotype_fill) +
  theme(legend.position = c(0.85,0.87), 
        axis.title.x = element_blank())+
  guides(fill = "none")

Se_VJD_plot = cases$Se_VJD %>% 
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  facet_grid(~ arm) +
  ylab("Severe cases") + xlab("Month") +
  scale_fill_manual(values = age_fill) +
  theme(legend.position = "none")+
  scale_x_discrete(labels  = c("13", 
                               "24",
                               "36", 
                               "60"))

Se_BVKJ_plot = cases$Se_BVKJ %>% 
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = serotype), position = "stack", stat = "identity") +
  facet_grid(age ~ serostatus) +
  ylab(" ") +
  scale_fill_manual(values = serotype_fill) +
  theme(legend.position = c("none"), 
        axis.title.x = element_blank())


# grid plot of all case  figures ------------------------------------------


legend_plot = guide_area()  # Creates a space for the legend at the top


final_plot = (
  legend_plot / (Sy_BVJ_plot | Sy_VK_plot) /
    (Ho_BVJD_plot | Ho_BVKJ_plot) /
    (Se_VJD_plot |  Se_BVKJ_plot)
) +
  plot_layout(guides = "collect", heights = c(0.1, 1, 1, 1)) +
  plot_annotation(tag_levels = 'a') &
  theme(legend.position = "top")

    
ggsave(final_plot, file = "output/figures/SupFig2.jpg",
       height = 19, width = 18, unit = "cm" )

