# Descriptive figures of the Qdenga data 
# set up -----------------------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(ggpubr)
library(patchwork)

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

# read data 
cases = readRDS("data/processed_Q/case_data.RDS")

serotype_fill = c("#8EAFBF", "#D55E00", "#CC79A7", "#016C59", "darkgrey")
age_fill = c("#5F9E9C", "#9A92B2", "#B07B4C")


# plot VCD data by serostatus and serotype/ age group --------------------------
serotype_plot = cases$VCD_BVKD %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = serotype),  
           position = "stack", stat = "identity",width = 0.5) +
  ylab("Symptomatic cases") +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140,50)) +
  scale_fill_manual(values = serotype_fill) +
  scale_x_discrete(labels  =  c("12", 
                                "18", 
                                "24",
                                "36", 
                                "48",
                                "54")) +
  facet_grid(serostatus~trial) + 
  theme(legend.position = c(0.11,0.86), 
        axis.title.x = element_blank()) 

age_plot =  cases$VCD_BVJA %>% 
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),   
           position = "stack", stat = "identity", width = 0.8) +
  ylab("Symptomatic cases") + 
  scale_y_continuous(limits = c(0, 150), breaks = seq(0,150,50)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("12", 
                               "18",
                               "24", 
                               "36")) +
  facet_grid(serostatus~trial) +
  theme(legend.position = c(0.11,0.86), 
        axis.title.x = element_blank())

# plot hosp data by serostatus and serotype/ age group -------------------------

hosp_serotype_plot = cases$N_hosp_BVK4 %>%  
  arrange(month) %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = serotype),  
           position = "stack", stat = "identity", width = 0.8) +
  ylab("Hospitalised cases") + 
  scale_y_continuous(limits = c(0,70), breaks = seq(0,70,20)) +
  scale_fill_manual(values = serotype_fill) +
  scale_x_discrete(labels  = c("24", 
                               "36", 
                               "48", 
                               "54")) +
  facet_grid(serostatus~trial) +
  theme(legend.position = "none", 
        axis.title.x = element_blank())

hosp_age_plot = cases$HOSP_BVJA %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),  
           position = "stack", stat = "identity", width = 0.8) +
  ylab("Hospitalised cases") + 
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,10)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("12", 
                               "18",
                               "24", 
                               "36")) +
  facet_grid(serostatus~trial) +
  theme(legend.position = "none", 
        axis.title.x = element_blank())

# plot serotype age data -------------------------------------------------------
plot_VCD_age_serotype = cases$VCD_KJ2 %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),  
           position = "stack", stat = "identity", width = 0.8) +
  ylab("Symptomatic cases") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,50)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("12", 
                               "24")) +
  facet_wrap(~serotype, ncol =2) +
  theme(legend.position = "none", 
        axis.title.x = element_blank())

plot_hosp_age_serotype = cases$HOSP_KJ2 %>%
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),
           position = "stack", stat = "identity", width = 0.8) +
  ylab("Hospitalised cases") +
  scale_y_continuous(limits = c(0,50), breaks = seq(0,50,20)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("12", 
                               "18",
                               "24", 
                               "36")) +
  facet_wrap(~serotype,  ncol =2) +
  theme(legend.position = "none", 
        axis.title.x = element_blank())

# grid plot of all case  figures ------------------------------------------

legend_plot = guide_area()  # Creates a space for the legend at the top

# Combine everything with explicit row ordering
final_plot = (
  legend_plot /  # Legend at the very top
    (serotype_plot | age_plot) / 
    (hosp_serotype_plot | hosp_age_plot) / 
    (plot_VCD_age_serotype | plot_hosp_age_serotype)
) +
  plot_layout(guides = "collect", heights = c(0.1, 1, 1, 1)) +  
  plot_annotation(tag_levels = 'a', caption = "Month") &  
  theme(
    legend.position = "top",
    plot.caption = element_text(hjust = 0.5, size = 14)
    
  )

ggsave(final_plot, file = "output/figures/SupFig1.jpg",
       height = 19, width = 18, unit = "cm" )
