# Descriptive figures of the Qdenga data 
# set up -----------------------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(ggpubr)
library(cowplot)

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 16),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) 
  
) 

# functions

source("TAK/R/factor_data.R")

# read data 

cases = readRDS("TAK/data/processed/case_data.RDS")
GMT  = read.csv(file = "TAK/data/processed/antibody_data.csv")[,-1]

age_fill = scales::brewer_pal(palette = "Blues")(4)[2:4]
serotype_fill = c(scales::brewer_pal(palette = "RdPu")(6)[2:5]) 
trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]

# plot VCD data by serostatus and serotype/ age group --------------------------
serotype_plot = cases$VCD_BVKD %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = serotype),  
           position = "stack", stat = "identity",width = 0.5) +
  ylab("Symptomatic cases") + xlab(" ") +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0,140,50)) +
  scale_fill_manual(values = serotype_fill) +
  scale_x_discrete(labels  =  c("1-12", 
                                "13-18", 
                                "19-24",
                                "25-36", 
                                "37-48",
                                "49-54")) +
  facet_grid(serostatus~trial) + 
  theme(legend.position = c(0.11,0.86)) 

age_plot =  cases$VCD_BVJA %>% 
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),   
           position = "stack", stat = "identity", width = 0.8) +
  ylab(" ") + xlab(" ") +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0,150,50)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("1-12", 
                               "13-18",
                               "19-24", 
                               "25-36")) +
  facet_grid(serostatus~trial) +
  theme(legend.position = c(0.11,0.86))

# plot hosp data by serostatus and serotype/ age group -------------------------

hosp_serotype_plot = cases$N_hosp_BVK4 %>%  
  arrange(month) %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = serotype),  
           position = "stack", stat = "identity", width = 0.8) +
  ylab("Hospitalised cases") + xlab(" ") +
  scale_y_continuous(limits = c(0,70), breaks = seq(0,70,20)) +
  scale_fill_manual(values = serotype_fill) +
  scale_x_discrete(labels  = c("1-24", 
                               "25-36", 
                               "37-48", 
                               "49-54")) +
  facet_grid(serostatus~trial) +
  theme(legend.position = "none")

hosp_age_plot = cases$HOSP_BVJA %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),  
           position = "stack", stat = "identity", width = 0.8) +
  ylab(" ") + xlab(" ") +
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,10)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("1-12", 
                               "13-18",
                               "19-24", 
                               "25-36")) +
  facet_grid(serostatus~trial) +
  theme(legend.position = "none")

# plot serotype age data -------------------------------------------------------
plot_VCD_age_serotype = cases$VCD_KJ2 %>%  
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),  
           position = "stack", stat = "identity", width = 0.8) +
  ylab(" ") + xlab(" ") +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,50)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("1-12", 
                               "13-24")) +
  facet_wrap(~serotype, ncol =1) +
  theme(legend.position = "none")

plot_hosp_age_serotype = cases$HOSP_KJ2 %>%
  ggplot(aes(x = as.factor(month), y = Y)) +
  geom_bar(aes(fill = age),
           position = "stack", stat = "identity", width = 0.8) +
  ylab(" ") + xlab(" ") +
  scale_y_continuous(limits = c(0,50), breaks = seq(0,50,20)) +
  scale_fill_manual(values = age_fill) +
  scale_x_discrete(labels  = c("1-12", 
                               "13-18",
                               "19-24", 
                               "25-36")) +
  facet_wrap(~serotype,  ncol =1) +
  theme(legend.position = "none")

# grid plot of all case / pop figures ------------------------------------------
vcd_plot = plot_grid(serotype_plot, age_plot, plot_VCD_age_serotype,
                     hosp_serotype_plot, hosp_age_plot, plot_hosp_age_serotype,
                     ncol = 3, align = "vh", axis = "rlbt",
                     rel_widths = c(1,1,0.4),
                     labels = c("a", "b", "c", "d", "e", "f"))

vcd_plot_annot = annotate_figure(vcd_plot, bottom = text_grob("Month", size = 28))

ggsave(vcd_plot_annot, file = "TAK/output/figures/all_case.jpg",
       height = 25, width = 45, unit = "cm" )



# plot titres ------------------------------------------------------------------

# plot observed titres by trial, serostatus and serotype 
plot_titres = GMT %>%
  mutate(serotype=factor(serotype, 
                         labels = c("DENV1", "DENV2", "DENV3", "DENV4"))) %>% 
  mutate(trial = factor(trial, 
                        levels = c("Placebo", "TAK"), 
                        labels = c("placebo", "vaccine"))) %>% 
  mutate(serostatus = factor(serostatus,
                             levels = c("SN", "SP"),
                             labels = c("seronegative",
                                        "seropositive"))) %>%  
  filter(month > 0, trial == "vaccine") %>% 
  ggplot(aes(x = month, y = (mean))) +
  geom_line(aes(color = serostatus)) +
  geom_point(aes(color = serostatus))+
  # geom_errorbar(aes(ymin =(lower), ymax = (upper),
  #                   color = serostatus)) +
  ylab("Neutralising titre \ninduced by Qdenga") + xlab("Month PD2")+
  scale_x_continuous(limits = c(0,54), breaks = seq(0,54,12)) +
  facet_wrap(~ serotype, ncol=4) +
  theme(legend.position = c(0.88,0.8))+
  scale_colour_manual(values = trial_fill)


ggsave(plot_titres, file = "TAK/output/figures/plot_titres.jpg",
       height = 7, width = 24, unit = "cm" )
