# Script to plot the main CYD figure including model fit and the vaccine  
# efficacy estimate 
rm(list = ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(readxl)
library(patchwork)

# source files 
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
path = "output/MS1/"

# Data 
VE_HOSP = readRDS("output/M17/VE_De.RDS")
VE = readRDS(paste0(path, "VE_De.RDS"))
AR = readRDS(paste0(path, "AR_De.RDS"))
VCD =  readRDS("data/processed_De/cases_stan_format.RDS")


# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59", "#111111")
trial_fill = c("#C51B8A", "#99CC99")
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")


theme_set(
  theme_bw()+
    theme(
      text = element_text(size = 12),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.title = element_blank())
  )



# add aggregated populations to data and calculate attack rates
AR_data = lapply(VCD, calc_attack_rates_De)
AR_model = extract_model_results_De(AR)  

AR_data$Ho_BVJD = AR_data$Ho_BVJD  %>%  
  mutate(time =  factor(time, labels = c("13", "24", "36", "60")))

AR_data$Se_VJD = AR_data$Se_VJD %>%  
  mutate(time =  factor(time, labels = c("13", "24", "36", "60")))


# plot attack rates ------------------------------------------------------------
# plot symp attack rate by serotype and trial arm
Sy_AR_plot_VK = AR_model %>%
  filter(group == "V_AR_VK_De") %>%
  separate(name, into = c("arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4"))) %>% 
  bind_rows(AR_data$Sy_VK) %>%
  mutate(arm = factor(arm, labels=c("P", "V"))) %>% 
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = serotype,group = interaction(type, serotype)),
    position = position_dodge(width = 0.5),size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
      linetype = type, color = serotype),
    position = position_dodge(width =  0.5), width =  0.4, linewidth = 0.8) +
  scale_color_manual(values = serotype_fill) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  guides(color  = "none") 

# plot symp attack rate by age, serostatus and trial arm -----------------------
Sy_AR_plot_BVJ = AR_model %>%
  filter(group == "V_AR_BVJ_De") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Sy_BVJ) %>%
  mutate(serostatus = factor(serostatus, labels=c("SN", "SP")),
         arm = factor(arm, labels=c("P", "V"))) %>% 
  unite(c(arm, serostatus), col = "x", sep = " ") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 0.8) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  scale_color_manual(values = age_fill) 
  

# plot hosp attack rate by age, serostatus, serotype and trial arm -------------
H_AR_plot_BVKJ = AR_model %>%
  filter(group == "H_AR_BVKJ_De") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Ho_BVKJ) %>%
  mutate(serostatus = factor(serostatus, labels=c("SN", "SP"))) %>% 
  ggplot(aes(x = serostatus, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.7),size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.7), width =  0.4, linewidth = 0.8) +
  scale_color_manual(values = serotype_fill) +
  facet_grid(arm~ age) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  guides(color  = "none") 


# plot hosp attack rate by age, serostatus, trial arm and time -----------------
H_AR_plot_BVJD =  AR_model %>%
  filter(group == "H_AR_BVJD_De") %>%
  separate(name, into = c("serostatus", "arm", "age", "time")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs")),
         time = factor(time, labels = c("13", "24", "36", "60"))) %>%
  bind_rows(AR_data$Ho_BVJD ) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.7),size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.7), width =  0.4, linewidth = 0.8) +
  facet_grid(serostatus ~ arm) +
  labs(x = "Month", y = "Hospitalisation \nattack rate (%)") +
  scale_color_manual(values = age_fill) 

# plot severe attack rates -----------------------------------------------------

# get severe attack rate by age if seropositive and vaccine 
Se_AR_SP_V = AR_model %>%
  filter(group == "S_AR_spvJ_De") %>% 
  separate(name, into = c("age")) %>% 
  mutate(age = factor(age, labels = c("2-8yrs", "9-16yrs")),
         arm = factor("vaccine"),
         serotype = factor("all"),
         serostatus = factor("seropositive")) %>%
  bind_rows(AR_data$Se_BVKJ) %>%  
  filter(arm == "vaccine", serostatus == "seropositive") 

# plot severe attack rate by serotype, trial arm, serotype and age
Se_AR_BVKJ = AR_model %>%
  filter(group == "S_AR_BVKJ_De") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(
    arm = factor(arm, labels=c("placebo", "vaccine")),
    serotype = factor(serotype, labels =c("DENV1","DENV2","DENV3", "DENV4")),
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")), 
    age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>%  
  bind_rows(AR_data$Se_BVKJ) %>%  
  filter(!(arm == "vaccine" & serostatus == "seropositive")) %>%  
  bind_rows(Se_AR_SP_V) %>%  
  mutate(serostatus = factor(serostatus, labels=c("SN", "SP"))) %>% 
  ggplot(aes(x = serostatus, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.7),size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.7), width =  0.4, linewidth = 0.8) +
  scale_color_manual(values = serotype_fill) +
facet_grid(age ~ arm)  +
  guides(shape = "none", linetype = "none") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

# plot severe attack rate by age, arm, time in SN   ----------------------------
Se_AR_VJD = AR_model %>%
  filter(group == "S_AR_snVJD_De") %>%
  separate(name, into = c("arm", "age", "time")) %>%
  mutate(
    time = factor(time, labels = c("13", "24", "36", "60")),
    arm = factor(arm, labels = c("placebo", "vaccine")),
    age = factor(age, labels = c("2-8yrs", "9-16yrs"))
  ) %>%  
  bind_rows(AR_data$Se_VJD) %>% 
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 0.8) +
  labs(x = "Month", y = "Severe \nattack rate (%)") +
  facet_wrap(~ arm, ncol=1) +
  scale_color_manual(values = age_fill) 

# save AR plot -----------------------------------------------------------------

plots = (Sy_AR_plot_BVJ | Sy_AR_plot_VK) /  (H_AR_plot_BVJD | H_AR_plot_BVKJ) / (Se_AR_VJD | Se_AR_BVKJ) + 
  plot_layout(guides = "collect", heights = c(1,2,1.5)) &  
  theme(legend.position = "top") 


ggsave(plots, file = "output/figures/SupFig6.jpg",
       height = 18, width = 20, units="cm")

# plot VE  ---------------------------------------------------------------------
VE_model = extract_model_results_De(VE)
VE_model_HOSP = extract_model_results_De(VE_HOSP)

VE_severe = VE_model %>% 
  filter(group == "VE_De") %>% 
  separate(name, into = c("serostatus", "serotype", "age", "outcome", "month")) %>%
  mutate(
    age = factor(age, labels = c("2-8yrs", "9-16yrs")), 
    outcome = factor(outcome, labels = c("symptomatic", "severe")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic")))

  
VE_plot = VE_severe %>%  
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = serostatus)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = serostatus), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(serotype ~ outcome + age) +
  theme(legend.position = c(0.09,0.082)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
      scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

# ggsave(
#   plot = VE_plot,
#   filename =  "output/figures/SupFig8.png",
#   height = 14,
#   width = 18,
#   units = "cm",
#   dpi = 600
# )

# severe vs hosp ---------------------------------------------------------------

VE_hosp = VE_model_HOSP %>%  
  filter(group == "VE_De") %>% 
  separate(name, into = c("serostatus", "serotype", "age", "outcome", "month")) %>%
  mutate(
    age = factor(age, labels = c("2-8yrs", "9-16yrs")), 
    outcome = factor(outcome, labels = c("symp", "hospitalised")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  filter(outcome == "hospitalised")

hosp_vs_severe = VE_severe %>%  
  filter(outcome == "severe") %>%  
  bind_rows(VE_hosp) %>%  
  filter(age == "2-8yrs") %>%  # VE is same for both
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = outcome)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = outcome), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(serostatus ~ serotype, scales = "free_y") +
  theme(legend.position = c(0.09,0.082)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

ggsave(
  plot = hosp_vs_severe,
  filename =  "output/figures/SupFig7.jpg",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600
)


