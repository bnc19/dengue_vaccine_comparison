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
path = "CYD/output/severe/M9/" # (M10 in plots but started fitting from M0)

# Data 
VE_HOSP = readRDS("CYD/output/M12/VE.RDS")
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
VCD =  readRDS("CYD/data/processed/cases_stan_format.RDS")


theme_set(
  theme_light() +
    theme(
      text = element_text(size = 16),
      legend.position = c(0.04,0.74),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))

age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = c(scales::brewer_pal(palette = "RdPu")(6)[2:5], "#CCCCCC") 
trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]
cols = c("#BCBDDC", "#FC9272")


# add aggregated populations to data and calculate attack rates
AR_data = lapply(VCD, calc_CYD_attack_rates)
AR_model = extract_CYD_model_results(AR)  

# plot attack rates ------------------------------------------------------------
# plot symp attack rate by serotype and trial arm
Sy_AR_plot_VK = AR_model %>%
  filter(group == "V_AR_VK") %>%
  separate(name, into = c("arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4"))) %>% 
  bind_rows(AR_data$Sy_VK) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = serotype,group = interaction(type, serotype)),
    position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
      linetype = type, color = serotype),
    position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "") +
  scale_color_manual(values = serotype_fill) +
  theme(legend.position = "none")

# plot symp attack rate by age, serostatus and trial arm -----------------------
Sy_AR_plot_BVJ = AR_model %>%
  filter(group == "V_AR_BVJ") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Sy_BVJ) %>%
  unite(c(arm, serostatus), col = "x", sep = "\n") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  scale_color_manual(values = age_fill) +
  theme(legend.position = "none")

# plot hosp attack rate by age, serostatus, serotype and trial arm -------------
H_AR_plot_BVKJ = AR_model %>%
  filter(group == "H_AR_BVKJ") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Ho_BVKJ) %>%
  unite(c(arm, serostatus), col = "x", sep = "\n") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "Hospitalisation \nattack rate (%)") +
  scale_color_manual(values = serotype_fill) +
  theme(legend.position = "none") +
  facet_wrap(~ age) 

# plot hosp attack rate by age, serostatus, trial arm and time -----------------
H_AR_plot_BVJD =  AR_model %>%
  filter(group == "H_AR_BVJD") %>%
  separate(name, into = c("serostatus", "arm", "age", "time")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs")),
         time = factor(time, labels = c("1-13", "14-24", "25-36", "37-60"))) %>%
  bind_rows(AR_data$Ho_BVJD ) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  facet_grid(serostatus ~ arm) +
  labs(x = "Month", y = "Hospitalisation \nattack rate (%)") +
  scale_color_manual(values = age_fill) +
  theme(legend.position = "none")  

# plot severe attack rates -----------------------------------------------------

# get severe attack rate by age if seropositive and vaccine 
Se_AR_SP_V = AR_model %>%
  filter(group == "S_AR_spvJ") %>% 
  separate(name, into = c("age")) %>% 
  mutate(age = factor(age, labels = c("2-8yrs", "9-16yrs")),
         arm = factor("vaccine"),
         serotype = factor("all"),
         serostatus = factor("seropositive")) %>%
  bind_rows(AR_data$Se_BVKJ) %>%  
  filter(arm == "vaccine", serostatus == "seropositive") 

# plot severe attack rate by serotype, trial arm, serotype and age
Se_AR_BVKJ = AR_model %>%
  filter(group == "S_AR_BVKJ") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(
    arm = factor(arm, labels=c("placebo", "vaccine")),
    serotype = factor(serotype, labels =c("DENV1","DENV2","DENV3", "DENV4")),
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")), 
    age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>%  
  bind_rows(AR_data$Se_BVKJ) %>%  
  filter(!(arm == "vaccine" & serostatus == "seropositive")) %>%  
  bind_rows(Se_AR_SP_V) %>%  
  ggplot(aes(x = serostatus, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "Severe \nattack rate (%)") +
  scale_color_manual(values = serotype_fill) +
  guides(shape = "none", linetype = "none")  +
facet_grid(age ~ arm)

# plot severe attack rate by age, arm, time in SN   ----------------------------
Se_AR_VJD = AR_model %>%
  filter(group == "S_AR_snVJD") %>%
  separate(name, into = c("arm", "age", "time")) %>%
  mutate(
    time = factor(time, labels = c("1-13", "14-24", "25-36", "37-60")),
    arm = factor(arm, labels = c("placebo", "vaccine")),
    age = factor(age, labels = c("2-8yrs", "9-16yrs"))
  ) %>%  
  bind_rows(AR_data$Se_VJD) %>%  
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = "Month", y = "Severe \nattack rate (%)") +
  facet_wrap(~ arm) +
  scale_color_manual(values = age_fill)

# save AR plot -----------------------------------------------------------------
AR_grid = cowplot::plot_grid(
  Se_AR_VJD,  
  Se_AR_BVKJ,
  cowplot::plot_grid( 
    Sy_AR_plot_BVJ,
    Sy_AR_plot_VK, ncol =2,
    labels = c("c", "d")), 
  H_AR_plot_BVJD,
  H_AR_plot_BVKJ,
  ncol = 1,
  rel_heights = c(1,1.3,1), 
  labels = c("a", "b", " ", "e", "f")
)

ggsave(
  plot = AR_grid,
  filename =  "CYD/output/figures/severe_fit_C.png",
  height = 50,
  width = 40,
  units = "cm",
  dpi = 600,
  scale = 0.9
)

# plot VE  ---------------------------------------------------------------------
VE_model = extract_CYD_model_results(VE)
VE_model_HOSP = extract_CYD_model_results(VE_HOSP)

VE_severe = VE_model %>% 
  filter(group == "VE") %>% 
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
  geom_line(aes(color = outcome)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = outcome), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(serostatus ~ age ~ serotype, scales = "free") + theme_light() + 
  theme(legend.position = c(0.88,0.9),
        text = element_text(size = 18),
        legend.title = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
      scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

ggsave(
  plot = VE_plot,
  filename =  "CYD/output/figures/severe_VE_plot_C.png",
  height = 30,
  width = 40,
  units = "cm",
  dpi = 600,
  scale = 0.9
)


# plot VE severe vs. VE hosp ---------------------------------------------------
VE_hosp = VE_model_HOSP %>%  
  filter(group == "VE_BKRT") %>% 
  separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>%
  mutate(
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  filter(outcome == "hospitalised")

hosp_vs_sev = VE_severe %>%
  filter(outcome == "severe") %>% 
  bind_rows(VE_hosp) %>% 
  pivot_wider(names_from = outcome, values_from = lower:upper) %>% 
  ggplot(aes(x = mean_severe, y = mean_hospitalised)) +
  geom_line(aes(color = serotype), linewidth = 2) +
  #geom_smooth(aes(color = serotype), method = "lm") + 
  facet_wrap(~ serostatus, scales = "free") +
  geom_abline() + 
  scale_color_manual(values = serotype_fill) + 
  labs(y = "Vaccine efficacy \nagainst hospitalsiation", 
       x = "Vaccine efficacy \nagainst severe disease") + 
  theme(legend.position = c(0.24,0.21)) + 
  stat_cor(method="pearson") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(-300, 100)),
      scale_x_continuous(limits = c(50, 100)),
      scale_x_continuous(limits = c(90, 100))
    ),
    y = list(
      scale_y_continuous(limits = c(-300, 100)),
      scale_y_continuous(limits = c(50, 100)),
      scale_y_continuous(limits = c(90, 100))
    )
  )

ggsave(
  plot = hosp_vs_sev,
  filename =  "CYD/output/figures/hosp_vs_sev_C.png",
  height = 10,
  width = 30,
  units = "cm",
  dpi = 600,
  scale = 0.9
)

  