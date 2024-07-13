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


# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59", "#111111")
trial_fill = c("#C51B8A", "#99CC99")
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")


theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)))

# plot attack rates 
AR_data = lapply(VCD, calc_CYD_attack_rates)
AR_model = extract_CYD_model_results(AR)  

# plot symp attack rate by serotype and trial arm ------------------------------
Sy_AR_plot_VK = AR_model %>%
  filter(group == "V_AR_VK") %>%
  separate(name, into = c("arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4"))) %>% 
  bind_rows(AR_data$Sy_VK) %>%
  mutate(serostatus = ifelse(serostatus == "seronegative", "SN", "SP")) %>% 
  mutate(arm = ifelse(arm == "placebo", "\nP", "\nV")) %>% 
  ggplot(aes(x = arm, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width = 0.5),
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serotype),
      linetype = type,
      color = serotype
    ),
    position = position_dodge(width =  0.5),
    width =  0.4,
    linewidth = 1
  ) +
  labs(x = " ", y = "") +
  scale_color_manual(values = serotype_fill) +
  guides(shape = "none", linetype = "none")+
  theme(legend.position = c(0.85,0.75))

# plot symp attack rate by age, serostatus and trial arm -----------------------
Sy_AR_plot_BVJ = AR_model %>%
  filter(group == "V_AR_BVJ") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Sy_BVJ) %>%
  mutate(serostatus = ifelse(serostatus == "seronegative", "SN", "SP")) %>% 
  mutate(arm = ifelse(arm == "placebo", "P", "V")) %>% 
  unite(c(arm, serostatus), col = "x", sep = "\n") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = age,
      group = interaction(type, age)),
    position = position_dodge(width = 0.5),
    size = 3) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, age),
      linetype = type,
      color = age),
    position = position_dodge(width =  0.5),
    width =  0.4,
    linewidth = 1) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  scale_color_manual(values = age_fill) +
  theme(legend.position = c(0.85,0.72))

# plot hosp attack rate by age, serostatus, serotype and trial arm -------------

H_AR_plot_BVKJ = AR_model %>%
  filter(group == "H_AR_BVKJ") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Ho_BVKJ) %>%
  mutate(serostatus = ifelse(serostatus == "seronegative", "SN", "SP")) %>% 
  mutate(arm = ifelse(arm == "placebo", "P", "V")) %>% 
  unite(c(arm, serostatus), col = "x", sep = "\n") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width = .8),
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serotype),
      linetype = type,
      color = serotype
    ),
    position = position_dodge(width =  .8),
    width =  0.4,
    linewidth = 1
  ) +
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
  geom_point(
    aes(
      shape = type,
      color = age,
      group = interaction(type, age)
    ),
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, age),
      linetype = type,
      color = age
    ),
    position = position_dodge(width =  0.6),
    width =  0.4,
    linewidth = 1
  ) +
  facet_grid(serostatus ~ arm) +
  labs(x = "Month", y = "Hospitalisation \nattack rate (%)") +
  scale_color_manual(values = age_fill) +
  theme(legend.position = "none")  

AR_plot = cowplot::plot_grid(
  cowplot::plot_grid( 
    Sy_AR_plot_BVJ,
    Sy_AR_plot_VK, ncol =2,
    labels = c("a", "b")), 
  H_AR_plot_BVJD,
  H_AR_plot_BVKJ,
  ncol = 1,
  rel_heights = c(1,1.3,1), 
  labels = c("", "c", "d")
)

ggsave(
  plot = AR_plot,
  filename =  "CYD/output/figures/main_fit_C.png",
  height = 22,
  width = 18,
  units = "cm",
  dpi = 600
)


# plot VE ----------------------------------------------------------------------
VE_model = extract_CYD_model_results(VE)

VE_f =  VE_model %>%
  filter(group == "VE_BKRT") %>% 
  separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>%
  mutate(
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) 
VE_plot = VE_f %>% 
ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = outcome)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = outcome), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(serostatus ~ serotype, scales = "free") +
  theme(legend.position = c(0.9,0.77),
        legend.title = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
      scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

ggsave(
  plot = VE_plot,
  filename =  "CYD/output/figures/VE_plot_C.png",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600
)


# plot fit compared to SA 
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")


# Data 
SA = readRDS("CYD/output/SA1/VE.RDS")
SA_model = extract_CYD_model_results(SA)

SA_f = SA_model %>%
  filter(group == "VE_BKRT") %>% 
  separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>%
  mutate(
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  mutate(cutoff = "48 months")


VE_f = VE_f %>% 
  mutate(cutoff = "60 months")


out = VE_f %>% 
  bind_rows(SA_f) %>% 
  ggplot(aes(x = month, y  = mean)) + 
  geom_line(aes(color = cutoff)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cutoff), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(outcome + serostatus ~ serotype, scales = "free") +
  theme(legend.position = c(0.9,0.89),
        legend.title = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

ggsave(
  plot = out,
  filename =  "CYD/output/figures/SA_cutoff.png",
  height = 16,
  width = 18,
  units = "cm",
  dpi = 600
)
