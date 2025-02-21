# Script to plot the main CYD figure including model fit and the vaccine  
# efficacy estimate 
rm(list = ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(readxl)
library(cowplot)

# source files 
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
path = "output/M1/"

# Data 
VE = readRDS(paste0(path, "VE_De.RDS"))
AR = readRDS(paste0(path, "AR_De.RDS"))
VCD =  readRDS("data/processed_De/cases_stan_format.RDS")

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
AR_data = lapply(VCD, calc_attack_rates_De)
AR_model = extract_model_results_De(AR)  

# plot symp attack rate by serotype and trial arm ------------------------------
Sy_AR_plot_VK = AR_model %>%
  filter(group == paste0("V_AR_VK_De")) %>%
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
  filter(group == "V_AR_BVJ_De") %>%
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
  filter(group == "H_AR_BVKJ_De") %>%
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
  filter(group == "H_AR_BVJD_De") %>%
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
  filename =  paste0(path, "/AR_plot_De.png"),
  height = 22,
  width = 18,
  units = "cm",
  dpi = 600
)


# plot VE ----------------------------------------------------------------------
VE_model = extract_model_results_De(VE)

VE_f =  VE_model %>%
  filter(group == "VE_BKRT_De") %>% 
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
  filename =  paste0(path, "/VE_plot_De.png"),
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600
)

