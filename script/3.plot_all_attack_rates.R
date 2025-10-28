# Script to plot the model fits of all three vaccines 
rm(list = ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(readxl)
library(patchwork)
library(grid)

# source files 
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
input = "output/M21/"
output = "output/figures"
# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
age_fill_Q = scales::brewer_pal(palette = "Blues")(4)[c(2:4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59", "#111111")
trial_fill = c("#C51B8A", "#99CC99")
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")


theme_set(
  theme_bw()+
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.subtitle = element_text(size = 8),
      plot.caption = element_text(size = 8),
      legend.text = element_text(size = 6.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0,-5,0,-5),
      legend.title = element_blank(),
  plot.margin = margin(t = 0, r = 0.2, b = -0.2, l = 0, "cm")
  )
)

# Data 

AR_Q = readRDS(paste0(input, "AR_Q.RDS"))
VCD_Q = read.csv("data/processed_Q/vcd_data.csv")
hosp_Q = read.csv("data/processed_Q/hosp_data.csv")

AR_D = readRDS(paste0(input, "AR_De.RDS"))
VCD_D =  readRDS("data/processed_De/cases_stan_format.RDS")

AR_B = readRDS(paste0(input, "AR_Bu.RDS"))
cases_B = readRDS("data/processed_Bu/cases_stan_format.RDS")


########################## PLOT QDENGA ATTACK RATES ############################
# Tidy data 
VCD_Q = factor_VCD_Q(VCD_Q)
hosp_Q = factor_VCD_Q(hosp_Q)

AR_model_Q = extract_model_results_Q(AR_Q)  

# calculate attack rates from data 
AR_data_Q = calc_attack_rates_Q(VCD_Q)
HR_data_Q = calc_hosp_rates_Q(VCD=VCD_Q, hosp=hosp_Q)
data_Q = bind_rows(AR_data_Q, HR_data_Q)


# Plot AR across time 
AR_VRD_Q = data_Q %>%
  filter(serostatus == "both",
         age == "all",
         trial != "both",
         serotype == "all") %>%
  mutate(trial = ifelse(trial == "TAK-003", "vaccine", trial)) %>%
  select(-serotype,-age,-serostatus,-X, -year) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  ))

time_plot_Q =  AR_model_Q %>%
  filter(group == "AR_VRD_Q") %>%
  separate(name, into = c("trial", "outcome", "month")) %>%
  mutate(
    trial = factor(trial, labels = c("placebo", "vaccine")),
    month = ifelse(month == 1, 12,
                   ifelse(
                     month == 2, 18,
                     ifelse(month == 3, 24,
                            ifelse(month == 4, 36,
                                   ifelse(month == 5, 48, 54))))),
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))) %>%
  bind_rows(AR_VRD_Q) %>%
  mutate(month = as.numeric(month)) %>%
  ggplot(aes(x = month, y = mean)) +
  geom_point(aes(color = trial, shape = type), position = position_dodge(width = 5), size = 1) +
  geom_errorbar( aes(ymin = lower ,ymax = upper ,color = trial, linetype = type),
                 position = position_dodge(width = 5),width =  0.4,linewidth = 0.5) +
  labs(x = "Month", y = "Qdenga attack rate (%)") +
  scale_x_continuous(breaks = c(12,18,24,36,48,54)) +
  scale_color_manual(values = trial_fill) +
  facet_wrap(~ outcome, ncol = 1, scales = "free_y") +
  theme(legend.position = c(0.85,0.86))


# Plot AR by age, trial, and serostatus 

AR_BVJR_Q = data_Q %>%
  filter(serotype == "all" |
           age != "all", serostatus != "both", month == 1000) %>%
  mutate(trial = ifelse(trial == "TAK-003", "vaccine", trial)) %>%
  select(-month,-X, -year,-serotype) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  ))

age_plot_Q =  AR_model_Q %>%
  filter(group == "AR_BVJR_Q") %>%
  separate(name, into = c("serostatus", "trial", "age", "outcome")) %>% 
  mutate(outcome = factor(outcome,labels = c("symptomatic", "hospitalised"))) %>% 
  mutate(serostatus = factor(serostatus, labels=c("seronegative", "seropositive")),
         trial = factor(trial, labels=c("placebo", "vaccine")),
         age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs"),
                      levels = c(1,2,3))) %>% 
  bind_rows(AR_BVJR_Q) %>%
  mutate(serostatus = factor(serostatus, labels=c("SN", "SP")),
         trial = factor(trial, labels=c("P", "V"))) %>% 
  unite(c(trial, serostatus), sep = " ", col = "x") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(
    aes( shape = type,  color = age, group = interaction(type, age)),
    position = position_dodge(width = 0.5), size = 1 ) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper, group = interaction(type, age),
        linetype = type,color = age),
    position = position_dodge(width =  0.5),
    width =  0.4,
    linewidth = 0.5
  ) +
  scale_color_manual(values = age_fill_Q) +
  guides(shape = "none",
         linetype = "none")+
  facet_wrap(~ outcome, ncol=1, scales = "free_y") +
  theme(legend.position = c(0.85,0.9),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())


# Plot AR by serotype, trial, and serostatus 

AR_BVKR_Q = data_Q %>%
  filter(serotype != "all" |
           age == "all", serostatus != "both", month == 1000) %>%
  mutate(trial = ifelse(trial == "TAK-003", "vaccine", trial)) %>%
  select(-month,-X, -year,-age) %>%
  mutate(outcome = factor(
    outcome,
    levels = c("VCD", "hosp"),
    labels = c("symptomatic", "hospitalised")
  )) %>%
  mutate(serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4"))) 

serotype_plot_Q = AR_model_Q %>%
  filter(group == "AR_BVKR_Q") %>%
  separate(name, into = c("serostatus", "trial", "serotype", "outcome")) %>%
  mutate(outcome = factor(outcome, labels = c("symptomatic", "hospitalised"))) %>%
  mutate(
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
    trial = factor(trial, labels = c("placebo", "vaccine")),
    serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4"))) %>%
  bind_rows(AR_BVKR_Q) %>% 
  mutate(serostatus = factor(serostatus, labels=c("SN", "SP")),
         trial = factor(trial, labels=c("P", "V"))) %>% 
  unite(c(trial, serostatus), sep = " ", col = "x") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.7), size = 1) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serotype),
      linetype = type,
      color = serotype
    ),
    position = position_dodge(width =  0.7),
    width =  0.4,
    linewidth = 0.5
  ) +
  scale_color_manual(values = serotype_fill) +
  facet_wrap(~ outcome, ncol=1, scales = "free_y") +
  theme(legend.position = c(0.85,0.86),
        axis.title.y = element_blank(),
        axis.title.x = element_blank() )

AR_plot_Q = time_plot_Q  + age_plot_Q + serotype_plot_Q  + 
  plot_layout(ncol = 3, guides = "collect") & 
  theme(legend.position = "right")


########################## PLOT DENGVAXIA ATTACK RATES #########################

# plot attack rates 
AR_data_D = lapply(VCD_D, calc_attack_rates_De)
AR_model_D = extract_model_results_De(AR_D)  

# plot symp attack rate by serotype and trial arm 
Sy_AR_plot_VK_D = AR_model_D %>%
  filter(group == paste0("V_AR_VK_De")) %>%
  separate(name, into = c("arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4"))) %>% 
  bind_rows(AR_data_D$Sy_VK %>%  mutate(serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4")))) %>%
  mutate( arm = factor(arm, labels=c("P", "V"))) %>% 
  ggplot(aes(x = arm, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width = 0.5),
    size = 1
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
    linewidth = 0.5
  ) +
  labs(x = " ", y = "Dengvaxia attack rate (%)") +
  scale_color_manual(values = serotype_fill) +
  guides(shape = "none", linetype = "none")+
  theme(legend.position = c(0.85,0.72),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) 

# plot symp attack rate by age, serostatus and trial arm
Sy_AR_plot_BVJ_D = AR_model_D %>%
  filter(group == "V_AR_BVJ_De") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data_D$Sy_BVJ) %>%
  mutate(serostatus = factor(serostatus, labels=c("SN", "SP")),
         arm = factor(arm, labels=c("P", "V"))) %>% 
  unite(c(arm, serostatus), sep = " ", col = "x") %>%
  mutate(outcome = "symptomatic") %>% 
  ggplot(aes(x = x, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = age,
      group = interaction(type, age)),
    position = position_dodge(width = 0.5),
    size = 1) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, age),
      linetype = type,
      color = age),
    position = position_dodge(width =  0.5),
    width =  0.4,
    linewidth = 0.5) +
  labs(x = " ", y = "Dengvaxia attack rate (%)") +
  scale_color_manual(values = age_fill) +
  theme(legend.position = c(0.85,0.72),
        axis.title.y = element_blank(),
        axis.title.x = element_blank())

# plot hosp attack rate by age, serostatus, serotype and trial arm
H_AR_plot_BVKJ_D = AR_model_D %>%
  filter(group == "H_AR_BVKJ_De") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data_D$Ho_BVKJ%>%  mutate(serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4")))) %>%
  mutate(serostatus = factor(serostatus, labels=c("SN", "SP")),
         arm = factor(arm, labels=c("P", "V"))) %>% 
  unite(c(arm, serostatus), sep = " ", col = "x") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width = .8),
    size = 1
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
    linewidth = 0.5
  ) +
  scale_color_manual(values = serotype_fill) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  facet_wrap(~ age, ncol = 1,  scales = "free_y")

# plot hosp attack rate by age, serostatus, trial arm and time

AR_data_D$Ho_BVJD = AR_data_D$Ho_BVJD  %>%  
  mutate(time =  factor(time, labels = c("13", "24", "36", "60")))
 
H_AR_plot_BVJD_D =  AR_model_D %>%
  filter(group == "H_AR_BVJD_De") %>%
  separate(name, into = c("serostatus", "arm", "age", "time")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs")),
         time = factor(time, labels = c("13", "24", "36", "60"))) %>%
  bind_rows(AR_data_D$Ho_BVJD ) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = age,
      group = interaction(type, age)
    ),
    position = position_dodge(width = 0.6),
    size = 1
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
    linewidth = 0.5
  ) +
  facet_grid(serostatus ~ arm,  scales = "free_y") +
  labs(x = "Month", y = "Dengvaxia hospitalisation\nattack rate (%)") +
  scale_color_manual(values = age_fill) +
  theme(legend.position = "none")  

# Combine all plots
shared_y_label = wrap_elements(
  grid::textGrob("Dengvaxia symptomatic\nattack rate (%)", rot = 90, gp = gpar(fontsize = 8))
)

Sy_plots = Sy_AR_plot_BVJ_D / Sy_AR_plot_VK_D

AR_plot_D = (shared_y_label | (Sy_plots | H_AR_plot_BVJD_D | H_AR_plot_BVKJ_D) + plot_layout(widths = c(0.65, 1, 0.65))) +  
  plot_layout(ncol = 2, widths = c(0.055, 1), guides = "collect") &  
  theme(legend.position = "right") 

######################## PLOT BUTANTAN-DV ATTACK RATES #########################

BVK_cases_B = cases_B$Sy_BVK
BVJ_cases_B = cases_B$Sy_BVJ

# plot attack rates 
# add aggregated populations to data and calculate attack rates

AR_age_data_B = calc_BUT_attack_rates(BVJ_cases_B)
AR_serotype_data_B = calc_BUT_attack_rates(BVK_cases_B)
AR_model_B = extract_BUT_model_results(AR_B)

# plot serotype serostatus attack rate 
AR_plot_BVK_B = AR_model_B %>%
  filter(group == "AR_BVK_Bu") %>%
  separate(name, into = c("serostatus", "arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         serotype = factor(serotype, labels = c(paste0("D", 1:2)))) %>% 
  bind_rows(AR_serotype_data_B %>%  mutate(serotype = factor(serotype, labels = c("D1", "D2")))) %>%
  mutate(arm = factor(arm, labels=c("P", "V"))) %>% 
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = serotype,
                 group = interaction(type, serotype)),
             position = position_dodge(width = 0.5),size = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.5),width =  0.4,     linewidth = 0.5) +
  labs(x = " ", y = "Butantan-DV \nsymptomatic \nattack rate (%)") +
  facet_wrap(~ serostatus,  scales = "free_y") + 
  scale_color_manual(values = serotype_fill) 

# plot symp attack rate by age and trial arm
AR_plot_BVJ_B = AR_model_B %>%
  filter(group == "AR_BVJ_Bu") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"))) %>% 
  bind_rows(AR_age_data_B) %>%
  mutate(arm = factor(arm, labels=c("P", "V"))) %>% 
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = age,
                 group = interaction(type, age)),
             position = position_dodge(width = 0.5), size = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4,    linewidth = 0.5) +
  guides(shape = "none",linetype = "none") +
  facet_wrap(~ serostatus,  scales = "free_y")  + 
  scale_color_manual(values = age_fill_Q) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank())

AR_plot_B = AR_plot_BVK_B + AR_plot_BVJ_B  + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "right")

############################ COMBINE ALL PLOTS #################################
grid_plot = (wrap_elements(AR_plot_Q) / wrap_elements(AR_plot_D) / wrap_elements(AR_plot_B)) +
  plot_annotation(tag_levels = list(c("a", "b", "c"))) +  
  plot_layout(heights = c(2, 2, 1)) & 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        plot.tag = element_text(size = 12))


ggsave(
  plot = grid_plot,
  filename =  paste0(output, "/Fig1.pdf"),
  height = 20,
  width = 18,
  units = "cm",
  dpi = 600,
)

