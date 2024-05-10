# Script to plot the main TAK figure including model fit and the vaccine  
# efficacy estimate 
rm(list=ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(wesanderson)

theme_set(
  theme_light() +
    theme(
      text = element_text(size = 22),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      legend.title = element_blank()
    ))


# colours 
age_fill= scales::brewer_pal(palette = "Blues")(4)[2:4]
serotype_fill = scales::brewer_pal(palette = "RdPu")(6)[2:5] 
trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]


# source files 
file.sources = paste0("TAK/R/", list.files(path = "TAK/R/"))
sapply(file.sources, source)
path = "TAK/output/M2/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
cases =  readRDS("TAK/data/processed/case_data.RDS")

# summarise over iterations 
AR_model = extract_TAK_model_results(AR)  

# Calculate data attack rates 
calc_AR = function(data){
  data %>%  
  mutate(outcome = factor(outcome,  levels = c("symp", "hosp"), 
                          labels = c("symptomatic", "hospitalised"))) %>% 
  mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
         lower = binconf(Y,N, method = "exact")[,2] * 100, 
         upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
  mutate(type = "data")
}

cases_AR = lapply(cases, calc_AR)

# Plot serotype attack rates 

AR_plot_BVKRD =  AR_model %>%
  filter(group == "AR_BVKRD") %>%
  separate(name, into = c("serostatus", "trial", "serotype", "outcome", "month")) %>% 
  filter(outcome == 1) %>% # symphas each time point 
  mutate(month = ifelse(month == 1, 12, 
                        ifelse(month == 2, 18, 
                               ifelse(month == 3,24,
                                      ifelse(month == 4, 36, 
                                             ifelse( month == 5, 48, 54)))))) %>%
  bind_rows(separate(filter(AR_model, group == "AR_BVKHD"), # add hosp which has fewer time points 
                     name, into = c("serostatus", "trial", "serotype", "time"))) %>%  
  mutate(serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
         trial = ifelse(trial == 1, "placebo", "vaccine"),
         outcome = ifelse(outcome == 1, "symptomatic", "hospitalised"),
         serotype = factor(serotype, labels= c("DENV1", "DENV2", "DENV3", "DENV4"))) %>% 
  mutate(month = ifelse(is.na(time), month,
                        ifelse(time == 1, 25,
                               ifelse(time == 2, 36,
                                      ifelse(time == 3, 48, 54))))) %>%
  mutate(outcome = ifelse(is.na(outcome), "hospitalised", outcome)) %>%  
  bind_rows(cases_AR$N_hosp_BVK4, cases_AR$VCD_BVKD) %>%
  mutate(month = ifelse(type == "data" & outcome == "hospitalised" & month == 24, 25, month)) %>% 
  mutate(month = factor(month, levels = c(12, 18, 24, 25, 36, 48, 54),
                        labels = c("1-12", "13-18", "19-24", "1-24", "25-36", "37-48",
                                   "49-54")))  %>% 
  mutate(outcome = factor(outcome, levels = c("symptomatic", "hospitalised"))) %>% 
  ggplot(aes(x = month, y = mean)) +
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
  facet_grid(trial  ~ serostatus +outcome, scales = "free_x" ) +
  labs(x = " ", y = "Attack rate (%)") +
  theme(legend.position = c(0.96,0.8))+
  scale_colour_manual(values = serotype_fill)


# Plot attack rates by age
AR_plot_BVJRD =  AR_model %>%
  filter(group == "AR_BVJRD") %>%
  separate(name, into = c("serostatus", "trial", "age", "outcome", "time")) %>% 
  mutate(serostatus = ifelse(serostatus== 1, "seronegative", "seropositive"),
         trial = ifelse(trial == 1, "placebo", "vaccine"),
         outcome = factor(outcome, 
                          levels = c(1,2),
                          labels = c("symptomatic", "hospitalised")),
         age = factor(age,
                      labels = c("4-5yrs", "6-11yrs", "12-16yrs"),
                      levels = c(1,2,3))) %>% 
  mutate(month = ifelse(time == 1, 12, ifelse(time == 2, 18, 
                                              ifelse(time == 3,24,
                                                     ifelse(time == 4, 36,
                                                            ifelse( time == 5, 48, 54)))))) %>%
  filter(month <= 36) %>% 
  bind_rows(cases_AR$VCD_BVJA, cases_AR$HOSP_BVJA) %>%
  mutate(month = factor(month, levels = c(12, 18, 24, 36),
                        labels = c("1-12", "13-18", "19-24", "25-36"))) %>% 
  ggplot(aes(x = month, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = age,
      group = interaction(type, age)
    ),
    position = position_dodge(width = .5),
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
    position = position_dodge(width = .5),
    width =  0.4,
    linewidth = 1
  ) +
  facet_grid(trial~ serostatus + outcome, scales = "free_x") +
  labs(x = "", y = "Attack rate (%)") +
  theme(legend.position = c(0.96,0.88)) +
  scale_colour_manual(values = age_fill) +
  guides(shape = "none", 
         linetype = "none")

 

# Plot by age and serotype --------------

AR_plot_KJRD =  AR_model %>%
  filter(group == "AR_KJRD") %>%
  separate(name, into = c("serotype", "age", "outcome", "month")) %>% 
  mutate(month = ifelse(month == 1, 12, 24),
         outcome = factor(outcome, levels = c(1,2), 
         labels = c("symptomatic", "hospitalised")),
         serotype = factor(serotype, labels= c("DENV1", "DENV2", "DENV3", "DENV4")),
         age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs"))) %>% 
  bind_rows(cases_AR$VCD_KJ2, cases_AR$HOSP_KJ2) %>%
  mutate(month = factor(month, levels = c(12,24),
                        labels = c("1-12", "13-24"))) %>% 
  ggplot(aes(x = month, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width =.5),
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
    position = position_dodge(width = .5),
    width =  0.4,
    linewidth = 1
  ) +
  facet_grid(outcome~age, scales="free_x") +
  labs(x = "Month", y = "Attack rate (%)") +
  theme(legend.position = "none")+
  scale_colour_manual(values = serotype_fill) 

  
  
# plot grid sep
  g3 = cowplot::plot_grid(
    AR_plot_BVKRD,
    AR_plot_BVJRD,
    AR_plot_KJRD,
    ncol = 1,
    axis = "tblr",
    align = "h",
    label_size = 26,
    rel_heights = c(2,2,1.7),
    labels = c("a", "b", "c")
  ) 
  
  
# output 
ggsave(
  plot = g3,
  filename =  "TAK/output/figures/full_fit_T.png",
  height = 40,
  width = 48,
  units = "cm",
  dpi = 600,
  scale = 0.98
)

