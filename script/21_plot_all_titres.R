
rm(list = ls())

library(readxl)
library(tidyverse)
source("compare_vaccines/R/plotting_functions.R")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))

B_titres_uncert = read_excel("BUT/data/raw/titres_uncert.xlsx")
C_raw_titres = read_excel("CYD/data/raw/titres.xlsx")
T_titres  = read.csv(file = "TAK/data/processed/antibody_data.csv")[,-1]



trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]


# plot B titres

B_plot_titre = B_titres_uncert %>%  
  mutate(Serostatus = factor(Serostatus, # make sure SN is first 
                             levels = c("SN", "SP"),
                             labels = c("seronegative", "seropositive"))) %>% 
  mutate(Serotype = factor(Serotype, levels = c("D1", "D2", "D3", "D4"),
                           labels = c("DENV1", "DENV2","DENV3", "DENV4"))) %>% 
  ggplot(aes(x = Serotype, y = Titre, group = Serostatus)) +
  geom_col(aes(fill = Serostatus), position = position_dodge(1)) +
  theme(legend.position = c("none")) +
  scale_fill_manual(values = trial_fill) +
  ylab("Observed \nneutralising \nantibody titre") +
  ggtitle("Butantan-DV")



# plot C titres


# age-group populations are balanced in CYD14 and all over 9 in CYD15 
# so don't need to take a weighted average over age 

## plot titres

pop_sero_age = data.frame(
  trial = c("CYD14", "CYD14", 
            "CYD14", "CYD14",
            "CYD15", "CYD15"),
  serostatus = c("SN", "SP",
                 "SN", "SP",
                 "SN", "SP"),
  N = c(round(423/2),round(423/2),
        round(900/2),round(900/2),
        251, 1048),
  age = c("2-8yrs", "2-8yrs", 
          "9-16yrs", "9-16yrs",
          "9-16yrs", "9-16yrs")
)


C_plot_titre = C_raw_titres %>% 
  mutate(age = ifelse(age == "9-14yrs", "9-16yrs", age)) %>% 
  pivot_longer(cols = t0:t4, names_to = "time", values_to = "titre") %>% 
  select(-N) %>% 
  left_join(pop_sero_age) %>% 
  mutate(serostatus = factor(serostatus, # make sure SN is first 
                             levels = c("SN", "SP"),
                             labels = c("seronegative", "seropositive"))) %>% 
  mutate(serotype = factor(serotype, levels = c("D1", "D2", "D3", "D4"),
                           labels = c("DENV1", "DENV2","DENV3", "DENV4"))) %>% 
  mutate(time = ifelse(time == "t0", 0, 
                       ifelse(time == "t1", 12,
                              ifelse(time == "t3", 24, 36)))) %>% 
  group_by(serotype, serostatus, age, time) %>% 
  summarise(titre = weighted.mean(titre, N)) %>%  # weighted mean by trial size 
  ggplot(aes(x = time, y = titre)) + 
  geom_point(aes(color = serostatus, shape = age)) +
  geom_line(aes(color = serostatus, linetype = age)) +
  facet_grid(~serotype) +
  ylab("Observed \nneutralising \nantibody titre") + xlab("Month") +
  scale_color_manual(values = trial_fill) +
  theme(legend.position = "top")+
  ggtitle("Dengvaxia")


# T plot titres 

# plot observed titres by trial, serostatus and serotype 
T_plot_titres = T_titres %>%
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
  ylab("Observed \nneutralising \nantibody titre")  + xlab("Month")+
  scale_x_continuous(limits = c(0,54), breaks = seq(0,54,12)) +
  facet_wrap(~ serotype, ncol=4) +
  theme(legend.position = "none")+
  scale_colour_manual(values = trial_fill) +
  ggtitle("Qdenga")

# combine titres

titre_grid = cowplot::plot_grid(
  C_plot_titre,
  T_plot_titres,
  rel_heights = c(1.2, 1, 1),
  B_plot_titre,
  ncol = 1,
  labels = c("a", "b", "c")
)


ggsave(
  titre_grid,
  file = "compare_vaccines/output/plot_titres.jpg",
  height = 16,
  width = 18,
  unit = "cm"
)



# Plot imputed titres ----------------------------------------------------------

# colours 

BUT_n = readRDS("BUT/output/M7/n.RDS")
TAK_n = readRDS("TAK/output/M32/n.RDS")
CYD_n  = readRDS("CYD/output/M12/n.RDS")

CYD_titres = CYD_n %>%   
  as.data.frame() %>%
  mutate(ni = row_number()) %>%
  pivot_longer(cols = -ni) %>%
  group_by(name) %>%
  mutate(value = (value)) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>%  
  filter(grepl("n\\[", name)) %>% 
  separate(name, into = c(NA, "name"), sep = "\\[") %>% 
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name, into = c("serostatus", "serotype", "time")) %>%  
  mutate(time= as.numeric(time),
         serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
         serotype = paste0("DENV", serotype),
         vaccine = "Dengvaxia")


BUT_titres = BUT_n %>%   
  as.data.frame() %>%
  mutate(ni = row_number()) %>%
  pivot_longer(cols = -ni) %>%
  group_by(name) %>%
  mutate(value = (value)) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>%  
  filter(grepl("n\\[", name)) %>% 
  separate(name, into = c(NA, "name"), sep = "\\[") %>% 
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name, into = c("serostatus", "serotype", "time")) %>%  
  mutate(time= as.numeric(time),
         serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
         serotype = paste0("DENV", serotype),
         vaccine = "Butantan-DV")


TAK_titres = TAK_n %>%   
  as.data.frame() %>%
  mutate(ni = row_number()) %>%
  pivot_longer(cols = -ni) %>%
  group_by(name) %>%
  mutate(value = (value)) %>%
  summarise(
    lower = quantile(value, 0.025),
    mean = mean(value),
    upper = quantile(value, 0.975)
  ) %>%  
  filter(grepl("n\\[", name)) %>% 
  separate(name, into = c(NA, "name"), sep = "\\[") %>% 
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name, into = c("serostatus", "serotype", "time")) %>%  
  mutate(time= as.numeric(time),
         serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
         serotype = paste0("DENV", serotype),
         vaccine = "Qdenga") 

TAK_titres_plot = TAK_titres %>%  
  ggplot(aes(x = time, y = mean)) +
  geom_line(aes(color = serostatus)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = serostatus), alpha = 0.5) +
  labs(x = "Month", y = "Imputed \nneutralising \nantibody titre") +
  scale_colour_manual(values = trial_fill) +
  scale_fill_manual(values = trial_fill) +
  facet_wrap(~serotype,  ncol = 4) +
  theme(legend.position = "none")+
  ggtitle("Qdenga")+
  scale_x_continuous(lim = c(0,54), breaks = seq(0,54,12))

CYD_titres_plot = CYD_titres %>%  
  ggplot(aes(x = time, y = mean)) +
  geom_line(aes(color = serostatus)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = serostatus), alpha = 0.5) +
  labs(x = "Month",  y = "Imputed \nneutralising \nantibody titre") +
  scale_colour_manual(values = trial_fill) +
  scale_fill_manual(values = trial_fill) +
  facet_wrap(~serotype,  ncol = 4)  +
  theme(legend.position = c(0.9,0.72))+
  ggtitle("Dengvaxia")+
  scale_x_continuous(lim = c(0,60), breaks = seq(0,60,12))

 BUT_titres_plot = BUT_titres %>%  
  ggplot(aes(x = time, y = mean)) +
  geom_line(aes(color = serostatus)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = serostatus), alpha = 0.5) +
  labs(x = "Month",  y = "Imputed \nneutralising \nantibody titre") +
  scale_colour_manual(values = trial_fill) +
  scale_fill_manual(values = trial_fill) +
  facet_wrap(~serotype,  ncol = 4) +
  theme(legend.position = "none") +
  ggtitle("Butantan-DV") +
  scale_x_continuous(lim = c(0,24), breaks = c(0,12,24))

  # combine titres
 
 
 imp_titre_grid = cowplot::plot_grid(
   CYD_titres_plot,
   TAK_titres_plot,
   cowplot::plot_grid(BUT_titres_plot, NULL, rel_widths = c(1, .85)),
   ncol = 1 ,
   labels = c("a", "b", "c")
 )
 
 
 ggsave(
   imp_titre_grid,
   file = "compare_vaccines/output/imputed_plot_titres.jpg",
   height = 16,
   width = 18,
   unit = "cm"
 )
 
  