rm(list = ls())

library(readxl)
library(tidyverse)
source("BUT/R/factor_data.R")

serotype_serostatus_cases = read_excel("BUT/data/raw/cases.xlsx")
age_cases = read_excel("BUT/data/raw/cases.xlsx", sheet = 2)
baseline_seropos = read_excel("BUT/data/raw/baseline_SP.xlsx")
titres = read_excel("BUT/data/raw/titres.xlsx")
titres_uncert = read_excel("BUT/data/raw/titres_uncert.xlsx")
trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]

# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[2:4]
serotype_fill = c(scales::brewer_pal(palette = "RdPu")(6)[2:5], "#CCCCCC") 

# create single tibble with all cases 
raw_cases = bind_rows(serotype_serostatus_cases, age_cases)

# factor all variables 
cases = factor_BUT_cases(raw_cases)

# plot case data ---------------------------------------------------------------
theme_set(
  theme_light() +
    theme(
      text = element_text(size = 16),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))

# Sy BVJ 
Sy_BVJ = cases %>% 
  filter(serotype == "all") 

Sy_BVJ_plot  = Sy_BVJ %>%   
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab(" ") + xlab("Trial arm") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_fill_manual(values = age_fill) +
  theme(legend.position = c(0.88,0.8)) +
  facet_grid(~serostatus) 

# Sy VK 
Sy_BVK = cases %>% 
  filter(serotype != "all")

Sy_BVK_plot = Sy_BVK %>%   
  ggplot(aes(x = arm, y = Y)) +
  geom_bar(aes(fill = serotype), position = "stack", stat = "identity") +
  ylab("Symptomatic cases") + xlab("Trial arm") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_fill_manual(values = serotype_fill) +
  theme(legend.position = c(0.89,0.82)) +
  facet_grid(~serostatus) 


grid_plot_data = plot_grid(Sy_BVK_plot,
                           Sy_BVJ_plot,
                           ncol = 2, 
                           labels = c("a", "b"))

ggsave(grid_plot_data, file = "BUT/output/figures/case_data.jpg",
       height = 10, width = 30, scale = .9, unit = "cm" )

# Save factorised data in list for Stan format 
out = list(
  Sy_BVJ = Sy_BVJ,
  Sy_BVK = Sy_BVK
)

saveRDS(out, file = "BUT/data/processed/cases_stan_format.RDS")


# tidy SP data -----------------------------------------------------------------

SP_data = baseline_seropos %>%  
  group_by(age) %>% 
  summarise(SP = sum(SP),
            N = sum(total)) %>% 
  mutate(age = factor(age, levels = c("2-6yrs", "7-17yrs", "18-59yrs"))) %>% 
  arrange(age)
 
saveRDS(SP_data, file = "BUT/data/processed/baseline_SP.RDS") 

# tidy titre data --------------------------------------------------------------
tidy_titres = titres %>%  
  mutate(Serostatus = factor(Serostatus, # make sure SN is first 
                              levels = c("SN", "SP"),
                              labels = c("seronegative", "seropositive"))) %>% 
  filter(Serotype %in% c("D1", "D2")) %>%  # only need D1 and D2 
  pivot_wider(names_from = Serotype, values_from = Titre) %>% 
  select(- c(Arm, Serostatus)) %>% 
  as.matrix()
  
saveRDS(tidy_titres, file = "BUT/data/processed/mu.RDS")   
  
# average titres for lc50 prior 
titres %>% 
  group_by(Serostatus) %>% 
  summarise(mean = mean(Titre)) %>% 
  mutate(log_titres = log(mean)) 
  

# plot initial titres

plot_titre = titres_uncert %>%  
  mutate(Serostatus = factor(Serostatus, # make sure SN is first 
                             levels = c("SN", "SP"),
                             labels = c("seronegative", "seropositive"))) %>% 
  mutate(Serotype = factor(Serotype, levels = c("D1", "D2", "D3", "D4"),
                           labels = c("DENV1", "DENV2","DENV3", "DENV4"))) %>% 
  ggplot(aes(x = Serotype, y = Titre, group = Serostatus)) +
  geom_point(aes(color = Serostatus), position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Serostatus), 
                position = position_dodge(0.5)) +
  theme(legend.position = c(0.85,0.9)) +
  scale_color_manual(values = trial_fill) +
  xlab("Neutralised serotype") +
  ylab("Neutralising titre \ninduced by Butantan-DV")
  

ggsave(plot_titre, file = "BUT/output/figures/plot_titre_B.jpg",
       height = 10, width = 15, unit = "cm" )
