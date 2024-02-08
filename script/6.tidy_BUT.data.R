# Script to tidy the BUT data 
rm(list = ls())
library(tidyverse)
library(cowplot)
library(readxl)
library(wesanderson)
source("BUt/R/factor_data.R")


serotype_serostatus_cases = read_excel("BUT/data/raw/cases.xlsx")
age_cases = read_excel("BUT/data/raw/cases.xlsx", sheet = 2)
baseline_seropos = read_excel("BUT/data/raw/baseline_SP.xlsx")
titres = read_excel("BUT/data/raw/titres.xlsx")



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
  ggplot(aes(x = serostatus, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab("Symptomatic cases") + xlab("Serostatus") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_fill_brewer(palette = "Accent") +
  facet_grid(~arm) +
  theme(legend.position = "none")

# Sy VK 
Sy_BVK = cases %>% 
  filter(serotype != "all")

Sy_BVK_plot = Sy_BVK %>%   
  ggplot(aes(x = serostatus, y = Y)) +
  geom_bar(aes(fill = serotype), position = "stack", stat = "identity") +
  ylab(" ") + xlab("Serotype") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = c(0.9,0.8)) +
  facet_grid(~arm) 


grid_plot_data = plot_grid(Sy_BVJ_plot,
                           Sy_BVK_plot, 
                           ncol = 2, 
                           labels = c("a", "b"))

ggsave(grid_plot_data, file = "BUT/output/figures/case_data.jpg",
       height = 50, width = 50, scale =0.68, unit = "cm" )

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
  mutate(age = factor(age, levels = c("2-6yrs", "7-17yrs", "18-59yrs")))
 
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
  filter(Serotype %in% c("D1", "D2")) %>%  # only need D1 and D2 
  group_by(Serostatus) %>% 
  summarise(mean = mean(Titre)) %>% 
  mutate(log_titres = log(mean)) 

