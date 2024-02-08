rm(list = ls())

library(readxl)
library(tidyverse)

pop  = read_excel("BUT/data/baseline_seropositive.xlsx", sheet = 2)
serotype_serostatus_cases = read_excel("BUT/data/cases.xlsx")
age_cases = read_excel("BUT/data/cases.xlsx", sheet = 2)
baseline_seropos = read_excel("BUT/data/baseline_seropositive.xlsx")
titres = read_excel("BUT/data/titres.xlsx")

# average titres for lc50 prior 
titres %>% 
  group_by(Serostatus) %>% 
  summarise(mean = mean(Titre)) %>% 
  mutate(log_titres = log(mean)) 
  
