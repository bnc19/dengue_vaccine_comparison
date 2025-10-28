
rm(list = ls())

library(readxl)
library(tidyverse)
source("R/plotting_functions.R")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.subtitle = element_text(size = 8),
      plot.caption = element_text(size = 8),
      legend.text = element_text(size = 6.5),
      legend.title = element_blank(),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)
    ))

cols = c("#8EAFBF", "#D0B0E0")  

# Plot imputed titres ----------------------------------------------------------

# colours 
model = "M21"
B_n = readRDS(paste0("output/", model, "/n_Bu.RDS"))
Q_n = readRDS(paste0("output/", model, "/n_Q.RDS"))
D_n  = readRDS(paste0("output/", model, "/n_De.RDS"))

D_titres = D_n %>%   
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
  separate(name, into = c(NA, "name"), sep = "\\[") %>% 
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name, into = c("serostatus", "serotype", "time")) %>%  
  mutate(time= as.numeric(time),
         serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
         serotype = paste0("DENV", serotype),
         vaccine = "Dengvaxia")


B_titres = B_n %>%   
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
  separate(name, into = c(NA, "name"), sep = "\\[") %>% 
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name, into = c("serostatus", "serotype", "time")) %>%  
  mutate(time= as.numeric(time),
         serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
         serotype = paste0("DENV", serotype),
         vaccine = "Butantan-DV")


Q_titres = Q_n %>%   
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
  separate(name, into = c(NA, "name"), sep = "\\[") %>% 
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name, into = c("serostatus", "serotype", "time")) %>%  
  mutate(time= as.numeric(time),
         serostatus = ifelse(serostatus == 1, "seronegative", "seropositive"),
         serotype = paste0("DENV", serotype),
         vaccine = "Qdenga") 

titres = bind_rows(Q_titres, B_titres, D_titres) %>%
  mutate(vaccine = factor(vaccine, levels = c("Qdenga", "Dengvaxia", "Butantan-DV"))) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_line(aes(color = serostatus)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = serostatus), alpha = 0.5) +
  labs(x = "Month", y = "Imputed neutralising antibody titre") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols) +
  facet_grid(vaccine ~ serotype, scales = "free_y") +
  theme(legend.position = "0.9,0.72") +
  scale_x_continuous(lim = c(0, 60), breaks = seq(0, 60, 12))

 ggsave(
   titres,
   file = "output/figures/SupFig4.jpg",
   height = 10,
   width = 18,
   unit = "cm"
 )
 
  