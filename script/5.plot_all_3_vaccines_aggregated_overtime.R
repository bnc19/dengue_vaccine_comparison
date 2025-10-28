# Plot 2 year VE for all vaccines by serostatus and serotype (but not age)
rm(list = ls())

source("R/plotting_functions.R")
library(tidyverse)

n = 1000

set.seed(42)

# colours 

palette_3 = c( "#D55E00", "#CC79A7", "#8EAFBF", "#A9D6C6")  

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
      legend.title = element_blank()
    ))

input = "output/M21/"

VE_Q = select(readRDS(paste0(input, "VE_Q.RDS")),contains("VE_Q"))  # without age 
VE_D = select(readRDS(paste0(input, "VE_De.RDS")),contains("VE_D")) 
VE_B = readRDS(paste0(input, "VE_Bu.RDS"))

# sample 1000 times from posterior distribution 
VE_Q = as.data.frame(sapply(VE_Q, sample, n))
VE_D = as.data.frame(sapply(VE_D, sample, n))
VE_B = as.data.frame(sapply(VE_B, sample, n))

# Tidy Qdenga data ----
tidy_VE_Q =  VE_Q %>%
  as.data.frame() %>%
  mutate(ni = row_number()) %>%
  pivot_longer(cols = -ni) %>%   # Pivot iterations to rows
  group_by(name) %>% 
  mutate(value = value * 100) %>%
  separate(name, into = c("group", "name"), sep = "\\[") %>%  # rename
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name,
           into = c("serostatus", "serotype", "age", "outcome", "month"))  
  
VE_Q_post = tidy_VE_Q  %>%  
  filter(month <= 24) %>% 
  group_by(serostatus, serotype, outcome, age,ni) %>%
  summarise(value = mean(value)) %>%
  group_by(serostatus, serotype, age, outcome) %>% 
  summarise( # mean and uncertainty of outcome posterior 
    lower = quantile(value, 0.025, na.rm = T),
    mean = mean(value, na.rm = T),
    upper = quantile(value, 0.975, na.rm = T)
  )

VE_Q_time = tidy_VE_Q  %>%  
  filter(month <= 24) %>% 
  group_by(serostatus, serotype, age, outcome) %>% 
  summarise( # mean and uncertainty of the posterior 
    lower_t = quantile(value, 0.025, na.rm = T),
    mean_t = mean(value, na.rm = T),
    upper_t = quantile(value, 0.975, na.rm = T)
  )

VE_Q_plot = VE_Q_post %>% 
  left_join(VE_Q_time) %>% 
  filter(age!=3) %>%   # age 2 and 3 are same 
  mutate( # factor for plotting 
    serotype = factor(paste0("DENV", serotype)),
    age = factor(age, levels = 1:2, labels = c("4-5yrs", "6-16yrs")),
    serostatus = factor(
      serostatus, levels = 1:3, 
      labels = c("seronegative", "monotypic", "multitypic")
    )
  ) %>% 
  filter(outcome == 1)  

# Tidy Dengvaxia data ----

tidy_VE_D =  VE_D %>%
  as.data.frame() %>%
  mutate(ni = row_number()) %>%
  pivot_longer(cols = -ni) %>%   # Pivot iterations to rows
  group_by(name) %>% 
  mutate(value = value * 100) %>%
  separate(name, into = c("group", "name"), sep = "\\[") %>%  # rename
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name,
           into = c("serostatus", "serotype", "age", "outcome", "month"))

VE_D_post = tidy_VE_D  %>%  
  filter(month <= 24, age == 1) %>%  # no age efficacy
  group_by(serostatus, serotype, outcome, age,ni) %>%
  summarise(value = mean(value)) %>%
  group_by(serostatus, serotype, age, outcome) %>% 
  summarise( # mean and uncertainty of outcome posterior 
    lower = quantile(value, 0.025, na.rm = T),
    mean = mean(value, na.rm = T),
    upper = quantile(value, 0.975, na.rm = T)
  )

VE_D_time = tidy_VE_D  %>%  
  filter(month <= 24, age == 1) %>% 
  group_by(serostatus, serotype, age, outcome) %>% 
  summarise( # mean and uncertainty of the posterior 
    lower_t = quantile(value, 0.025, na.rm = T),
    mean_t = mean(value, na.rm = T),
    upper_t = quantile(value, 0.975, na.rm = T)
  )

VE_D_plot = VE_D_post %>% 
  left_join(VE_D_time) %>% 
  mutate( # factor for plotting 
    serotype = factor(paste0("DENV", serotype)),
    serostatus = factor(
      serostatus, levels = 1:3, 
      labels = c("seronegative", "monotypic", "multitypic")
    )
  ) %>% 
  filter(outcome ==1) 


# Tidy Butantan-DV data ----

tidy_VE_B =  VE_B %>%
  as.data.frame() %>%
  mutate(ni = row_number()) %>%
  pivot_longer(cols = -ni) %>%   # Pivot iterations to rows
  group_by(name) %>% 
  mutate(value = value * 100) %>%
  separate(name, into = c("group", "name"), sep = "\\[") %>%  # rename
  separate(name, into = c("name", NA), sep = "\\]") %>% 
  separate(name,
           into = c("serostatus", "serotype", "age", "month"))
  
VE_B_post = tidy_VE_B  %>%  
  group_by(serostatus, serotype, ni) %>%
  summarise(value = mean(value)) %>%
  group_by(serostatus, serotype) %>% 
  summarise( # mean and uncertainty of the posterior 
    lower = quantile(value, 0.025, na.rm = T),
    mean = mean(value, na.rm = T),
    upper = quantile(value, 0.975, na.rm = T)
  )

VE_B_time = tidy_VE_B  %>%  
  group_by(serostatus, serotype) %>% 
  summarise( # mean and uncertainty of the posterior 
    lower_t = quantile(value, 0.025, na.rm = T),
    mean_t = mean(value, na.rm = T),
    upper_t = quantile(value, 0.975, na.rm = T)
  )


VE_B_plot = VE_B_post %>% 
  left_join(VE_B_time) %>% 
  mutate( # factor for plotting 
    serotype = factor(paste0("DENV", serotype)),
    serostatus = factor(
      serostatus, levels = 1:3, 
      labels = c("seronegative", "monotypic", "multitypic")
    )
  ) 


VE_Q_plot$vaccine = "Qdenga"
VE_D_plot$vaccine = "Dengvaxia"
VE_B_plot$vaccine = "Butantan-DV"
VE_B_plot$age = "2-59yrs"
VE_D_plot$age = "2-16yrs"

year2_comp_data = VE_Q_plot %>%
  filter(age != 3) %>% # Qdenga VE is same in age groups 2 and 3
  bind_rows(VE_D_plot, VE_B_plot) %>%
  mutate(vaccine = factor(vaccine, levels = c("Qdenga", "Dengvaxia", "Butantan-DV"))) %>% 
  mutate(age = paste0(vaccine, "\n(", age, ")"))

year2_comp = year2_comp_data %>%
  ggplot(aes(x =  vaccine, y = mean, group = age)) +
  geom_errorbar(aes(ymin = lower_t, ymax = upper_t, color = age), 
                size = 2, width = 0, position = position_dodge(0.6), alpha = 0.5) +
  geom_errorbar(data = subset(year2_comp_data, grepl("Qde", age)), 
                aes(ymin = lower, ymax = upper, color = age), 
                width = 0.4, size =0.2, position = position_dodge(0.6)) +
  geom_errorbar(data = subset(year2_comp_data, age == "Butantan-DV\n(2-59yrs)"), 
                aes(ymin = lower, ymax = upper), 
                color =  "#D55E00", width = 0.2, size =0.2, position = position_dodge(0.6)) + 
  geom_errorbar(data = subset(year2_comp_data, age ==  "Dengvaxia\n(2-16yrs)"), 
                aes(ymin = lower, ymax = upper), 
                color =  "#CC79A7", width = 0.2, size =0.2, position = position_dodge(0.6)) + 
  geom_point(aes(color = age), size = 0.2, position = position_dodge(0.6)) +
  facet_grid(serostatus ~ serotype) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(colour = "age group", x = "", y = "Vaccine efficacy (%)") +
  scale_color_manual(values = palette_3) +
  scale_fill_manual(values = palette_3) +
  ggtitle("Cumulative vaccine efficacy across months 1 to 24")

year2_comp

ggsave(
  plot = year2_comp,
  filename = "output/figures/Fig2.pdf",
  height = 9,
  width = 10,
  units = "cm",
  dpi = 600
)
