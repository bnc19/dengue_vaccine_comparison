  # Create a figure comparing the parameters of the best fitting models for each vaccine
  library(tidyverse)

  
  theme_set(
    theme_bw() +
      theme(
        text = element_text(size = 12),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0)))
  
  
  model_dep = data.frame(
  vaccine = c("Dengvaxia", "Qdenga", "Butantan-DV"), 
  L = c("serostatus & serotype", "serostatus", NA),
  delta = c("age", "serotype", NA),
  n50 = c("serostatus & serotype_MO\nMU", 
         "serostatus & serotype_MO\nMU",
         "serostatus & serotype_MO\nMU"), 
  beta = c(NA, "age_youngest", NA),
  p = c("age", "age", "age"),
  tau = c("global", "serotype", NA ), 
  eta = c(NA, NA, "age_youngest"),
  rho = c("global", "global", "global"),
  w = c("global", "global", "global"),
  alpha = c("global", "global", NA)
  )

# format data ready to plot 
df = model_dep %>%
  pivot_longer(cols = L:alpha,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  separate(depends, into = c("depends", "extra"), sep = "_") %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "serotype",
        "age",
        "serostatus",
        "serostatus & serotype",
        "global",
        "not included"
      ),
      labels = c(
        "serotype",
        "age group",
        "serostatus",
        "serostatus & serotype",
        "global",
        "not included"  
      ))) %>% 
  mutate(group = 
           ifelse(Parameter %in% c("L", "tau"), "Vaccine efficacy \n(enhancement)",
           ifelse(Parameter %in% c("n50", "beta", "alpha", "w"), "Vaccine efficacy \n(titres)",                  
           ifelse(Parameter %in% c("p", "eta"), "Dengue exposure",
           ifelse(Parameter %in% c("rho", "delta"), "Disease risk",NA))))) %>% 
  mutate(
    Parameter = factor(
      Parameter,
      levels = c(
        "L", "tau", "w", "n50", "beta", "alpha", "p", "eta","rho", "delta"
      ))) 
  

  

mycols = c("#666699","#1C9099",  "#FEEBE2",
           "#FBB4B9", "#CCCCCC", "#FFFFFF")

# plot parameter dependencies 
p1 = df %>%
  ggplot(aes(x = Parameter, y = vaccine)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_fill_manual(values=mycols) +
  geom_text(aes(label=extra),size=3) + # add extra dependency info as letters 
  scale_x_discrete(labels = c(
                              "beta" = expression(beta),
                              "delta" = expression(delta),
                              "alpha" = expression(alpha),
                              "tau" = expression(tau),
                              "rho" = expression(rho),
                              "eta" = expression(eta))) +
  ylab("Vaccine") +
  facet_grid(~group, scales = "free_x", space = "free_x")




ggsave(p1, file = "compare_vaccines/output/compare_model_param.jpg",
       height = 8, width = 18, units="cm")

