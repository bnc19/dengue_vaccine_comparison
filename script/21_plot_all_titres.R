
rm(list = ls())

library(readxl)
library(tidyverse)

theme_set(
  theme_light() +
    theme(
      text = element_text(size = 16),
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
  # geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Serostatus), 
  #               position = position_dodge(0.5)) +
  theme(legend.position = c("none")) +
  scale_fill_manual(values = trial_fill) +
  xlab("Neutralised serotype") +
  ylab("Neutralising titre \ninduced by Butantan-DV")



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
  ylab("Neutralising titre \ninduced by Dengvaxia") + xlab("Month PD3") +
  scale_color_manual(values = trial_fill) +
  theme(legend.position = c(0.92,0.7))




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
  # geom_errorbar(aes(ymin =(lower), ymax = (upper),
  #                   color = serostatus)) +
  ylab("Neutralising titre \ninduced by Qdenga") + xlab("Month PD2")+
  scale_x_continuous(limits = c(0,54), breaks = seq(0,54,12)) +
  facet_wrap(~ serotype, ncol=4) +
  theme(legend.position = "none")+
  scale_colour_manual(values = trial_fill)

# combine titres

titre_grid = cowplot::plot_grid(C_plot_titre, T_plot_titres,
                   B_plot_titre, ncol = 1, labels = c("a", "b", "c"))


ggsave(titre_grid, file = "compare_vaccines/output/plot_titres.jpg",
       height = 20, width = 30, unit = "cm" )
