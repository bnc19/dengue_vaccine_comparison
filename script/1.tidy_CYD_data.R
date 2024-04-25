# Script to tidy the CYD case data 
rm(list = ls())
library(tidyverse)
library(cowplot)
library(readxl)
library(wesanderson)
source("CYD/R/factor_data.R")

# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = c(scales::brewer_pal(palette = "RdPu")(6)[2:5], "#CCCCCC") 
trial_fill = scales::brewer_pal(palette = "PuBuGn")(3)[2:3]

# read in titres 
raw_titres = read_excel("CYD/data/raw/titres.xlsx")

# read in population
pop = read_excel("CYD/data/raw/population.xlsx")

# read in all sheets of raw case data
# each sheet corresponds to a published table of data 
sheets = excel_sheets("CYD/data/raw/cases_raw.xlsx")

raw_case_list = lapply(sheets, function(X) 
  read_excel("CYD/data/raw/cases_raw.xlsx", sheet = X))

names(raw_case_list) = sheets

# select sheets used for model calibration 
raw_case_list2 = raw_case_list[c(1:10)]

# save each table in global env 
list2env(raw_case_list2, globalenv())


# create single tibble with all cases 
raw_cases = bind_rows(raw_case_list2)  %>%  # collapse age into two categories 
  mutate(age = ifelse(age == "9-14yrs" | age == "9-16yrs" | age == "9-11yrs" | 
                        age == "12-16yrs" | age =="12-14yrs", "9-16yrs",
                      ifelse(age == "2-8yrs" | age == "4-8yrs" | age == "2-5yrs" |
                               age == "6-8yrs", "2-8yrs",
                             ifelse(age=="all", "all", NA)))) %>% 
  # sum over trial arm and new age groups
  group_by(time, outcome, arm, age, serotype, serostatus) %>% 
  summarise(Y = sum(Y),
            N = sum(N)) %>% 
  # remove cumulative Se, placebo, SP because have disaggregated by serotye
  # but don't have Se, vaccine, SP disaggregated by serotype so leave in 
  filter(!(outcome == "Se" & arm == "C" & serotype == "all" & time == "13-72"))

# factor all variables 
cases = factor_CYD_cases(raw_cases, raw = T)

# Add scaled population to cases -----------------------------------------------

# format scaled pop
pop_f = select(factor_CYD_cases(pop, raw = T), -serotype, -time)

# round to integer 
pop_f$N = round(pop_f$N)

# aggregate over trial to get total pop by b/v/j
pop_all = pop_f %>% 
  group_by(arm, age, serostatus) %>% 
  summarise(N = sum(N))

# total pop only in CYD14/15
pop_14_15 = pop_f %>% 
  filter(trial != "CYD23") %>% 
  group_by(arm, age, serostatus) %>% 
  summarise(N = sum(N)) 

pop_symp = pop_14_15 %>% 
  group_by(arm) %>% 
  summarise(N = sum(N)) %>% 
  mutate(age = factor("all"), serostatus = factor("both")) %>% 
  bind_rows(pop_14_15) 

# add CYD14/15 pop to Symp cases
Sy_cases = cases %>% 
  filter(outcome == "Sy") %>% 
  rename(N_cohort = N) %>%  # current N is size of subcohort  
  left_join(pop_symp, by = c("arm", "age", "serostatus"))

# add new pop to cases
cases2 = cases %>% 
  # if Sy then only CYD14 and CYD14
  filter(outcome != "Sy") %>%  
  rename(N_cohort = N) %>%  # current N is size of subcohort  
  left_join(pop_all, by = c("arm", "age", "serostatus")) %>% 
  bind_rows(Sy_cases) # add Symp back in 

# check
check = cases2 %>% 
  mutate(check = round(N_cohort / N * 100,1)) 

range(check$check, na.rm = T)
# all about 10%
 
# censor cases -----------------------------------------------------------------

# only first (hospital) case reported so censor cases after each time interval 

cases_cens_D = cases2 %>%  
  select(- N_cohort) %>% 
  filter(time != "all", outcome == "Ho") %>% # select Ho / Se over time 
  pivot_wider(names_from = time, values_from = c(Y,N)) %>%  # pivot N and Y wider 
  mutate(`N_14-24` =  `N_14-24` -  `Y_1-13`,  # censor previous time interval cases
         `N_25-36` =  `N_25-36` - (`Y_1-13` + `Y_14-24`),
         `N_37-60` =  `N_37-60` - (`Y_1-13` + `Y_14-24` + `Y_25-36`)) %>% 
  pivot_longer(cols = `Y_1-13`:`N_37-60`) %>% # return data to original format 
  separate(name, into = c("X", "time"), sep = "_")  %>% 
  pivot_wider(id_cols = c(outcome, arm, age, serotype, serostatus, time),
              names_from = X, values_from = value) %>%  
  ungroup()


# for cases only reported over full trial duration, take mean population 
mean_pop_over_time = cases_cens_D %>%  
  group_by(arm, age, serostatus) %>% 
  summarise(N = mean(N))

# left join new population to cases across full trial
cases_cens_all = cases2 %>%  
  select(- N_cohort) %>% 
  filter(time == "all") %>% 
  left_join(mean_pop_over_time, 
            by = c("arm", "age", "serostatus")) %>% 
  select(-N.x) %>%  rename("N" = N.y) # after comparing N, remove old N


# for severe cases across time, left join censored population over time 
cases_cens_sev = cases2 %>%  
  select(- N_cohort) %>% 
  filter(time != "all", outcome == "Se") %>% 
  left_join(select(cases_cens_D, arm,age,serostatus, time,N),
            by = c("arm", "age", "serostatus", "time")) %>% 
  select(-N.x) %>%  rename("N" = N.y) # after comparing N, remove old N

# symptomatic cases are only during the first time interval so no censoring
cases_symp = cases2 %>%  
  select(- N_cohort) %>% 
  filter(outcome == "Sy")

# add all cases together 
cases_final = bind_rows(cases_cens_D, cases_cens_all, cases_cens_sev, cases_symp)

# save as RDS
saveRDS(cases_final, file = "CYD/data/processed/cases.RDS")

# check censored cases AR similar to uncensored 
original = cases2 %>% 
  arrange(time, outcome, arm, age, serotype, serostatus) %>% 
  mutate("x" = Y/N * 100) %>%  pull(x)

censored = cases_final %>% 
  arrange(time, outcome, arm, age, serotype, serostatus) %>% 
  mutate("x" = Y/N * 100) %>%  pull(x)

min(original - censored)

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

# Ho BVJD 
Ho_BVJD = cases_final %>% 
  filter(outcome == "Ho" & serotype == "all") 

Ho_BVJD_plot = Ho_BVJD %>% 
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab("Hospitalised cases") + xlab("Month") +
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20)) +
  scale_fill_manual(values = age_fill) +
  facet_grid(serostatus~arm) +
  theme(legend.position = c(0.1,0.8))

# Ho BVKJ
Ho_BVKJ = cases_final %>% 
  filter(outcome == "Ho" & serotype != "all") 
  
Ho_BVKJ_plot = Ho_BVKJ %>% 
  ggplot(aes(x = serotype, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab(" ") + xlab("Serotype") +
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20)) +
  scale_fill_manual(values = age_fill) +
  facet_grid(serostatus~arm) +
  theme(legend.position = "none") 

# Sy BVJ 
Sy_BVJ = cases_final %>% 
  filter(outcome == "Sy", serotype == "all") 

Sy_BVJ_plot  = Sy_BVJ %>%   
  ggplot(aes(x = serostatus, y = Y)) +
  geom_bar(aes(fill = age), position = "stack", stat = "identity") +
  ylab("Symptomatic cases") + xlab("Serostatus") +
  scale_y_continuous(limits = c(0,280), breaks = seq(0,280,50)) +
  scale_fill_manual(values = age_fill) +
  facet_grid(~arm) +
  theme(legend.position = "none")

# Sy VK 
Sy_VK = cases_final %>% 
  filter(outcome == "Sy", serotype != "all")
  
Sy_VK_plot = Sy_VK %>%   
  ggplot(aes(x = serotype, y = Y)) +
  geom_bar(aes(fill = arm), position = "stack", stat = "identity") +
  ylab(" ") + xlab("Serotype") +
  scale_y_continuous(limits = c(0,250), breaks = seq(0,250,50)) +
  scale_fill_manual(values = trial_fill) +
  theme(legend.position = c(0.9,0.8))

# Se VJD 
Se_VJD = cases_final %>% 
  filter(outcome == "Se", time != "all")

Se_VJD_plot = Se_VJD %>% 
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = arm), position = "stack", stat = "identity") +
  facet_grid(~ age) +
  ylab("Severe cases") + xlab("Month") +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10,2)) +
  scale_fill_manual(values = trial_fill) +
  theme(legend.position = "none")

# Se BVKJ 

Se_BVKJ = cases_final %>% 
  filter(outcome == "Se", time == "all") 

Se_BVKJ_plot = Se_BVKJ %>% 
  ggplot(aes(x = serostatus, y = Y)) +
  geom_bar(aes(fill = serotype), position = "stack", stat = "identity") +
  facet_grid(arm ~ age) +
  ylab(" ") + xlab("Month") +
  scale_y_continuous(limits = c(0,60), breaks = seq(0,60,10)) +
  scale_fill_manual(values = serotype_fill) +
  theme(legend.position = c(0.1,0.8))


grid_plot_data = plot_grid(Ho_BVJD_plot,
                           Ho_BVKJ_plot,
                           Sy_BVJ_plot,
                           Sy_VK_plot, 
                           Se_VJD_plot,
                           Se_BVKJ_plot,
                           ncol = 2, 
                           labels = c("a", "b", "c", "d", "e", "f"))

ggsave(grid_plot_data, file = "CYD/output/figures/case_data.jpg",
       height = 50, width = 50, scale =0.68, unit = "cm" )

# Save factorised data in list for Stan format ---------------------------------
out = list(
  Ho_BVJD = Ho_BVJD,
  Ho_BVKJ = Ho_BVKJ,
  Sy_BVJ = Sy_BVJ,
  Sy_VK = Sy_VK,
  Se_VJD = Se_VJD,
  Se_BVKJ = Se_BVKJ
)

out2 = lapply(out, factor_CYD_cases)

saveRDS(out2, file = "CYD/data/processed/cases_stan_format.RDS")

# summarise titres by serotype and serostatus ----------------------------------

# age-group populations are balanced in CYD14 and all over 9 in CYD15 
# so don't need to take a weighted average over age 

# define immuno subset population by trial and serostatus 
pop_sero = data.frame(
trial = c("CYD14", "CYD14", "CYD15", "CYD15"),
serostatus = c("SN", "SP","SN", "SP"),
N = c(423, 900, 251, 1048)
)

tidy_titres = raw_titres %>% 
  pivot_longer(cols = t0:t4, names_to = "time", values_to = "titre") %>%  
  group_by(serotype, serostatus, trial, time) %>% 
  summarise(titre = mean(titre))  %>% # mean over age 
  left_join(pop_sero) %>% # add population 
  group_by(serotype, serostatus, time) %>% 
  summarise(titre = weighted.mean(titre, N)) %>%  # weighted mean by trial size 
  pivot_wider(names_from = serotype, values_from = titre) # pivot wider for output


write.csv(tidy_titres, "CYD/data/processed/tidy_titres.csv")  

# calculate mean titre by serostatus for lc50 prior

tidy_titres %>% 
  pivot_longer(cols = D1:D4, names_to = "serotype", values_to = "titre") %>%  
group_by(serostatus) %>% 
  summarise(mean = log(mean(titre)))


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


plot_titre = raw_titres %>% 
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


ggsave(plot_titre, file = "CYD/output/figures/plot_titre_C.jpg",
       height = 7, width = 24, unit = "cm" )
