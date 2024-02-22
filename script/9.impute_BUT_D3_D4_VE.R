# Script to impute D3 and D3 efficacy using M7 parameters (serotype / 
# serostatus lc50), as well as M8 (single lc50 parameter). 
rm(list=ls())
library(tidyverse)

M7_post = read.csv("BUT/output/M7/posterior_chains.csv")
M8_post = read.csv("BUT/output/M8/posterior_chains.csv")
titres = readxl::read_excel("BUT/data/raw/titres.xlsx")

# colours 
serotype_fill = scales::brewer_pal(palette = "RdPu")(6)[2:5] 

mu = titres %>% 
  pivot_wider(names_from = Serotype, values_from = Titre) %>% 
  select(-Serostatus, - Arm) %>% 
  as.matrix()

time = 1:24
B = V  = 2
K = 4
C = 3
T = length(time)
I = 1000
ts = c(-2.1,0.31)
hs = c(1.93,4.34)
hl = 72.14
pi_1 = -log(2) / hs
pi_2 = -log(2) / hl

# Titres
n = array(NA, dim = c(B,K,T))
n_C = array(NA, dim = c(C,K,T))

for(b in 1:B){
  for(k in 1:K){
    for(t in 1:T){
      n[b,k,t] = mu[b,k] * (exp(pi_1[b] * time[t] + pi_2 * ts[b]) + exp(pi_2 * time[t] + pi_1[b] * ts[b])) / (exp(pi_1[b] * ts[b]) + exp(pi_2 * ts[b])) 
    }}}

for(k in 1:K){
  for(t in 1:T){
    n_C[1,k,t] =  n[1,k,t];
    n_C[2,k,t] =  n[2,k,t];
    n_C[3,k,t] =  n[2,k,t];
  }}

# sample posterior chains
set.seed(10)

# M7 --------------------------------------------------------------------------- 
M7_post_sample = as.data.frame(sapply(M7_post, sample, I))
lc3 = rnorm(I, 5.09, 0.5)
lc4 = rnorm(I, 4.32, 0.5)

M7_lc_SN = exp(replicate(K, M7_post_sample$lc.1.1.))
M7_lc_MO = exp(cbind(M7_post_sample$lc.2.1., M7_post_sample$lc.2.2., lc3, lc4))
M7_lc_MU = M7_lc_MO * exp(-M7_post_sample$omega) 
M7_nc50 = array(NA, dim = c(C,K,I))
M7_nc50[1,,] = t(M7_lc_SN)
M7_nc50[2,,] = t(M7_lc_MO)
M7_nc50[3,,] = t(M7_lc_MU)

# RR
M7_RR_symp = array(NA, dim=c(C,K,T, I))

for(c in 1:C){
  for(k in 1:K){
      for(t in 1:T){
        for(i in 1:I){
          M7_RR_symp[c,k,t,i] =  1 / (1 +  (n_C[c,k,t] / (M7_nc50[c,k,i]))^M7_post_sample$w.1.[i]) 
        }}}}                                



M7_VacE = array(NA, dim=c(C,K,T,I))

for(c in 1:C){
  for(k in 1:K){
    for(t in 1:T){
      for(i in 1:I){
        M7_VacE[c,k,t,i] = 1 - M7_RR_symp[c,k,t,i]
      }}}}    


M7_VE= M7_VacE %>%  
  reshape2::melt() %>% 
  rename("serostatus" = Var1, "serotype"= Var2, "time"= Var3) %>% 
  mutate(value = value * 100) %>%  
  group_by(serostatus, serotype, time) %>% 
  summarise(
    lower = quantile(value, 0.025, na.rm = T),
    mean = mean(value, na.rm = T),
    upper = quantile(value, 0.975, na.rm = T)
  ) 

M7_VE_plot =  M7_VE %>% 
  ungroup() %>%  
  mutate(serotype = factor(serotype,
                           labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
         serostatus = factor(serostatus, 
                             labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = time , y = mean)) +
  geom_line(aes(color = serotype)) +
  geom_ribbon( aes( ymin = lower, ymax = upper, fill = serotype), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  facet_grid(serostatus~serotype ) + theme_light() + 
  theme(text = element_text(size=16), legend.position = "none") +
  scale_fill_manual(values = serotype_fill) +
  scale_color_manual(values = serotype_fill) +
  scale_x_continuous(limits = c(0,24), breaks = seq(0,24,6)) +
  scale_y_continuous(limits = c(0,100))
  
ggsave(
  plot = M7_VE_plot,
  filename =  "BUT/output/figures/M7_imputed_VE_plot.png",
  height = 30,
  width = 40,
  units = "cm",
  dpi = 600,
  scale = 0.8
)

  
# M8 --------------------------------------------------------------------------- 
M8_post_sample = as.data.frame(sapply(M8_post, sample, I))
  
M8_lc = exp(M8_post_sample$lc.1.1.)
  
# RR
M8_RR_symp = array(NA, dim=c(C,K,T, I))
  
for(c in 1:C){
  for(k in 1:K){
    for(t in 1:T){
      for(i in 1:I){
        M8_RR_symp[c,k,t,i] =  1 / (1 +  (n_C[c,k,t] / (M8_lc[i]))^M8_post_sample$w.1.[i]) 
        }}}}                                
  
M8_VacE = array(NA, dim=c(C,K,T,I))
  
  for(c in 1:C){
    for(k in 1:K){
      for(t in 1:T){
        for(i in 1:I){
          M8_VacE[c,k,t,i] = 1 - M8_RR_symp[c,k,t,i]
        }}}}    
  
  
  M8_VE = M8_VacE %>%  
    reshape2::melt() %>% 
    rename("serostatus" = Var1, "serotype"= Var2, "time"= Var3) %>% 
    mutate(value = value * 100) %>%  
    group_by(serostatus, serotype, time) %>% 
    summarise(
      lower = quantile(value, 0.025, na.rm = T),
      mean = mean(value, na.rm = T),
      upper = quantile(value, 0.975, na.rm = T)
    ) 
  
  M8_VE_plot =  M8_VE %>% 
    ungroup() %>%  
    mutate(serostatus = ifelse(serostatus == 1, 1, 2)) %>% 
    mutate(serotype = factor(serotype,
                             labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
           serostatus = factor(serostatus, 
                               labels = c("seronegative", "seropositive"))) %>% 
    ggplot(aes(x = time , y = mean)) +
    geom_line(aes(color = serotype)) +
    geom_ribbon( aes( ymin = lower, ymax = upper, fill = serotype), alpha = 0.5) +
    labs(x = "Month", y = "Vaccine efficacy (%)") +
    facet_grid(serostatus~serotype ) + theme_light() + 
    theme(text = element_text(size=16), legend.position = "none") +
    scale_fill_manual(values = serotype_fill) +
    scale_color_manual(values = serotype_fill) +
    scale_x_continuous(limits = c(0,24), breaks = seq(0,24,6))  +
    scale_y_continuous(limits = c(0,100))
  
  ggsave(
    plot = M8_VE_plot,
    filename =  "BUT/output/figures/M8_imputed_VE_plot.png",
    height = 30,
    width = 40,
    units = "cm",
    dpi = 600,
    scale = 0.8
  )
