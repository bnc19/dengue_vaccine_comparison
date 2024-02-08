# Script to impute D3 and D3 efficacy using M4 parameters
rm(list=ls())

library(readxl)
library(tidyverse)
post = read.csv("BUT/output/M4/posterior_chains.csv")
titres = read_excel("BUT/data/titres.xlsx")

mu = titres %>%  
  pivot_wider(values_from = "Titre", names_from = "Serotype") %>% 
  select(- Serostatus, - Arm) %>% 
  as.matrix()

time = 1:24
B = V  = 2
K = 4
C = 3
T = length(time)
I = 1000
ts = c(-2.12,0.31)
hs = c(1.9,4.3)
hl = 73.5

# sample posterior chains
set.seed(10)
post_samp = as.data.frame(sapply(post, sample, I))

# Select parameters
pi_1 = -log(2) / hs
pi_2 = -log(2) / hl

lc3 = rnorm(I, 5.1150, 0.5)
lc4 = rnorm(I, 4.4794, 0.5)

omega = post_samp$omega
w = post_samp$w.1.
lc_SN = exp(replicate(K, post_samp$lc.1.1.))
MO = cbind(post_samp$lc.2.1., post_samp$lc.2.2., lc3, lc4)
lc_MO = exp(MO)
lc_MU = exp(MO) * exp(-omega) 
nc50 = array(NA, dim = c(C,K,I))
nc50[1,,] = t(lc_SN)
nc50[2,,] = t(lc_MO)
nc50[3,,] = t(lc_MU)


# Arrays
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

# RR
RR_symp = array(NA, dim=c(C,K,T, I))

for(c in 1:C){
  for(k in 1:K){
      for(t in 1:T){
        for(i in 1:I){
          RR_symp[c,k,t,i] =  1 / (1 +  (n_C[c,k,t] / (nc50[c,k,i]))^w[i]) 
        }}}}                                



VacE = array(NA, dim=c(C,K,T,I))

for(c in 1:C){
  for(k in 1:K){
    for(t in 1:T){
      for(i in 1:I){
            VacE[c,k,t,i] = 1 - RR_symp[c,k,t,i]
      }}}}    


VE= VacE %>%  
  reshape2::melt() %>% 
  rename("serostatus" = Var1, "serotype"= Var2, "time"= Var3) %>% 
  mutate(value = value * 100) %>%  
  group_by(serostatus, serotype, time) %>% 
  summarise(
    lower = quantile(value, 0.025, na.rm = T),
    mean = mean(value, na.rm = T),
    upper = quantile(value, 0.975, na.rm = T)
  ) 

VE_plot =  VE %>% 
  ungroup() %>%  
  mutate(serotype = factor(serotype,labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
         serostatus = factor(serostatus, labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = time , y = mean)) +
  geom_line(aes(color = serostatus)) +
  geom_ribbon( aes( ymin = lower, ymax = upper, fill = serostatus), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  facet_grid(serostatus~serotype ) + theme_light() + 
  theme(text = element_text(size=16), legend.position = "none") +
  scale_x_continuous(limits = c(0,24), breaks = seq(0,24,6))

ggsave(
  plot = VE_plot,
  filename =  "BUT/output/figures/imputed_VE_plot.png",
  height = 30,
  width = 40,
  units = "cm",
  dpi = 600,
  scale = 0.8
)
