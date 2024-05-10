# Script to run model comparison 
rm(list=ls())

# set up 
library(tidyverse)
library(cowplot)

n_param = c(15,14,13, 13, 13,
            14, 14, 11)
# read waic
index_files = which(grepl("M", list.files(path = "BUT/output/")))
source_files = paste0("BUT/output/", list.files(path = "BUT/output/")[index_files], "/WAIC.RDS")
WAIC = (lapply(source_files, readRDS))
n1 = gsub("/WAIC.RDS", "", gsub("BUT/output/M", "", source_files))
names(WAIC) = n1
WAIC2 = WAIC[as.character(sort(as.numeric(n1)))]

# compare WAIC
comp_WAIC = loo::loo_compare(WAIC2)

# M7 is best 
best_model = as.numeric(rownames(comp_WAIC)[1])


# extract waic and elpd 
waic_df = comp_WAIC %>%  
  as.data.frame() %>% 
  rownames_to_column(var = "model") %>% 
  select(model, elpd_diff, waic) %>% 
  mutate(model = as.numeric(model))

# get log lik

ll_source = paste0("BUT/output/", list.files(path = "BUT/output/")[index_files], "/posterior.csv")
post = bind_rows(lapply(ll_source, read.csv))

ll = post %>%  
  filter(variable == "ll") %>% 
  mutate(model = as.numeric(n1)) %>% 
  select(mean, q5, q95, model) %>% 
  arrange(model) %>% 
  mutate(param = n_param) %>% 
  left_join(waic_df)

# models 
models = as.numeric(paste0(1:length(n_param)))
parameters= c("L", "beta", "lc", "omega", "eta")

# create empty matrix 
m = matrix(nrow = length(models), ncol = length(parameters))

# format 
d = m %>%
  as.data.frame() 

rownames(d) = models
colnames(d) = parameters

# L
d[1,1] = "serostatus"

# beta
d[c(1:2,6),2] = "age_2-6"

# lc 
d[1:3,3] = "serostatus & serotype_MO"
d[4:7,3] = "serostatus & serotype_MO\nMU"
d[5,3] = "serostatus & serotype"
d[8,3] = "global"

# omega
d[4:7,4] = "serostatus_MU" 
d[5,4] = "serostatus_SN\nMU" 

# eta
d[7:8,5] = "age_2-6"

# add log lik to dependency matrix 
d$model = models

# format data ready to plot 
df = left_join(d, ll) %>%
  pivot_longer(cols = L:eta,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  separate(depends, into = c("depends", "extra"), sep = "_") %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "age",
        "serostatus", 
        "serostatus & serotype",
        "global",
        "not included"
      ) )) %>%
  mutate(color = ifelse(model == best_model, "red", "black"))

mycols = c("#1C9099", "#FEEBE2",
           "#FBB4B9",
           "#CCCCCC", "#FFFFFF")

# plot parameter dependencies 
p1 = df %>%
  ggplot(aes(x = model, y = Parameter)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_x_continuous(limits = c(0.5,8.5),breaks = 1:8) +
  scale_fill_manual(values=mycols) +
  geom_text(aes(label=extra),size=3) + # add extra dependency info as letters 
  scale_y_discrete(labels = c('omega' = expression(omega),
                              "beta" = expression(beta),
                              "eta" = expression(eta))) +
  xlab("Model")

# plot log -lilelihood 
p2 = df %>% 
  ggplot(aes(x = model-0.5, y = mean)) +
  geom_point(aes(color=color, size = color)) + geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha =0.2) + 
  scale_x_continuous(limits = c(0,8), breaks = 1:8) + 
  ylab("Log-likelihood") +
  xlab(" ")  + theme_bw() +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

# plot number of parameters 
p3 =df %>% 
  ggplot(aes(x = model-0.5, y = param)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab("Number of \nparameters") +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,8), breaks = 1:8) +
  scale_y_continuous( breaks = 11:15) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 


# plot ELPD difference  
p4 =df %>% 
  ggplot(aes(x = model-0.5, y = elpd_diff)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab(paste("Difference in \nELPD compared \nto model", best_model)) +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,8), breaks = 1:8) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

g1 = cowplot::plot_grid(p2, NULL, p3, NULL,p4, NULL, p1, ncol=1, 
                        rel_heights = c(1,-.4,1,-.4,1,-.4, 2), 
                        axis = "tblr", align = "hv")


ggsave(g1, file = "BUT/output/figures/model_variants_B.jpg",
       height = 25, width = 30, units="cm", scale = 0.8)


