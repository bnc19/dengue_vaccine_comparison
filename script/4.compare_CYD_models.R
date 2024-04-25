# Script to run model comparison 

# set up 
rm(list = ls())
library(tidyverse)
library(cowplot)

n_param = c(37,40,34,34,35,
            32,34,34,34,31,
            34,34,31)

# read waic
index_files = which(grepl("M", list.files(path = "CYD/output/")))
source_files = paste0("CYD/output/", list.files(path = "CYD/output/")[index_files], "/WAIC.RDS")
WAIC = (lapply(source_files, readRDS))
n1 = gsub("/WAIC.RDS", "", gsub("CYD/output/M", "", source_files))
names(WAIC) = n1
WAIC2 = WAIC[as.character(sort(as.numeric(n1)))]

# compare WAIC
comp_WAIC = loo::loo_compare(WAIC2) # M12 

best_model = as.numeric(rownames(comp_WAIC)[1])

# extract waic and elpd 
waic_df = comp_WAIC %>%  
  as.data.frame() %>% 
  rownames_to_column(var = "model") %>% 
  select(model, elpd_diff, waic) %>% 
  mutate(model = as.numeric(model))

# get log lik
ll_source = paste0("CYD/output/", list.files(path = "CYD/output/")[index_files], "/posterior.csv")
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
parameters= c("L", "delta", "lc", "beta", "p", "tau", "omega", "kappa")

# create empty matrix 
m = matrix(nrow = length(models), ncol = length(parameters))

# format 
d = m %>%
  as.data.frame() 

rownames(d) = models
colnames(d) = parameters

# L
d[ ,1] = "serostatus & serotype"
d[11:12,1] = "serotype"

# delta
d[ ,2] = "serotype"
d[4,2] = "global"
d[5:13,2] = "age"

# lc 
d[ ,3] = "serostatus & serotype (serostatus = MO)"
d[c(8,12,13),3] = "serostatus & serotype (serostatus = MO, MU)"
d[9:11,3] = "serostatus & serotype"

# beta
d[ ,4] = "age"
d[7:13,4] = NA

# p 
d[ ,5] = "age"
d[2,5] = "serotype & age"

# tau 
d[ ,6] = "serotype"
d[c(3,6,10:13) ,6] = "global"

# omega
d[8:13,7] = "global"

# kappa
d[9:11, 8] = "global"

# add log lik to dependency matrix 
d$model = models

# format data ready to plot 
df = left_join(d, ll) %>%
  pivot_longer(cols = L:kappa,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "age",
        "serotype & age",
        "serotype",
        "serostatus & serotype",
        "serostatus & serotype (serostatus = MO)",
        "serostatus & serotype (serostatus = MO, MU)",
        "global",
        "not included"
      ) )) %>%
  mutate(color = ifelse(model == best_model, "red", "black"))

mycols = c("#A6BDDB", "#1C9099", "#6BAED6",
           "#666699", "#FBB4B9", 
           "#FEEBE2", "#CCCCCC", "#FFFFFF")

# plot parameter dependencies 
p1 = df %>%
  ggplot(aes(x = model, y = Parameter)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_x_continuous(limits = c(0.5,13.5),breaks = 1:13) +
  scale_fill_manual(values=mycols) +
  scale_y_discrete(labels = c('omega' = expression(omega),
                              'tau'   = expression(tau),
                              "kappa" = expression(kappa),
                              "beta" = expression(beta),
                              "delta" = expression(delta))) +
  ylab("Model")


# plot log -likelihood 
p2 = df %>% 
  ggplot(aes(x = model-0.5, y = mean)) +
  geom_point(aes(color=color, size = color)) + geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha =0.2) + 
  scale_x_continuous(limits = c(0,13), breaks = 1:13) + ylab("Log-likelihood") +
  xlab(" ")  + theme_bw() +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

# plot number of parameters 
p3 = df %>% 
  ggplot(aes(x = model-0.5, y = param)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab("Number of parameters") +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,13), breaks = 1:13) +
  scale_y_continuous(limits = c(30,40), breaks = seq(30,40,by=2)) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 


# plot ELPD difference  
p4 =df %>% 
  ggplot(aes(x = model-0.5, y = elpd_diff)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab(paste("Difference in ELPD \ncompared to model", best_model)) +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,13), breaks = 1:13) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

g1 = cowplot::plot_grid(p2, NULL, p3, NULL,p4, NULL, p1, ncol=1, 
                        rel_heights = c(1,0,1,0,1,0, 1.8), 
                        axis = "tblr", align = "hv")


ggsave(g1, file = "CYD/output/figures/model_variants_C.jpg",
       height = 30, width = 30, units="cm", scale = 0.7)


################################################################################
##################          Compare severe model              ##################       
################################################################################

index_files_sev = which(grepl("M", list.files(path = "CYD/output/severe/")))
source_files_sev = paste0("CYD/output/severe/", 
                          list.files(path = "CYD/output/severe/")[index_files_sev], "/WAIC.RDS")

WAIC_sev = (lapply(source_files_sev, readRDS))
n2 = gsub("/WAIC.RDS", "", gsub("CYD/output/severe/M", "", source_files_sev))
names(WAIC_sev) = n2
WAIC_sev2 = WAIC_sev[as.character(sort(as.numeric(n2)))]

n_param_sev = c(35, 32, 33, 33, 34, 
                35, 36, 34, 35, 39, 
                39, 38, 37, 38, 38,
                38)
# compare WAIC
comp_WAIC_sev = loo::loo_compare(WAIC_sev2) # models start from M0 
best_model_sev = as.numeric(rownames(comp_WAIC_sev)[1])

# M9 is the best fitting model

waic_df_sev = comp_WAIC_sev %>%  
  as.data.frame() %>% 
  rownames_to_column(var = "model") %>% 
  select(model, elpd_diff, waic) %>% 
  mutate(model = as.numeric(model) + 1) # add 1 so on the same scale as other models 

best_model_sev = best_model_sev + 1 # add 1 so on the same scale as other models 

ll_source_sev = paste0("CYD/output/severe/", 
                       list.files(path = "CYD/output/severe/")[index_files_sev], "/posterior.csv")

post_sev = bind_rows(lapply(ll_source_sev, read.csv))

ll_sev = post_sev %>%  
  filter(variable == "ll") %>% 
  mutate(model = as.numeric(n2) +1) %>% 
  select(mean, q5, q95, model) %>% 
  arrange(model) %>% 
  mutate(param = n_param_sev) %>% 
  left_join(waic_df_sev)

# models 
models_sev = as.numeric(paste0(1:length(n_param_sev)))
parameters_sev = c("L", "tau",  "beta", "epsilon", "psi", "delta")

# create empty matrix 
m_sev = matrix(nrow = length(models_sev), ncol = length(parameters_sev))

# format 
d_sev = m_sev %>%
  as.data.frame() 

rownames(d_sev) = models_sev
colnames(d_sev) = parameters_sev

# TO DO: CHECK VARIANT FIGURE AGAINST MODEL RUNS

# L
d_sev[ ,1] = "global"
d_sev[c(1,10,12:16), 1] = "serotype"

# tau
d_sev[11,2] = "serotype"

# beta
d_sev[c(3,6,8,14) ,3] = "age"
d_sev[c(7,9,10,11,12,15,16) ,3] = "age & outcome"

# epsilon
d_sev[c(5,6,7,10,11,13:16) ,4] = "global"

# psi 
d_sev[ ,5] = "global"
d_sev[c(4:14,16) ,5] = "age"


# delta
d_sev[ ,6] = "age"
d_sev[16 ,6] = "global"

# add log lik to dependency matrix 
d_sev$model = models_sev

# format data ready to plot 
df_sev = left_join(d_sev, ll_sev) %>%
  pivot_longer(cols = L:delta,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "age",
        "serotype",
        "age & outcome",
        "global",
        "not included"
      ) )) %>%
  mutate(color = ifelse(model == best_model_sev, "red", "black"))

mycols_sev = c("#A6BDDB", "#1C9099", "#FBB4B9", 
           "#CCCCCC", "#FFFFFF")

# plot parameter dependencies 
p1_sev = df_sev %>%
  ggplot(aes(x = model, y = Parameter)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_x_continuous(limits = c(0.5,16.5),breaks = 1:16) +
  scale_fill_manual(values=mycols_sev) +
  scale_y_discrete(labels = c('psi' = expression(psi),
                              'tau'   = expression(tau),
                              "beta" = expression(beta),
                              "epsilon" = expression(epsilon),
                              "delta" = expression(delta))) +
  ylab("Model")


# plot log -likelihood 
p2_sev = df_sev %>% 
  ggplot(aes(x = model-0.5, y = mean)) +
  geom_point(aes(color=color, size = color)) + geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha =0.2) + 
  scale_x_continuous(limits = c(0,16), breaks = 1:16) + ylab("Log-likelihood") +
  xlab(" ")  + theme_bw() +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

# plot number of parameters 
p3_sev = df_sev %>% 
  ggplot(aes(x = model-0.5, y = param)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab("Number of parameters") +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,16), breaks = 1:16) +
  scale_y_continuous(limits = c(30,40), breaks = seq(30,40,by=2)) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 


# plot ELPD difference  
p4_sev = df_sev %>% 
  ggplot(aes(x = model-0.5, y = elpd_diff)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab(paste("Difference in ELPD \ncompared to model", best_model_sev)) +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,16), breaks = 1:16) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

g1_sev = cowplot::plot_grid(p2_sev, NULL, p3_sev, NULL,p4_sev, NULL, p1_sev, ncol=1, 
                        rel_heights = c(1,-0.3,1,-0.3,1,-0.3, 1.8), 
                        axis = "tblr", align = "hv")


ggsave(g1_sev, file = "CYD/output/figures/severe_model_variants_C.jpg",
       height = 25, width = 30, units="cm", scale = 0.9)

