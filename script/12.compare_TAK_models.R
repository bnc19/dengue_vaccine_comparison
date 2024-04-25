
rm(list=ls())
# set up 
library(tidyverse)
library(cowplot)
library(grid)

# read waic
index_files = which(grepl("M", list.files(path = "TAK/output/")))
source_files = paste0("TAK/output/", list.files(path = "TAK/output/")[index_files], "/WAIC.RDS")
WAIC = (lapply(source_files, readRDS))
n1 = gsub("/WAIC.RDS", "", gsub("TAK/output/M", "", source_files))
names(WAIC) = n1
WAIC2 = WAIC[as.character(sort(as.numeric(n1)))]

# compare WAIC
comp_WAIC = loo::loo_compare(WAIC2)
best_model = as.numeric(rownames(comp_WAIC)[1])

# Model 2 
# extract waic and elpd 
waic_df = comp_WAIC %>%  
  as.data.frame() %>% 
  rownames_to_column(var = "model") %>% 
  select(model, elpd_diff, waic) %>% 
  mutate(model = as.numeric(model))

# get log lik
ll_source = paste0("TAK/output/", list.files(path = "TAK/output/")[index_files], "/posterior.csv")
post = bind_rows(lapply(ll_source, read.csv))

ll = post %>%  
  filter(variable == "ll") %>% 
  mutate(model = as.numeric(n1)) %>% 
select(mean, q5, q95, model) %>% 
  arrange(model) %>% 
  left_join(waic_df)
  
models = as.numeric(paste0(1:4))
parameters= c("lc",  "omega", "kappa")

# create empty matrix 
m = matrix(nrow = length(models), ncol = length(parameters))

# format 
d = m %>%
  as.data.frame() 

rownames(d) = models
colnames(d) = parameters

# fill in matrix with dependencies of each model variant 

# lc
d[1,1] = "serostatus & serotype (serostatus = MO)"
d[2,1] = "serostatus & serotype (serostatus = MO, MU)"
d[3,1] = "serostatus & serotype (serostatus = MO, SN)"
d[4,1] = "serostatus & serotype"

# omega 
d[c(2,4), 2] = "global"

# kappa 
d[3:4, 3] = "global"

# add log lik to dependency matrix 
d$model = models

# format data ready to plot 
df = left_join(d, ll) %>%
  pivot_longer(cols = lc:kappa,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "serostatus & serotype",
        "serostatus & serotype (serostatus = MO)",
        "serostatus & serotype (serostatus = MO, MU)", 
        "serostatus & serotype (serostatus = MO, SN)",
        "global",
        "not included"
      ) )) %>%
  mutate(color = ifelse(model == best_model, "red", "black"))

mycols = c("#A6BDDB", "#1C9099",
           "#666699", "#FBB4B9", 
           "#CCCCCC", "#FFFFFF")

# plot parameter dependencies 
p1 = df %>%
  ggplot(aes(x = model, y = Parameter)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_x_continuous(limits = c(0.5,4.5),breaks = 1:4) +
  scale_fill_manual(values=mycols) +
  scale_color_manual(values=rep("black",8)) + # make all letters black 
  guides(col = guide_legend(order = 1)) + # add extra legend 
  scale_y_discrete(labels = c('omega' = expression(omega),
                              "kappa" = expression(kappa))) +
                     ylab("Model")

# plot log -lielihood 
p2 = df %>% 
  ggplot(aes(x = model-0.5, y = mean)) +
  geom_point(aes(color=color, size = color)) + geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha =0.2) + 
  scale_x_continuous(limits = c(0,4), breaks = 1:4) + ylab("Log-likelihood") +
  xlab(" ") + scale_y_continuous(limits = c(-610,-560)) + theme_bw() +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 


# plot ELPD difference  
p3 =df %>% 
  ggplot(aes(x = model-0.5, y = elpd_diff)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab(paste("Difference in ELPD \ncompared to model", best_model)) +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,4), breaks = 1:4) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

g1 = cowplot::plot_grid(p2, NULL, p3, NULL, p1, ncol=1, 
                        rel_heights = c(1,0,1,0, 1.8), 
                        axis = "tblr", align = "hv")


ggsave(g1, file = "TAK/output/figures/model_variants_T.jpg",
       height = 20, width = 30, units="cm", scale = 0.7)


