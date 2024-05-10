  
  
  # set up 
  rm(list=ls())
  # set up 
  library(tidyverse)
  library(cowplot)
  library(grid)
  
  
  n_param = c(65,68,62,65,62,65,59,62,68,65,
              70,71,72,73,64,65,66,67,64,64,
              67,66,60,61,64,67,58,61,59,55,
              58, 55, 55, 55)
  
  # read waic
  index_files = which(grepl("M", list.files(path = "TAK/output/")))
  source_files = paste0("TAK/output/", list.files(path = "TAK/output/")[index_files], "/WAIC.RDS")
  WAIC = (lapply(source_files, readRDS))
  n1 = gsub("/WAIC.RDS", "", gsub("TAK/output/M", "", source_files))
  n1[n1 == "32_FINAL"] = "32"
  names(WAIC) = n1
  WAIC2 = WAIC[as.character(sort(as.numeric(n1)))]
  
  # compare WAIC
  comp_WAIC = loo::loo_compare(WAIC2)
  
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
    mutate(param = n_param) %>% 
    left_join(waic_df)
    
  
  models = as.numeric(paste0(1:34))
  parameters= c("L", "w", "lc", "alpha", "beta", "p", "epsilon", "rho", "tau", "omega")
  depend = c("global", "serotype", "serostatus", "age", "outcome")
  
  # create empty matrix 
  m = matrix(nrow = length(models), ncol = length(parameters))
  
  # format 
  d = m %>%
    as.data.frame() 
  
  rownames(d) = models
  colnames(d) = parameters
  
  # fill in matrix with dependencies of each model variant 
  
  # L
  d[ ,1] = "serostatus"
  d[c(28), 1] = paste0("serostatus", " & ", "serotype")
  
  # w
  d[  ,2] = "global"
  d[11:18, 2] = c("serostatus", "serotype")
  d[21:22, 2] = c("serotype", "serostatus")
  
  # lc
  d[ ,3] = "serostatus & serotype"
  d[c(3,4,26),3] = "serostatus & serotype_SN\nMO"
  d[c(5:6),3] = "serostatus & serotype_MO\nMU"
  d[c(7:8,15:25,27:31) ,3] = "serostatus & serotype_MO"
  d[32 ,3] = "serostatus & serotype_MO\nMU"
  d[33 ,3] = "serostatus & serotype_SN\nMO"
  d[34 ,3] = "serostatus & serotype"
  
  # alpha
  d[ ,4] = "global" 
  d[c(13:14,17:22),4] = "serostatus" 
  d[25:26,4] = "serotype"
  
  # beta
  d[ ,5] = "age"
  d[c(9:10,20:22) ,5] = paste0("age", " & ", "outcome",  "_4-5")
  d[23,5] = NA
  d[24:34, 5] = paste0("age", "_4-5")
  
  # p
  d[ ,6] = paste0("serotype" , " & " , "age")
  d[30:34,6] = "age"
  
  # epsilon 
  d[29, 7] = "global"
  
  # rho 
  d[,8] = "serotype"
  d[27:34,8] = "global"
  
  # tau 
  d[ ,9] = "serotype"
  d[c(1,3,5,7,10), 9] = "global"
  
  # omega
  d[32,10] =  "serostatus_MU" 
  d[33,10] =  "serostatus_SN" 
  d[34,10] =  "serostatus_SN\nMU" 
    

# add log lik to dependency matrix 
d$model = models

# format data ready to plot 
df = left_join(d, ll) %>%
  pivot_longer(cols = L:omega,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  separate(depends, into = c("depends", "extra"), sep = "_") %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "age & outcome",
        "age",
        "serotype & age",
        "serotype",
        "serostatus & serotype",
        "serostatus",
        "global",
        "not included"
      ) )) %>%
  mutate(extra = ifelse(is.na(extra), " ", extra)) %>% 
  mutate(extra = factor(extra, levels =  c("4-5", 
                                           "MO",
                                           "MO\nMU",
                                           "SN\nMO",
                                           "MU",
                                           "SN",
                                           "SN\nMU",
                                           " "))) %>% 
  # mutate(extra2 = factor(extra,  labels = c(letters[1:7], " "))) %>%
  mutate(color = ifelse(model == 32, "red", "black"))


mycols = c("#A6BDDB", "#1C9099", "#6BAED6",
           "#666699", "#FBB4B9", 
           "#FEEBE2", "#CCCCCC", "#FFFFFF")

# plot parameter dependencies 
p1 = df %>%
  ggplot(aes(x = model, y = Parameter, group = extra)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_x_continuous(limits = c(0.5,34.5),breaks = 1:34) +
  scale_fill_manual(values = mycols) +
  geom_text(aes(label=extra),size=3) + # add extra dependency info as letters 
  scale_y_discrete(labels = c('omega' = expression(omega),
                              'tau'   = expression(tau),
                              "rho" = expression(rho),
                              "alpha" = expression(alpha),
                              "beta" = expression(beta),
                              "epsilon" = expression(epsilon))) +
  ylab("Parameter") + xlab("Model")
                          
# plot log -lielihood 
p2 = df %>% 
  ggplot(aes(x = model-0.5, y = mean)) +
  geom_point(aes(color=color, size = color)) + geom_line() +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha =0.2) + 
  scale_x_continuous(limits = c(0,34), breaks = 1:34) + ylab("Log-likelihood") +
  xlab(" ") + scale_y_continuous(limits = c(-610,-560)) + theme_bw() +
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
  scale_x_continuous(limits = c(0,34), breaks = 1:34) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 


# plot ELPD difference  
p4 =df %>% 
  ggplot(aes(x = model-0.5, y = elpd_diff)) + # - 0.5 so aligns with other plots 
  geom_point(aes(color=color, size = color)) + geom_line() +
  theme_bw() + ylab("Difference \nin ELPD compared \nto model 32") +
  scale_color_manual(values=c("#000000", "#CC0033"))+
  scale_size_manual(values=c(1,2.5)) +
  scale_x_continuous(limits = c(0,34), breaks = 1:34) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) 

g1 = cowplot::plot_grid(p2, NULL, p3, NULL,p4, NULL, p1, ncol=1, 
                        rel_heights = c(1,-0.45,1,-0.45,1,-0.45, 1.8), 
                        axis = "tblr", align = "hv")


ggsave(g1, file = "TAK/output/figures/full_model_variants_T.jpg",
       height = 20, width = 25, units="cm")


