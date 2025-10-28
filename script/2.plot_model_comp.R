# Script to run model comparison 

# set up 
rm(list = ls())
library(tidyverse)
library(cowplot)
library(posterior)

mycols = c("#B39EB5", "#B2D3A8", "#89A497", "#AEC6CF", "#B8B8AC", "#FFFFFF")

parameters = c("eta", "delta", "n50", "beta", "L", "tau", "rho")

my_theme = theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 8, hjust = 0.5),
      plot.subtitle = element_text(size = 8),
      plot.caption = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_blank(),
      strip.text = element_text(size = 8),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)))

# read waic 
index_files = which(grepl("M", list.files(path = "output/")))
folder_names = list.files(path = "output/")[index_files]

# WAIC files
waic_files = file.path("output", folder_names, "WAIC.RDS")
exists_waic = file.exists(waic_files)
waic_files = waic_files[exists_waic]

waic_list = lapply(waic_files, readRDS)
names(waic_list) = folder_names[exists_waic]

# read loglik 
loglik_files = file.path("output", folder_names, "loglik.RDS")
exists_ll = file.exists(loglik_files)
loglik_files = loglik_files[exists_ll]

loglik_list = lapply(loglik_files, readRDS)
names(loglik_list) = folder_names[exists_ll]

# compare WAIC 
# WAIC = WAIC[!grepl("S", names(WAIC))]
# WAIC = WAIC[!grepl("rho", names(WAIC))]

comp_WAIC = loo::loo_compare(waic_list)
best_model = as.numeric(str_replace(rownames(comp_WAIC)[1], "M", ""))

# Plot lok likelihood 

results = lapply(names(loglik_list), function(name) {
  loglik_array = loglik_list[[name]]
  
  # Convert draws_array to draws_matrix
  loglik_mat = as_draws_matrix(loglik_array)
  
  # Total log-likelihood for each posterior draw
  loglik_sum = rowSums(loglik_mat)
  
  # Summarise
  tibble(
    model = name,
    mean = mean(loglik_sum),
    low = quantile(loglik_sum, 0.025),
    high = quantile(loglik_sum, 0.975)
  )
}) %>%
  bind_rows() %>% 
  arrange(desc(mean))

best_model_ll = as.numeric(str_replace(results[1,]$model, "M", ""))

# Plot number of parameters 

param = data.frame(
  model = 1:nrow(results),
  n = c(96,96,94,96,95,
        96,97,99,99,100,
        102,102,100,103,101,
        103,99,99,99, 102,
        96)
)

plot_data = results %>% 
  mutate(model = str_remove(model, "M"),
         model = as.numeric(model)) %>%
  left_join(param, by = "model")

# Rescale n to the range of log likelihoods
range_ll = range(plot_data$mean, na.rm = TRUE)
range_n  = range(plot_data$n, na.rm = TRUE)

plot_data = plot_data %>%
  mutate(n_rescaled = scales::rescale(n, to = range_ll, from = range_n))

# create empty matrix 
n_models = nrow(results)
m = array(dim = c(3, n_models, ncol = length(parameters)))

# eta
m[3,13:n_models,1] = "age group_2-6yrs"

# delta
m[1, ,2] = "serotype"
m[2, ,2] = "age group"
m[2,1:2,2] = "serotype"
m[3, ,2] = NA

# n50
m[, ,3] = "serostatus & serotype"
m[,1,3] = "serostatus & serotype_MO"

# beta
m[1,4,4] = "age group"
m[1,5:n_models,4] = "age group_4-5yrs"
m[1,c(17),4] = NA

m[2,6:n_models,4] = "age group"
m[2,18:n_models,4] = NA

m[3,7,4] = "age group"
m[3,15,4] = "age group_2-6yrs"

# L
m[ , ,5] = "serostatus"
m[1,c(8, 16),5] = "serostatus & serotype"
m[2,9:n_models,5] = "serostatus & serotype"
m[2,19:n_models,5] = "serostatus"

m[3,10,5] = "serostatus & serotype"

# tau
m[ , ,6]  = "serostatus"
m[1,c(11,20),6] = "serostatus & serotype"
m[2,c(12,19,20),6] = "serostatus & serotype"
m[3, ,6] = NA

# rho
m[ , ,7]  = "global"
m[ ,14,7] = "serotype"

Q = as.data.frame(m[1,,])
D = as.data.frame(m[2,,])
B = as.data.frame(m[3,,])

# models 

colnames(Q) = parameters
Q$model = 1:n_models
colnames(D) = parameters
D$model = 1:n_models
colnames(B) = parameters
B$model = 1:n_models

# format df ready to plot 
df = bind_rows(Q,D,B, .id = "vaccine") %>%  
  pivot_longer(cols = eta:rho,
               names_to = "Parameter",
               values_to = "depends") %>%
  mutate(depends = ifelse(is.na(depends), "not included", depends)) %>%
  separate(depends, into = c("depends", "extra"), sep = "_") %>%
  mutate(
    depends = factor(
      depends,
      levels = c(
        "age group",
        "serotype",
        "serostatus & serotype",
        "serostatus",
        "global",
        "not included"
      ) )) %>%
  mutate(color = ifelse(model == best_model, "red", "black")) %>% 
  mutate(vaccine = factor(vaccine, levels = 1:3, labels = c("Qdenga", "Dengvaxia", "Butantan-DV")))

# plot parameter dependencies 
plot_param = df %>%
  ggplot(aes(x = model, y = Parameter)) +
  geom_tile(aes(fill = depends), color = "black", alpha =0.6) +
  theme_classic() +
  geom_text(aes(label=extra),size=2) + # add extra dependency info as letters 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_x_continuous(limits = c(0.5,n_models+0.5),breaks = 1:n_models) +
  scale_fill_manual(values=mycols) +
  scale_y_discrete(labels = c('rho' = expression(rho),
                              'tau'   = expression(tau),
                              "beta" = expression(beta),
                              "eta" = expression(eta),
                              "delta" = expression(delta))) +
  xlab("Model") + ylab("Parameters") +
  facet_wrap(~vaccine, ncol=1)

# plot ELPD difference  
plot_elpd = comp_WAIC %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "model") %>% 
  mutate(model2 =  as.numeric(str_replace(model, "M", ""))) %>% 
  mutate(color = ifelse(model2 == 21, "red", "black")) %>% 
  ggplot(aes(x = model2, #reorder(model, elpd_diff)
             y = elpd_diff, 
             ymin = elpd_diff - se_diff, 
             ymax = elpd_diff + se_diff,
             color = color)) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "ELPD difference\n(relative to best)") +
  scale_color_manual(values=c("#000000", "#CC0033")) +
  my_theme +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank()
  ) +
  scale_x_continuous(limits = c(0.5,n_models+0.5)) 
  


plot_ll = ggplot(plot_data, aes(x = model)) +
  geom_point(aes(y = mean)) +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  geom_hline(yintercept = results$mean[1], linetype = "dashed") +
  geom_line(aes(y = n_rescaled), colour = "red") +   # add number of parameters
  labs(y = "Log likelihood") +
  my_theme +
  scale_x_continuous(limits = c(0.5,n_models+0.5)) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank()
  ) +
  scale_y_continuous(
    name = "Log likelihood",
    sec.axis = sec_axis(
      ~ scales::rescale(., to = range_n, from = range_ll),
      name = "Number of\nparameters",
      breaks = seq(min(range_n), max(range_n), by = 2),  # every 2
      labels = scales::number_format(accuracy = 1)       # integers only
    )
  )


# output -------------
g1 = cowplot::plot_grid(plot_elpd, NULL, plot_ll,NULL,plot_param, 
                        ncol=1, 
                        rel_heights = c(1,-.45, 1,-.4, 2), 
                        labels = c("a", "", "b", "", "c"),
                        axis = "tblr", align = "hv")


ggsave(g1, file = "output/figures/SupFig5.jpg",
       height = 20, width = 18, units="cm")

