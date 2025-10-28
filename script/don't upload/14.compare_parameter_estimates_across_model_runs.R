
# set up 
library(tidyverse)
library(cowplot)
library(grid)

# read waic
index_files = which(grepl("M", list.files(path = "output/")))
source_files = paste0("output/", list.files(path = "output/")[index_files], "/posterior.csv")
posterior = (lapply(source_files, read.csv))
n1 = gsub("/posterior.csv", "", gsub("output/M", "", source_files))
names(posterior) = n1
posterior2 = posterior[as.character(sort(as.numeric(n1)))]


param = c("rho[1]", "gamma", "phi", "L[1,1]","L[3,1]", "L[2,1]")

# Model 4 
# extract waic and elpd 
post_df = posterior %>%  
  bind_rows(.id = "model") %>% 
  select(variable,	mean, q5, q95, model) %>% 
  mutate(model = as.numeric(model)) %>%  
  filter(variable  %in% param) %>% 
  mutate_if(is.numeric, round, 2) %>%  
  arrange(variable, model)

print(post_df)
