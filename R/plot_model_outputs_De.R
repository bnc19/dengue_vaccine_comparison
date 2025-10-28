# function to add extra populations to VCD data and calculate AR --------------- 

calc_attack_rates_De = function(VCD, titre = F) {
  out = VCD %>% 
    mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
           lower = binconf(Y,N, method = "exact")[,2] * 100, 
           upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
    select(- c(Y,N)) %>% 
    mutate(type = "data") 
  
  return(out)
}

# function to extract model results --------------------------------------------
extract_model_results_De = function(fit_ext){
  
  out =  fit_ext %>%  
    as.data.frame() %>%  
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value * 100) %>%  
    summarise(
      lower = quantile(value, 0.025, na.rm = T),
      mean = mean(value, na.rm = T),
      upper = quantile(value, 0.975, na.rm = T)
    )
  
  out2 = out %>%  
    separate(name, into = c("group", "name"), sep = "\\[") %>% 
    separate(name, into = c("name", NA), sep = "\\]") %>% 
    mutate(type = "model") 
  
  return(out2)
}

