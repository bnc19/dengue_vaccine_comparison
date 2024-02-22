
# function to extract VE at 2 years  -------------------------------------------
format_2_year_VE = function(fit, vaccine){

  out1 =  fit %>%
    as.data.frame() %>%
    mutate(ni = row_number()) %>%
    pivot_longer(cols = -ni) %>%   # Pivot iterations to rows
    group_by(name) %>%
    mutate(value = value * 100) %>%
    separate(name, into = c("group", "name"), sep = "\\[") %>%  # rename
    separate(name, into = c("name", NA), sep = "\\]")
  
  if(vaccine != "BUT"){
    out2 =  out1 %>%  
    separate(name,
             into = c("Serostatus", "Serotype", "Outcome", "Month")) %>% # each variable a column
    filter(Outcome == 1) 
  } else{
    out2 =  out1 %>%  
      separate(name,
               into = c("Serostatus", "Serotype", "Age", "Month"))  %>%  # each variable a column
    group_by(Serostatus, Serotype, ni, Month) %>% # average over age 
      summarise(value = mean(value, na.rm = T))  
    
  }
  
  out3 = out2 %>%  
    mutate(Month = as.numeric(Month)) %>% 
    filter(Month <= 24) %>% # Butantan data only goes to month 24 
    group_by(Serostatus, Serotype) %>%
    summarise( # mean and uncertainty of the posterior and over time
      lower = quantile(value, 0.025, na.rm = T),
      mean = mean(value, na.rm = T),
      upper = quantile(value, 0.975, na.rm = T)
    ) %>%
    mutate( # factor for plotting 
      Serotype = factor(paste0("DENV", Serotype)),
      Serostatus = factor(
        Serostatus, levels = 1:3, 
        labels = c("seronegative", "monotypic", "multitypic")
      )
    )
  
  return(out3)
}



# function to extract VE at 4.5 years for Q and D  -----------------------------
format_4.5_year_VE = function(fit_ext){
  
  out1 =  fit_ext %>%
    as.data.frame() %>%
    mutate(ni = row_number()) %>%
    pivot_longer(cols = -ni) %>%   # Pivot iterations to rows
    group_by(name) %>%
    mutate(value = value * 100) %>%
    separate(name, into = c("group", "name"), sep = "\\[") %>%  # rename
    separate(name, into = c("name", NA), sep = "\\]") %>%  
      separate(name,
               into = c("Serostatus", "Serotype", "Outcome", "Month")) %>% # each variable a column
    mutate(Month = as.numeric(Month)) %>% 
    filter(Month <= 54) %>% # Takeda data only goes to month 54  
    group_by(Serostatus, Serotype,Outcome) %>% 
    summarise( # uncertainity over time and across iterations 
      lower = quantile(value, 0.025, na.rm = T),
      mean = mean(value, na.rm = T),
      upper = quantile(value, 0.975, na.rm = T)
    ) %>%
    mutate( # factor for plotting 
      Serotype = factor(paste0("DENV", Serotype)),
      Outcome = factor(Outcome, labels = c("Symptomatic", "Hospitalised")),
      Serostatus = factor(
        Serostatus, levels = 1:3, 
        labels = c("seronegative", "monotypic", "multitypic")
      )
    )
  
  return(out1)
}

# function to extract model results --------------------------------------------

extract_model_results = function(fit_ext){
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
