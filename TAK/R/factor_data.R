# Functions relating to trial data formatting and plotting ---------------------


# 1: function to factorise VCD data when importing it at the start of any script
factor_TAK_VCD = function(VCD){
  VCD = VCD  %>%  
    mutate(age = factor(age, levels = c("4-5yrs", "6-11yrs", "12-16yrs", "all")), 
           serostatus = factor(serostatus,
                               levels = c("SN", "SP", "both"),
                               labels = c("seronegative", "seropositive", "both")),
           serotype = factor(serotype,levels = c("D1", "D2", "D3", "D4", "all"),
                             labels = c("DENV-1", "DENV-2", "DENV-3", "DENV-4", "all")),
           trial = factor(trial, levels = c("P", "V","both"), 
                          labels = c("placebo", "TAK-003", "both")))
  
  return(VCD)
}



# 2: function to factorise serology data when importing it at the start of any script
factor_TAK_serology= function(data){
  data = data  %>%  
    mutate(serostatus = factor(serostatus,
                               levels = c("SN", "SP"),
                               labels = c("seronegative", "seropositive")),
           serotype = factor(serotype),
           trial = factor(trial, levels = c("Placebo", "TAK"), 
                          labels = c("placebo", "TAK-003")))
  
  return(data)
}

