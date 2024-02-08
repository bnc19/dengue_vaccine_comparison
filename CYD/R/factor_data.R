# Functions relating to trial data formatting  ---------------------------------


# 1: function to factorise VCD data when importing it at the start of any script
factor_CYD_cases= function(cases, raw = F){
  
  if(is.null(cases$serotype)) cases$serotype = "both"
  if(is.null(cases$serostatus)) cases$serostatus = "all"
  if(is.null(cases$age)) cases$age = "all"
  if(is.null(cases$arm)) cases$arm = "both"
  if(is.null(cases$time)) cases$time = "all"
  
if(raw == TRUE) { 
  cases = cases  %>%  
    mutate(age = factor(age, 
                        levels = c("2-8yrs", "9-16yrs", "all"),
                        labels = c("2-8yrs", "9-16yrs", "all")), 
           serostatus = factor(serostatus, 
                               levels = c("SN", "SP", "both"), 
                               labels = c("seronegative", "seropositive", "both")),
           serotype = factor(serotype,
                             levels = c("D1", "D2", "D3", "D4", "all"), 
                             labels = c("DENV1", "DENV2", "DENV3", "DENV4", "all")),
           arm = factor(arm, 
                        levels = c("C", "V", "both"), 
                        labels = c("placebo", "vaccine", "both")),
           time = factor(time,
                         levels = c("13-25", "26-36", "37-48", "49-72", "13-72"),
                         labels = c("1-13", "14-24", "25-36", "37-60", "all")))
} else {
  cases = cases  %>%  
    mutate(age = factor(age, 
                        levels = c("2-8yrs", "9-16yrs", "all"),
                        labels = c("2-8yrs", "9-16yrs", "all")), 
           serostatus = factor(serostatus, 
                               levels = c("seronegative", "seropositive", "both"), 
                               labels = c("seronegative", "seropositive", "both")),
           serotype = factor(serotype,
                             levels = c("DENV1", "DENV2", "DENV3", "DENV4", "all"), 
                             labels = c("DENV1", "DENV2", "DENV3", "DENV4", "all")),
           arm = factor(arm, 
                        levels = c("placebo", "vaccine", "both"), 
                        labels = c("placebo", "vaccine", "both")),
           time = factor(time,
                         levels = c("1-13", "14-24", "25-36", "37-60", "all"),
                         labels = c("1-13", "14-24", "25-36", "37-60", "all")))
  
  }
  return(cases)
}



# 2: function to factorise serology data when importing it at the start of any script
factor_CYD_serology= function(data){
  data = data  %>%  
    mutate(serostatus = factor(serostatus,
                               levels = c("SN", "SP"),
                               labels = c("seronegative", "seropositive")),
           serotype = factor(serotype),
           trial = factor(trial, levels = c("Placebo", "TAK"), 
                          labels = c("placebo", "TAK-003")))
  
  return(data)
}

