# Functions relating to trial data formatting ----------------------------------

# 1: function to factorise VCD data when importing it at the start of any script
factor_BUT_cases = function(cases){
  
  if(is.null(cases$serostatus)) cases$serostatus = "both"
  if(is.null(cases$age)) cases$age = "all"
  if(is.null(cases$serotype)) cases$serotype = "all"
  if(is.null(cases$arm)) cases$arm = "both"
  
  cases = cases  %>%  
    mutate(age = factor(age,
                        levels = c("2-6yrs", "7-17yrs", "18-59yrs", "all")), 
           serostatus = factor(serostatus,
                               levels = c("SN", "SP"),
                               labels = c("seronegative", "seropositive")),
           serotype = factor(serotype,
                             levels = c("D1", "D2", "all"),
                             labels = c("DENV1", "DENV2", "all")),
           arm = factor(arm, 
                        levels = c("C", "V"),
                        labels = c("placebo", "vaccine")),
           time = factor("1-24"))
  
  return(cases)
}
