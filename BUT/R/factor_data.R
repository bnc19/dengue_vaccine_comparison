# Functions relating to trial data formatting ----------------------------------

# 1: function to factorise VCD data when importing it at the start of any script
factor_BUT_cases = function(cases){
  
  if(is.null(cases$Serostatus)) cases$Serostatus = "both"
  if(is.null(cases$Age)) cases$Age = "all"
  if(is.null(cases$Serotype)) cases$Serotype = "all"
  if(is.null(cases$Arm)) cases$Arm = "both"
  
  cases = cases  %>%  
    mutate(Age = factor(Age,
                        levels = c("2-6yrs", "7-17yrs", "18-59yrs", "all")), 
           Serostatus = factor(Serostatus,
                               levels = c("SN", "SP", "both"),
                               labels = c("seronegative", "seropositive", "both")),
           Serotype = factor(Serotype,
                             levels = c("D1", "D2", "all"),
                             labels = c("DENV1", "DENV2", "all")),
           Arm = factor(Arm, 
                        levels = c("placebo", "vaccine", "both")))
  
  return(cases)
}
