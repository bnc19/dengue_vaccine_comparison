
# Function to calculate the DIC ------------------------------------------------

calculate_DIC = function(
  stan_fit, 
  stan_fit_fix 
){
  library(loo)
  
  log_lik = extract_log_lik(stan_fit)
  log_lik_mean = extract_log_lik(stan_fit_fix)
  
  # calculate total model log likelihood 
  
  sum_log_lik = rowSums(log_lik) # total LL by iteration 
  
  avg_log_lik = mean(sum_log_lik) # mean LL 
  CrI_log_lik = quantile(sum_log_lik, probs = c(0.025,0.975))
  
  sum_log_lik_mean = sum(log_lik_mean) # total LL at mean of posterior 
  
  # calculate deviance 
  
  deviance = - 2 * sum_log_lik 
  deviance_at_mean = - 2 * sum_log_lik_mean
  mean_of_deviance = mean(deviance)
  
  # calculate pD (effective number of parameters)
  
  pD = mean_of_deviance - deviance_at_mean
  
  # calculate dic
  dic = deviance_at_mean + 2 * pD
  
  out = data.frame(
    DIC = dic,
    mean_LL = avg_log_lik,
    lower_LL = CrI_log_lik[1],
    upper_LL = CrI_log_lik[2]
  )
  rownames(out) = ""
  return(out)
}

# Function to run diagnostics on stan model ------------------------------------

diagnose_stan_fit = function(
  stan_fit,
  file_path,
  pars = c("lambda[1,1]", "rho[1]", "L[1]", "w[1]", "p[1]", "sens", "lc50[1]", "delta[1]", "gamma")
){
  
  # required package 
  library(bayesplot)
  library(loo)

  # # log lik 
  # 
  # log_lik = extract_log_lik(stan_fit)
  # sum_log_lik = rowSums(log_lik) # total LL by iteration 
  # out = data.frame(
  #   mean_log_lik = mean(sum_log_lik), # mean LL 
  #   lower_log_lik = quantile(sum_log_lik, probs = 0.025),
  #   upper_log_lik = quantile(sum_log_lik, probs = 0.975)
  # )
  # 
  # write.csv(out,  paste0(file_path,"/log_lik.csv"))
  # 
  # 
  # write.csv(log_lik,  paste0(file_path, "/log_lik_full.csv"))
  # 
  # 
  # # loo   
  # loo = stan_fit$loo()
  # # saveRDS(loo, paste0(file_path,"/loo.RDS"))
  # 
  # 
  # 
  # 
  # # pars 
  pars1 = pars[1:round(length(pars)/2)]
  pars2 = pars[(round(length(pars)/2) + 1):length(pars)]


  # rhat 
  
  if(any(rhat(stan_fit) > 1.01, na.rm = T)){
    print("Rhat above 1.01")
  } else{
    print("Rhat good")}
  

  # divergence #
  # rstan::check_divergences(stan_fit)
  diag_sum <- stan_fit$diagnostic_summary()
  message(paste0("# divergences = ",diag_sum$num_divergent))
  
  # get posterior 
  
  stan_fit_post= stan_fit$draws()
  
  # plot hist by chain
  hist_chain1 =  mcmc_hist_by_chain(stan_fit_post, pars = pars1)
  hist_chain2 =  mcmc_hist_by_chain(stan_fit_post, pars = pars2)
 
  ggsave(
   hist_chain1,
   file =  paste0(file_path, "/hist_chain1.png"),
   height = 20,
   width = 50,
   unit = "cm",
   dpi = 720
 )
  
 ggsave(
   hist_chain2,
   file =  paste0(file_path, "/hist_chain2.png"),
   height = 20,
   width = 50,
   unit = "cm",
   dpi = 720
 )
 
  # markov chain trace plots   
  
  markov_trace = mcmc_trace(stan_fit_post, pars = pars)
  
  ggsave(
    markov_trace,
    file =  paste0(file_path, "/trace_plot.png"),
    height = 20,
    width = 50,
    unit = "cm",
    dpi = 720
  )
  
  # pairs plot 
  np_cp = nuts_params(stan_fit)
  pairs_plot1 = mcmc_pairs(stan_fit_post,
                           np = np_cp,
                           pars = pars1)
  
  pairs_plot2 = mcmc_pairs(stan_fit_post, 
                           np = np_cp,
                           pars = pars2)
  
  ggsave(
    pairs_plot1,
    file =  paste0(file_path, "/pairs_plot1.png"),
    height = 50,
    width = 50,
    unit = "cm",
    dpi = 300
  )
  
  ggsave(
    pairs_plot2,
    file = paste0(file_path, "/pairs_plot2.png"),
    height = 50,
    width = 50,
    unit = "cm",
    dpi = 720
  )
  
  # density 
  dens_plot1 = mcmc_dens(stan_fit_post, pars = pars[1:round(length(pars)/2)])
  dens_plot2 = mcmc_dens(stan_fit_post, pars = pars[(round(length(pars)/2) + 1):length(pars)])
  
  ggsave(
    dens_plot1,
    file =  paste0(file_path, "/dens_plot1.png"),
    height = 50,
    width = 50,
    unit = "cm",
    dpi = 300
  )
  
  ggsave(
    dens_plot2,
    file = paste0(file_path, "/dens_plot2.png"),
    height = 50,
    width = 50,
    unit = "cm",
    dpi = 720
  )
  
 
  
  write.csv(loo$estimates, file = paste0(file_path, "/loo.csv"))
  write.csv(data.frame(pareto_k = sum(loo$diagnostics$pareto_k > 0.7)), 
            paste0(file_path, "/pk.csv"))
}