calculate_CYD_bayes = function(stan_fit,
                           L_mean,
                           L_sd,
                           file_path,
                           lower){
  
  # Bayes factor
  library(truncnorm)
  library(logspline)
  
  # Calculate prior density at 0 
  
  prior_density = dtruncnorm(0,a=lower, mean =L_mean, sd = L_sd)
  
  # Calculate posterior density at 0 
  
  fit_ext = stan_fit$draws(format = "df")
  post_samples = fit_ext$`L[1]`
  
  fit.posterior = logspline(post_samples)
  
  posterior_density = dlogspline(0, fit.posterior)
  BF = prior_density / posterior_density
  
  out = data.frame(
    L_mean = L_mean,
    L_sd = L_sd,
    BF = BF,
    post_dens = posterior_density,
    prior_dens = prior_density
  )
  
  write.csv(out, paste0(file_path, "/BF.csv"))  
}