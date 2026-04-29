
# Calculate Savage-Dickey Bayes factors for the waning parameter L
#
# Tests the hypothesis that the antibody waning rate (L) equals zero — i.e.
# that there is no waning — separately for each of the three vaccines.
# Uses the Savage-Dickey density ratio: BF = p(L=0 | prior) / p(L=0 | posterior).
# A BF > 1 indicates the data favour no waning over waning; BF < 1 indicates
# evidence for waning. Results are written to BF.csv in the specified directory.

calculate_bayes = function(stan_fit,
                           L_mean,
                           L_sd,
                           file_path,
                           lower){
  
  # Bayes factor
  library(truncnorm)
  library(logspline)
  
  # ── Prior density at L = 0 ────────────────────────────────────────────────
  
  prior_density = dtruncnorm(0, a=lower, mean = L_mean, sd = L_sd)
  
  # ── Extract posterior draws for L ────────────────────────────────────────
  
  fit_ext = stan_fit$draws(format = "df")
  
  post_samples_Q = fit_ext$`L[1,1]`
  post_samples_D = fit_ext$`L[2,1]`
  post_samples_B = fit_ext$`L[3,1]`
  
  # ── Fit logspline density to posterior draws ──────────────────────────────
  
  fit.posterior_Q = logspline(post_samples_Q)
  fit.posterior_D = logspline(post_samples_D)
  fit.posterior_B = logspline(post_samples_B)
  
  # Evaluate the fitted posterior density at L = 0 (denominator of Savage-Dickey)
  
  posterior_density_Q = dlogspline(0, fit.posterior_Q)
  posterior_density_D = dlogspline(0, fit.posterior_D)
  posterior_density_B = dlogspline(0, fit.posterior_B)
  
  # ── Compute Bayes factors ─────────────────────────────────────────────────
  
  
  BF_Q = prior_density / posterior_density_Q
  BF_D = prior_density / posterior_density_D
  BF_B = prior_density / posterior_density_B
  
  out = data.frame(
    L_mean = L_mean,
    L_sd = L_sd,
    BF_Q = BF_Q,
    BF_D = BF_D,
    BF_B = BF_B,
    post_dens_Q = posterior_density_Q,
    post_dens_D = posterior_density_D,
    post_dens_B = posterior_density_B,
    prior_dens = prior_density
  )
  
  write.csv(out, paste0(file_path, "/BF.csv"))  
}