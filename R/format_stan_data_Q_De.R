format_stan_data_Q_De = function(baseline_SP_Q,
                            VCD_Q,
                            hosp_Q,
                            mu_Q,
                            VCD_years_Q,
                            time_Q = time_Q,   # Time points for Q model
                            baseline_SP_De, 
                            cases_De,
                            mu_De,
                            VCD_years_De, 
                            time_De = time_De,  # Time points for D model
                            C = C,   # No. current immune status
                            B = B,   # No. baseline serostatus
                            K = K,   # No. serotypes
                            V = V,   # No. trial arms
                            J_Q = J_Q, # No. age groups (Qdenga)
                            J_De = J_De, # No. age groups (Dengvaxia)
                            R = R,   # No. outcomes
                            TR = TR,   # No. vaccine trials (Q and D)
                            HI = HI,       # Period of heterotypic immunity
                            uniform = uniform, # Use uniform priors on probabilities of reporting? (0/1)
                            include_pK3 = include_pK3, # Include serotype-specific probability? (0/1)
                            include_eps = include_eps, # Include enhanced secondary hospitalization? (0/1)
                            L_mean = L_mean, # Mean of L prior (real)
                            L_sd = L_sd, # SD of L prior (real)
                            enhancement = enhancement, # Vaccine enhancement of SN? (0/1)
                            mono_lc_SN = mono_lc_SN, # Serotype-specific lc SN? (0/1/2)
                            mono_lc_MU = mono_lc_MU, # Serotype-specific lc MO? (0/1/2)
                            rho_K = rho_K, # Mono/serotype-specific rho (0/1)
                            w_CK = w_CK, # Mono/serostatus/serotype w (0/1/2)
                            alpha_CK = alpha_CK, # Mono/serostatus/serotype alpha (0/1/2)
                            tau_K = tau_K, # Mono/serotype-specific tau (0/1)
                            include_beta = include_beta, # Include age-specific nc50? (array[TR], 0-3)
                            L_K = L_K, # Mono/serotype-specific L (array[TR], 0/1)
                            delta_KJ = delta_KJ, # Single/serotype-specific/age delta (array[TR], 0-2),
                            psi_J = psi_J # Age specific prob symp is severe (not currently used)
                            ) {

list_Q_data = format_stan_data_Q(
    baseline_SP = baseline_SP_Q,
    VCD = VCD_Q,
    hosp = hosp_Q,
    B = B,
    R = R,
    K = K,
    V = V,
    J = J_Q,
    C = C,
    HI = HI,
    include_pK3 = include_pK3,
    include_eps = include_eps,
    include_beta = include_beta[1],
    mono_lc_SN = mono_lc_SN,
    mono_lc_MU = mono_lc_MU,
    rho_K = rho_K,
    L_K = L_K[1],
    w_CK = w_CK,
    delta_KJ = delta_KJ[1],
    alpha_CK = alpha_CK,
    tau_K = tau_K,
    L_sd = L_sd,
    L_mean = L_mean,
    enhancement = enhancement,
    VCD_years = VCD_years_Q,
    time = time_Q,
    mu = mu_Q)

list_De_data = format_stan_data_De(
  baseline_SP = baseline_SP_De,
  cases = cases_De,
  VCD_years = VCD_years_De,
  time = time_De,
  mu = mu_De,
  B = B,
  R = R,
  K = K,
  V = V,
  J = J_De,
  C = C,
  HI = HI,
  include_pK2 = include_pK3,
  include_eps = include_eps,
  include_beta = include_beta[2],
  mono_lc_SN = mono_lc_SN,
  mono_lc_MU = mono_lc_MU,
  rho_K = rho_K,
  L_K = L_K[2],
  w_CK = w_CK,
  alpha_CK = alpha_CK,
  tau_K = tau_K,
  delta_KJ = delta_KJ[2],
  L_mean = L_mean,
  L_sd = L_sd,
  enhancement = enhancement,
  psi_J = psi_J
)

# combine MU 
mu = array(NA, dim = c(TR, B, K)) # Create an empty array of correct size
mu[1,,] = list_Q_data$mu  # Assign Qdenga data to TR=1
mu[2,,] = list_De_data$mu   # Assign Dengvaxia data to TR=2 
# data -------------------------------------------------------------------------

data_list = list(
  C = C,   # No. current immune status
  B = B,   # No. baseline serostatus
  K = K,   # No. serotypes
  D_Q = list_Q_data$D, # No. Qdenga VCD data
  D_De = list_De_data$D, # No. Dengvaxia VCD data
  V = V,   # No. trial arms
  J_Q = list_Q_data$J, # No. age groups (Qdenga)
  J_De = list_De_data$J, # No. age groups (Dengvaxia)
  A = list_Q_data$A,   # No. age-specific VCD data
  T_Q = list_Q_data$T, # No. Q model time points
  T_De = list_De_data$T, # No. D model time points
  R = R,   # No. outcomes
  TR = TR,   # No. vaccine trials (Q and D)
  time_Q = list_Q_data$time,   # Time points for Q model
  time_De = list_De_data$time,  # Time points for D model
  HI = HI,       # Period of heterotypic immunity
  mu = mu,       # Titres at T0 (array[TR,B,K])
  
  # Dengvaxia data
  SP_J_De = list_De_data$SP_J,  # Baseline seropositive (array[J_De])
  pop_J_De = list_De_data$pop_J, # Baseline population (array[J_De])
  VCD_De = list_De_data$VCD,   # Total VCD (int)
  HOSP_D_De = list_De_data$HOSP_D, # Total hospitalizations (array[D_De])
  VCD_VK_De = list_De_data$VCD_VK, # VCD (array[V*K])
  VCD_BVJ_De = list_De_data$VCD_BVJ, # VCD (array[B*V*J_De])
  HOSP_BVJD_De = list_De_data$HOSP_BVJD, # Hospitalizations (array[B*V*J_De,D_De])
  HOSP_BVKJ_De = list_De_data$HOSP_BVKJ, # Hospitalizations (array[B*V*K*J_De])
  pop_VCD_BVJ_De = list_De_data$pop_VCD_BVJ, # Population (array[B,V,J_De])
  pop_HOSP_BVJD_De = list_De_data$pop_HOSP_BVJD, # Population (array[B,V,J_De,D_De])
  pop_VCD_De = list_De_data$pop_VCD, # Population for VCD (int)
  pop_HOSP_De = list_De_data$pop_HOSP, # Mean hospitalization population across trial (int)
  pop_HOSP_D_De = list_De_data$pop_HOSP_D, # Hospitalization population over time (array[D_De])
  
  # Qdenga data
  SP_J_Q = list_Q_data$SP_J,  # Baseline seropositive (array[J_Q])
  pop_J_Q = list_Q_data$pop_J, # Baseline population (array[J_Q])
  VCD_D_Q = list_Q_data$VCD_D, # Total VCD (array[D_Q])
  HOSP_D_Q = list_Q_data$HOSP_D, # Total hospitalizations (array[D_Q])
  VCD_BVKD_Q = list_Q_data$VCD_BVKD, # VCD (array[B*V*K,D_Q])
  VCD_BVJA_Q = list_Q_data$VCD_BVJA, # VCD (array[B*V*J_Q,A])
  VCD_KJ2_Q = list_Q_data$VCD_KJ2, # VCD at 12 and 24 months (array[K*J_Q,2])
  HOSP_BVJA_Q = list_Q_data$HOSP_BVJA, # Hospitalizations (array[B*V*J_Q,A])
  HOSP_BVK4_Q = list_Q_data$HOSP_BVK4, # Hospitalizations (array[B*V*K,4])
  HOSP_KJ2_Q = list_Q_data$HOSP_KJ2, # Hospitalizations at 12 and 24 months (array[K*J_Q,2])
  pop_Q = list_Q_data$pop, # Population (array[B,V,J_Q,D_Q])
  pop_BVD_Q = list_Q_data$pop_BVD, # Population (array[B,V,D_Q])
  N_VCD_BV5_Q = list_Q_data$N_VCD_BV5, # Cases in t_5 (not age-specific) (array[B,V])
  
  # Model flags
  uniform = uniform, # Use uniform priors on probabilities of reporting? (0/1)
  include_pK3 = include_pK3, # Include serotype-specific probability? (0/1)
  include_eps = include_eps, # Include enhanced secondary hospitalization? (0/1)
  L_mean = L_mean, # Mean of L prior (real)
  L_sd = L_sd, # SD of L prior (real)
  enhancement = enhancement, # Vaccine enhancement of SN? (0/1)
  mono_lc_SN = mono_lc_SN, # Serotype-specific lc SN? (0/1/2)
  mono_lc_MU = mono_lc_MU, # Serotype-specific lc MO? (0/1/2)
  rho_K = rho_K, # Mono/serotype-specific rho (0/1)
  w_CK = w_CK, # Mono/serostatus/serotype w (0/1/2)
  alpha_CK = alpha_CK, # Mono/serostatus/serotype alpha (0/1/2)
  tau_K = tau_K, # Mono/serotype-specific tau (0/1)
  
  # Trial specific flags 
  include_beta = include_beta, # Include age-specific nc50? (array[TR], 0-3)
  L_K = L_K, # Mono/serotype-specific L (array[TR], 0/1)
  delta_KJ = delta_KJ # Single/serotype-specific/age delta (array[TR], 0-2)
)
  return(data_list)
}
