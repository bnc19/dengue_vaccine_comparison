format_CYD_stan_data = function(baseline_SP,
                            cases,
                            mu,
                            B,
                            K,
                            R,
                            C,
                            V,
                            J,
                            HI,
                            VCD_years,
                            time,
                            include_pK2,
                            include_eps,
                            include_beta,
                            mono_lc_SN,
                            mono_lc_MU,
                            rho_K,
                            L_K,
                            w_CK,
                            alpha_CK,
                            tau_K,
                            delta_KJ,
                            psi_J,
                            FOI_J,
                            L_mean,
                            L_sd,
                            lower_bound_L,
                            MU_test_SN,
                            MU_symp,
                            enhancement) {
  
T = length(time)
D = length(VCD_years)
  
# SP at baseline 
SP_J = baseline_SP %>%  # make sure youngest first 
  mutate(age = factor(age, levels = c("2-8yrs", "9-16yrs"))) %>% 
  arrange(age)
  
# calculate total hosp over time
HOSP_D = cases$Ho_BVJD %>% 
  group_by(time) %>%  
  summarise(Y=sum(Y))

# calculate total severe over time 
SE_SN_D = cases$Se_VJD %>% 
  group_by(time) %>% 
  summarise(Y=sum(Y))

# Hosp BVJD 
HOSP_BVJD = cases$Ho_BVJD %>% 
  arrange(time, serostatus, arm , age)

m_HOSP_BVJD = array(HOSP_BVJD$Y, dim = c(B * V * J, D))

# Hosp population over time 
Pop_HOSP_BVJD = HOSP_BVJD %>%
  arrange(time, age, arm, serostatus)

m_Pop_HOSP_BVJD = array(Pop_HOSP_BVJD$N, dim = c(B, V, J, D))

# Total hosp population over time 
Pop_HOSP_D = HOSP_BVJD %>%
  group_by(time) %>%  summarise(N = sum(N))

# Hosp BVJK   
HOSP_BVKJ = cases$Ho_BVKJ %>%
  arrange(serostatus, arm , age, serotype) 

m_HOSP_BVKJ = array(HOSP_BVKJ$Y, dim = c(B * V * J * K))

# Symp BVJ 
VCD_BVJ = cases$Sy_BVJ %>%  
arrange(serostatus, arm , age)

m_VCD_BVJ = array(VCD_BVJ$Y, dim = c(B * V * J))
  
# Symp VK (CYD14 and 15 only, PP pop)  
VCD_VK = cases$Sy_VK %>%
  arrange(arm, serotype)

m_VCD_VK = array(VCD_VK$Y, dim = c(V*K))

# Symp population  (CYD14 and 15 only, ITT pop, single time point )  
Pop_VCD_BVJ = VCD_BVJ %>%
  arrange(age, arm, serostatus)

m_Pop_VCD_BVJ = array(Pop_VCD_BVJ$N, dim = c(B, V, J))

# Severe SN VJD 
SE_SN_VJD = cases$Se_VJD %>%  
  arrange(time, arm, age)

m_SE_SN_VJD = array(SE_SN_VJD$Y,dim = c(V * J, D))

# Severe BVJK 
SE_BVKJ = cases$Se_BVKJ %>%   # Vaccine SP is missing (final 8)
  filter(serotype!= "all") %>%  
  arrange(serostatus, arm , age, serotype) 

m_SE_BVKJ = array(SE_BVKJ$Y, dim = c(B * V * J * K - 8))

# Severe SPvJ (seropos vaccine by age)
SE_SPvJ = cases$Se_BVKJ %>%  # Vaccine SP is missing (final 8)
  filter(serotype == "all") %>%  
  arrange(age) 

m_SE_SPvJ = array(SE_SPvJ$Y)

# Initial antibody titres by serostatus and serotype 
mu_t = mu %>%
  select(-X) %>% 
  as.matrix()

# data -------------------------------------------------------------------------
  stan_data = list(
    time = time,
    T = T,
    J = J,
    K = K,
    B = B,
    V = V,
    D = D, 
    C = C,
    HI = HI, 
    R = R, 
    include_pK2=include_pK2,
    include_eps=include_eps,
    include_beta=include_beta,
    mono_lc_SN=mono_lc_SN,
    mono_lc_MU=mono_lc_MU,
    rho_K=rho_K,
    L_K=L_K,
    w_CK=w_CK,
    alpha_CK=alpha_CK,
    tau_K=tau_K,
    delta_KJ=delta_KJ,
    psi_J=psi_J,
    FOI_J=FOI_J,
    L_mean=L_mean,
    L_sd=L_sd,
    lower_bound_L=lower_bound_L,
    MU_test_SN=MU_test_SN,
    MU_symp=MU_symp,
    enhancement=enhancement,
    SP_J = SP_J$SP,
    pop_J =SP_J$N, 
    VCD = sum(cases$Sy_VK$Y),
    HOSP_D = as.numeric(HOSP_D$Y),
    SE_SN_D = as.numeric(SE_SN_D$Y), # SN severe over time 
    SE = sum(cases$Se_BVKJ$Y), # total severe 
    VCD_BVJ = m_VCD_BVJ,  
    VCD_VK = m_VCD_VK,
    HOSP_BVJD = m_HOSP_BVJD,
    HOSP_BVKJ = m_HOSP_BVKJ,
    SE_BVKJ = m_SE_BVKJ,
    SE_SN_VJD = m_SE_SN_VJD,
    SE_SPvJ = m_SE_SPvJ,
    pop_VCD_BVJ = m_Pop_VCD_BVJ,
    pop_HOSP_BVJD = m_Pop_HOSP_BVJD,
    pop_VCD = sum(m_Pop_VCD_BVJ),
    pop_HOSP_D = as.numeric(Pop_HOSP_D$N),
    pop_HOSP = mean(Pop_HOSP_D$N),
    mu = mu_t
  )
  
  return(stan_data)
}
