data{
   int<lower = 1> C; // No. current immune status
   int<lower = 1> B; // No. baseline serostatus
   int<lower = 1> K; // No. serotypes
   int<lower = 1> D_Q; // No. Qdenga vcd data
   int<lower = 1> D_De; // No. Deng vcd data
   int<lower = 1> V; // No. trial arms 
   int<lower = 1> J_Q; // No. age groups
   int<lower = 1> J_De; // No. age groups
   int<lower = 1> J_Bu; // No. age groups
   int<lower = 1> A; // No. age vcd data 
   int<lower = 1> T_Q; // No. Q model time points 
   int<lower = 1> T_De; // No. D model time points
   int<lower = 1> T_Bu; // No. D model time points
   int<lower = 1> R; // No. outcomes 
   int<lower = 1> TR; // Number of vaccine trials (Q, D, B)

// model time 
array[T_Q] int<lower=0> time_Q;
array[T_De] int<lower=0> time_De;
array[T_Bu] int<lower=0> time_Bu;

int<lower = 0> HI; // period of heterotypic immunity
   
array[TR, B, K] real<lower=0> mu; // trial specific titres at T0
  
   // Dengvaxia data 
   array[J_De] int<lower=0> SP_J_De;     // baseline seropositive
   array[J_De] int<lower=0> pop_J_De;    // baseline pop
   int<lower=0> VCD_De;                  // total VCD 
   array[D_De] int<lower=0> HOSP_D_De;   // total hosp
   
   array[V*K]     int<lower=0> VCD_VK_De;        
   array[B*V*J_De]   int<lower=0> VCD_BVJ_De;    
   array[B*V*J_De,D_De] int<lower=0> HOSP_BVJD_De;    
   array[B*V*K*J_De] int<lower=0> HOSP_BVKJ_De;       
   array[B,V,J_De]   int<lower=0> pop_VCD_BVJ_De;    
   array[B,V,J_De,D_De] int<lower=0> pop_HOSP_BVJD_De; 
   
   int<lower=0> pop_VCD_De ; 
   int<lower=0> pop_HOSP_De; // mean hosp pop across trial 
   array[D_De] int<lower=0> pop_HOSP_D_De ; // hosp pop over time 
   
   // Qdenga data 
   array[J_Q] int<lower=0> SP_J_Q;     // baseline seropositive
   array[J_Q] int<lower=0> pop_J_Q;    // baseline pop
   array[D_Q] int<lower=0> VCD_D_Q;    // total VCD 
   array[D_Q] int<lower=0> HOSP_D_Q;   // total hosp  
   
   array[B*V*K,D_Q] int<lower=0> VCD_BVKD_Q;  
   array[B*V*J_Q,A] int<lower=0> VCD_BVJA_Q;  
   array[K*J_Q,2] int<lower=0> VCD_KJ2_Q;     //  12 and 24 months
   array[B*V*J_Q,A] int<lower=0> HOSP_BVJA_Q; 
   array[B*V*K,4] int<lower=0> HOSP_BVK4_Q; 
   array[K*J_Q,2] int<lower=0> HOSP_KJ2_Q;    //  12 and 24 months
   array[B,V,J_Q,D_Q] int<lower=0> pop_Q;       
   array[B,V,D_Q] int<lower=0> pop_BVD_Q;    
   array[B,V] int<lower=0> N_VCD_BV5_Q;     // cases in t_5 (not age-specific)
   
   // Butantan-DV data
   
   array[J_Bu] int<lower=0> SP_J_Bu; // baseline seropositive 
   array[J_Bu] int<lower=0> pop_J_Bu; // baseline pop
   int<lower=0> VCD_Bu; // total VCD 
   array[B*V*2] int<lower=0> VCD_BVK_Bu; 
   array[B*V*J_Bu] int<lower=0> VCD_BVJ_Bu; 
   array[B,V,J_Bu] int<lower=0> pop_BVJ_Bu;  
   int<lower=0> pop_Bu; // total pop  

  
  // FLAGS
   
   // shared flags 
   int<lower = 0, upper = 1> share_n_param;     // share hs, hl, ts parameters across trials (T/F)
   int<lower = 0, upper = 1> uniform;          // uniform priors on probabilities of reporting (T/F)
   int<lower = 0, upper = 1> include_pK3;      // include serotype specific p? (T/F)
   int<lower = 0, upper = 1> include_eps;      // include enhanced secondary hosp? (T/F)
   real<lower = 0> L_mean;                     // mean (trunc normal) L prior
   real<lower = 0> L_sd;                       // sd (trunc normal) L prior 
   int<lower = 0, upper = 1> enhancement;      // 0 for no vac enhancement of SN, 1 for enhancement 
   int<lower = 0, upper = 2> mono_lc_SN;       // 0 for serotype specific lc SN / 1 for single / 2 for offset from MO
   int<lower = 0, upper = 2> mono_lc_MU;       // 0 for serotype specific lc MO / 1 for single / 2 for offset from MO
   int<lower = 0, upper = 1> rho_K;           // 0 for mono rho, 1 for serotype rho 
   int<lower = 0, upper = 2> w_CK;            // 0 for mono w, 1 for serostatus w, 2 for serotype w 
   int<lower = 0, upper = 2> alpha_CK;        // 0 for mono alpha, 1 for serostatus alpha, 2 for serotype alpha 
   int<lower = 0, upper = 1> tau_K;           // 0 for mono tau, 1 for serotype tau

// vaccine specific flags (age and hospitalisation)
array[TR] int<lower = 0, upper = 3> include_beta;   // include age-specific nc50?: 1= change age groups 1 & 2 for both outcomes / 2 = only age grp 1, sep for each outcome / 3 = only age grp 1 for both outcomes 
array[TR] int<lower = 0, upper = 1> L_K;            // 0 for mono L, 1 for serotype K 
array[TR] int<lower = 0, upper = 2> delta_KJ;       // 0 for single delta, 1 for serotype-specific delta, 2 for age delta
array[TR]  int<lower = 0, upper = 1> inc_FOIJ;     // 0 for no age-specific FOI, 1 for youngest age-specific 

}

transformed data {

  // Define Qdenga data 
  array[B,J_Q,2] real pop_BJD_Q; 
  array[J_Q,2] real pop_JD_Q; 
  array[B,V,D_Q] real pop_BVD_real_Q;     
  array[B,V] real pop_BV_Q;  
  array[B,V] real pop_BV_K_Q;  
  array[V,D_Q] int pop_VD_Q;
  array[B,D_Q] int pop_BD_Q;
  array[D_Q] int pop_D_Q;
  vector<lower=0>[D_Q] v_pop_D_Q;
  vector[J_Q] v_SP_J_Q = to_vector(SP_J_Q);     // baseline seropositive
  vector[J_Q] v_pop_J_Q = to_vector(pop_J_Q);    // baseline pop
  vector[D_Q] v_VCD_D_Q = to_vector(VCD_D_Q);    // total VCD
  vector[D_Q] v_HOSP_D_Q = to_vector(HOSP_D_Q);   // total hosp
  matrix[B*V*K,D_Q] v_VCD_BVKD_Q=to_matrix(VCD_BVKD_Q);  // VCD
  matrix[B*V*J_Q,A] v_VCD_BVJA_Q=to_matrix(VCD_BVJA_Q);  // VCD
  matrix[K*J_Q,2] v_VCD_KJ2_Q=to_matrix(VCD_KJ2_Q);     // VCD at 12 and 24 months
  matrix[B*V*J_Q,A] v_HOSP_BVJA_Q=to_matrix(HOSP_BVJA_Q); // hosp
  matrix[B*V*K,4] v_HOSP_BVK4_Q=to_matrix(HOSP_BVK4_Q); // hosp
  matrix[K*J_Q,2] v_HOSP_KJ2_Q =to_matrix(HOSP_KJ2_Q);    // hosp at 12 and 24 months

  // Define Dengvaxia data 
 array [V,J_De,D_De] int pop_HOSP_VJD_De;
 array [V,J_De] int pop_VCD_VJ_De;
 

// QDENGA 
// POPULATION AGGREGATION 
for(b in 1:B)
 for(j in 1:J_Q){
 pop_BJD_Q[b,j,1] = sum(pop_Q[b, ,j,1]);
 pop_BJD_Q[b,j,2] = sum(pop_Q[b, ,j,3]);
 }
 
for(j in 1:J_Q)
 for(d in 1:2)
 pop_JD_Q[j,d] = sum(pop_BJD_Q[ ,j,d]); 

for(b in 1:B)
 for(v in 1:V)
  for(d in 1:D_Q)
   pop_BVD_real_Q[b,v,d] = pop_BVD_Q[b,v,d]; 

for(b in 1:B)
 for(v in 1:V)
 pop_BV_Q[b,v] = mean(pop_BVD_real_Q[b,v, ]);
 
for(b in 1:B)
 for(v in 1:V)
 pop_BV_K_Q[b,v] = mean(pop_BVD_real_Q[b,v,4:6]); // serotype hosp data from 36months  
 
for(v in 1:V)
 for(d in 1:D_Q)
 pop_VD_Q[v,d] = sum(pop_BVD_Q[,v,d]);

for(b in 1:B)
 for(d in 1:D_Q)
  pop_BD_Q[b,d] = sum(pop_BVD_Q[b, ,d]); 

for(d in 1:D_Q){
 pop_D_Q[d] = sum(pop_BD_Q[,d]);
 }
  
  v_pop_D_Q = to_vector(pop_D_Q);
  
  // DENGVAXIA 
  // pop not by serostatus 
  for (v in 1:V)
    for (j in 1:J_De) {
      pop_VCD_VJ_De[v,j] = sum(pop_VCD_BVJ_De[ ,v,j]); // Moved outside d-loop
      for (d in 1:D_De)
        pop_HOSP_VJD_De[v,j,d] = sum(pop_HOSP_BVJD_De[ ,v,j,d]); 
    }
}

parameters{
// FOI 
array[2] real<lower = 0> lambda_K_Bu;            // FOR B  (not time specific, only 2 serotypes)
array[K, D_Q] real<lower = 0> lambda_D_Q;        // FOI Q
array[D_De] real<lower = 0> lambda_D_De;         // FOI De not serotype specific 
simplex[K] theta;                                // proportion of each serotype circulating in the Dengvaxia trial               
  
array[TR, J_Q] real<lower = 0, upper = 1> p;     // probability of exposure before T0 (J = 3 for Dengvaxia ignored)
array[TR, K] real<lower = 0, upper = 1> pK3;     // serotype-specific p in age group 3

array[TR, K] real<lower = 0, upper = 1> delta;   // prob hosp K is serotype or age
array[TR, K] real<lower = 0> L;                  // enhancement parameter 
array[TR, K] real<lower = 0> tau;                // multiply L for hosp
array[TR, K] real<lower = 0> w;                  // shape parameter 
array[TR, C, K] real<lower = 0> lc;              // 50% symptomatic protection
array[TR, K] real<lower = 0> alpha;              // reduction in titre for protection against hosp
array[TR, (J_Q - 1)] real beta;                  // increase in lc50 for protection in younger

array[TR] real<lower = 0, upper = 1> sens;       // baseline test sensitivity 
array[TR] real<lower = 0, upper = 1> spec;       // baseline test specificity 
array[TR] real<lower = 0> omega;                 // reduction in titre for MU vs. MO
array[TR] real kappa;                            // increase/decrease in titre for SN vs. MO
array[TR] real<lower = 0>  FOI_J1;               // age specific FOI 

// JOINT
real<lower = 1> epsilon;                         // prob hosp 2' rel to 1'/3'/4'
real<lower = 0, upper = 1> gamma;                // prob symptomatic 2' 
array[K] real<lower = 0> rho;                    // 1/risk symptomatic 1' rel to 2'
real<lower = 0, upper = 1> phi;                  // prob symptomatic 3'/4' rel to 1' 
  
// Testing whether antibody shape parameters can be shared 
array[TR, B] real<lower = 0> hs;                     // short term decay 
array[TR]    real<lower = 12> hl;                    // long term decay 
array[TR, B] real<lower = 0> ts;                     // time switch decay rate

}

transformed parameters{
  
// #########################  JOINT PARAMETERS // ######################### 
vector<lower = 0>[K] rhoT;   // Ensure rhoT is positive
vector<lower = 0, upper = 1>[K] gammaT;  // Ensure gammaT is between 0 and 1
vector<lower = 0, upper = 1>[K] phiT;    // Ensure phiT is between 0 and 1

  // model selection 
if(rho_K==1) {
  for(k in 1:K) {
    rhoT[k] = rho[k];
    gammaT[k] = gamma;
    phiT[k] = phi;
  }
} else {
  for(k in 1:K) {
    rhoT[k] = rho[1];
    gammaT[k] = gamma; 
    phiT[k] = phi;
  }
}
// rho is reduction in primary symp compared to sec
// gamma is prob sec symp but gammaT is prob primary (divided by rho)
    for(k in 1:K) gammaT[k] /= rhoT[k]; 
  
  
// ######################  QDENGA PARAMETERS TO ACCESS // ######################
  array[J_Q] real pSP_Q;; // prob SP (baseline for likelihood)
  real ll_Q; // log-likelihood passed to model block and then added to target
  array[K, J_Q, R, D_Q] real<lower=0> C_KJRD_Q;
  array[B, V, K, J_Q, R, D_Q] real<lower=0> C_BVKJRD_Q;
  array[B, K, T_Q] real<lower=0> n_Q;  // titres
  array[C, V, K, J_Q, T_Q] real<lower=0> RR_symp_Q;
  array[C, V, K, J_Q, T_Q] real<lower=0> RR_hosp_Q;

  
  // censor population for last time interval 
  array[B, V, J_Q, D_Q] real<lower=0> pop_BVJD_Q;
  // Matrix distribution of cases
  matrix[K*J_Q, 2]   pD_KJ2_Q;    
  matrix[B*V*K, D_Q] pD_BVKD_Q;   
  matrix[B*V*J_Q, A] pD_BVJA_Q;   
  matrix[K*J_Q, 2]   pH_KJ2_Q;    
  matrix[B*V*K, 4]   pH_BVK4_Q;   
  matrix[B*V*J_Q, A] pH_BVJA_Q;
  real<lower=0> shape1_Q; // pk3 prior
  real<lower=0> shape2_Q; 
  //  distribution of cases
  array[R, D_Q] real pC_RD_Q;
  
  
// ######################  DENGVAXIA PARAMETERS TO ACCESS ######################
  array[J_De] real<lower=0, upper=1> pSP_De; // prob SP (baseline for likelihood)
  real ll_De; // log-likelihood passed to model block and then added to target
  
  array[V,K,J_De]  real<lower = 0> Sy_VKJ_De;
  array[B,V,K,J_De]   real<lower = 0> Sy_BVKJ_De;          // all cases 
  array[B,V,K,J_De,D_De] real<lower = 0> Ho_BVKJD_De;  

  array[B,K,T_De]  real<lower = 0> n_De;                // titres
  array[C,V,K,J_De,T_De] real<lower = 0> RR_symp_De;                  
  array[C,V,K,J_De,T_De] real<lower = 0> RR_hosp_De;

  // matrix distribution of cases
  vector [B*V*J_De] mD_BVJ_De ;    
  vector [V*K]      mD_VK_De ;  
  matrix [B*V*J_De,D_De] mH_BVJD_De ;    
  vector [B*V*K*J_De] mH_BVKJ_De ;   
  
  real<lower=0> shape1_De; // pk3 prior
  real<lower=0> shape2_De; 
  
  //  distribution of cases
  real<lower=0, upper=1> pC_De;
  array[D_De] real<lower=0, upper=1> pH_De;

// ######################  BUTANTAN PARAMETERS TO ACCESS #######################

  array[J_Bu] real<lower=0, upper=1> pSP_Bu; // prob SP (baseline for likelihood)
  real ll_Bu; // log-likelihood passed to model block and then added to target
  
  array[B,V,2,J_Bu] real<lower = 0> Sy_BVKJ_Bu;
  array[V,2,J_Bu]   real<lower = 0> Sy_VKJ_Bu;

  array[B,2,T_Bu] real<lower = 0> n_Bu; // titres
  array[C,V,2,J_Bu,T_Bu] real<lower = 0> RR_symp_Bu;                  

  // matrix distribution of cases
  vector[B*V*2] mD_BVK_Bu;    
  vector[B*V*J_Bu] mD_BVJ_Bu;  
  
  real<lower=0> shape1_Bu; // pk2 prior
  real<lower=0> shape2_Bu; 
  
  // distribution of cases
  real<lower=0, upper=1> pC_Bu;

  
// ################ Biphasic titre decay rates shared? #########################

array[TR, B] real pi_1; // decay rate 1 
array[TR]    real pi_2;          // decay rate 2
array[TR,B ] real ts2; 

if(share_n_param == 1){
  for(t in 1:TR){
    for(b in 1:B) pi_1[t, b] = -log(2) / hs[1, b];
    pi_2[t] = -log(2) / hl[1] ; 
    for(b in 1:B) ts2[t, b] = ts[1, b];
  }
  } else {
    for(t in 1:TR){
      for(b in 1:B) pi_1[t, b] = -log(2) / hs[t, b];
      pi_2[t] = -log(2) / hl[t] ; 
      for(b in 1:B) ts2[t, b] = ts[t, b];  
    }
    }
// ######################## BUTANTAN BLOCK  ####################################
  { // this { defines a block within which variables declared are local (can't have lower or upper)
  
  real lambda_m;
  real pm; 
  
  // p by lambda 
  lambda_m = mean(lambda_K_Bu); // mean across D1 and D1
  pm = 1-exp(-38.5*12*lambda_m); // always prob exposure to a single serotype 
 
  if(pm<1e-3) pm=1e-3;
  if(pm>0.99) pm=0.99;
  real vari = 0.02*pm *(1-pm) ; 
  shape1_Bu = (((1-pm) / vari) - (1/pm)) * (pm^2) ;
  shape2_Bu = shape1_Bu * (1/pm - 1);

  array [K,J_Bu] real p_KJ; // serotype age prob exposure 
  
  array [J_Bu] real TpSN;              // true SN
  array [K,J_Bu] real TpMO;            // true monotypic  
  array [(K+2),J_Bu] real TpMU2;       // true multitypic
  array [K,J_Bu] real TpMU3;           // true multitypic  
  array [J_Bu]   real TpMU4;           // true multitypic  
  
  array [B,V,J_Bu,(T_Bu+1)]  real pSN;       // prob SN time plus baseline 
  array [B,V,K,J_Bu,(T_Bu+1)] real pMO ;     // prob monotypic time plus baseline 
  array [B,V,(K+2),J_Bu,(T_Bu+1)] real pMU2; // prob multitypic time plus baseline 
  array [B,V,K,J_Bu,(T_Bu+1)] real pMU3;     // prob multitypic time plus baseline 
  
  array [C,2,T_Bu] real n_C;                   // allow for multitypic titres 
  array [C,2]      real L_C;                   // MO and MU L = 1 
  array [C,2,J_Bu] real nc50;                  // nc50
  array [2,J_Bu,T_Bu] real lambda;                // FOI 

  array [B,V,2,J_Bu,T_Bu] real Inc1;
  array [B,V,2,J_Bu,T_Bu] real Inc2;
  array [B,V,2,J_Bu,T_Bu] real Inc3;
  array [B,V,2,J_Bu,T_Bu] real Inc4;
     
  array [B,V,2,J_Bu,T_Bu] real D1;
  array [B,V,2,J_Bu,T_Bu] real D2;
  array [B,V,2,J_Bu,T_Bu] real D34;
  array [B,V,2,J_Bu]   real Di;
  
  // cases 
  array [2,J_Bu]   real Sy_KJ;
  array [J_Bu]     real Sy_J;
  real Sy;
  
  // prop cases
  array [B,V,2]   real pD_BVK; 
  array [B,V,J_Bu]   real pD_BVJ; 

// BIPHASIC TITRES 

for(b in 1:B)
 for(k in 1:2)
  for(t in 1:T_Bu)
   n_Bu[b,k,t] = mu[3,b,k] * (exp(pi_1[3,b] * time_Bu[t] + pi_2[3] * ts2[3,b]) + exp(pi_2[3] * time_Bu[t] + pi_1[3,b] * ts2[3,b])) / (exp(pi_1[3,b] * ts2[3,b]) + exp(pi_2[3] * ts2[3,b])) ; 

// MO and MU have SP titres 
 for(k in 1:2)
  for(t in 1:T_Bu){
    n_C[1,k,t] =  n_Bu[1,k,t];
    n_C[2,k,t] =  n_Bu[2,k,t]; 
    n_C[3,k,t] =  n_Bu[2,k,t];
  }
  
      
// INITIAL CONDITIONS
for(m in 1:K)
 for(j in 1:J_Bu) p_KJ[m,j] = p[3,j] ; // EVEN IF pk3 = 1, Butantan-DV always just p

// this defines the true serostatus populations 
 for(j in 1:J_Bu){
    TpSN[j] = (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1
    TpMO[1,j] = p_KJ[1,j]   * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 2
    TpMO[2,j] = p_KJ[2,j]   * (1-p_KJ[1,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 3
    TpMO[3,j] = p_KJ[3,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 4
    TpMO[4,j] = p_KJ[4,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 1,2
    TpMU2[1,j] = p_KJ[1,j] * p_KJ[2,j] * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1,3
    TpMU2[2,j] = p_KJ[1,j] * p_KJ[3,j] * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 1,4
    TpMU2[3,j] = p_KJ[1,j] * p_KJ[4,j] * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 2,3
    TpMU2[4,j] = p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[1,j]) * (1-p_KJ[4,j]);
    // 2,4
    TpMU2[5,j] = p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[3,j]);
    // 3,4 
    TpMU2[6,j] = p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[2,j]);
    // not 1 
    TpMU3[1,j] = p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) ;
    // not 2
    TpMU3[2,j] = p_KJ[1,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[2,j]) ;
    // not 3
    TpMU3[3,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[3,j]) ;
    // not 4 
    TpMU3[4,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[4,j]) ;
    // all 
    TpMU4[j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j];
 }

// this accounts for imperfect test performance 
for(v in 1:V)
 for(j in 1:J_Bu)
   for(k in 1:K){
    pSN[1,v,j,1] = spec[3] * TpSN[j] ;
    pSN[2,v,j,1] = (1-spec[3]) * TpSN[j] ;
    pMO[1,v,k,j,1] = (1- sens[3]) * TpMO[k,j] ;
    pMO[2,v,k,j,1] = sens[3] * TpMO[k,j] ;
   }

  for(v in 1:V)
   for(j in 1:J_Bu) {
    for(k in 1:6){  // there are 6 MU_2 combinations  
        pMU2[1,v,k,j,1] = 0;
        pMU2[2,v,k,j,1] = TpMU2[k,j] ; 
    }
    for(k in 1:K){  // there are 4 MU_3 combinations 
        pMU3[1,v,k,j,1] = 0; 
        pMU3[2,v,k,j,1] = TpMU3[k,j] ;
    } }


// probabilities of testing seropositive 

 for(j in 1:J_Bu)  pSP_Bu[j] = (1-spec[3]) * TpSN[j] + sens[3] * sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j];


// VACCINE RISK RATIO 

// age-group titre offsets 
array [2] real e_beta = {exp(beta[3,1]), exp(beta[3,2])};
      
if(mono_lc_SN == 1){ // SN, oldest
  for(k in 1:2) nc50[1,k,3] = exp(lc[3,1,1]); 
   } else if (mono_lc_SN == 0) {
  for(k in 1:2) nc50[1,k,3] = exp(lc[3,1,k]); 
} else if(mono_lc_SN == 2){
  for(k in 1:2)  nc50[1,k,3] = exp(lc[3,2,k]) * exp(kappa[3]); 
}

if(mono_lc_MU == 1){ // MU, oldest
 for(k in 1:2)  nc50[3,k,3] = exp(lc[3,3,1]); 
   } else if(mono_lc_MU == 0) {
 for(k in 1:2) nc50[3,k,3] = exp(lc[3,3,k]); 
} else if (mono_lc_MU == 2){
   for(k in 1:2)  nc50[3,k,3] = exp(lc[3,2,k]) * exp(-omega[3]); 
}

// MO, oldest
for(k in 1:2)  nc50[2,k,3] = exp(lc[3,2,k]); // MO, serotype, oldest

// age group offsets 
if(include_beta[3] == 1){
  for(c in 1:C)
   for(k in 1:2){
   nc50[c,k,1] =   nc50[c,k,3] * e_beta[1]; // youngest
   nc50[c,k,2] =   nc50[c,k,3] * e_beta[2]; // middle 
 } 
} else if(include_beta[3] == 2) {
  for(c in 1:C)
   for(k in 1:2){
   nc50[c,k,1] =   nc50[c,k,3] * e_beta[1]; // youngest
   nc50[c,k,2] =   nc50[c,k,3]; // middle 
   }
 } else {
  for(c in 1:C)
   for(k in 1:2){ // no beta offset 
   nc50[c,k,1] = nc50[c,k,3];
   nc50[c,k,2] = nc50[c,k,3];   
  }
}



// Enhancement 
for(k in 1:2){
  if(enhancement == 1){ 
   if(L_K[3] == 0){
   L_C[1,k] = 1 + L[3,1]; // symp 
  } else if (L_K[3] == 1){
   L_C[1,k] = 1 + L[3,k]; // symp 
 } 
 } else if (enhancement == 0){ // no SN enhancement 
  L_C[1,k] = 1 ; 
 } 
// no enhancement if seropositive 
        L_C[2,k] = 1; // MO
        L_C[3,k] = 1; // MU
   }
 
// Risk ratios - w can be mono, serostatus or serotype dependent  
for(c in 1:C)
 for(k in 1:2)
  for(j in 1:J_Bu)
   for(t in 1:T_Bu){
     if(w_CK == 0){
     RR_symp_Bu[c,2,k,j,t] = L_C[c,k] / (1 +  (n_C[c,k,t] / nc50[c,k,j])^w[3,1]) ;
     } else if(w_CK == 1){
     RR_symp_Bu[c,2,k,j,t] = L_C[c,k] / (1 +  (n_C[c,k,t] / nc50[c,k,j])^w[3,c]) ;
     } else{
     RR_symp_Bu[c,2,k,j,t] = L_C[c,k] / (1 +  (n_C[c,k,t] / nc50[c,k,j])^w[3,k]) ;
     }
     RR_symp_Bu[c,1,k,j,t] =  1 ;
    }

 // FOI 
 
array [J_Bu] real FOI_J = {FOI_J1[3], 1, 1}  ; // scale youngest group only 

if(inc_FOIJ[3] == 0) {
for(k in 1:2)
 for(j in 1:J_Bu)
  for(t in 1:T_Bu)  lambda[k,j,t] = exp(-lambda_K_Bu[k]) ; 
} else{
for(k in 1:2)
 for(j in 1:J_Bu)
  for(t in 1:T_Bu)  lambda[k,j,t] = exp(-lambda_K_Bu[k] * FOI_J[j]) ; 
}

// SURVIVAL MODEL - prob of surviving each time point without infection - only D1 and D2 lambda 
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_Bu)
   for(t in 1:T_Bu) pSN[b,v,j,(t+1)] = lambda[1,j,t] * lambda[2,j,t] * pSN[b,v,j,t];
   
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J_Bu) {
    for(t in 1:T_Bu) pMO[b,v,1,j,(t+1)] = lambda[2,j,t] * pMO[b,v,1,j,t] ; // escape 2 
    for(t in 1:T_Bu) pMO[b,v,2,j,(t+1)] = lambda[1,j,t] * pMO[b,v,2,j,t] ; // escape 1 
    for(t in 1:T_Bu) pMO[b,v,3,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t] * pMO[b,v,3,j,t] ; // escape 1 and 2 
    for(t in 1:T_Bu) pMO[b,v,4,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t] * pMO[b,v,4,j,t] ; // escape 1 and 2 
    
    for(t in (HI+1):T_Bu) 
     for(k in 1:2) pMO[b,v,k,j,(t+1)] += (1 - gammaT[k] * RR_symp_Bu[1,v,k,j,(t-HI)]) * (1 - lambda[k,j,(t - HI)]) * pSN[b,v,j,(t-HI)] ; // new D1 and D2 infections during trial
    }
    
    
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J_Bu) {
      for(t in 1:T_Bu){
          // MU_12
          pMU2[b,v,1,j,(t+1)] =  pMU2[b,v,1,j,t] ; 
           // MU_13
          pMU2[b,v,2,j,(t+1)] = lambda[2,j,t] *  pMU2[b,v,2,j,t] ; 
           // MU_14
          pMU2[b,v,3,j,(t+1)] = lambda[2,j,t] *  pMU2[b,v,3,j,t] ; 
          // MU_23
          pMU2[b,v,4,j,(t+1)] = lambda[1,j,t] *  pMU2[b,v,4,j,t] ; 
          // MU_24
          pMU2[b,v,5,j,(t+1)] = lambda[1,j,t] *  pMU2[b,v,5,j,t] ; 
          // MU_34
          pMU2[b,v,6,j,(t+1)] = lambda[1,j,t] * lambda[2,j,t] *  pMU2[b,v,6,j,t] ; 
      }
       for(t in (HI+1):T_Bu) { // new D1 and D2 infections 
          // MU_12
          pMU2[b,v,1,j,(t+1)] += (1 - gammaT[1] * rhoT[1] * RR_symp_Bu[2,v,1,j,(t-HI)]) * (1 - lambda[1,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] + (1 - gammaT[2] * rhoT[2] * RR_symp_Bu[2,v,2,j,(t-HI)]) *  (1-lambda[2,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_13
          pMU2[b,v,2,j,(t+1)] += (1 - gammaT[1] * rhoT[1] * RR_symp_Bu[2,v,1,j,(t-HI)]) * (1 - lambda[1,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] ;
          // MU_14
          pMU2[b,v,3,j,(t+1)] += (1 - gammaT[1] * rhoT[1] * RR_symp_Bu[2,v,1,j,(t-HI)]) * (1 - lambda[1,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] ;
          // MU_23
          pMU2[b,v,4,j,(t+1)] += (1 - gammaT[2] * rhoT[2] * RR_symp_Bu[2,v,2,j,(t-HI)]) * (1 - lambda[2,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)]  ;
          // MU_24
          pMU2[b,v,5,j,(t+1)] += (1 - gammaT[2] * rhoT[2] * RR_symp_Bu[2,v,2,j,(t-HI)]) * (1 - lambda[2,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] ;
        }}

for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J_Bu){
     for(t in 1:T_Bu){ 
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] = lambda[1,j,t]  *  pMU3[b,v,1,j,t] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] = lambda[2,j,t]  *  pMU3[b,v,2,j,t] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] = pMU3[b,v,3,j,t] ; 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] = pMU3[b,v,4,j,t] ; 
     }
    for(t in (HI+1):T_Bu) {
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] +=  (1 - phiT[2] * gammaT[2] * RR_symp_Bu[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; // 34 then 2 infection 
        // MU3 -2           
        pMU3[b,v,2,j,(t+1)] +=  (1 - phiT[1] * gammaT[1] * RR_symp_Bu[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; // 34 then 1 infection 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] +=  (1 - phiT[2] * gammaT[2] * RR_symp_Bu[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] * RR_symp_Bu[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] ;  // 14 then 2 or 24 then 1 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] +=  (1 - phiT[2] * gammaT[2] * RR_symp_Bu[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[1] * gammaT[1] * RR_symp_Bu[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] ;  // 13 then 2 or 23 then 1 
      }
    }  

// calculate infection incidence, symptomatic, hospitalised 
for(b in 1:B)
  for(v in 1:V)
    for(j in 1:J_Bu) {
      for(t in 1:T_Bu){
         Inc3[b,v,1,j,t] = (1 - lambda[1,j,t]) * (pMU2[b,v,4,j,t] + pMU2[b,v,5,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,2,j,t] = (1 - lambda[2,j,t]) * (pMU2[b,v,2,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,6,j,t]);  
         
         for(k in 1:2){
           Inc1[b,v,k,j,t] = (1 - lambda[k,j,t]) * pSN[b,v,j,t];  
           Inc2[b,v,k,j,t] = (1 - lambda[k,j,t]) * (sum(pMO[b,v, ,j,t]) - pMO[b,v,k,j,t]);
           Inc4[b,v,k,j,t] = (1 - lambda[k,j,t]) * pMU3[b,v,k,j,t];  
           D1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * RR_symp_Bu[1,v,k,j,t] ;
           D2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k] * RR_symp_Bu[2,v,k,j,t];
           D34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * RR_symp_Bu[3,v,k,j,t];  
         }}}

// Aggregate to match published time points e.g. 1-24 months 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:2)
   for(j in 1:J_Bu){
   Di[b,v,k,j] = sum(D1[b,v,k,j,1:T_Bu])  + sum(D2[b,v,k,j,1:T_Bu]) + sum(D34[b,v,k,j,1:T_Bu]); 
   }

// Cases 
for(b in 1:B) 
 for(v in 1:V) 
  for(k in 1:2)
   for(j in 1:J_Bu) 
          Sy_BVKJ_Bu[b,v,k,j] = Di[b,v,k,j] * sum(pop_BVJ_Bu[ ,v,j]); 

// sum over symp
for(v in 1:V) 
 for(k in 1:2) 
  for(j in 1:J_Bu) 
  Sy_VKJ_Bu[v,k,j] = sum(Sy_BVKJ_Bu[ ,v,k,j]) ; 

    
for(k in 1:2) 
 for(j in 1:J_Bu) 
 Sy_KJ[k,j] = sum(Sy_VKJ_Bu[ ,k,j]) ;  
   

for(j in 1:J_Bu) 
 Sy_J[j] = sum(Sy_KJ[ ,j]) ; 

 
 Sy = sum(Sy_J) ; 

for(b in 1:B)
 for(v in 1:V)
   for(k in 1:2)
    for(j in 1:J_Bu){
  pD_BVK[b,v,k] = sum(Sy_BVKJ_Bu[b,v,k, ]) / Sy ; // serotype serostatus symp
  pD_BVJ[b,v,j] = sum(Sy_BVKJ_Bu[b,v, ,j]) / Sy ;   // age serostatus symp 
  }

// matrix for likelihood function 

for(k in 1:2){ // 2 serotypes 
  mD_BVK_Bu[k]   =   pD_BVK[1,1,k] ; // SN C
  mD_BVK_Bu[k+2] =   pD_BVK[1,2,k] ; // SN V
  mD_BVK_Bu[k+4] =   pD_BVK[2,1,k] ; // SP C
  mD_BVK_Bu[k+6] =   pD_BVK[2,2,k] ; // SP V
}

for(j in 1:J_Bu){ // 3 ages 
   mD_BVJ_Bu[j]     =  pD_BVJ[1,1,j]; //  SN C 
   mD_BVJ_Bu[j+3]   =  pD_BVJ[1,2,j]; //  SN V
   mD_BVJ_Bu[j+6]   =  pD_BVJ[2,1,j]; //  SP C 
   mD_BVJ_Bu[j+9]   =  pD_BVJ[2,2,j]; //  SP V
  }
  
 pC_Bu = Sy / pop_Bu; 


// log likelihood 
 ll_Bu=0;
 ll_Bu += binomial_lpmf(SP_J_Bu | pop_J_Bu, pSP_Bu);
 ll_Bu += binomial_lpmf(VCD_Bu  | pop_Bu, pC_Bu);
 ll_Bu += multinomial_lpmf(VCD_BVK_Bu  | mD_BVK_Bu) ;
 ll_Bu += multinomial_lpmf(VCD_BVJ_Bu  | mD_BVJ_Bu) ;
 
}
// ######################### QDENGA BLOCK  #####################################

  { // { defines a block within which variables declared are local (can't have lower or upper)
  
  array[D_Q] real lambda_mD;
  real lambda_m;
  real pm; 
  
  // contrain pk3 / p3 by lambda 
  for(d in 1:D_Q) lambda_mD[d] = mean(lambda_D_Q[ ,d]); // mean FOI at time D 
  lambda_m = (lambda_mD[1] * 12 + lambda_mD[2] * 6 + lambda_mD[3] * 6 + lambda_mD[4] * 12 + lambda_mD[5] * 12+ lambda_mD[6] * 6) / 54 ; // weighted mean across time 
  pm = 1-exp(-14*12*lambda_m); // always prob exposure to a single serotype 
  if(pm<1e-3) pm=1e-3;
  if(pm>0.99) pm=0.99;
  shape1_Q = 49 * pm ;
  shape2_Q = shape1_Q * (1/pm - 1);

// Age-related hazard and serotype proportion
array[K] real hK3; // age3 historic hazard
array[K] real q;   // historic serotype proportion 
array[K] real hK1; // age1 historic hazard 
array[K] real hK2; // age2 historic hazard 

// Serotype exposure probabilities
array[K, J_Q] real p_KJ; // serotype age prob exposure

// True serostatus probabilities at baseline 
array[J_Q] real TpSN;       // true SN
array[K, J_Q] real TpMO;    // true monotypic  
array[K+2, J_Q] real TpMU2; // true multitypic
array[K, J_Q] real TpMU3;   // true multitypic  
array[J_Q] real TpMU4;      // true multitypic  

// Time-dependent probabilities
array[B, V, J_Q, T_Q+1] real pSN;      // prob SN time plus baseline 
array[B, V, K, J_Q, T_Q+1] real pMO;   // prob monotypic time plus baseline 
array[B, V, K+2, J_Q, T_Q+1] real pMU2; // prob multitypic time plus baseline 
array[B, V, K, J_Q, T_Q+1] real pMU3;  // prob multitypic time plus baseline 

// RR parameters 
array[C, K, T_Q] real n_C;  // allow for multitypic titres 
array[C, K, R] real L_C;    // MO and MU L = 1 
array[C, K, J_Q, R] real nc50; // nc50

// FOI 
array[K, T_Q] real lambda; // FOI 

// Incidence and disease outcomes
array[B, V, K, J_Q, T_Q] real Inc1;
array[B, V, K, J_Q, T_Q] real Inc2;
array[B, V, K, J_Q, T_Q] real Inc3;
array[B, V, K, J_Q, T_Q] real Inc4;
array[B, V, K, J_Q, T_Q] real D1;
array[B, V, K, J_Q, T_Q] real D2;
array[B, V, K, J_Q, T_Q] real D34;
array[B, V, K, J_Q, T_Q] real H1;
array[B, V, K, J_Q, T_Q] real H2;
array[B, V, K, J_Q, T_Q] real H34;

// Disease cases and hospitalisations
array[B, V, K, J_Q, D_Q] real Di;
array[B, V, K, J_Q, D_Q] real H;

// Population structure
array[V, J_Q, D_Q] real pop_VJD;

// Case distributions
array[B, V, K, R, D_Q] real pC_BVKRD;
array[B, V, J_Q, R, A] real pC_BVJRA;
array[B, V, K, 4] real pC_BVK4;
array[K, J_Q, R, 2] real pC_KJR2;
array[J_Q] real pJ_VCD; // age dist of cases for month 37-48

// Case counts for different strata
array[B, V, J_Q] real C_BVJ5; // serostatus trial age cases for month 37-48
array[B, J_Q] real C_BJ5;     // serostatus age cases for month 37-48
array[J_Q] real C_J5;         // age cases for month 37-48
array[B, V, K, J_Q, R, 4] real C_BVKJR4;
array[V, K, J_Q, R, D_Q] real C_VKJRD;
array[J_Q, R, D_Q] real C_JRD;
array[R, D_Q] real C_RD;

for(b in 1:B)
 for(k in 1:K)
  for(t in 1:T_Q)
   n_Q[b,k,t] = mu[1,b,k] * (exp(pi_1[1,b] * time_Q[t] + pi_2[1] * ts2[1,b]) + exp(pi_2[1] * time_Q[t] + pi_1[1,b] * ts2[1,b])) / (exp(pi_1[1,b] * ts2[1,b]) + exp(pi_2[1] * ts2[1,b])) ; 

 for(k in 1:K)
  for(t in 1:T_Q){
    n_C[1,k,t] =  n_Q[1,k,t];
    n_C[2,k,t] =  n_Q[2,k,t]; // MO and MU have SP titres 
    n_C[3,k,t] =  n_Q[2,k,t];
  }
  

// INITIAL CONDITIONS
if(include_pK3 == 0) {  // p common to all serotypes
  for(k in 1:K){
    for(j in 1:J_Q) p_KJ[k,j] = p[1,j]; // if include_pK3 == 0, all p are prob exposure to a single serotype 
    hK3[k] = 0;
    q[k] = 0;
    hK1[k] = 0;
    hK2[k] = 0;
    }
 } else { // serotype-specific p
  for(k in 1:K) hK3[k] = -log(1-pK3[1,k]) ; // pK3 prob exposure to serotype k 
  real shK3;
  shK3=sum(hK3);
  for(k in 1:K) {
    q[k] = hK3[k] / shK3 ;
    hK1[k] = q[k] * -log(1-p[1,1]) ; // if include_pK3 == 1, p1 and p2 is cumulative prob exposure 
    hK2[k] = q[k] * -log(1-p[1,2]) ;
    p_KJ[k,1] = 1 - exp(-hK1[k]) ;
    p_KJ[k,2] = 1 - exp(-hK2[k]) ;
    p_KJ[k,3] = pK3[1,k];
  }
}

// this defines the true serostatus populations 
 for(j in 1:J_Q){
    TpSN[j] = (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1
    TpMO[1,j] = p_KJ[1,j]   * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 2
    TpMO[2,j] = p_KJ[2,j]   * (1-p_KJ[1,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 3
    TpMO[3,j] = p_KJ[3,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 4
    TpMO[4,j] = p_KJ[4,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 1,2
    TpMU2[1,j] = p_KJ[1,j] * p_KJ[2,j] * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1,3
    TpMU2[2,j] = p_KJ[1,j] * p_KJ[3,j] * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 1,4
    TpMU2[3,j] = p_KJ[1,j] * p_KJ[4,j] * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 2,3
    TpMU2[4,j] = p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[1,j]) * (1-p_KJ[4,j]);
    // 2,4
    TpMU2[5,j] = p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[3,j]);
    // 3,4 
    TpMU2[6,j] = p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[2,j]);
    // not 1 
    TpMU3[1,j] = p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) ;
    // not 2
    TpMU3[2,j] = p_KJ[1,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[2,j]) ;
    // not 3
    TpMU3[3,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[3,j]) ;
    // not 4 
    TpMU3[4,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[4,j]) ;
    // all 
    TpMU4[j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j];
 }

// this accounts for imperfect test performance 
for(v in 1:V)
 for(j in 1:J_Q)
   for(k in 1:K){
    pSN[1,v,j,1] = spec[1] * TpSN[j] ;
    pSN[2,v,j,1] = (1-spec[1]) * TpSN[j] ;
    pMO[1,v,k,j,1] = (1- sens[1]) * TpMO[k,j] ;
    pMO[2,v,k,j,1] = sens[1] * TpMO[k,j] ;
   }

// Multitypic individuals are all classified as seropositive
for (v in 1:V) {
  for (j in 1:J_Q) {
    // MU2 combinations (6 possibilities)
    for (k in 1:6) {
      pMU2[1, v, k, j, 1] = 0;
      pMU2[2, v, k, j, 1] = TpMU2[k, j];
    }
    // MU3 combinations (K possibilities)
    for (k in 1:K) {
      pMU3[1, v, k, j, 1] = 0;
      pMU3[2, v, k, j, 1] = TpMU3[k, j];
    }
  }
}

// Probabilities of testing seropositive 
 for(j in 1:J_Q)  pSP_Q[j] = (1-spec[1]) * TpSN[j] + sens[1] * sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j];

// VACCINE RISK RATIO 

// outcome and age-group titre offsets 
array[2] real e_beta = {exp(beta[1,1]), exp(beta[1,2])};
array[C,K] real e_alpha; 


for(k in 1:K)
 for(c in 1:C){
   if(alpha_CK == 0){ 
     e_alpha[c,k] = exp(-alpha[1,1]);  // mono outcome offset 
     } else if(alpha_CK == 1) {
       e_alpha[c,k] = exp(-alpha[1,c]);  // serostatus outcome offset 
       } else {
         e_alpha[c,k] = exp(-alpha[1,k]);  // serotype outcome offset 
         }}
        
if(mono_lc_SN == 1){ // SN, oldest, symp 
  for(k in 1:K) nc50[1,k,3,1] = exp(lc[1,1,1]); 
   } else if (mono_lc_SN == 0) {
  for(k in 1:K) nc50[1,k,3,1] = exp(lc[1,1,k]); 
} else if(mono_lc_SN == 2){
  for(k in 1:K)  nc50[1,k,3,1] = exp(lc[1,2,k]) * exp(kappa[1]); 
}

if(mono_lc_MU == 1){ // MU, oldest, symp
 for(k in 1:K)  nc50[3,k,3,1] = exp(lc[1,3,1]); 
   } else if(mono_lc_MU == 0) {
 for(k in 1:K) nc50[3,k,3,1] = exp(lc[1,3,k]); 
} else if (mono_lc_MU == 2){
   for(k in 1:K)  nc50[3,k,3,1] = exp(lc[1,2,k]) * exp(-omega[1]); 
}

// MO, oldest, symp
for(k in 1:K) nc50[2,k,3,1] = exp(lc[1,2,k]); 

// age group offsets 
if(include_beta[1] == 1){
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] =   nc50[c,k,3,1] * e_beta[1]; // youngest
   nc50[c,k,2,1] =   nc50[c,k,3,1] * e_beta[2]; // middle 
 } 
} else if(include_beta[1] == 3) { // age differences 
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] =   nc50[c,k,3,1] * e_beta[1]; // youngest
   nc50[c,k,2,1] =   nc50[c,k,3,1]; // middle 
   }
 } else {
  for(c in 1:C)
   for(k in 1:K){ // no age differences
   nc50[c,k,1,1] = nc50[c,k,3,1];
   nc50[c,k,2,1] = nc50[c,k,3,1];   
  }
}

// hosp 
for(c in 1:C)
 for(k in 1:K) 
  for(j in 1:J_Q) // outcome differences 
    nc50[c,k,j,2] =   nc50[c,k,j,1] * e_alpha[c,k];

if(include_beta[1] == 2){ // age and outcome 
  for(c in 1:C)
   for(k in 1:K){
     nc50[c,k,1,1] =   nc50[c,k,1,1] * e_beta[1]; 
     nc50[c,k,1,2] =   nc50[c,k,1,2] * e_beta[2]; 
   }
}

// Enhancement 
for(k in 1:K){
  if(enhancement == 1){ 
   if(L_K[1] == 0){
   L_C[1,k,1] = 1 + L[1,1]; // symp 
   L_C[1,k,2] = 1 + L[1,1] * ((tau_K==0)?tau[1,1]:tau[1,k]);
  } else if (L_K[1] == 1){
   L_C[1,k,1] = 1 + L[1,k]; // symp 
   L_C[1,k,2] = 1 + L[1,k] * ((tau_K==0)?tau[1,1]:tau[1,k]);
 } 
 } else if (enhancement == 0){ // no SN enhancement 
    for(r in 1:R)  L_C[1,k,r] = 1 ; 
 } 

 for(r in 1:R) { // no enhancement if seropositive 
        L_C[2,k,r] = 1; // MO
        L_C[3,k,r] = 1; // MU
      }
   }
 
// Risk ratios - w can be mono, serostatus or serotype dependent  
for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J_Q)
   for(t in 1:T_Q){
     if(w_CK == 0){
     RR_symp_Q[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[1,1]) ;
     RR_hosp_Q[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[1,1]) ;
     } else if(w_CK == 1){
     RR_symp_Q[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[1,c]) ;
     RR_hosp_Q[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[1,c]) ;  
     } else{
     RR_symp_Q[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[1,k]) ;
     RR_hosp_Q[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[1,k]) ;  
     }
     RR_symp_Q[c,1,k,j,t] =  1 ;
     RR_hosp_Q[c,1,k,j,t] =  1 ;
    }

// FOI  
for(k in 1:K){
  for(t in 1:12)  lambda[k,t] = exp(-lambda_D_Q[k,1]);
  for(t in 13:18) lambda[k,t] = exp(-lambda_D_Q[k,2]);
  for(t in 19:24) lambda[k,t] = exp(-lambda_D_Q[k,3]);
  for(t in 25:36) lambda[k,t] = exp(-lambda_D_Q[k,4]);
  for(t in 37:48) lambda[k,t] = exp(-lambda_D_Q[k,5]);
  for(t in 49:54) lambda[k,t] = exp(-lambda_D_Q[k,6]);
}

// SURVIVAL MODEL - prob of surviving each time point without infection 
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_Q)
   for(t in 1:T_Q) pSN[b,v,j,(t+1)] = lambda[1,t]*lambda[2,t]*lambda[3,t]*lambda[4,t] * pSN[b,v,j,t];

  
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J_Q) {
    for(t in 1:T_Q) pMO[b,v,k,j,(t+1)] = lambda[1,t]*lambda[2,t]*lambda[3,t]*lambda[4,t]/lambda[k,t]* pMO[b,v,k,j,t] ; 
    for(t in (HI+1):T_Q) pMO[b,v,k,j,(t+1)] += (1 - gammaT[k] * RR_symp_Q[1,v,k,j,(t-HI)]) * (1- lambda[k,(t-HI)]) * pSN[b,v,j,(t-HI)] ;
    }  
   
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J_Q) {
      for(t in 1:T_Q){
          // MU_12
          pMU2[b,v,1,j,(t+1)] = lambda[3,t] * lambda[4,t] *  pMU2[b,v,1,j,t] ; 
           // MU_13
          pMU2[b,v,2,j,(t+1)] = lambda[2,t] * lambda[4,t] *  pMU2[b,v,2,j,t] ; 
           // MU_14
          pMU2[b,v,3,j,(t+1)] = lambda[2,t] * lambda[3,t] *  pMU2[b,v,3,j,t] ; 
          // MU_23
          pMU2[b,v,4,j,(t+1)] = lambda[1,t] * lambda[4,t] *  pMU2[b,v,4,j,t] ; 
          // MU_24
          pMU2[b,v,5,j,(t+1)] = lambda[1,t] * lambda[3,t] *  pMU2[b,v,5,j,t] ; 
          // MU_34
          pMU2[b,v,6,j,(t+1)] = lambda[1,t] * lambda[2,t] *  pMU2[b,v,6,j,t] ; 
      }
      for(t in (HI+1):T_Q) {
          // MU_12
          pMU2[b,v,1,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp_Q[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,2,j,(t-HI)] + (1 - gammaT[2] * rhoT[2] * RR_symp_Q[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_13
          pMU2[b,v,2,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp_Q[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp_Q[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_14
          pMU2[b,v,3,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp_Q[2,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp_Q[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_23
          pMU2[b,v,4,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp_Q[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp_Q[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_24
          pMU2[b,v,5,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp_Q[2,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp_Q[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_34
          pMU2[b,v,6,j,(t+1)] +=  (1 - gammaT[3] * rhoT[3] * RR_symp_Q[2,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp_Q[2,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMO[b,v,3,j,(t-HI)] ;
        }}

for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J_Q){
     for(t in 1:T_Q){ 
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] = lambda[1,t]  *  pMU3[b,v,1,j,t] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] = lambda[2,t]  *  pMU3[b,v,2,j,t] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] = lambda[3,t]  *  pMU3[b,v,3,j,t] ; 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] = lambda[4,t]  *  pMU3[b,v,4,j,t] ; 
     }
     for(t in (HI+1):T_Q) {
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp_Q[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp_Q[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp_Q[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp_Q[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp_Q[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp_Q[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp_Q[3,v,4,j,(t-HI)]) * (1-lambda[4,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp_Q[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp_Q[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] ;  
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] += (1 - phiT[3] * gammaT[3] *  RR_symp_Q[3,v,3,j,(t-HI)]) * (1-lambda[3,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp_Q[3,v,2,j,(t-HI)]) * (1-lambda[2,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp_Q[3,v,1,j,(t-HI)]) * (1-lambda[1,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] ;  
      }
    }  

// calculate infection incidence, symptomatic, hospitalised 
for(b in 1:B)
  for(v in 1:V)
    for(j in 1:J_Q) {
      for(t in 1:T_Q){
         Inc3[b,v,1,j,t] = (1 - lambda[1,t]) * (pMU2[b,v,4,j,t] + pMU2[b,v,5,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,2,j,t] = (1 - lambda[2,t]) * (pMU2[b,v,2,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,3,j,t] = (1 - lambda[3,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,5,j,t]);  
         Inc3[b,v,4,j,t] = (1 - lambda[4,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,2,j,t] + pMU2[b,v,4,j,t]);
         
         for(k in 1:K){
           Inc1[b,v,k,j,t] = (1 - lambda[k,t]) * pSN[b,v,j,t]; 
           Inc2[b,v,k,j,t] = (1 - lambda[k,t]) * (sum(pMO[b,v, ,j,t]) - pMO[b,v,k,j,t]);
           Inc4[b,v,k,j,t] = (1 - lambda[k,t]) * pMU3[b,v,k,j,t];  
           D1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * RR_symp_Q[1,v,k,j,t] ;
           D2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k] * RR_symp_Q[2,v,k,j,t];
           D34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * RR_symp_Q[3,v,k,j,t];  
           H1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * delta[1,k] * RR_hosp_Q[1,v,k,j,t];
           if(include_eps == 0){
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * delta[1,k]  * RR_hosp_Q[2,v,k,j,t];
           } else{
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * delta[1,k] * epsilon * RR_hosp_Q[2,v,k,j,t]; 
           }
           H34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * delta[1,k] * RR_hosp_Q[3,v,k,j,t] ;   
         }}}

// Aggregate to match published time points e.g. 1-12 months 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J_Q){
   Di[b,v,k,j,1] = sum(D1[b,v,k,j,1:12])  + sum(D2[b,v,k,j,1:12])   ;
   H[b,v,k,j,1]  = sum(H1[b,v,k,j,1:12])  + sum(H2[b,v,k,j,1:12])   ;
   Di[b,v,k,j,2] = sum(D1[b,v,k,j,13:18]) + sum(D2[b,v,k,j,13:18])  ;
   H[b,v,k,j,2]  = sum(H1[b,v,k,j,13:18]) + sum(H2[b,v,k,j,13:18])  ;
   Di[b,v,k,j,3] = sum(D1[b,v,k,j,19:24]) + sum(D2[b,v,k,j,19:24])  ;
   H[b,v,k,j,3]  = sum(H1[b,v,k,j,19:24]) + sum(H2[b,v,k,j,19:24])  ;
   Di[b,v,k,j,4] = sum(D1[b,v,k,j,25:36]) + sum(D2[b,v,k,j,25:36])  ;
   H[b,v,k,j,4]  = sum(H1[b,v,k,j,25:36]) + sum(H2[b,v,k,j,25:36])  ;
   Di[b,v,k,j,5] = sum(D1[b,v,k,j,37:48]) + sum(D2[b,v,k,j,37:48])  ;
   H[b,v,k,j,5]  = sum(H1[b,v,k,j,37:48]) + sum(H2[b,v,k,j,37:48])  ;
   Di[b,v,k,j,6] = sum(D1[b,v,k,j,49:54]) + sum(D2[b,v,k,j,49:54])  ;
   H[b,v,k,j,6]  = sum(H1[b,v,k,j,49:54]) + sum(H2[b,v,k,j,49:54])  ;
   
   Di[b,v,k,j,1] += sum(D34[b,v,k,j,1:12])  ;
   H[b,v,k,j,1]  += sum(H34[b,v,k,j,1:12])  ;
   Di[b,v,k,j,2] += sum(D34[b,v,k,j,13:18])  ;
   H[b,v,k,j,2]  += sum(H34[b,v,k,j,13:18])  ;
   Di[b,v,k,j,3] += sum(D34[b,v,k,j,19:24])  ;
   H[b,v,k,j,3]  += sum(H34[b,v,k,j,19:24])  ;
   Di[b,v,k,j,4] += sum(D34[b,v,k,j,25:36])  ;
   H[b,v,k,j,4]  += sum(H34[b,v,k,j,25:36])  ;
   Di[b,v,k,j,5] += sum(D34[b,v,k,j,37:48])  ;
   H[b,v,k,j,5]  += sum(H34[b,v,k,j,37:48])  ;
   Di[b,v,k,j,6] += sum(D34[b,v,k,j,49:54])  ;
   H[b,v,k,j,6]  += sum(H34[b,v,k,j,49:54])  ;
   }


for(v in 1:V)
 for(j in 1:J_Q)
  for(d in 1:(D_Q-1)) 
    pop_VJD[v,j,d] = sum(pop_Q[ ,v,j,d]) ; // pop not by serostatus 

// Cases 
for(j in 1:J_Q) {
   for(b in 1:B) {
    for(v in 1:V) {
     for(k in 1:K) {
       for(d in 1:(D_Q-1)) {
          C_BVKJRD_Q[b,v,k,j,1,d] = Di[b,v,k,j,d] * pop_VJD[v,j,d]; // VCD up to month 48
        }
      }
          C_BVJ5[b,v,j] = sum(C_BVKJRD_Q[b,v, ,j,1,5]) ;
    }
    C_BJ5[b,j] = sum(C_BVJ5[b, ,j]);
  }
  C_J5[j] = sum(C_BJ5[ ,j]);
}

real s_C_J5 = sum(C_J5);

for(j in 1:J_Q) {
 pJ_VCD[j] = C_J5[j] / s_C_J5 ;
 
 for(v in 1:V) {
   for(b in 1:B) {
     for(d in 1:(D_Q-1)) {
       pop_BVJD_Q[b,v,j,d] = pop_Q[b,v,j,d]  ; // editable population 
     }
     pop_BVJD_Q[b,v,j,6] = pop_Q[b,v,j,6] - (N_VCD_BV5_Q[b,v] * pJ_VCD[j])  ; // censor population 
   }
   
   pop_VJD[v,j,6] = sum(pop_BVJD_Q[ ,v,j,6]) ; // pop not by serostatus
   
   for(b in 1:B) {
      for(k in 1:K) {
        C_BVKJRD_Q[b,v,k,j,1,6] = Di[b,v,k,j,6] * pop_VJD[v,j,6];  // VCD at 54 months
      }
    }
  }
}

for (j in 1:J_Q) {
  for (k in 1:K) {
    for (v in 1:V) {
      for (b in 1:B) {
        for (d in 1:D_Q) {
          C_BVKJRD_Q[b, v, k, j, 2, d] = H[b, v, k, j, d] * pop_VJD[v, j, d]; // hosp
        }

        for (r in 1:R) {
          C_BVKJR4[b, v, k, j, r, 1] = sum(C_BVKJRD_Q[b, v, k, j, r, 1:3]); // aggregate 1:24 months for hosp

          for (d in 2:4) {
            C_BVKJR4[b, v, k, j, r, d] = C_BVKJRD_Q[b, v, k, j, r, (d + 2)]; 
            // 1 is 1:24, 2 is 36 (4), 3 is 48 (5), 4 is 54 (6)
          }
        }
      }
    }
  
  
  for(r in 1:R)
   for(d in 1:D_Q) {
     for(v in 1:V) C_VKJRD[v,k,j,r,d] = sum(C_BVKJRD_Q[ ,v,k,j,r,d]); 
     C_KJRD_Q[k,j,r,d] = sum(C_VKJRD[ ,k,j,r,d]) ; 
   }
 }
 
 for(r in 1:R)
  for(d in 1:D_Q)
    C_JRD[j,r,d] = sum(C_KJRD_Q[ ,j,r,d]);
}


for(r in 1:R)
 for(d in 1:D_Q)
  C_RD[r,d] = sum(C_JRD[ ,r,d]) ; 


for(b in 1:B)
 for(v in 1:V)
  for(r in 1:R) {
   for(k in 1:K)
    for(d in 1:D_Q)
     pC_BVKRD[b,v,k,r,d] = (sum(C_BVKJRD_Q[b,v,k, ,r,d])) / C_RD[r,d];
   for(j in 1:J_Q)
    for(a in 1:A)
     pC_BVJRA[b,v,j,r,a] = sum(C_BVKJRD_Q[b,v, ,j,r,a]) / C_RD[r,a];
  }

// sum over year 2 for serotype and age cases
 for(k in 1:K)
  for(j in 1:J_Q)
   for(r in 1:R) {
     pC_KJR2[k,j,r,1] = C_KJRD_Q[k,j,r,1] / C_RD[r,1] ;
     pC_KJR2[k,j,r,2] = sum(C_KJRD_Q[k,j,r,2:3]) / sum(C_RD[r,2:3]) ;
  }

// serotype hosp combined for year 1 and 2 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K) {
   pC_BVK4[b,v,k,1] = sum(C_BVKJR4[b,v,k, ,2,1]) / sum(C_RD[2,1:3]) ; // total hosp 1:24 months
   for(d in 2:4)
     pC_BVK4[b,v,k,d] = sum(C_BVKJR4[b,v,k, ,2,d]) / C_RD[2,(d+2)]; // 4, 5, 6
  }

// matrix for likelihood function 
for(j in 1:J_Q)
  for(d in 1:2) {
   pD_KJ2_Q[j,d]   = pC_KJR2[1,j,1,d];
   pD_KJ2_Q[j+3,d] = pC_KJR2[2,j,1,d];
   pD_KJ2_Q[j+6,d] = pC_KJR2[3,j,1,d];
   pD_KJ2_Q[j+9,d] = pC_KJR2[4,j,1,d];
   pH_KJ2_Q[j,d]   = pC_KJR2[1,j,2,d];
   pH_KJ2_Q[j+3,d] = pC_KJR2[2,j,2,d];
   pH_KJ2_Q[j+6,d] = pC_KJR2[3,j,2,d];
   pH_KJ2_Q[j+9,d] = pC_KJR2[4,j,2,d];
 }

for(k in 1:K)
  for(d in 1:D_Q){
   pD_BVKD_Q[k,d]    =  pC_BVKRD[1,1,k,1,d];
   pD_BVKD_Q[k+4,d]  =  pC_BVKRD[1,2,k,1,d]; 
   pD_BVKD_Q[k+8,d]  =  pC_BVKRD[2,1,k,1,d]; 
   pD_BVKD_Q[k+12,d] =  pC_BVKRD[2,2,k,1,d]; 
  }

// serotype hosp published for 4 time intervals
for(k in 1:K)
  for(d in 1:4){
   pH_BVK4_Q[k,d]    =  pC_BVK4[1,1,k,d];
   pH_BVK4_Q[k+4,d]  =  pC_BVK4[1,2,k,d]; 
   pH_BVK4_Q[k+8,d]  =  pC_BVK4[2,1,k,d]; 
   pH_BVK4_Q[k+12,d] =  pC_BVK4[2,2,k,d]; 
 }

for(j in 1:J_Q)
  for(a in 1:A){
   pD_BVJA_Q[j,a]     =  pC_BVJRA[1,1,j,1,a];
   pD_BVJA_Q[j+3,a]   =  pC_BVJRA[1,2,j,1,a];
   pD_BVJA_Q[j+6,a]   =  pC_BVJRA[2,1,j,1,a];
   pD_BVJA_Q[j+9,a]   =  pC_BVJRA[2,2,j,1,a];
   pH_BVJA_Q[j,a]     =  pC_BVJRA[1,1,j,2,a];
   pH_BVJA_Q[j+3,a]   =  pC_BVJRA[1,2,j,2,a];
   pH_BVJA_Q[j+6,a]   =  pC_BVJRA[2,1,j,2,a];
   pH_BVJA_Q[j+9,a]   =  pC_BVJRA[2,2,j,2,a];
 }

for(d in 1:D_Q)
 for(r in 1:R)
  pC_RD_Q[r,d] = C_RD[r,d] / pop_D_Q[d]; 

// log likelihood 
 ll_Q=0;
 
 ll_Q += binomial_lpmf(SP_J_Q   | pop_J_Q, pSP_Q);
 ll_Q += binomial_lpmf(VCD_D_Q  | pop_D_Q, pC_RD_Q[1,]);
 ll_Q += binomial_lpmf(HOSP_D_Q | pop_D_Q, pC_RD_Q[2,]);


   ll_Q +=  multinomial_lpmf(VCD_BVKD_Q[ ,1] | pD_BVKD_Q[ ,1]);
  for(d in 2:D_Q)   ll_Q += multinomial_lpmf(VCD_BVKD_Q[ ,d] | pD_BVKD_Q[ ,d]);


 for(a in 1:A)    ll_Q += multinomial_lpmf(VCD_BVJA_Q[ ,a] | pD_BVJA_Q[ ,a]) ;
 for(d in 1:2)    ll_Q += multinomial_lpmf(VCD_KJ2_Q[ ,d]  | pD_KJ2_Q[ ,d]) ;
 for(d in 1:4)    ll_Q += multinomial_lpmf(HOSP_BVK4_Q[ ,d] | pH_BVK4_Q[ ,d]);
 for(a in 1:A)    ll_Q += multinomial_lpmf(HOSP_BVJA_Q[ ,a] | pH_BVJA_Q[ ,a]) ;
 for(d in 1:2)    ll_Q += multinomial_lpmf(HOSP_KJ2_Q[ ,d]  | pH_KJ2_Q[ ,d]) ;

}

// ######################## DENGVAXIA BLOCK  ###################################
  { // this { defines a block within which variables declared are local (can't have lower or upper)
  
  array[D_De] real lambda_k;
  real lambda_m;
  real pm; 
  
  // contrain pk3 / p by lambda 
  for(d in 1:D_De) lambda_k[d] = lambda_D_De[d]  / 4; // lambda_D is cumulative FOI so divide by 4 to make it mean serotype FOI
  lambda_m =(lambda_k[1] * 13 + lambda_k[2] * 11 + lambda_k[3] * 12 + lambda_k[4] * ( T_De-13 - 11 - 12)) / T_De ; // weighted FOI over trial 
  pm = 1-exp(-12.5*12*lambda_m); // always prob exposure to a single serotype 
 
  if(pm<1e-3) pm=1e-3;
  if(pm>0.99) pm=0.99;
  real vari = 0.02*pm *(1-pm) ; 
  shape1_De = (((1-pm) / vari) - (1/pm)) * (pm^2) ;
  shape2_De = shape1_De * (1/pm - 1);

vector[K] deltaT;

// Age-specific historic hazard and serotype proportions
array[K] real hK2;       // Age 3 historic hazard
array[K] real q;         // Historic serotype proportion
array[K] real hK1;       // Age 1 historic hazard
array[K, J_De] real p_KJ; // Serotype age probability of exposure 

// True infection states
array[J_De] real TpSN;               // True seronegative (SN)
array[K, J_De] real TpMO;            // True monotypic (MO)
array[K + 2, J_De] real TpMU2;       // True multitypic (MU2)
array[K, J_De] real TpMU3;           // True multitypic (MU3)
array[J_De] real TpMU4;              // True multitypic (MU4) 
  
// Probability of infection states over time + baseline
array[B, V, J_De, T_De + 1] real pSN;       // Probability seronegative (SN)
array[B, V, K, J_De, T_De + 1] real pMO;    // Probability monotypic (MO)
array[B, V, K + 2, J_De, T_De + 1] real pMU2; // Probability multitypic (MU2)
array[B, V, K, J_De, T_De + 1] real pMU3;    // Probability multitypic (MU3)
  
// Other model variables
array[C, K, T_De] real n_C;        // Allow for multitypic titres
array[C, K, R] real L_C;           // MO and MU: L = 1
array[C, K, J_De, R] real nc50;    // nc50
array[K, J_De, T_De] real lambda;  // Force of infection (FOI)
array[J_De] real scale_FOI_J;      // Scale FOI by J

// Incidence variables
array[B, V, K, J_De, T_De] real Inc1;
array[B, V, K, J_De, T_De] real Inc2;
array[B, V, K, J_De, T_De] real Inc3;
array[B, V, K, J_De, T_De] real Inc4;

// Disease progression variables
array[B, V, K, J_De, T_De] real D1;
array[B, V, K, J_De, T_De] real D2;
array[B, V, K, J_De, T_De] real D34;
array[B, V, K, J_De, T_De] real H1;
array[B, V, K, J_De, T_De] real H2;
array[B, V, K, J_De, T_De] real H34;
  
// Aggregate disease metrics
array[B, V, K, J_De] real Di;
array[B, V, K, J_De, D_De] real H;
  
// Case-related variables
array[V, K, J_De, D_De] real Ho_VKJD;
array[K, J_De, D_De] real Ho_KJD;
array[K, J_De] real Sy_KJ;
array[J_De, D_De] real Ho_JD;
array[J_De] real Sy_J;
array[D_De] real Ho_D;
real Sy; // Total symptomatic cases

// Proportion of cases
array[B, V, J_De] real pD_BVJ;
array[V, K] real pD_VK;
array[B, V, J_De, D_De] real pH_BVJD;
array[B, V, K, J_De] real pH_BVKJ;

// BIPHASIC TITRES 

for(b in 1:B)
 for(k in 1:K)
  for(t in 1:T_De)
   n_De[b,k,t] = mu[2,b,k] * (exp(pi_1[2,b] * time_De[t] + pi_2[2] * ts2[2,b]) + exp(pi_2[2] * time_De[t] + pi_1[2,b] * ts2[2,b])) / (exp(pi_1[2,b] * ts2[2,b]) + exp(pi_2[2] * ts2[2,b])) ; 

// MO and MU have SP titres 
 for(k in 1:K)
  for(t in 1:T_De){
    n_C[1,k,t] =  n_De[1,k,t];
    n_C[2,k,t] =  n_De[2,k,t]; 
    n_C[3,k,t] =  n_De[2,k,t];
  }
  
      
if(delta_KJ[2]==1) {
  for(k in 1:K) deltaT[k] = delta[2,k];
} else if(delta_KJ[2] == 2) {
    for(j in 1:J_De) deltaT[j] = delta[2,j];
    deltaT[3] = 1;
    deltaT[4] = 1;
} else {
  for(k in 1:K) deltaT[k] = delta[2,1];
}      
      
// INITIAL CONDITIONS
if(include_pK3 == 0) {  // p common to all serotypes
  for(k in 1:K){
    for(j in 1:J_De) p_KJ[k,j] = p[2,j]; // if include_pK3 == 0, all p are prob exposure to a single serotype 
    hK2[k] = 0;
    q[k] = 0;
    hK1[k] = 0;
    }
 } else { // serotype-specific p
  for(k in 1:K) hK2[k] = -log(1-pK3[2,k]) ; // pK3 prob exposure to serotype k 
  real shK2;
  shK2=sum(hK2);
  for(k in 1:K) {
    q[k] = hK2[k] / shK2 ;
    hK1[k] = q[k] * -log(1-p[2,1]) ; // if include_pK3 == 1, p1 is cumulative prob exposure 
    p_KJ[k,1] = 1 - exp(-hK1[k]) ;
    p_KJ[k,2] = pK3[2,k];
  }
}

// this defines the true serostatus populations 
 for(j in 1:J_De){
    TpSN[j] = (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1
    TpMO[1,j] = p_KJ[1,j]   * (1-p_KJ[2,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 2
    TpMO[2,j] = p_KJ[2,j]   * (1-p_KJ[1,j]) * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 3
    TpMO[3,j] = p_KJ[3,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 4
    TpMO[4,j] = p_KJ[4,j]   * (1-p_KJ[1,j]) * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 1,2
    TpMU2[1,j] = p_KJ[1,j] * p_KJ[2,j] * (1-p_KJ[3,j]) * (1-p_KJ[4,j]);
    // 1,3
    TpMU2[2,j] = p_KJ[1,j] * p_KJ[3,j] * (1-p_KJ[2,j]) * (1-p_KJ[4,j]);
    // 1,4
    TpMU2[3,j] = p_KJ[1,j] * p_KJ[4,j] * (1-p_KJ[2,j]) * (1-p_KJ[3,j]);
    // 2,3
    TpMU2[4,j] = p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[1,j]) * (1-p_KJ[4,j]);
    // 2,4
    TpMU2[5,j] = p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[3,j]);
    // 3,4 
    TpMU2[6,j] = p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) * (1-p_KJ[2,j]);
    // not 1 
    TpMU3[1,j] = p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[1,j]) ;
    // not 2
    TpMU3[2,j] = p_KJ[1,j] * p_KJ[3,j] * p_KJ[4,j] * (1-p_KJ[2,j]) ;
    // not 3
    TpMU3[3,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[4,j] * (1-p_KJ[3,j]) ;
    // not 4 
    TpMU3[4,j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * (1-p_KJ[4,j]) ;
    // all 
    TpMU4[j] = p_KJ[1,j] * p_KJ[2,j] * p_KJ[3,j] * p_KJ[4,j];
 }

// this accounts for imperfect test performance 
for(v in 1:V)
 for(j in 1:J_De)
   for(k in 1:K){
    pSN[1,v,j,1] = spec[2] * TpSN[j] ;
    pSN[2,v,j,1] = (1-spec[2]) * TpSN[j] ;
    pMO[1,v,k,j,1] = (1- sens[2]) * TpMO[k,j] ;
    pMO[2,v,k,j,1] = sens[2] * TpMO[k,j] ;
   }

// Assume all multitypics test seropositive 
  for(v in 1:V)
   for(j in 1:J_De) {
    for(k in 1:6){  // there are 6 MU_2 combinations  
        pMU2[1,v,k,j,1] = 0;
        pMU2[2,v,k,j,1] = TpMU2[k,j] ; 
    }
    for(k in 1:K){  // there are 4 MU_3 combinations 
        pMU3[1,v,k,j,1] = 0; 
        pMU3[2,v,k,j,1] = TpMU3[k,j] ;
    } }
   

// probabilities of testing seropositive 
for(j in 1:J_De)  pSP_De[j] = (1-spec[2]) * TpSN[j] + sens[2] * sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j];


// VACCINE RISK RATIO 

// outcome and age-group titre offsets 
array [K] real e_beta = {exp(beta[2,1]), exp(beta[2,2])};
array [C,K] real e_alpha; 

for(k in 1:K)
 for(c in 1:C){
   if(alpha_CK == 0){ 
     e_alpha[c,k] = exp(-alpha[2,1]);  // mono outcome offset 
     } else if(alpha_CK == 1) {
       e_alpha[c,k] = exp(-alpha[2,c]);  // serostatus outcome offset 
       } else {
         e_alpha[c,k] = exp(-alpha[2,k]);  // serotype outcome offset 
         }}
        
if(mono_lc_SN == 1){ // SN, oldest, symp 
  for(k in 1:K) nc50[1,k,2,1] = exp(lc[2,1,1]); // single 
   } else if (mono_lc_SN == 0) {
  for(k in 1:K) nc50[1,k,2,1] = exp(lc[2,1,k]); // serotype 
} else if(mono_lc_SN == 2){
  for(k in 1:K)  nc50[1,k,2,1] = exp(lc[2,2,k]) * exp(kappa[2]); // offset MO 
}

if(mono_lc_MU == 1){ // MU, oldest, symp
 for(k in 1:K)  nc50[3,k,2,1] = exp(lc[2,3,1]); // single 
   } else if(mono_lc_MU == 0) {
 for(k in 1:K) nc50[3,k,2,1] = exp(lc[2,3,k]); // serotype 
} else if (mono_lc_MU == 2){
   for(k in 1:K)  nc50[3,k,2,1] = exp(lc[2,2,k]) * exp(-omega[2]); // offset MO 
}

// MO, oldest, symp
for(k in 1:K)  nc50[2,k,2,1] = exp(lc[2,2,k]); // MO, serotype, oldest, symp

// age group offsets 
if(include_beta[2] == 1){ // youngest same for both outcomes 
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] =   nc50[c,k,2,1] * e_beta[1]; // youngest
 } 
} else { // no beta 
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1,1] = nc50[c,k,2,1];
  }
}

// hosp 
for(c in 1:C)
 for(k in 1:K) 
  for(j in 1:J_De)
    nc50[c,k,j,2] = nc50[c,k,j,1] * e_alpha[c,k];

if(include_beta[2] == 2){ // youngest by outcome 
  for(c in 1:C)
   for(k in 1:K){
     nc50[c,k,1,1] =   nc50[c,k,1,1] * e_beta[1]; 
     nc50[c,k,1,2] =   nc50[c,k,1,2] * e_beta[2]; 
   }
}

// Enhancement 
for(k in 1:K){
  if(enhancement == 1){ 
   if(L_K[2] == 0){
   L_C[1,k,1] = 1 + L[2,1]; // symp 
   L_C[1,k,2] = 1 + L[2,1] * ((tau_K==0)?tau[2,1]:tau[2,k]);
  } else if (L_K[2] == 1){
   L_C[1,k,1] = 1 + L[2,k]; // symp 
   L_C[1,k,2] = 1 + L[2,k] * ((tau_K==0)?tau[2,1]:tau[2,k]);
 } 
 } else if (enhancement == 0){ // no SN enhancement 
    for(r in 1:R)  L_C[1,k,r] = 1 ; 
 } 

 for(r in 1:R) { // no enhancement if seropositive 
        L_C[2,k,r] = 1; // MO
        L_C[3,k,r] = 1; // MU
      }
   }
 
// Risk ratios - w can be mono, serostatus or serotype dependent  
for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J_De)
   for(t in 1:T_De){
     if(w_CK == 0){
     RR_symp_De[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[2,1]) ;
     RR_hosp_De[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[2,1]) ;
     } else if(w_CK == 1){
     RR_symp_De[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[2,c]) ;
     RR_hosp_De[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[2,c]) ;  
     } else{
     RR_symp_De[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[2,k]) ;
     RR_hosp_De[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[2,k]) ;  
     }
     RR_symp_De[c,1,k,j,t] =  1 ;
     RR_hosp_De[c,1,k,j,t] =  1 ;
    }

// FOI  
for(k in 1:K)
 for(j in 1:J_De){
  for(t in 1:13)    lambda[k,j,t] = exp(-lambda_D_De[1] *  theta[k]) ; // 13-25
  for(t in 14:24)   lambda[k,j,t] = exp(-lambda_D_De[2] *  theta[k]) ; // 26-36
  for(t in 25:36)   lambda[k,j,t] = exp(-lambda_D_De[3] *  theta[k]) ; // 37-48
  for(t in 37:T_De) lambda[k,j,t] = exp(-lambda_D_De[4] *  theta[k]) ; // 49- endpoint 
}


// SURVIVAL MODEL - prob of surviving each time point without infection 
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_De)
   for(t in 1:T_De) pSN[b,v,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t]*lambda[3,j,t]*lambda[4,j,t] * pSN[b,v,j,t];
   
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J_De) {
    for(t in 1:T_De) pMO[b,v,k,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t]*lambda[3,j,t]*lambda[4,j,t]/lambda[k,j,t] * pMO[b,v,k,j,t] ;
    for(t in (HI+1):T_De) pMO[b,v,k,j,(t+1)] += (1 - gammaT[k] * RR_symp_De[1,v,k,j,(t-HI)]) * (1- lambda[k,j,(t-HI)]) * pSN[b,v,j,(t-HI)] ;
    }  
   
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J_De) {
      for(t in 1:T_De){
          // MU_12
          pMU2[b,v,1,j,(t+1)] = lambda[3,j,t] * lambda[4,j,t] *  pMU2[b,v,1,j,t] ; 
           // MU_13
          pMU2[b,v,2,j,(t+1)] = lambda[2,j,t] * lambda[4,j,t] *  pMU2[b,v,2,j,t] ; 
           // MU_14
          pMU2[b,v,3,j,(t+1)] = lambda[2,j,t] * lambda[3,j,t] *  pMU2[b,v,3,j,t] ; 
          // MU_23
          pMU2[b,v,4,j,(t+1)] = lambda[1,j,t] * lambda[4,j,t] *  pMU2[b,v,4,j,t] ; 
          // MU_24
          pMU2[b,v,5,j,(t+1)] = lambda[1,j,t] * lambda[3,j,t] *  pMU2[b,v,5,j,t] ; 
          // MU_34
          pMU2[b,v,6,j,(t+1)] = lambda[1,j,t] * lambda[2,j,t] *  pMU2[b,v,6,j,t] ; 
      }
      for(t in (HI+1):T_De) {
          // MU_12
          pMU2[b,v,1,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp_De[2,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] + (1 - gammaT[2] * rhoT[2] * RR_symp_De[2,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_13
          pMU2[b,v,2,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp_De[2,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp_De[2,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_14
          pMU2[b,v,3,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp_De[2,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp_De[2,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_23
          pMU2[b,v,4,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp_De[2,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp_De[2,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_24
          pMU2[b,v,5,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp_De[2,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp_De[2,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_34
          pMU2[b,v,6,j,(t+1)] +=  (1 - gammaT[3] * rhoT[3] * RR_symp_De[2,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp_De[2,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] ;
        }}

for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J_De){
     for(t in 1:T_De){ 
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] = lambda[1,j,t]  *  pMU3[b,v,1,j,t] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] = lambda[2,j,t]  *  pMU3[b,v,2,j,t] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] = lambda[3,j,t]  *  pMU3[b,v,3,j,t] ; 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] = lambda[4,j,t]  *  pMU3[b,v,4,j,t] ; 
     }
       for(t in (HI+1):T_De) {
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp_De[3,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp_De[3,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp_De[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp_De[3,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp_De[3,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp_De[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp_De[3,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp_De[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp_De[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] ;  
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] += (1 - phiT[3] * gammaT[3] *  RR_symp_De[3,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp_De[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp_De[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] ;  
      }
    }  

// calculate infection incidence, symptomatic, hospitalised 
for(b in 1:B)
  for(v in 1:V)
    for(j in 1:J_De) {
      for(t in 1:T_De){
         Inc3[b,v,1,j,t] = (1 - lambda[1,j,t]) * (pMU2[b,v,4,j,t] + pMU2[b,v,5,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,2,j,t] = (1 - lambda[2,j,t]) * (pMU2[b,v,2,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,3,j,t] = (1 - lambda[3,j,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,5,j,t]);  
         Inc3[b,v,4,j,t] = (1 - lambda[4,j,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,2,j,t] + pMU2[b,v,4,j,t]);
         
         for(k in 1:K){
           Inc1[b,v,k,j,t] = (1 - lambda[k,j,t]) * pSN[b,v,j,t]; 
           Inc2[b,v,k,j,t] = (1 - lambda[k,j,t]) * (sum(pMO[b,v, ,j,t]) - pMO[b,v,k,j,t]);
           Inc4[b,v,k,j,t] = (1 - lambda[k,j,t]) * pMU3[b,v,k,j,t];  
           D1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * RR_symp_De[1,v,k,j,t] ;
           D2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k] * RR_symp_De[2,v,k,j,t];
           D34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * RR_symp_De[3,v,k,j,t];  
  
           if(delta_KJ[2] != 2){
           H1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * deltaT[k] * RR_hosp_De[1,v,k,j,t];
           if(include_eps == 0){
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[k]  * RR_hosp_De[2,v,k,j,t];
           } else{
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[k] * epsilon * RR_hosp_De[2,v,k,j,t]; 
           }
           H34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * deltaT[k] * RR_hosp_De[3,v,k,j,t] ;   
         } else {
           H1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * deltaT[j] * RR_hosp_De[1,v,k,j,t];
           if(include_eps == 0){
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[j]  * RR_hosp_De[2,v,k,j,t];
           } else{
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[j] * epsilon * RR_hosp_De[2,v,k,j,t]; 
           }
           H34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * deltaT[j] * RR_hosp_De[3,v,k,j,t] ;   
         }
         
         }}}

// Aggregate to match published time points e.g. 1-12 months 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J_De){ 
   Di[b,v,k,j] = sum(D1[b,v,k,j,1:13])   + sum(D2[b,v,k,j,1:13]) ; 
   Di[b,v,k,j]   +=  sum(D34[b,v,k,j,1:13]) ;
   
   H[b,v,k,j,1]  = sum(H1[b,v,k,j,1:13])  + sum(H2[b,v,k,j,1:13]) ;
   H[b,v,k,j,2]  = sum(H1[b,v,k,j,14:24]) + sum(H2[b,v,k,j,14:24]) ;
   H[b,v,k,j,3]  = sum(H1[b,v,k,j,25:36]) + sum(H2[b,v,k,j,25:36]) ;
   H[b,v,k,j,4]  = sum(H1[b,v,k,j,37:T_De])  + sum(H2[b,v,k,j,37:T_De]) ;

   H[b,v,k,j,1]  += sum(H34[b,v,k,j,1:13])  ;
   H[b,v,k,j,2]  += sum(H34[b,v,k,j,14:24]) ;
   H[b,v,k,j,3]  += sum(H34[b,v,k,j,25:36]) ;
   H[b,v,k,j,4]  += sum(H34[b,v,k,j,37:T_De])  ;
   }



// Cases 
for(b in 1:B) 
 for(v in 1:V) 
  for(k in 1:K)
   for(j in 1:J_De) 
    for(d in 1:D_De){
          Sy_BVKJ_De[b,v,k,j] = Di[b,v,k,j] * pop_VCD_VJ_De[v,j]; // year 1 symp 
          Ho_BVKJD_De[b,v,k,j,d] = H[b,v,k,j,d] * pop_HOSP_VJD_De[v,j,d]; // hosp across trial 
            }

// Sum over hosp 
for(v in 1:V) 
 for(k in 1:K) 
  for(j in 1:J_De) 
   for(d in 1:D_De)
   Ho_VKJD[v,k,j,d] = sum(Ho_BVKJD_De[ ,v,k,j,d]); 
      
for(k in 1:K) 
 for(j in 1:J_De) 
  for(d in 1:D_De)
  Ho_KJD[k,j,d] = sum(Ho_VKJD[ ,k,j,d]); 

for(j in 1:J_De) 
 for(d in 1:D_De)
 Ho_JD[j,d] = sum(Ho_KJD[ ,j,d]);    

for(d in 1:D_De)
Ho_D[d] = sum(Ho_JD[ ,d]); 


// sum over symp

for(v in 1:V) 
 for(k in 1:K) 
  for(j in 1:J_De) 
  Sy_VKJ_De[v,k,j] = sum(Sy_BVKJ_De[ ,v,k,j]) ; 
    
for(k in 1:K) 
 for(j in 1:J_De) 
 Sy_KJ[k,j] = sum(Sy_VKJ_De[ ,k,j]) ;  

for(j in 1:J_De) 
 Sy_J[j] = sum(Sy_KJ[ ,j]) ;  
 Sy = sum(Sy_J) ; 

for(b in 1:B)
 for(v in 1:V)
   for(k in 1:K)
    for(j in 1:J_De)
     for(d in 1:D_De){
  pD_VK[v,k] = (sum(Sy_VKJ_De[v,k, ])) / Sy ;  // serotype symp y1 
  pD_BVJ[b,v,j] = (sum(Sy_BVKJ_De[b,v, ,j])) / Sy ; // symp year 1 
  pH_BVJD[b,v,j,d] = (sum(Ho_BVKJD_De[b,v, ,j,d])) / Ho_D[d]; // hosp not k
  pH_BVKJ[b,v,k,j] = (sum(Ho_BVKJD_De[b,v,k,j, ])) / (sum(Ho_D)); // hosp not d
  }


// matrix for likelihood function 

for(k in 1:K){
  mD_VK_De[k] =   pD_VK[1,k] ; 
  mD_VK_De[k+4] = pD_VK[2,k] ; 
}

for(j in 1:J_De){
   mD_BVJ_De[j]     =  pD_BVJ[1,1,j]; // SN C 
   mD_BVJ_De[j+2]   =  pD_BVJ[1,2,j]; // SN V
   mD_BVJ_De[j+4]   =  pD_BVJ[2,1,j]; // SP C
   mD_BVJ_De[j+6]   =  pD_BVJ[2,2,j]; // SP V 
  }
  
for(j in 1:J_De)
  for(d in 1:D_De){   
   mH_BVJD_De[j,d]     =  pH_BVJD[1,1,j,d]; // SN C 
   mH_BVJD_De[j+2,d]   =  pH_BVJD[1,2,j,d]; // SN V
   mH_BVJD_De[j+4,d]   =  pH_BVJD[2,1,j,d]; // SP C
   mH_BVJD_De[j+6,d]   =  pH_BVJD[2,2,j,d]; // SP V
    }
 
 for(k in 1:K){
   mH_BVKJ_De[k] =       pH_BVKJ[1,1,k,1]; // SN C < 9
   mH_BVKJ_De[k + 4] =   pH_BVKJ[1,1,k,2]; // SN C > 9
   mH_BVKJ_De[k + 8] =   pH_BVKJ[1,2,k,1]; // SN V < 9
   mH_BVKJ_De[k + 12] =  pH_BVKJ[1,2,k,2]; // SN V > 9
   mH_BVKJ_De[k + 16] =  pH_BVKJ[2,1,k,1]; // SP C < 9
   mH_BVKJ_De[k + 20] =  pH_BVKJ[2,1,k,2]; // SP C > 9
   mH_BVKJ_De[k + 24] =  pH_BVKJ[2,2,k,1]; // SP V < 9
   mH_BVKJ_De[k + 28] =  pH_BVKJ[2,2,k,2]; // SP V > 9 
 }

 pC_De = Sy / pop_VCD_De; 
for(d in 1:D_De)
  pH_De[d] = Ho_D[d] / pop_HOSP_D_De[d]; 
  
// log likelihood 
 ll_De=0;
 

    ll_De += binomial_lpmf(HOSP_D_De| pop_HOSP_D_De, pH_De);
 for(d in 1:D_De)    ll_De += multinomial_lpmf(HOSP_BVJD_De[ ,d]  | mH_BVJD_De[ ,d]) ;
                  ll_De += multinomial_lpmf(HOSP_BVKJ_De | mH_BVKJ_De);

 
 ll_De += binomial_lpmf(SP_J_De  | pop_J_De, pSP_De);
 ll_De += binomial_lpmf(VCD_De   | pop_VCD_De, pC_De);
 ll_De += multinomial_lpmf(VCD_BVJ_De | mD_BVJ_De) ;
 ll_De += multinomial_lpmf(VCD_VK_De  | mD_VK_De) ;
}

}

model {
  target += ll_Q;
  target += ll_De;
  target += ll_Bu;
   
 
// SHARED PRIORS 
    if(uniform == 0){ // use science posterior 
  gamma ~ normal(0.85,0.2);  // prob secondary
  rho ~ normal(1.97,0.40);   // RR secondary
  phi ~ normal(0.25,0.05);
  } else {
  gamma ~ uniform(0,1);  // prob secondary
  rho ~ normal(2,2);   // RR secondary
  phi ~ uniform(0,1);
  }
  
    epsilon ~ normal(1,1);

  for(t in 1:TR){
    delta[t, ] ~ normal(0.25,0.10);
    L[t, ] ~ normal(L_mean,L_sd);
    tau[t, ] ~ normal(1,1);
    w[t, ] ~ normal(1,2);
    omega[t] ~ normal(0,2);
    alpha[t, ] ~ normal(0,2);
    kappa[t] ~  normal(0,2);
    beta[t, ] ~ normal(0,2);
    FOI_J1[t] ~ normal(1,1);
  } 
    

  
for(j in 1:2){ // age groups 1 and 2 
  p[1,j] ~ beta(3,5);
  p[2,j] ~ beta(3,5);
  p[3,j] ~ beta(3,5);
  }
  
  if(include_pK3==1) {
  pK3[1] ~ beta(shape1_Q,shape2_Q) ; // prob exposure to serotype k 
  pK3[2] ~ beta(shape1_De,shape2_De) ; // prob exposure to serotype k 
  pK3[3] ~ beta(shape1_Bu,shape2_Bu) ; // prob exposure to serotype k 
  p[1,3] ~ beta(3,5) ;
  p[2,3] ~ beta(3,5) ;
  p[3,3] ~ beta(3,5) ;
  } else {
  pK3[1] ~ beta(3,5) ;
  pK3[2] ~ beta(3,5) ;
  pK3[3] ~ beta(3,5) ;

  p[1,3] ~ beta(shape1_Q,shape2_Q) ; // prob exposure to a serotype   
  p[2,3] ~ beta(shape1_De,shape2_De) ; // prob exposure to a serotype  
  p[3,3] ~ beta(shape1_Bu,shape2_Bu) ; // prob exposure to a serotype  
  }
  
// Qdenga priors   
  for(k in 1:K) lambda_D_Q[k, ] ~ lognormal(-7,2);
  lc[1,1, ] ~ normal(4.5,1);
  for(c in 2:C)  lc[1,c, ] ~ normal(6.5,1);
  sens[1] ~ normal(0.9,0.05) ;
  spec[1] ~ normal(0.995, 0.01) ;
  hs[1,1] ~ normal(1.65,0.5); // SN
  hs[1,2] ~ normal(4.20,0.5);  // SP
  hl[1] ~ normal(84,12); // fits plus https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7557381/
  ts[1,1] ~ normal(-2.21,0.5) ;
  ts[1,2] ~ normal(0.15,0.5) ;
  
// Dengvaxia priors 
  lambda_D_De ~ lognormal(-7,2);
  theta ~ dirichlet(rep_vector(1.0, K));
  lc[2,1, ] ~ normal(3.85,1);
  for(c in 2:C)  lc[2,c, ] ~ normal(5.98,1);
  sens[2] ~ normal(0.94,0.05) ;
  spec[2] ~ normal(0.74,0.05) ;
  hs[2,1] ~ normal(2.93,0.5); // SN
  hs[2,2] ~ normal(20.0,0.5);  // SP
  hl[2] ~ normal(72.14,12); // fits plus https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7557381/
  ts[2,1] ~ normal(1.34,0.5) ;
  ts[2,2] ~ normal(1.34,0.5) ;
  
// Butantan-DV
  lc[3,1,] ~ normal(5.18,1);
  for(c in 2:C)  lc[3,c,] ~ normal(6.27,1);
  lambda_K_Bu ~ lognormal(-7,2);
  sens[3] ~ normal(0.9,0.05) ;
  spec[3] ~ normal(0.995, 0.01) ;
  hs[3,1] ~ normal(1.65,0.5); // SN
  hs[3,2] ~ normal(4.20,0.5);  // SP
  hl[3] ~ normal(84,12); // fits plus https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7557381/
  ts[3,1] ~ normal(-2.21,0.5) ;
  ts[3,2] ~ normal(0.15,0.5) ;

}

generated quantities{
 
// QDENGA  
// Attack Rates (AR)
array[B, V, K, R, D_Q] real<lower=0> C_BVKRD_Q;    // Case count by B, V, K, R, D_Q
array[B, V, J_Q, R, D_Q] real<lower=0> C_BVJRD_Q;  // Case count by B, V, J_Q, R, D_Q
array[B, V, R, D_Q] real<lower=0> C_BVRD_Q;        // Case count by B, V, R, D_Q
array[B, V, J_Q, R, D_Q] real<lower=0> AR_BVJRD_Q; // Attack rate by B, V, J_Q, R, D_Q
array[B, V, K, R, D_Q] real<lower=0> AR_BVKRD_Q;   // Attack rate by B, V, K, R, D_Q
array[B, V, K, 4] real<lower=0> AR_BVKHD_Q;        // Attack rate by B, V, K, 4
array[B, V, K, R] real<lower=0> AR_BVKR_Q;         // Attack rate by B, V, K, R
array[B, V, J_Q, R] real<lower=0> AR_BVJR_Q;       // Attack rate by B, V, J_Q, R
array[K, J_Q, R, 2] real<lower=0> AR_KJRD_Q;       // Attack rate by K, J_Q, R, 2
array[B, R, D_Q] real<lower=0> AR_BRD_Q;           // Attack rate by B, R, D_Q
array[V, R, D_Q] real<lower=0> AR_VRD_Q;           // Attack rate by V, R, D_Q

// Vaccine Efficacy (VE)
array[C, K, J_Q, R, T_Q] real<upper=1> VE_Q;       // Vaccine efficacy by C, K, J_Q, R, T_Q
array[C, K, R, T_Q] real<upper=1> VE_BKRT_Q;       // Vaccine efficacy by C, K, R, T_Q
array[C, J_Q, R, T_Q] real<upper=1> VE_BJRT_Q;     // Vaccine efficacy by C, J_Q, R, T_Q

// DENGVAXIA  

// Attack Rates (AR)
array[B, V, J_De, D_De] real<lower=0> H_AR_BVJD_De;  // Hospitalized AR by B, V, J_De, D_De
array[B, V, K, J_De] real<lower=0> H_AR_BVKJ_De;     // Hospitalized AR by B, V, K, J_De
array[B, V, J_De] real<lower=0> V_AR_BVJ_De;         // Vaccine-associated AR by B, V, J_De
array[V, K] real<lower=0> V_AR_VK_De;               // Vaccine-associated AR by V, K

// Vaccine Efficacy (VE)
array[C, K, J_De, R, T_De] real<upper=1> VE_De;     // Vaccine efficacy by C, K, J_De, R, T_De
array[C, K, R, T_De] real<upper=1> VE_BKRT_De;      // Vaccine efficacy by C, K, R, T_De
array[C, J_De, R, T_De] real<upper=1> VE_BJRT_De;   // Vaccine efficacy by C, J_De, R, T_De

// BUTANTAN-DV

array[B,V,2] real<lower = 0> AR_BVK_Bu;
array[B,V,J_Bu] real<lower = 0> AR_BVJ_Bu;

// VE
array[C,2,J_Bu,T_Bu] real<upper = 1> VE_Bu;

// JOINT LOG LIKELIHOOD
vector[J_Q+ 2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A+K*J_Q*4+B*V*K*4 + 2*D_De + 3 + 3 + 6] log_lik; // plus 6 for Butantan-DV

 {
// LL
vector[J_Q] ll1; // pSP
vector[D_Q] ll2; // C 
vector[D_Q] ll3; // H
vector[B*V*K*D_Q] ll4; // VCD_BVKD
vector[B*V*J_Q*A] ll5; // VCD_BVJA
vector[B*V*J_Q*A] ll6; // HOSP_BVJA
vector[K*J_Q*2] ll7;   // VCD_KJ2
vector[K*J_Q*2] ll8;   // HOSP_KJ2
vector[B*V*K*4] ll9;   // HOSP_BVK4
vector[J_Q] p_tmp1;    
vector[D_Q] p_tmp2;

p_tmp1 = to_vector(pSP_Q);
ll1= v_SP_J_Q .* log(p_tmp1) + (v_pop_J_Q-v_SP_J_Q) .* log(1-p_tmp1);
p_tmp2 = to_vector(pC_RD_Q[1,]);
ll2= v_VCD_D_Q .* log(p_tmp2) + (v_pop_D_Q-v_VCD_D_Q) .* log(1-p_tmp2);
p_tmp2 = to_vector(pC_RD_Q[2,]);
ll3= v_HOSP_D_Q .* log(p_tmp2) + (v_pop_D_Q-v_HOSP_D_Q) .* log(1-p_tmp2);

ll4 = to_vector(v_VCD_BVKD_Q .* log(pD_BVKD_Q));
ll5 = to_vector(v_VCD_BVJA_Q .* log(pD_BVJA_Q));
ll6 = to_vector(v_HOSP_BVJA_Q .* log(pH_BVJA_Q));
ll7 = to_vector(v_VCD_KJ2_Q .* log(pD_KJ2_Q));
ll8 = to_vector(v_HOSP_KJ2_Q .* log(pH_KJ2_Q));
ll9 = to_vector(v_HOSP_BVK4_Q .* log(pH_BVK4_Q));

log_lik[1:J_Q] = ll1;
log_lik[(J_Q+1):(J_Q+D_Q)] = ll2;
log_lik[(J_Q+D_Q+1):(J_Q+2*D_Q)] = ll3;
log_lik[(J_Q+2*D_Q+1):(J_Q+2*D_Q+B*V*K*D_Q)] = ll4;
log_lik[(J_Q+2*D_Q+B*V*K*D_Q+1):(J_Q+2*D_Q+B*V*K*D_Q+B*V*J_Q*A)] = ll5;
log_lik[(J_Q+2*D_Q+B*V*K*D_Q+B*V*J_Q*A+1):(J_Q+2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A)] = ll6;
log_lik[(J_Q+2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A+1):(J_Q+2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A+K*J_Q*2)] = ll7;
log_lik[(J_Q+2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A+K*J_Q*2+1):(J_Q+2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A+K*J_Q*4)] = ll8;
log_lik[(J_Q+2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A+K*J_Q*4+1):(J_Q+2*D_Q+B*V*K*D_Q+2*B*V*J_Q*A+K*J_Q*4+B*V*K*4)] = ll9;


for(d in 1:D_De) log_lik[319 + d] = binomial_lpmf(HOSP_D_De[d] | pop_HOSP_D_De[d] , pH_De[d]);
for(d in 1:D_De) log_lik[d + D_De + 319]  = multinomial_lpmf(HOSP_BVJD_De[ ,d]  | mH_BVJD_De[ ,d]);
log_lik[2*D_De + 1 + 319] = multinomial_lpmf(HOSP_BVKJ_De | mH_BVKJ_De);
for(j in 1:J_De) log_lik[j + 2*D_De + 1 + 319] = binomial_lpmf(SP_J_De[j] | pop_J_De[j], pSP_De[j]);
log_lik[2*D_De + 3 + 1 + 319] = binomial_lpmf(VCD_De | pop_VCD_De, pC_De);
log_lik[2*D_De + 3 + 2 + 319] = multinomial_lpmf(VCD_BVJ_De | mD_BVJ_De);
log_lik[2*D_De + 3 + 3 + 319] = multinomial_lpmf(VCD_VK_De | mD_VK_De);
}


for(j in 1:J_Bu) log_lik[2*D_De + 3 + 3 + 319 + j] = binomial_lpmf(SP_J_Bu[j] | pop_J_Bu[j], pSP_Bu[j]);
log_lik[2*D_De + 3 + 3 + 319 + 4] = binomial_lpmf(VCD_Bu  | pop_Bu, pC_Bu);
log_lik[2*D_De + 3 + 3 + 319 + 5] = multinomial_lpmf(VCD_BVK_Bu  | mD_BVK_Bu) ;
log_lik[2*D_De + 3 + 3 + 319 + 6] = multinomial_lpmf(VCD_BVJ_Bu  | mD_BVJ_Bu) ;

// ################################# QDENGA  ################################### 
// Attack rates
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(r in 1:R)
    for(d in 1:D_Q)
    C_BVKRD_Q[b,v,k,r,d] = sum(C_BVKJRD_Q[b,v,k, ,r,d]) ;

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_Q)
   for(r in 1:R)
    for(d in 1:D_Q)
    C_BVJRD_Q[b,v,j,r,d] = sum(C_BVKJRD_Q[b,v, ,j,r,d]);

for(b in 1:B)
 for(v in 1:V)
  for(r in 1:R)
   for(d in 1:D_Q)
    C_BVRD_Q[b,v,r,d] = sum(C_BVKRD_Q[b,v, ,r,d]) ;

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_Q)
   for(r in 1:R)
    for(d in 1:D_Q)
     AR_BVJRD_Q[b,v,j,r,d] = C_BVJRD_Q[b,v,j,r,d] / pop_BVJD_Q[b,v,j,d];

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(r in 1:R)
    for(d in 1:D_Q)
     AR_BVKRD_Q[b,v,k,r,d] = C_BVKRD_Q[b,v,k,r,d] / pop_BVD_Q[b,v,d];

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(r in 1:R)
   AR_BVKR_Q[b,v,k,r] = sum(C_BVKRD_Q[b,v,k,r, ]) / pop_BV_Q[b,v] ;

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_Q)
   for(r in 1:R)
    AR_BVJR_Q[b,v,j,r] = sum(C_BVJRD_Q[b,v,j,r,1:4]) / mean(pop_BVJD_Q[b,v,j,1:4]) ; // age data up to 36 months

for(v in 1:V)
 for(r in 1:R)
  for(d in 1:D_Q)
   AR_VRD_Q[v,r,d] = sum(C_BVRD_Q[ ,v,r,d]) / pop_VD_Q[v,d];

for(b in 1:B)
 for(r in 1:R)
  for(d in 1:D_Q)
   AR_BRD_Q[b,r,d] = sum(C_BVRD_Q[b, ,r,d]) / pop_BD_Q[b,d];
  
  for(k in 1:K)
   for(j in 1:J_Q)
    for(r in 1:R){
   AR_KJRD_Q[k,j,r,1] = C_KJRD_Q[k,j,r,1] / pop_JD_Q[j,1]; 
   AR_KJRD_Q[k,j,r,2] = sum(C_KJRD_Q[k,j,r,2:3]) / pop_JD_Q[j,2]; 
}

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K){
    AR_BVKHD_Q[b,v,k,1] = sum(C_BVKRD_Q[b,v,k,2,1:3]) / pop_BVD_Q[b,v,2];
    AR_BVKHD_Q[b,v,k,2] = C_BVKRD_Q[b,v,k,2,4] / pop_BVD_Q[b,v,4];
    AR_BVKHD_Q[b,v,k,3] = C_BVKRD_Q[b,v,k,2,5] / pop_BVD_Q[b,v,5];
    AR_BVKHD_Q[b,v,k,4] = C_BVKRD_Q[b,v,k,2,5] / pop_BVD_Q[b,v,6];
    }
    
// VE
for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J_Q)
    for(t in 1:T_Q) VE_Q[c,k,j,1,t] = 1 - RR_symp_Q[c,2,k,j,t] ;

for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J_Q)
    for(t in 1:T_Q) VE_Q[c,k,j,2,t] = 1 - RR_hosp_Q[c,2,k,j,t] ; 

for(c in 1:C)
 for(k in 1:K)
  for(r in 1:R)
   for(t in 1:T_Q) VE_BKRT_Q[c,k,r,t] = mean(VE_Q[c,k, ,r,t]);

for(c in 1:C)
 for(j in 1:J_Q)
  for(r in 1:R)
   for(t in 1:T_Q) VE_BJRT_Q[c,j,r,t] = mean(VE_Q[c, ,j,r,t]);
   
// ################################ DENGVAXIA // ###############################

// Attack rates

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_De)
   for(d in 1:D_De)
     H_AR_BVJD_De[b,v,j,d] = sum(Ho_BVKJD_De[b,v, ,j,d]) / pop_HOSP_BVJD_De[b,v,j,d];

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J_De)
    H_AR_BVKJ_De[b,v,k,j] = sum(Ho_BVKJD_De[b,v,k,j, ]) / mean(pop_HOSP_BVJD_De[b,v,j, ]);
   
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_De)
   V_AR_BVJ_De[b,v,j] = sum(Sy_BVKJ_De[b,v, ,j]) / pop_VCD_BVJ_De[b,v,j]; 
   
 for(v in 1:V)
  for(k in 1:K)
   V_AR_VK_De[v,k] = sum(Sy_VKJ_De[v,k, ]) / sum(pop_VCD_VJ_De[v, ]) ; 


// VE

for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J_De)
    for(t in 1:T_De){
VE_De[c,k,j,1,t] = 1 - RR_symp_De[c,2,k,j,t] ;
VE_De[c,k,j,2,t] = 1 - RR_hosp_De[c,2,k,j,t] ; 

}

for(c in 1:C)
 for(k in 1:K)
  for(r in 1:R)
   for(t in 1:T_De)
    VE_BKRT_De[c,k,r,t] = mean(VE_De[c,k, ,r,t]);

for(c in 1:C)
 for(j in 1:J_De)
  for(r in 1:R)
   for(t in 1:T_De)
   VE_BJRT_De[c,j,r,t] =mean(VE_De[c, ,j,r,t]);
   
// ################################ BUTANTAN-DV ################################

// Attack rates
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:2)
   AR_BVK_Bu[b,v,k] = sum(Sy_BVKJ_Bu[b,v,k, ]) / sum(pop_BVJ_Bu[b,v, ]); 
   
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J_Bu)
   AR_BVJ_Bu[b,v,j] = sum(Sy_BVKJ_Bu[b,v, ,j]) / pop_BVJ_Bu[b,v,j] ; 
   
// VE
for(c in 1:C)
 for(k in 1:2)
  for(j in 1:J_Bu)
    for(t in 1:T_Bu)
VE_Bu[c,k,j,t] = (1 - RR_symp_Bu[c,2,k,j,t]) ;
}
