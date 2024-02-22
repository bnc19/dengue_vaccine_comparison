data{
   int<lower = 1> C; // No. current immune status
   int<lower = 1> B; // No. baseline serostatus
   int<lower = 1> K; // No. serotypes in trial (2)
   int<lower = 1> M; // No. serotypes (4)
   int<lower = 1> V; // No. trial arms 
   int<lower = 1> J; // No. age groups
   int<lower = 1> T; // No. model time points (0:53)
   int<lower = 0> time[T]; // time 
   real<lower = 0> mu[B,K]; // titres at T0
   array [B] real<lower = 0> hs;               // short term decay 
            real<lower = 12> hl;               // long term decay 
   array [B] real<lower = -12> ts;             // time switch decay rat
   real<lower = 0, upper = 1> sens ;           // baseline test sensitivity 
   real<lower = 0, upper = 1> spec ;           // baseline test specificity 
   int<lower = 0> HI; // period of heterotypic immunity
   
   array[J] int<lower=0> SP_J; // baseline seropositive 
   array[J] int<lower=0> pop_J;    // baseline pop
   int<lower=0> VCD;    // total VCD 
   array[B*V*K] int<lower=0> VCD_BVK;  // VCD 
   array[B*V*J] int<lower=0> VCD_BVJ;   // VCD 
   array[B,V,J] int<lower=0> pop_BVJ;   // pop 
   int<lower=0> pop;   // total pop 
   
   // FLAGS
   int<lower = 0, upper = 1> single_lc;        // 0 for serostatus and or serotype lc50 based on other flags / 1 for single lc50 
   int<lower = 0, upper = 2> include_beta;     // 0 if no age-specific nc50, 1= change age groups 1 and 2, 2 = change age grp 1 
   int<lower = 0, upper = 2> mono_lc_SN;       // 0 for serotype specific lc SN / 1 for single / 2 for serotype-specific lc MO with kappa offset  
   int<lower = 0, upper = 2> mono_lc_MU;       // 0 for serotype specific lc MU / 1 for single / 2 for serotype-specific lc MO with omega offset 
   int<lower = 0, upper = 1> rho_K;            // 0 for mono rho, 1 for serotype rho 
   int<lower = 0, upper = 1> L_K;              // 0 for mono L, 1 for serotype K 
   int<lower = 0, upper = 2> w_CK;             // 0 for mono w, 1 for serostatus w, 2 for serotype w 
   real<lower = 0> L_mean;                     // mean (trunc normal) L prior
   real<lower = 0> L_sd;                       // sd (trunc normal) L prior 
   int<lower = -1, upper = 0> lower_bound_L;   // lower bound of L  (-1, 0)
   int<lower = 0, upper = 1> MU_test_SN;       // multitypic test SN at baseline (T/F)
   int<lower = 0, upper = 1>  MU_symp;         // 0 for only 1' and 2' symp infections, 1 for post-sec symp inf   
   int<lower = 0, upper = 1> enhancement;      // 0 for no vac enhancement of SN, 1 for enhancement 
   int<lower = 0, upper = 1> inc_FOIJ;          // 0 for no age-specific FOI, 1 for youngest age-specific 
}

parameters{
  array [K] real<lower = 0> lambda_K;         // foi  
  array [J] real<lower = 0, upper = 1> p;     // probability of exposure before T0
  real<lower = 0, upper = 1> gamma;           // prob symp 1' 
  array [K] real<lower = 0> rho;              // risk symp 2' rel to 1'
  real<lower = 0, upper = 1> phi;             // prob symp 3'/4' rel to 1'
  array [K] real<lower = lower_bound_L> L;    // enhancement  
  array [C] real<lower = 0> w;                // shape parameter 
  array [C,K] real<lower = 0> lc;             // 50% symp protection
  real<lower = 0> omega;                      // reduction in titre for protection against MU compared to MO 
  real<lower = 0> kappa;                      // increase in titre for protection against SN compared to MO 
  array [(J-1)] real  beta;                   // increase in lc50 for protection in younger 
  real<lower = 0>  FOI_J1;                    // age specific FOI 
}

transformed parameters{
  array[J] real<lower=0, upper=1> pSP; // prob SP (baseline for likelihood)
  real ll; // log-likelihood passed to model block and then added to target
  
  array[B,V,K,J] real<lower = 0> Sy_BVKJ;
  array [V,K,J]  real<lower = 0> Sy_VKJ;

  array[B,K,T]     real<lower = 0> n;   // titres
  array[C,V,K,J,T] real<lower = 0> RR_symp;                  

  // matrix distribution of cases
  vector [B*V*K]   mD_BVK ;    
  vector [B*V*J]   mD_BVJ  ;  
  
  real<lower=0> shape1; // pk2 prior
  real<lower=0> shape2; 
  
  //  distribution of cases
  real<lower=0, upper=1> pC;
  
  { // this { defines a block within which variables declared are local (can't have lower or upper)
  
  real lambda_m;
  real pm; 
  
  // p by lambda 
  lambda_m = mean(lambda_K); // mean across D1 and D1
  pm = 1-exp(-38.5*12*lambda_m); // always prob exposure to a single serotype 
 
  if(pm<1e-3) pm=1e-3;
  if(pm>0.99) pm=0.99;
  real vari = 0.02*pm *(1-pm) ; 
  shape1 = (((1-pm) / vari) - (1/pm)) * (pm^2) ;
  shape2 = shape1 * (1/pm - 1);

  vector[K] rhoT; // keep code simple by always referring to vector length K
  vector[K] gammaT;
  vector[K] phiT;
  
  real p_KJ[M,J]; // serotype age prob exposure 
  
  array [J] real TpSN;              // true SN
  array [M,J] real TpMO;            // true monotypic  
  array [(M+2),J] real TpMU2;       // true multitypic
  array [M,J] real TpMU3;           // true multitypic  
  array [J]   real TpMU4;           // true multitypic  
  
  array [B,V,J,(T+1)]  real pSN;       // prob SN time plus baseline 
  array [B,V,M,J,(T+1)] real pMO ;     // prob monotypic time plus baseline 
  array [B,V,(M+2),J,(T+1)] real pMU2; // prob multitypic time plus baseline 
  array [B,V,M,J,(T+1)] real pMU3;     // prob multitypic time plus baseline 
  
  array [C,K,T] real n_C;                   // allow for multitypic titres 
  array [C,K]   real L_C;                   // MO and MU L = 1 
  array [C,K,J] real nc50;                  // nc50
  array [K,J,T] real lambda;                // FOI 

  array [B,V,K,J,T] real Inc1;
  array [B,V,K,J,T] real Inc2;
  array [B,V,K,J,T] real Inc3;
  array [B,V,K,J,T] real Inc4;
     
  array [B,V,K,J,T] real D1;
  array [B,V,K,J,T] real D2;
  array [B,V,K,J,T] real D34;
  array [B,V,K,J]   real Di;
  
  // cases 
  array [K,J]   real Sy_KJ;
  array [J]     real Sy_J;
  real Sy;
  
  // prop cases
  array [B,V,K]   real pD_BVK; 
  array [B,V,J]   real pD_BVJ; 

// BIPHASIC TITRES 

  array [B] real pi_1;   // decay rate 1 
  real pi_2;      // decay rate 2

for(b in 1:B) pi_1[b] = -log(2) / hs[b];
pi_2 = -log(2) / hl ; 

for(b in 1:B)
 for(k in 1:K)
  for(t in 1:T)
   n[b,k,t] = mu[b,k] * (exp(pi_1[b] * time[t] + pi_2 * ts[b]) + exp(pi_2 * time[t] + pi_1[b] * ts[b])) / (exp(pi_1[b] * ts[b]) + exp(pi_2 * ts[b])) ; 

// MO and MU have SP titres 
 for(k in 1:K)
  for(t in 1:T){
    n_C[1,k,t] =  n[1,k,t];
    n_C[2,k,t] =  n[2,k,t]; 
    n_C[3,k,t] =  n[2,k,t];
  }
  
// MODEL SELECTION 
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
      
      
// INITIAL CONDITIONS
for(m in 1:M)
 for(j in 1:J) p_KJ[m,j] = p[j] ; 

// this defines the true serostatus populations 
 for(j in 1:J){
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
 for(j in 1:J)
   for(m in 1:M){
    pSN[1,v,j,1] = spec * TpSN[j] ;
    pSN[2,v,j,1] = (1-spec) * TpSN[j] ;
    pMO[1,v,m,j,1] = (1- sens) * TpMO[m,j] ;
    pMO[2,v,m,j,1] = sens * TpMO[m,j] ;
   }

if(MU_test_SN == 0){ // assume multitypic all classified SP
  for(v in 1:V)
   for(j in 1:J) {
    for(m in 1:6){  // there are 6 MU_2 combinations  
        pMU2[1,v,m,j,1] = 0;
        pMU2[2,v,m,j,1] = TpMU2[m,j] ; 
    }
    for(m in 1:M){  // there are 4 MU_3 combinations 
        pMU3[1,v,m,j,1] = 0; 
        pMU3[2,v,m,j,1] = TpMU3[m,j] ;
    } }
     } else { // allow mulitypic misclassification 
    for(v in 1:V)
     for(j in 1:J){
      for(m in 1:6){  // there are 6 MU_2 combinations  
          pMU2[1,v,m,j,1] = (1- sens) * TpMU2[m,j] ;
          pMU2[2,v,m,j,1] = sens * TpMU2[m,j] ; 
      }
      for(m in 1:M){  // there are 4 MU_3 combinations 
          pMU3[1,v,m,j,1] = (1- sens) * TpMU3[m,j] ; 
          pMU3[2,v,m,j,1] = sens * TpMU3[m,j] ;
      }}}


// probabilities of testing seropositive 
if(MU_test_SN == 0){  
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j];
} else{
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * (sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j]) ;
}

// VACCINE RISK RATIO 

// age-group titre offsets 
array [K] real e_beta = {exp(beta[1]), exp(beta[2])};
      
if(single_lc == 1){
  for(c in 1:C)
   for(k in 1:K) 
    for(j in 1:J)
    nc50[c,k,j] = exp(lc[1,1]); 
} else {

if(mono_lc_SN == 1){ // SN, oldest
  for(k in 1:K) nc50[1,k,3] = exp(lc[1,1]); 
   } else if (mono_lc_SN == 0) {
  for(k in 1:K) nc50[1,k,3] = exp(lc[1,k]); 
} else if(mono_lc_SN == 2){
  for(k in 1:K)  nc50[1,k,3] = exp(lc[2,k]) * exp(kappa); 
}

if(mono_lc_MU == 1){ // MU, oldest
 for(k in 1:K)  nc50[3,k,3] = exp(lc[3,1]); 
   } else if(mono_lc_MU == 0) {
 for(k in 1:K) nc50[3,k,3] = exp(lc[3,k]); 
} else if (mono_lc_MU == 2){
   for(k in 1:K)  nc50[3,k,3] = exp(lc[2,k]) * exp(-omega); 
}

// MO, oldest
for(k in 1:K)  nc50[2,k,3] = exp(lc[2,k]); // MO, serotype, oldest

// age group offsets 
if(include_beta == 1){
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1] =   nc50[c,k,3] * e_beta[1]; // youngest
   nc50[c,k,2] =   nc50[c,k,3] * e_beta[2]; // middle 
 } 
} else if(include_beta == 2) {
  for(c in 1:C)
   for(k in 1:K){
   nc50[c,k,1] =   nc50[c,k,3] * e_beta[1]; // youngest
   nc50[c,k,2] =   nc50[c,k,3]; // middle 
   }
 } else {
  for(c in 1:C)
   for(k in 1:K){ // no beta offset 
   nc50[c,k,1] = nc50[c,k,3];
   nc50[c,k,2] = nc50[c,k,3];   
  }
}

} 

// Enhancement 
for(k in 1:K){
  if(enhancement == 1){ 
   if(L_K == 0){
   L_C[1,k] = 1 + L[1]; // symp 
  } else if (L_K == 1){
   L_C[1,k] = 1 + L[k]; // symp 
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
 for(k in 1:K)
  for(j in 1:J)
   for(t in 1:T){
     if(w_CK == 0){
     RR_symp[c,2,k,j,t] = L_C[c,k] / (1 +  (n_C[c,k,t] / nc50[c,k,j])^w[1]) ;
     } else if(w_CK == 1){
     RR_symp[c,2,k,j,t] = L_C[c,k] / (1 +  (n_C[c,k,t] / nc50[c,k,j])^w[c]) ;
     } else{
     RR_symp[c,2,k,j,t] = L_C[c,k] / (1 +  (n_C[c,k,t] / nc50[c,k,j])^w[k]) ;
     }
     RR_symp[c,1,k,j,t] =  1 ;
    }

 // FOI 
 
real FOI_J[J] = {FOI_J1, 1, 1}  ; // scale youngest group only 

if(inc_FOIJ == 0) {
for(k in 1:K)
 for(j in 1:J)
  for(t in 1:T)  lambda[k,j,t] = exp(-lambda_K[k]) ; 
} else{
for(k in 1:K)
 for(j in 1:J)
  for(t in 1:T)  lambda[k,j,t] = exp(-lambda_K[k] * FOI_J[j]) ; 
}

// SURVIVAL MODEL - prob of surviving each time point without infection - only D1 and D2 lambda 
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(t in 1:T) pSN[b,v,j,(t+1)] = lambda[1,j,t] * lambda[2,j,t] * pSN[b,v,j,t];
   
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J) {
    for(t in 1:T) pMO[b,v,1,j,(t+1)] = lambda[2,j,t] * pMO[b,v,1,j,t] ; // escape 2 
    for(t in 1:T) pMO[b,v,2,j,(t+1)] = lambda[1,j,t] * pMO[b,v,2,j,t] ; // escape 1 
    for(t in 1:T) pMO[b,v,3,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t] * pMO[b,v,3,j,t] ; // escape 1 and 2 
    for(t in 1:T) pMO[b,v,4,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t] * pMO[b,v,4,j,t] ; // escape 1 and 2 
    
    for(t in (HI+1):T) 
     for(k in 1:K) pMO[b,v,k,j,(t+1)] += (1 - gammaT[k] * RR_symp[1,v,k,j,(t-HI)]) * (1 - lambda[k,j,(t - HI)]) * pSN[b,v,j,(t-HI)] ; // new D1 and D2 infections during trial
    }
    
    
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J) {
      for(t in 1:T){
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
       for(t in (HI+1):T) { // new D1 and D2 infections 
          // MU_12
          pMU2[b,v,1,j,(t+1)] += (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1 - lambda[1,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] + (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) *  (1-lambda[2,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_13
          pMU2[b,v,2,j,(t+1)] += (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1 - lambda[1,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] ;
          // MU_14
          pMU2[b,v,3,j,(t+1)] += (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1 - lambda[1,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] ;
          // MU_23
          pMU2[b,v,4,j,(t+1)] += (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1 - lambda[2,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)]  ;
          // MU_24
          pMU2[b,v,5,j,(t+1)] += (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1 - lambda[2,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] ;
        }}

for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J){
     for(t in 1:T){ 
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] = lambda[1,j,t]  *  pMU3[b,v,1,j,t] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] = lambda[2,j,t]  *  pMU3[b,v,2,j,t] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] = pMU3[b,v,3,j,t] ; 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] = pMU3[b,v,4,j,t] ; 
     }
    for(t in (HI+1):T) {
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] +=  (1 - phiT[2] * gammaT[2] * RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; // 34 then 2 infection 
        // MU3 -2           
        pMU3[b,v,2,j,(t+1)] +=  (1 - phiT[1] * gammaT[1] * RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; // 34 then 1 infection 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] +=  (1 - phiT[2] * gammaT[2] * RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] * RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] ;  // 14 then 2 or 24 then 1 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] +=  (1 - phiT[2] * gammaT[2] * RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[1] * gammaT[1] * RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] ;  // 13 then 2 or 23 then 1 
      }
    }  

// calculate infection incidence, symptomatic, hospitalised 
for(b in 1:B)
  for(v in 1:V)
    for(j in 1:J) {
      for(t in 1:T){
         Inc3[b,v,1,j,t] = (1 - lambda[1,j,t]) * (pMU2[b,v,4,j,t] + pMU2[b,v,5,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,2,j,t] = (1 - lambda[2,j,t]) * (pMU2[b,v,2,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,6,j,t]);  
         
         for(k in 1:K){
           Inc1[b,v,k,j,t] = (1 - lambda[k,j,t]) * pSN[b,v,j,t];  
           Inc2[b,v,k,j,t] = (1 - lambda[k,j,t]) * (sum(pMO[b,v, ,j,t]) - pMO[b,v,k,j,t]);
           Inc4[b,v,k,j,t] = (1 - lambda[k,j,t]) * pMU3[b,v,k,j,t];  
           D1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * RR_symp[1,v,k,j,t] ;
           D2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k] * RR_symp[2,v,k,j,t];
           D34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * RR_symp[3,v,k,j,t];  
         }}}

// Aggregate to match published time points e.g. 1-24 months 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J){
   Di[b,v,k,j] = sum(D1[b,v,k,j,1:T])  + sum(D2[b,v,k,j,1:T]) ; 
   if(MU_symp == 1)  Di[b,v,k,j]   +=  sum(D34[b,v,k,j,1:T]) ;// if post-secondary can be symptomatic 
   }

// Cases 
for(b in 1:B) 
 for(v in 1:V) 
  for(k in 1:K)
   for(j in 1:J) 
          Sy_BVKJ[b,v,k,j] = Di[b,v,k,j] * sum(pop_BVJ[ ,v,j]); 

// sum over symp
for(v in 1:V) 
 for(k in 1:K) 
  for(j in 1:J) 
  Sy_VKJ[v,k,j] = sum(Sy_BVKJ[ ,v,k,j]) ; 
    
for(k in 1:K) 
 for(j in 1:J) 
 Sy_KJ[k,j] = sum(Sy_VKJ[ ,k,j]) ;  

for(j in 1:J) 
 Sy_J[j] = sum(Sy_KJ[ ,j]) ;  
 Sy = sum(Sy_J) ; 

for(b in 1:B)
 for(v in 1:V)
   for(k in 1:K)
    for(j in 1:J){
  pD_BVK[b,v,k] = sum(Sy_BVKJ[b,v,k, ]) / Sy ; // serotype serostatus symp
  pD_BVJ[b,v,j] = sum(Sy_BVKJ[b,v, ,j]) / Sy ;   // age serostatus symp 
  }

// matrix for likelihood function 

for(k in 1:K){ // 2 serotypes 
  mD_BVK[k]   =   pD_BVK[1,1,k] ; // SN C
  mD_BVK[k+2] =   pD_BVK[1,2,k] ; // SN V
  mD_BVK[k+4] =   pD_BVK[2,1,k] ; // SP C
  mD_BVK[k+6] =   pD_BVK[2,2,k] ; // SP V
}

for(j in 1:J){ // 3 ages 
   mD_BVJ[j]     =  pD_BVJ[1,1,j]; //  SN C 
   mD_BVJ[j+3]   =  pD_BVJ[1,2,j]; //  SN V
   mD_BVJ[j+6]   =  pD_BVJ[2,1,j]; //  SP C 
   mD_BVJ[j+9]   =  pD_BVJ[2,2,j]; //  SP V
  }
  
 pC = Sy / pop; 

// log likelihood 
 ll=0;
 ll += binomial_lpmf(SP_J | pop_J, pSP);
 ll += binomial_lpmf(VCD  | pop, pC);
 ll += multinomial_lpmf(VCD_BVK  | mD_BVK) ;
 ll += multinomial_lpmf(VCD_BVJ  | mD_BVJ) ;
 
}
}

model {
 
 target += ll;

// priors
  p[1] ~ beta(3,5);
  p[2] ~ beta(3,5);
  p[3] ~ beta(shape1,shape2) ; // prob exposure to a serotype   
  lambda_K ~ lognormal(-7,2);
  gamma ~ normal(0.85,0.2);  // prob secondary
  rho ~ normal(1.97,0.40);   // RR secondary
  phi ~ normal(0.25,0.05);
  L ~ normal(L_mean,L_sd);
  w ~ normal(1,2);
  lc[1,] ~ normal(5.18,1);
  for(c in 2:C)  lc[c,] ~ normal(6.27,1);
  omega ~ normal(0,2);
  kappa ~  normal(0,2);
  beta ~ normal(0,2);
  FOI_J1 ~ normal(1,1);
}

generated quantities{
// AR
array [B,V,K] real<lower = 0 > AR_BVK ;
array [B,V,J] real<lower = 0 > AR_BVJ ;

// VE
array [C,K,J,T] real<upper =1 > VE  ;
  
// log lik
vector [6] log_lik;

for(j in 1:J) log_lik[j] = binomial_lpmf(SP_J[j] | pop_J[j], pSP[j]);
log_lik[4] = binomial_lpmf(VCD   | pop, pC);
log_lik[5] = multinomial_lpmf(VCD_BVK  | mD_BVK) ;
log_lik[6] = multinomial_lpmf(VCD_BVJ  | mD_BVJ) ;

// Attack rates
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   AR_BVK[b,v,k] = sum(Sy_BVKJ[b,v,k, ]) / sum(pop_BVJ[b,v, ]); 
   
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   AR_BVJ[b,v,j] = sum(Sy_BVKJ[b,v, ,j]) / pop_BVJ[b,v,j] ; 


// VE
for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J)
    for(t in 1:T)
VE[c,k,j,t] = (1 - RR_symp[c,2,k,j,t]) ;

}
