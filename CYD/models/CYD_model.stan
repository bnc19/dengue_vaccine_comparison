data{
   int<lower = 1> C; // No. current immune status
   int<lower = 1> B; // No. baseline serostatus
   int<lower = 1> K; // No. serotypes
   int<lower = 1> D; // No. vcd data
   int<lower = 1> V; // No. trial arms 
   int<lower = 1> J; // No. age groups
   int<lower = 1> T; // No. model time points (0:53)
   int<lower = 1> R; // No. outcomes 
   
   int<lower = 0> time[T]; // time 
   real<lower = 0> mu[B,K]; // titres at T0
   
   int<lower = 0> HI; // period of heterotypic immunity
   
   array[J] int<lower=0> SP_J;     // baseline seropositive
   array[J] int<lower=0> pop_J;    // baseline pop
   int<lower=0> VCD;               // total VCD 
   array[D] int<lower=0> HOSP_D;   // total hosp 
   
   array[V*K]     int<lower=0> VCD_VK;        // VCD 
   array[B*V*J]   int<lower=0> VCD_BVJ;       // VCD 
   array[B*V*J,D] int<lower=0> HOSP_BVJD;     // hosp 
   array[B*V*K*J] int<lower=0> HOSP_BVKJ;     // hosp  
   array[B,V,J]   int<lower=0> pop_VCD_BVJ;   // pop 
   array[B,V,J,D] int<lower=0> pop_HOSP_BVJD; // pop
   
   int<lower=0> pop_VCD ; 
   int<lower=0> pop_HOSP; // mean hosp pop across trial 
   array[D] int<lower=0> pop_HOSP_D ; // hosp pop over time 
   
   // FLAGS
   int<lower = 0, upper = 1> include_pK2;      // include serotype specific p? (T/F)
   int<lower = 0, upper = 1> include_eps;      // include enhanced secondary hosp? (T/F)
   int<lower = 0, upper = 2> include_beta;     // 0 if no age-specific nc50, 1 = change age groups 1 for both outcomes, 2 = age grp 1 sep for each outcome
   int<lower = 0, upper = 2> mono_lc_SN;       // 0 for serotype specific lc SN / 1 for single / 2 for serotype-specific lc MO with kappa offset  
   int<lower = 0, upper = 2> mono_lc_MU;       // 0 for serotype specific lc MU / 1 for single / 2 for serotype-specific lc MO with omega offset 
   int<lower = 0, upper = 1> rho_K;            // 0 for mono rho, 1 for serotype rho 
   int<lower = 0, upper = 1> L_K;              // 0 for mono L, 1 for serotype K 
   int<lower = 0, upper = 2> w_CK;             // 0 for mono w, 1 for serostatus w, 2 for serotype w 
   int<lower = 0, upper = 2> alpha_CK;         // 0 for mono alpha, 1 for serostatus alpha, 2 for serotype alpha 
   int<lower = 0, upper = 1> tau_K;            // 0 for mono tau, 1 for serotype tau
   int<lower = 0, upper = 2> delta_KJ;         // 0 for single delta, 1 for serotype-specific delta, 2 for age delta
   int<lower = 0, upper = 1> FOI_J;            // 0 if same FOI for all ages, 1 if FOI scaled for youngest age-group
   real<lower = 0> L_mean;                     // mean (trunc normal) L prior
   real<lower = 0> L_sd;                       // sd (trunc normal) L prior 
   int<lower = -1, upper = 0> lower_bound_L;   // lower bound of L  (-1, 0)
   int<lower = 0, upper = 1> MU_test_SN;       // multitypic test SN at baseline (T/F)
   int<lower = 0, upper = 1>  MU_symp;         // 0 for only 1' and 2' symp infections, 1 for post-sec symp inf   
   int<lower = 0, upper = 1> enhancement;      // 0 for no vac enhancement of SN, 1 for enhancement 
   
}

transformed data {
  array [V,J,D] int pop_HOSP_VJD;
  array [V,J] int pop_VCD_VJ;
// pop not by serostatus 
for(v in 1:V)
 for(j in 1:J)
  for(d in 1:D){
    pop_HOSP_VJD[v,j,d] = sum(pop_HOSP_BVJD[ ,v,j,d]) ; 
    pop_VCD_VJ[v,j] = sum(pop_VCD_BVJ[ ,v,j]) ;   
}
}

parameters{
  array [B] real<lower = 0> hs;               // short term decay 
            real<lower = 12> hl;              // long term decay 
  real<lower = -12> ts;                       // time switch decay rate
  array [D] real<lower = 0> lambda_D;         // foi  
  simplex[K] theta ;                          // proportion of each serotype circulating in the trial               
  array [J] real<lower = 0, upper = 1> p;     // probability of exposure before T0
  array [K] real<lower = 0, upper = 1> pK2;   // serotype-specific p in age group 2
  real<lower = 1> epsilon ;                   // prob hosp 2' rel to 1'/3'/4'
  real<lower = 0, upper = 1> gamma;           // prob symp 1' 
  array [K] real<lower = 0> rho;              // risk symp 2' rel to 1'
  real<lower = 0, upper = 1> phi;             // prob symp 3'/4' rel to 1'
  array [K] real<lower = 0, upper = 1> delta; // prob hosp
  array [K] real<lower = lower_bound_L> L;    // enhancement  
  array [K] real<lower = 0> tau;              // multiply (L-1) for hosp
  array [K] real<lower = 0> w;                // shape parameter 
  array [C,K] real<lower = 0> lc;             // 50% symp protection
  array [K] real<lower = 0> alpha;            // reduction in titre for protection against hosp
  real<lower = 0> omega;                      // reduction in titre for protection against MU compared to MO 
  real<lower = 0> kappa;                      // increase in titre for protection against SN compared to MO 
  array [J] real  beta;                       // increase in lc50 for protection in younger 
  real<lower = 0, upper = 1> sens ;           // baseline test sensitivity 
  real<lower = 0, upper = 1> spec ;           // baseline test specificity 
  real<lower = 0> scale_FOI;                  // scale youngest FOI 
}

transformed parameters{
  array[J] real<lower=0, upper=1> pSP; // prob SP (baseline for likelihood)
  real ll; // log-likelihood passed to model block and then added to target
  
  array [V,K,J]    real<lower = 0> Sy_VKJ;
  array[B,V,K,J]   real<lower = 0> Sy_BVKJ;          // all cases 
  array[B,V,K,J,D] real<lower = 0> Ho_BVKJD;  

  array[B,K,T]     real<lower = 0> n;                // titres
  array[C,V,K,J,T] real<lower = 0> RR_symp;                  
  array[C,V,K,J,T] real<lower = 0> RR_hosp;

  // matrix distribution of cases
  vector [B*V*J]   mD_BVJ ;    
  vector [V*K]     mD_VK ;  
  matrix [B*V*J,D] mH_BVJD ;    
  vector [B*V*K*J] mH_BVKJ ;   
  
  // so can access in generated quanitites 
  array [B,V,K,J,T] real<lower = 0> D1;
  array [B,V,K,J,T] real<lower = 0> D2;
  array [B,V,K,J,T] real<lower = 0> D34;
  array [B,V,K,J,T] real<lower = 0> H1;
  array [B,V,K,J,T] real<lower = 0> H2;
  array [B,V,K,J,T] real<lower = 0> H34;
  
  real<lower=0> shape1; // pk2 prior
  real<lower=0> shape2; 
  
  //  distribution of cases
  real<lower=0, upper=1> pC;
  array[D] real<lower=0, upper=1> pH;
  
  { // this { defines a block within which variables declared are local (can't have lower or upper)
  
  array[D] real lambda_k;
  real lambda_m;
  real pm; 
  
  // contrain pk2 / p by lambda 
  for(d in 1:D) lambda_k[d] = lambda_D[d]  / 4; // lambda_D is cumulative FOI so divide by 4 to make it mean serotype FOI
  lambda_m =(lambda_k[1] * 13 + lambda_k[2] * 11 + lambda_k[3] * 12 + lambda_k[4] * ( T-13 - 11 - 12)) / T ; // weighted FOI over trial 
  pm = 1-exp(-12.5*12*lambda_m); // always prob exposure to a single serotype 
 
  if(pm<1e-3) pm=1e-3;
  if(pm>0.99) pm=0.99;
  real vari = 0.02*pm *(1-pm) ; 
  shape1 = (((1-pm) / vari) - (1/pm)) * (pm^2) ;
  shape2 = shape1 * (1/pm - 1);

  vector[K] rhoT; // keep code simple by always referring to vector length K
  vector[K] gammaT;
  vector[K] phiT;
  vector[K] deltaT;

  real hK2[K]; // age3 historic hazard
  real q[K];   // historic serotype proportion 
  real hK1[K]; // age1 historic hazard 
  real p_KJ[K,J]; // serotype age prob exposure 
  
  array [J] real TpSN;              // true SN
  array [K,J] real TpMO;            // true monotypic  
  array [(K+2),J] real TpMU2;       // true multitypic
  array [K,J] real TpMU3;           // true multitypic  
  array [J]   real TpMU4;           // true multitypic  
  
  array [B,V,J,(T+1)]  real pSN;       // prob SN time plus baseline 
  array [B,V,K,J,(T+1)] real pMO ;     // prob monotypic time plus baseline 
  array [B,V,(K+2),J,(T+1)] real pMU2; // prob multitypic time plus baseline 
  array [B,V,K,J,(T+1)] real pMU3;     // prob multitypic time plus baseline 
  
  array [C,K,T] real n_C;                     // allow for multitypic titres 
  array [C,K,R] real L_C;                     // MO and MU L = 1 
  array [C,K,J,R] real nc50;                  // nc50
  array [K,J,T] real lambda;                  // FOI 
  array [J] real scale_FOI_J ;                // scale FOI_J 

  array [B,V,K,J,T] real Inc1;
  array [B,V,K,J,T] real Inc2;
  array [B,V,K,J,T] real Inc3;
  array [B,V,K,J,T] real Inc4;
  array [B,V,K,J]   real Di;
  array [B,V,K,J,D] real H;
  
  // cases 
  array [V,K,J,D] real Ho_VKJD;
  array [K,J,D]   real Ho_KJD;
  array [K,J]     real Sy_KJ;
  array [J,D]     real Ho_JD;
  array [J]       real Sy_J;
  array [D]       real Ho_D;
  real Sy;
  
  // prop cases
  array [B,V,J]   real pD_BVJ; 
  array [V,K]     real pD_VK; 
  array [B,V,J,D] real pH_BVJD;
  array [B,V,K,J] real pH_BVKJ;

// BIPHASIC TITRES 

  array [B] real pi_1;   // decay rate 1 
  real pi_2;      // decay rate 2

for(b in 1:B) pi_1[b] = -log(2) / hs[b];
pi_2 = -log(2) / hl ; 

for(b in 1:B)
 for(k in 1:K)
  for(t in 1:T)
   n[b,k,t] = mu[b,k] * (exp(pi_1[b] * time[t] + pi_2 * ts) + exp(pi_2 * time[t] + pi_1[b] * ts)) / (exp(pi_1[b] * ts) + exp(pi_2 * ts)) ; 

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
      
if(delta_KJ==1) {
  for(k in 1:K) deltaT[k] = delta[k];
} else if(delta_KJ == 2) {
    for(j in 1:J) deltaT[j] = delta[j];
    deltaT[3] = 1;
    deltaT[4] = 1;
} else {
  for(k in 1:K) deltaT[k] = delta[1];
}      
      
      
// INITIAL CONDITIONS
if(include_pK2 == 0) {  // p common to all serotypes
  for(k in 1:K){
    for(j in 1:J) p_KJ[k,j] = p[j]; // if include_pK2 == 0, all p are prob exposure to a single serotype 
    hK2[k] = 0;
    q[k] = 0;
    hK1[k] = 0;
    }
 } else { // serotype-specific p
  for(k in 1:K) hK2[k] = -log(1-pK2[k]) ; // pK2 prob exposure to serotype k 
  real shK2;
  shK2=sum(hK2);
  for(k in 1:K) {
    q[k] = hK2[k] / shK2 ;
    hK1[k] = q[k] * -log(1-p[1]) ; // if include_pK2 == 1, p1 is cumulative prob exposure 
    p_KJ[k,1] = 1 - exp(-hK1[k]) ;
    p_KJ[k,2] = pK2[k];
  }
}

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
   for(k in 1:K){
    pSN[1,v,j,1] = spec * TpSN[j] ;
    pSN[2,v,j,1] = (1-spec) * TpSN[j] ;
    pMO[1,v,k,j,1] = (1- sens) * TpMO[k,j] ;
    pMO[2,v,k,j,1] = sens * TpMO[k,j] ;
   }

if(MU_test_SN == 0){ // assume multitypic all classified SP
  for(v in 1:V)
   for(j in 1:J) {
    for(k in 1:6){  // there are 6 MU_2 combinations  
        pMU2[1,v,k,j,1] = 0;
        pMU2[2,v,k,j,1] = TpMU2[k,j] ; 
    }
    for(k in 1:K){  // there are 4 MU_3 combinations 
        pMU3[1,v,k,j,1] = 0; 
        pMU3[2,v,k,j,1] = TpMU3[k,j] ;
    } }
     } else { // allow mulitypic misclassification 
    for(v in 1:V)
     for(j in 1:J){
      for(k in 1:6){  // there are 6 MU_2 combinations  
          pMU2[1,v,k,j,1] = (1- sens) * TpMU2[k,j] ;
          pMU2[2,v,k,j,1] = sens * TpMU2[k,j] ; 
      }
      for(k in 1:K){  // there are 4 MU_3 combinations 
          pMU3[1,v,k,j,1] = (1- sens) * TpMU3[k,j] ; 
          pMU3[2,v,k,j,1] = sens * TpMU3[k,j] ;
      }}}


// probabilities of testing seropositive 
if(MU_test_SN == 0){  
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j];
} else{
 for(j in 1:J)  pSP[j] = (1-spec) * TpSN[j] + sens * (sum(TpMO[ ,j]) + sum(TpMU2[ ,j]) + sum(TpMU3[ ,j]) + TpMU4[j]) ;
}

// VACCINE RISK RATIO 

// outcome and age-group titre offsets 
array [K] real e_beta = {exp(beta[1]), exp(beta[2])};
array [C,K] real e_alpha; 

for(k in 1:K)
 for(c in 1:C){
   if(alpha_CK == 0){ 
     e_alpha[c,k] = exp(-alpha[1]);  // mono outcome offset 
     } else if(alpha_CK == 1) {
       e_alpha[c,k] = exp(-alpha[c]);  // serostatus outcome offset 
       } else {
         e_alpha[c,k] = exp(-alpha[k]);  // serotype outcome offset 
         }}
        
if(mono_lc_SN == 1){ // SN, oldest, symp 
  for(k in 1:K) nc50[1,k,2,1] = exp(lc[1,1]); // single 
   } else if (mono_lc_SN == 0) {
  for(k in 1:K) nc50[1,k,2,1] = exp(lc[1,k]); // serotype 
} else if(mono_lc_SN == 2){
  for(k in 1:K)  nc50[1,k,2,1] = exp(lc[2,k]) * exp(kappa); // offset MO 
}

if(mono_lc_MU == 1){ // MU, oldest, symp
 for(k in 1:K)  nc50[3,k,2,1] = exp(lc[3,1]); // single 
   } else if(mono_lc_MU == 0) {
 for(k in 1:K) nc50[3,k,2,1] = exp(lc[3,k]); // serotype 
} else if (mono_lc_MU == 2){
   for(k in 1:K)  nc50[3,k,2,1] = exp(lc[2,k]) * exp(-omega); // offset MO 
}


// MO, oldest, symp
for(k in 1:K)  nc50[2,k,2,1] = exp(lc[2,k]); // MO, serotype, oldest, symp

// age group offsets 
if(include_beta == 1){ // youngest same for both outcomes 
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
  for(j in 1:J)
    nc50[c,k,j,2] = nc50[c,k,j,1] * e_alpha[c,k];

if(include_beta == 2){ // youngest by outcome 
  for(c in 1:C)
   for(k in 1:K){
     nc50[c,k,1,1] =   nc50[c,k,1,1] * e_beta[1]; 
     nc50[c,k,1,2] =   nc50[c,k,1,2] * e_beta[2]; 
   }
}

// Enhancement 
for(k in 1:K){
  if(enhancement == 1){ 
   if(L_K == 0){
   L_C[1,k,1] = 1 + L[1]; // symp 
   L_C[1,k,2] = 1 + L[1] * ((tau_K==0)?tau[1]:tau[k]);
  } else if (L_K == 1){
   L_C[1,k,1] = 1 + L[k]; // symp 
   L_C[1,k,2] = 1 + L[k] * ((tau_K==0)?tau[1]:tau[k]);
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
  for(j in 1:J)
   for(t in 1:T){
     if(w_CK == 0){
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[1]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[1]) ;
     } else if(w_CK == 1){
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[c]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[c]) ;  
     } else{
     RR_symp[c,2,k,j,t] = L_C[c,k,1] / (1 +  (n_C[c,k,t] / nc50[c,k,j,1])^w[k]) ;
     RR_hosp[c,2,k,j,t] = L_C[c,k,2] / (1 +  (n_C[c,k,t] / nc50[c,k,j,2])^w[k]) ;  
     }
     RR_symp[c,1,k,j,t] =  1 ;
     RR_hosp[c,1,k,j,t] =  1 ;
    }

if(FOI_J == 0){
  scale_FOI_J = {1,1}; 
} else{
  scale_FOI_J = {scale_FOI, 1};
}

// FOI  
for(k in 1:K)
 for(j in 1:J){
  for(t in 1:13)  lambda[k,j,t] = exp(-lambda_D[1] * scale_FOI_J[j] * theta[k]) ; // 13-25
  for(t in 14:24) lambda[k,j,t] = exp(-lambda_D[2] * scale_FOI_J[j] * theta[k]) ; // 26-36
  for(t in 25:36) lambda[k,j,t] = exp(-lambda_D[3] * scale_FOI_J[j] * theta[k]) ; // 37-48
  for(t in 37:T)  lambda[k,j,t] = exp(-lambda_D[4] * scale_FOI_J[j] * theta[k]) ; // 49- endpoint 
}


// SURVIVAL MODEL - prob of surviving each time point without infection 
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(t in 1:T) pSN[b,v,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t]*lambda[3,j,t]*lambda[4,j,t] * pSN[b,v,j,t];
   
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J) {
    for(t in 1:T) pMO[b,v,k,j,(t+1)] = lambda[1,j,t]*lambda[2,j,t]*lambda[3,j,t]*lambda[4,j,t]/lambda[k,j,t] * pMO[b,v,k,j,t] ;
    for(t in (HI+1):T) pMO[b,v,k,j,(t+1)] += (1 - gammaT[k] * RR_symp[1,v,k,j,(t-HI)]) * (1- lambda[k,j,(t-HI)]) * pSN[b,v,j,(t-HI)] ;
    }  
   
for(b in 1:B)
 for(v in 1:V)
   for(j in 1:J) {
      for(t in 1:T){
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
      for(t in (HI+1):T) {
          // MU_12
          pMU2[b,v,1,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] + (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_13
          pMU2[b,v,2,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_14
          pMU2[b,v,3,j,(t+1)] +=  (1 - gammaT[1] * rhoT[1] * RR_symp[2,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMO[b,v,1,j,(t-HI)] ;
          // MU_23
          pMU2[b,v,4,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] + (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_24
          pMU2[b,v,5,j,(t+1)] +=  (1 - gammaT[2] * rhoT[2] * RR_symp[2,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMO[b,v,2,j,(t-HI)] ;
          // MU_34
          pMU2[b,v,6,j,(t+1)] +=  (1 - gammaT[3] * rhoT[3] * RR_symp[2,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMO[b,v,4,j,(t-HI)] + (1 - gammaT[4] * rhoT[4] * RR_symp[2,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMO[b,v,3,j,(t-HI)] ;
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
        pMU3[b,v,3,j,(t+1)] = lambda[3,j,t]  *  pMU3[b,v,3,j,t] ; 
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] = lambda[4,j,t]  *  pMU3[b,v,4,j,t] ; 
     }
       for(t in (HI+1):T) {
        // MU3 -1 
        pMU3[b,v,1,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] + (1 - phiT[3] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -2 
        pMU3[b,v,2,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,6,j,(t-HI)] ; 
        // MU3 -3
        pMU3[b,v,3,j,(t+1)] += (1 - phiT[4] * gammaT[4] *  RR_symp[3,v,4,j,(t-HI)]) * (1-lambda[4,j,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,3,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,5,j,(t-HI)] ;  
        // MU3 -4
        pMU3[b,v,4,j,(t+1)] += (1 - phiT[3] * gammaT[3] *  RR_symp[3,v,3,j,(t-HI)]) * (1-lambda[3,j,(t-HI)]) * pMU2[b,v,1,j,(t-HI)] + (1 - phiT[2] * gammaT[2] *  RR_symp[3,v,2,j,(t-HI)]) * (1-lambda[2,j,(t-HI)]) * pMU2[b,v,2,j,(t-HI)] + (1 - phiT[1] * gammaT[1] *  RR_symp[3,v,1,j,(t-HI)]) * (1-lambda[1,j,(t-HI)]) * pMU2[b,v,4,j,(t-HI)] ;  
      }
    }  

// calculate infection incidence, symptomatic, hospitalised 
for(b in 1:B)
  for(v in 1:V)
    for(j in 1:J) {
      for(t in 1:T){
         Inc3[b,v,1,j,t] = (1 - lambda[1,j,t]) * (pMU2[b,v,4,j,t] + pMU2[b,v,5,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,2,j,t] = (1 - lambda[2,j,t]) * (pMU2[b,v,2,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,6,j,t]);  
         Inc3[b,v,3,j,t] = (1 - lambda[3,j,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,3,j,t] + pMU2[b,v,5,j,t]);  
         Inc3[b,v,4,j,t] = (1 - lambda[4,j,t]) * (pMU2[b,v,1,j,t] + pMU2[b,v,2,j,t] + pMU2[b,v,4,j,t]);
         
         for(k in 1:K){
           Inc1[b,v,k,j,t] = (1 - lambda[k,j,t]) * pSN[b,v,j,t]; 
           Inc2[b,v,k,j,t] = (1 - lambda[k,j,t]) * (sum(pMO[b,v, ,j,t]) - pMO[b,v,k,j,t]);
           Inc4[b,v,k,j,t] = (1 - lambda[k,j,t]) * pMU3[b,v,k,j,t];  
           D1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * RR_symp[1,v,k,j,t] ;
           D2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k] * RR_symp[2,v,k,j,t];
           D34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * RR_symp[3,v,k,j,t];  
  
           if(delta_KJ != 2){
           H1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * deltaT[k] * RR_hosp[1,v,k,j,t];
           if(include_eps == 0){
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[k]  * RR_hosp[2,v,k,j,t];
           } else{
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[k] * epsilon * RR_hosp[2,v,k,j,t]; 
           }
           H34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * deltaT[k] * RR_hosp[3,v,k,j,t] ;   
         } else {
           H1[b,v,k,j,t]  = Inc1[b,v,k,j,t] * gammaT[k] * deltaT[j] * RR_hosp[1,v,k,j,t];
           if(include_eps == 0){
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[j]  * RR_hosp[2,v,k,j,t];
           } else{
           H2[b,v,k,j,t]  = Inc2[b,v,k,j,t] * gammaT[k] * rhoT[k]  * deltaT[j] * epsilon * RR_hosp[2,v,k,j,t]; 
           }
           H34[b,v,k,j,t] = (Inc3[b,v,k,j,t] + Inc4[b,v,k,j,t]) * gammaT[k] * phiT[k] * deltaT[j] * RR_hosp[3,v,k,j,t] ;   
         }
         
         }}}

// Aggregate to match published time points e.g. 1-12 months 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J){
   Di[b,v,k,j] = sum(D1[b,v,k,j,1:13])   + sum(D2[b,v,k,j,1:13]) ; 
   
   H[b,v,k,j,1]  = sum(H1[b,v,k,j,1:13])  + sum(H2[b,v,k,j,1:13]) ;
   H[b,v,k,j,2]  = sum(H1[b,v,k,j,14:24]) + sum(H2[b,v,k,j,14:24]) ;
   H[b,v,k,j,3]  = sum(H1[b,v,k,j,25:36]) + sum(H2[b,v,k,j,25:36]) ;
   H[b,v,k,j,4]  = sum(H1[b,v,k,j,37:T])  + sum(H2[b,v,k,j,37:T]) ;
   }

if(MU_symp == 1){ // if post-secondary can be symptomatic 
for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J){
   Di[b,v,k,j]   +=  sum(D34[b,v,k,j,1:13]) ;
   
   H[b,v,k,j,1]  += sum(H34[b,v,k,j,1:13])  ;
   H[b,v,k,j,2]  += sum(H34[b,v,k,j,14:24]) ;
   H[b,v,k,j,3]  += sum(H34[b,v,k,j,25:36]) ;
   H[b,v,k,j,4]  += sum(H34[b,v,k,j,37:T])  ;
   }
  }
   

// Cases 
for(b in 1:B) 
 for(v in 1:V) 
  for(k in 1:K)
   for(j in 1:J) 
    for(d in 1:D){
          Sy_BVKJ[b,v,k,j] = Di[b,v,k,j] * pop_VCD_VJ[v,j]; // year 1 symp 
          Ho_BVKJD[b,v,k,j,d] = H[b,v,k,j,d] * pop_HOSP_VJD[v,j,d]; // hosp across trial 
            }

// Sum over hosp 
for(v in 1:V) 
 for(k in 1:K) 
  for(j in 1:J) 
   for(d in 1:D)
   Ho_VKJD[v,k,j,d] = sum(Ho_BVKJD[ ,v,k,j,d]); 
      
for(k in 1:K) 
 for(j in 1:J) 
  for(d in 1:D)
  Ho_KJD[k,j,d] = sum(Ho_VKJD[ ,k,j,d]); 

for(j in 1:J) 
 for(d in 1:D)
 Ho_JD[j,d] = sum(Ho_KJD[ ,j,d]);    

for(d in 1:D)
Ho_D[d] = sum(Ho_JD[ ,d]); 


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
    for(j in 1:J)
     for(d in 1:D){
  pD_VK[v,k] = (sum(Sy_VKJ[v,k, ])) / Sy ;  // serotype symp y1 
  pD_BVJ[b,v,j] = (sum(Sy_BVKJ[b,v, ,j])) / Sy ; // symp year 1 
  pH_BVJD[b,v,j,d] = (sum(Ho_BVKJD[b,v, ,j,d])) / Ho_D[d]; // hosp not k
  pH_BVKJ[b,v,k,j] = (sum(Ho_BVKJD[b,v,k,j, ])) / (sum(Ho_D)); // hosp not d
  }


// matrix for likelihood function 

for(k in 1:K){
  mD_VK[k] =   pD_VK[1,k] ; 
  mD_VK[k+4] =   pD_VK[2,k] ; 
}

for(j in 1:J){
   mD_BVJ[j]     =  pD_BVJ[1,1,j]; // SN C 
   mD_BVJ[j+2]   =  pD_BVJ[1,2,j]; // SN V
   mD_BVJ[j+4]   =  pD_BVJ[2,1,j]; // SP C
   mD_BVJ[j+6]   =  pD_BVJ[2,2,j]; // SP V 
  }
  
for(j in 1:J)
  for(d in 1:D){   
   mH_BVJD[j,d]     =  pH_BVJD[1,1,j,d]; // SN C 
   mH_BVJD[j+2,d]   =  pH_BVJD[1,2,j,d]; // SN V
   mH_BVJD[j+4,d]   =  pH_BVJD[2,1,j,d]; // SP C
   mH_BVJD[j+6,d]   =  pH_BVJD[2,2,j,d]; // SP V
    }
 
 for(k in 1:K){
   mH_BVKJ[k] =       pH_BVKJ[1,1,k,1]; // SN C < 9
   mH_BVKJ[k + 4] =   pH_BVKJ[1,1,k,2]; // SN C > 9
   mH_BVKJ[k + 8] =   pH_BVKJ[1,2,k,1]; // SN V < 9
   mH_BVKJ[k + 12] =  pH_BVKJ[1,2,k,2]; // SN V > 9
   mH_BVKJ[k + 16] =  pH_BVKJ[2,1,k,1]; // SP C < 9
   mH_BVKJ[k + 20] =  pH_BVKJ[2,1,k,2]; // SP C > 9
   mH_BVKJ[k + 24] =  pH_BVKJ[2,2,k,1]; // SP V < 9
   mH_BVKJ[k + 28] =  pH_BVKJ[2,2,k,2]; // SP V > 9 
 }

 pC = Sy / pop_VCD; 
for(d in 1:D)
  pH[d] = Ho_D[d] / pop_HOSP_D[d]; 
  
// log likelihood 
 ll=0;
 

    ll += binomial_lpmf(HOSP_D| pop_HOSP_D, pH);
 for(d in 1:D)    ll += multinomial_lpmf(HOSP_BVJD[ ,d]  | mH_BVJD[ ,d]) ;
                  ll += multinomial_lpmf(HOSP_BVKJ | mH_BVKJ);

 
 ll += binomial_lpmf(SP_J  | pop_J, pSP);
 ll += binomial_lpmf(VCD   | pop_VCD, pC);
 ll += multinomial_lpmf(VCD_BVJ | mD_BVJ) ;
 ll += multinomial_lpmf(VCD_VK  | mD_VK) ;
 
}
}

model {
  target += ll;
// priors
  hs[1] ~ normal(2.93,0.5);  // SN
  hs[2] ~ normal(20.0,0.5);  // SP
  hl ~ normal(72.14,12); // TAK posterior 
  ts ~ normal(1.34,0.5) ;

  p[1] ~ beta(3,5);

  if(include_pK2==1) {
  pK2 ~ beta(shape1,shape2) ; // prob exposure to serotype k 
  p[2] ~ beta(3,5) ;
  } else {
  pK2 ~ beta(3,5) ;
  p[2] ~ beta(shape1,shape2) ; // prob exposure to a serotype   
  }
  
  lambda_D ~ lognormal(-7,2);
  scale_FOI ~ normal(1,2); 
  theta ~ dirichlet(rep_vector(1.0, K));
  
  gamma ~ normal(0.85,0.2);  // prob secondary
  rho ~ normal(1.97,0.40);   // RR secondary
  delta ~ normal(0.25,0.10);
  phi ~ normal(0.25,0.05);
  epsilon ~ normal(1,1);
  
  L ~ normal(L_mean,L_sd);
  tau ~ normal(1,1);
  w ~ normal(1,2);
  lc[1,] ~ normal(3.85,1);
  for(c in 2:C)  lc[c,] ~ normal(5.98,1);
  omega ~ normal(0,2);
  alpha ~ normal(0,2);
  kappa ~  normal(0,2);
  beta ~ normal(0,2);
  
  sens ~ normal(0.94,0.05) ;
  spec ~ normal(0.74,0.05) ;
}

generated quantities{
  
  // LL
  vector[2*D + 3 + 3] log_lik;
  
   // AR
  array [B,V,J,D] real<lower = 0 > H_AR_BVJD;
  array [B,V,K,J] real<lower = 0 > H_AR_BVKJ;
  array [B,V,J]   real<lower = 0 > V_AR_BVJ ;
  array [V,K]     real<lower = 0 > V_AR_VK ;
  
  // VE
  array [C,K,J,R,T]  real<upper =1 > VE  ;
  array [C,K,R,T]    real<upper =1 > VE_BKRT  ;
  array [C,J,R,T]    real<upper =1 > VE_BJRT  ;
  
// log lik
 for(d in 1:D) log_lik[d] = binomial_lpmf(HOSP_D[d] | pop_HOSP_D[d] , pH[d]);
 for(d in 1:D) log_lik[d+D]  = multinomial_lpmf(HOSP_BVJD[ ,d]  | mH_BVJD[ ,d]) ;
 log_lik[2*D+1] = multinomial_lpmf(HOSP_BVKJ | mH_BVKJ);
 for(j in 1:J) log_lik[j + 2*D + 1] = binomial_lpmf(SP_J[j] | pop_J[j], pSP[j]);
 log_lik[2*D + 3 + 1] = binomial_lpmf(VCD   | pop_VCD, pC);
 log_lik[2*D + 3 + 2] = multinomial_lpmf(VCD_BVJ   | mD_BVJ) ;
 log_lik[2*D + 3 + 3] = multinomial_lpmf(VCD_VK   | mD_VK) ;


// Attack rates

for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   for(d in 1:D)
     H_AR_BVJD[b,v,j,d] = sum(Ho_BVKJD[b,v, ,j,d]) / pop_HOSP_BVJD[b,v,j,d];

for(b in 1:B)
 for(v in 1:V)
  for(k in 1:K)
   for(j in 1:J)
    H_AR_BVKJ[b,v,k,j] = sum(Ho_BVKJD[b,v,k,j, ]) / mean(pop_HOSP_BVJD[b,v,j, ]);
   
for(b in 1:B)
 for(v in 1:V)
  for(j in 1:J)
   V_AR_BVJ[b,v,j] = sum(Sy_BVKJ[b,v, ,j]) / pop_VCD_BVJ[b,v,j]; 
   
 for(v in 1:V)
  for(k in 1:K)
   V_AR_VK[v,k] = sum(Sy_VKJ[v,k, ]) / sum(pop_VCD_VJ[v, ]) ; 


// VE

for(c in 1:C)
 for(k in 1:K)
  for(j in 1:J)
    for(t in 1:T){
VE[c,k,j,1,t] = 1 - RR_symp[c,2,k,j,t] ;
VE[c,k,j,2,t] = 1 - RR_hosp[c,2,k,j,t] ; 

}

for(c in 1:C)
 for(k in 1:K)
  for(r in 1:R)
   for(t in 1:T)
    VE_BKRT[c,k,r,t] = mean(VE[c,k, ,r,t]);

for(c in 1:C)
 for(j in 1:J)
  for(r in 1:R)
   for(t in 1:T)
   VE_BJRT[c,j,r,t] =mean(VE[c, ,j,r,t]);
}
