data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  int Intra[N];
  int Other[N];
  matrix[N,S] SpMatrix;
  //int<lower = 1> Plot;
  real<lower = 1>Phos[N];
  real<lower = 1> Shade[N]; // shade (canopy cover)  per Plot, phos (Colwell P per block)- being treated for now as continuous
}
parameters{
  real<lower = 0> lambda;
  real alpha_intra;
  real alpha_mean;
  vector[S] alpha_sp;
  real Beta_1[S]; // Phos effect
  real Beta_2[S]; // Shade effect
  real Beta_3[S]; // interaction effect
}
model{
  // create a vector of predictions, a baseline standard deviation for priors on alpha_sp, 
  //    an object to hold the sums of interactions, and the realized interspecific alphas
  vector[N] F_hat;
  real sigma_sp[S];
  vector[N] interaction_effects;
  vector[S] alpha_inter;

  // set priors
  alpha_intra ~ normal(0, 1000);
  alpha_mean ~ normal(0, 1000);
  alpha_sp ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);
  Beta_1 ~ normal(0,1000);
  Beta_2 ~ normal(0,1000);
  Beta_3 ~ normal(0,1000);

  // implement the biological model 
  for(i in 1:N){ // just not sure how to index this 
  lambda = Beta_1[1]*Phos[i] + Beta_2[1]*Shade[i] + Beta_3[1]*Phos[i]*Shade[i];
  alpha_intra = Beta_1[2]*Phos[i] + Beta_2[2]*Shade[i] + Beta_3[2]*Phos[i]*Shade[i];
  alpha_inter = alpha_mean + alpha_sp;
  interaction_effects = SpMatrix * alpha_inter;
  }
  for(i in 1:N){
    F_hat[i] = lambda * exp(alpha_intra * Intra[i] + interaction_effects[i] + alpha_mean * Other[i]);
  }
  Fecundity ~ poisson(F_hat);
}


  
  // set priors (super wide for now)
  //Beta_1 ~ normal(0,1000);
  //Beta_2 ~ normal(0,1000);
  //Beta_3 ~ normal(0,1000);
  
//  for(i in 1:N){
//  lambda[i] = Beta_1[1]*Phos[i] + Beta_2[1]*Shade[i] + Beta_3[1]*Phos[i]*Shade[i]; // not sure about this indexing
//   alpha_generic[i] = Beta_1[2]*Phos[i] + Beta_2[2]*Shade[i] + Beta_3[2]*Phos[i]*Shade[i]; 
//   beta[i] = Beta_1[3]*Phos[i] + Beta_2[3]*Shade[i] + Beta_3[3]*Phos[i]*Shade[i]; 
//   }
  
//}
