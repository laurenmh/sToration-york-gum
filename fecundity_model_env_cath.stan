// This Stan model is a fecundity model (Ricker version) for Arctotheca calendula
// seed count (already extrapolated)

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  int Intra[N];
  int Other[N];
  matrix[N,S] SpMatrix;
  real shade;
  real phos;
}
parameters{
  real<lower = 0> lambda;
  real alpha_intra;
  real alpha_mean;
  vector[S] alpha_sp;
  vector[4] b; // for alpha_mean enviro regression parameters 
  vector[4] c; // for alpha_intra enviro regression parameters
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

  // implement the biological model
  for(s in 1:shade){
    for(p in 1:phos){
  alpha_mean = b[1] + b[2]*shade[s] + b[3]*phos[p] + b[4]*shade[s]*phos[p]
  alpha_intra = c[1] + c[2]*shade[s] + c[3]*phos[p] + c[4]*shade[s]*phos[p]
   }
  }
  
  alpha_inter = alpha_mean + alpha_sp;
  interaction_effects = SpMatrix * alpha_inter;
  for(i in 1:N){
    F_hat[i] = lambda * exp(alpha_intra * Intra[i] + interaction_effects[i] + alpha_mean * Other[i]);
  }
  Fecundity ~ poisson(F_hat);
}






