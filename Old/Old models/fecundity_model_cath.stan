// This Stan model is a fecundity model (Ricker version) 
// seed count (already extrapolated)

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  int Intra[N];
  int Other[N];
  matrix[N,S] SpMatrix;
}
parameters{
  real<lower = 1> lambda;
  real <lower = 0> alpha_intra;
  real <lower = 0> alpha_mean;
  vector[S] alpha_sp;
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
  alpha_inter = alpha_mean + alpha_sp;
  interaction_effects = SpMatrix * alpha_inter;
  for(i in 1:N){
    F_hat[i] = lambda*exp(alpha_intra * Intra[i] + interaction_effects[i] + alpha_mean * Other[i]);
  }
  Fecundity ~ poisson(F_hat);
}






