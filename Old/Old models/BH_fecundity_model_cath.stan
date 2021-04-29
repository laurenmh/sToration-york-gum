// This Stan model is a fecundity model (Ricker version) for Arctotheca calendula
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
  real lambda_0;
  vector[S] alpha_sp;
}

transformed parameters{
  real< lower = 0> lambda;
  vector[N] alpha_mean;
  vector[N] alpha_intra;
  lambda = exp(lambda_0); //exp to restrict alphas to competitive 
}

model{
  // create a vector of predictions, a baseline standard deviation for priors on alpha_sp, 
  //    an object to hold the sums of interactions, and the realized interspecific alphas
  vector[N] F_hat;
  real sigma_sp[S];
  matrix[N,S] interaction_effects;
  matrix[N,S] alpha_inter;

  // set priors
  // alpha_intra ~ normal(0, 1000);
  // alpha_mean ~ normal(0, 1000);
  alpha_sp ~ normal(0, 1000);
  lambda_0 ~ normal(0, 1000);

  // implement the biological model
  for(i in 1:N){ 
    for(s in 1:S){
      alpha_inter[i,s] = exp(alpha_mean[i] + alpha_sp[s]);
      interaction_effects[i,s] = SpMatrix[i,s] * alpha_inter[i,s];
    }
    F_hat[i] = lambda/(1+alpha_intra[i] * Intra[i] + sum(interaction_effects[i,]) + exp(alpha_mean[i]) * Other[i]);
  }
  Fecundity ~ poisson(F_hat);
}








