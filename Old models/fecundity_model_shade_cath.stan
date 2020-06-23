// This Stan model is a fecundity model (Ricker version) for Arctotheca calendula
// seed count (already extrapolated)

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Fecundity[N];
  int Intra[N];
  int Other[N];
  matrix[N,S] SpMatrix;
  vector[N] shade; // make these in R from dataset
  vector[N] phos;
}

parameters{
  real <lower = 1> lambda;
  real<lower = 0> sigma;
  vector[S] alpha_sp;
  vector[2] b; // for alpha_mean enviro regression parameters 
  vector[2] c; // for alpha_intra enviro regression parameters
}

transformed parameters{
  vector[N] alpha_mean;
  vector[N] alpha_intra;
  for(i in 1:N){
    alpha_mean[i] = b[1] + b[2]*shade[i];
    alpha_intra[i] = exp(c[1] + c[2]*shade[i]);
  }
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
  lambda ~ gamma(0.001, 0.001);
  //lambda ~ gamma(0.001, 0.001);
  b ~ normal(0, 1000);
  c ~ normal(0, 1000);

  // implement the biological model
  for(i in 1:N){ 
    for(s in 1:S){
      alpha_inter[i,s] = exp(alpha_mean[i] + alpha_sp[s]);
      interaction_effects[i,s] = SpMatrix[i,s] * alpha_inter[i,s];
    }
    F_hat[i] = lambda*exp(alpha_intra[i] * Intra[i] + sum(interaction_effects[i,]) + alpha_mean[i] * Other[i]);
  }
  Fecundity ~ poisson(F_hat);
}








