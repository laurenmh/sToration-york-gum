// This script fits a Beverton-Holt generalized competition model after using a 
//   Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) to shrink species-specific deviations
//   in competitive effects to the average. This refits the model constraining the species-specific
//   parameters that did not pull away from zero in the original fit to exactly zero.

data{
  int<lower = 1> N;                // Number of plots
  int<lower = 1> S;                // Number of species (not including the focal)
  int<lower = 1> P_alpha;          // Number of covariates determining competition coefficients (including the intercept)
  int<lower = 1> P_alpha_d;        // Number of variables that vary with species id (varying intercepts, slopes,...)
  int<lower = 1> P_eta;            // Number of covariates determining strength of intra-specific competition
  int<lower = 1> P_lambda;         // Number of covariates determining LDGR
  int Fecundity[N];                // Fecundity of the focal species in each plot
  matrix[N, S] A;                  // Matrix of abundances for each species (not including non-focal conspecifics)
  vector[N] C;                     // vector of conspecific abundances in each plot
  matrix[N, P_alpha] X_alpha;      // model matrix for fixed effects
  matrix[N, P_alpha_d] Z_alpha;    // model matrix for species-level deviations from the mean
  matrix[P_alpha_d, P_alpha_d] M;  // mapping from b parameters to alpha_hat
  matrix[P_alpha_d, S] Q;          // matrix of 0's and 1's to constrain some alpha_hats to zero
  matrix[N, P_eta] X_eta;          // model matrix for modeling variability in the stength of intra-specific competition
  matrix[N, P_lambda] X_lambda;    // model matrix for modeling variability in the LDGR
}

transformed data{
  row_vector[S] ones = rep_row_vector(1, S);
  matrix[P_alpha_d, P_alpha_d] M_inv = inverse(M);
}

parameters{
  vector[P_lambda] beta_lambda_std;                // coefficients to model variability in LDGR
  vector[P_eta] beta_eta_std;                      // coefficients to model variability in strength of intra-specific competition
  vector[P_alpha] beta_alpha_std;                  // coefficients to model variability in strength of average inter-specific competition
  matrix[P_alpha_d, S] alpha_hat_std;              // species-specifc deviations from the average, concatenated into a matrix
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the standardized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  
  vector[P_alpha] beta_alpha;                     // scaled and centered version of beta_alpha_std
  vector[P_eta] beta_eta;                         // scaled and centered version of beta_eta_std
  vector[P_lambda] beta_lambda;                   // scaled and centered version of beta_lambda_std
  vector[N] mu;                                   // means based on BH model
  
  vector[N] eta;                                  // vector of intraspecific competition coefficients
  vector[N] lambda;                               // vector of LDGRs
  matrix[N, S] alpha_mat;                         // NxS matrix of cometition coefficents
    
  // apply constraints to deviations from the mean effects and map back to regression parameterization
  matrix[P_alpha_d, S] B = M_inv * (alpha_hat_std .* Q);
  
  // scale and center the regression coefficients
  beta_alpha = beta_alpha_std * 0.5;
  beta_alpha[1] = 3 * beta_alpha_std[1] - 6;
  beta_eta = beta_eta_std * 0.5;
  beta_eta[1] = 3 * beta_eta_std[1] - 6;
  beta_lambda = beta_lambda_std;
  
  // transform from log-linear models
  alpha_mat = exp((X_alpha * beta_alpha) * ones + Z_alpha * B);
  eta = exp(X_eta * beta_eta);
  lambda = exp(X_lambda * beta_lambda);
  
  // construct means based on BH model
  for(i in 1:N){
    mu[i] = lambda[i] / (1 + eta[i] * C[i] + alpha_mat[i,] * A[i, ]');
  }

}

model{
  
  // Priors
  beta_lambda_std ~ std_normal();
  beta_eta_std ~ std_normal();
  beta_alpha_std ~ std_normal();
  for(p in 1:P_alpha_d){
    alpha_hat_std[p, ] ~ std_normal();
  }

  // Likelihood
  Fecundity ~ poisson(mu);
  
}

generated quantities{
  
  int y_rep[N];
  
  // generate samples for posterior predictive checks
  y_rep = poisson_rng(mu);
  
}