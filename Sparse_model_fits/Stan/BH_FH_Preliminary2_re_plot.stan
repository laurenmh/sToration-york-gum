// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
//  This script futher includes structure to account for the nested experimental design
//  

data{
  int<lower = 1> N;                // Number of observations
  int<lower = 0> K_plot;           // Number of plots
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
  matrix[N, P_eta] X_eta;          // model matrix for modeling variability in the stength of intra-specific competition
  matrix[N, P_lambda] X_lambda;    // model matrix for modeling variability in the LDGR
  
  // Below defines the nested structure of the experiment
  int plot_id[N];                  // index for the plot
  
  // The below values define the regularized horseshoe priors used for species-specific parameters
  real tau0; 		                   // determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	               // scale for alpha_sp values that deviate from the mean
  real slab_df;		                 // effective degrees of freedom for significant alpha_sp values
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
  row_vector[S] ones = rep_row_vector(1, S);
  matrix[P_alpha_d, P_alpha_d] M_inv = inverse(M);
}

parameters{
  vector[P_lambda] beta_lambda_std;                // coefficients to model variability in LDGR
  vector[P_eta] beta_eta_std;                      // coefficients to model variability in strength of intra-specific competition
  vector[P_alpha] beta_alpha_std;                  // coefficients to model variability in strength of average inter-specific competition
  matrix[P_alpha_d, S] alpha_hat_std;              // species-specifc deviations from the average, concatenated into a matrix
  
  // parameters involved in experimental design
  vector[K_plot] u_plot_std;                       // unscaled varying means for the plot
  real<lower = 0> sigma_plot;                      // scale of variation among plots
  
  // parameters involved in shrinkage priors
  matrix<lower = 0>[P_alpha_d, S] local_scale;     // non-regularized local scale
  real<lower = 0> c2_std;                          // unscaled version of c2
  real<lower = 0> tau_std;                         // unscaled version of tau
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the standardized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  
  vector[P_alpha] beta_alpha;                     // scaled and centered version of beta_alpha_std
  vector[P_eta] beta_eta;                         // scaled and centered version of beta_eta_std
  // vector[P_lambda] beta_lambda;                  // scaled and centered version of beta_lambda_std
  
  vector[N] eta;                                  // vector of intraspecific competition coefficients
  vector[N] lambda;                               // vector of LDGRs
  matrix[N, S] alpha_mat;                         // NxS matrix of cometition coefficents
  
  // scale tau and c2
  real c2 = slab_scale2*c2_std;                   // c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)
  real tau = tau0*tau_std;                        // tau ~ cauchy(0, tau0)
  
  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  matrix[P_alpha_d, S] local_scale_tilde = 
    sqrt(c2 * square(local_scale) ./ (c2 + square(tau) * square(local_scale)));
    
  // scale the deviations from the mean effects and map back to regression parameterization
  matrix[P_alpha_d, S] B = M_inv * (tau * local_scale_tilde .* alpha_hat_std);
  
  // scale the effects of the plot and block
  vector[K_plot] u_plot = u_plot_std * sigma_plot;
  
  // scale and center the regression coefficients
  beta_alpha = beta_alpha_std * 0.5;
  beta_alpha[1] = 3 * beta_alpha_std[1] - 6;
  beta_eta = beta_eta_std * 0.5;
  beta_eta[1] = 3 * beta_eta_std[1] - 6;
  
  // transform from log-linear models
  alpha_mat = exp((X_alpha * beta_alpha) * ones + Z_alpha * B);
  eta = exp(X_eta * beta_eta);
  lambda = exp(X_lambda * beta_lambda_std);
}

model{
  // Declare vector of expected fecundity values (mu),
  vector[N] mu;

  // Priors
  beta_lambda_std ~ std_normal();
  beta_eta_std ~ std_normal();
  beta_alpha_std ~ std_normal();
  u_plot_std ~ std_normal();
  sigma_plot ~ student_t(2, 0, 1);
  
  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  for(p in 1:P_alpha_d){
    alpha_hat_std[p, ] ~ std_normal();
    local_scale[p, ] ~ cauchy(0, 1);
  }
  
  tau_std ~ cauchy(0, 1);
  c2_std ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  for(i in 1:N){
    mu[i] = (lambda[i] / (1 + eta[i] * C[i] + alpha_mat[i,] * A[i, ]')) * exp(u_plot[plot_id[i]]);
  }
  
  Fecundity ~ poisson(mu);
}
