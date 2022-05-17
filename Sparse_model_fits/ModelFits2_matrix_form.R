# This script will run the empirical model fits for each focal species with
#   any number of environmental covariates. A description of the model formulation
#   can be found in the file "alt_model_form_ideas.pdf". 
#   A separate script will then make the empirical figures for the manuscript.

  rm(list = ls())
  library(rstan)
  library(here)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
  
# load helper functions
  source(here("Sparse_model_fits/functions.R"))

  FocalLetter <- "T" # "W", "A", "T", "H"
  FocalPrefix <- "TRCY" # "WAAC", "ARCA", "HYGL", "TRCY"
  FocalSpecies <- "Trachymene.cyanopetala" # "Waitzia.acuminata", "Arctotheca.calendula", "Trachymene.cyanopetala", "Hypochaeris.glabra"
  mod <- "PSI-lambda_PS_alpha" # Phosphorous + shade + interaction for lambda, phosphorous + shade for alpha
  
# Load in the data and subset out the current focal species.
  SpData <- read.csv("water_full_env.csv")
  SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
  SpData <- na.omit(SpData) 
  FocalLetter
  SpDataFocal <- subset(SpData, Focal.sp.x == FocalLetter)

# Now calculate the total number of (non-focal) species to use for the model, discounting
  #   any species columns with 0 abundance. Save a vector of the species names
  #   corresponding to each column for easy matching later.
  AllSpAbunds <- SpDataFocal[,10:69]
  AllSpNames <- names(AllSpAbunds)
  SpTotals <- colSums(AllSpAbunds)
  SpToKeep <- AllSpNames[SpTotals > 0 & AllSpNames != FocalSpecies]
  
# create matrix A of heterospecific abundances
  A <- as.matrix(
    AllSpAbunds[, which(names(AllSpAbunds) %in% SpToKeep)]
  )
  
# create vector of conspecific abundances
  Consp_A <- AllSpAbunds[, which(names(AllSpAbunds) == FocalSpecies)]
  
# create model matrices, X_alpha, X_lambda, X_eta (see "alt_model_form_ideas.pdf")
  env_covs <- c("Reserve.x", "Colwell.P", "Canopy")
  env_data <- SpDataFocal[, which(names(SpDataFocal) %in% env_covs)]
  env_data$Reserve.x <- as.factor(env_data$Reserve.x)
# standardize continuous variables
  env_data$Colwell.P_std <- scale(env_data$Colwell.P)
  env_data$Canopy_std <- scale(env_data$Canopy)
  
# model matrix for "generic" effects
  X_alpha <- model.matrix(
    ~ Colwell.P_std + Canopy_std, 
    data = env_data
  )
  
# model matrices for strength of intraspecific competition (X_eta) and LDGR (X_lambda)
  # can differ from X_alpha if it makes sense
   X_lambda <- model.matrix(
     ~ Reserve.x * Colwell.P_std * Canopy_std, 
     data = env_data
   )
   X_eta <- X_alpha
  
# create model matrix for species-level deviations from the mean
  Z_alpha <- model.matrix(
    ~ Reserve.x*Colwell.P_std + Reserve.x*Canopy_std, 
    data = env_data
  )
  
# define mapping from parameter vectors b_1, b_2, ..., b_S to alpha_hat_1, alpha_hat_2,...
#  (see "alt_model_form_ideas.pdf" and "functions.R" for more details)
  M <- construct_M_dc(ncovs = 3)

# define hyperpriors for Finnish Horseshoe
  tau0 <- 1
  slab_scale <- sqrt(2)
  slab_df <- 4
  Fecundity <- as.integer(SpDataFocal$Number.flowers.total)
  
# compile data into list for model fitting
  data_prelim <- list(
    N = nrow(SpDataFocal),
    S = ncol(A),
    P_alpha = ncol(X_alpha),
    P_alpha_d = ncol(Z_alpha),
    P_eta = ncol(X_eta),
    P_lambda = ncol(X_lambda),
    Fecundity = Fecundity,
    A = A,
    C = Consp_A,
    X_alpha = X_alpha,
    Z_alpha = Z_alpha,
    M = M,
    X_lambda = X_lambda,
    X_eta = X_eta,
    tau0 = tau0,
    slab_scale = slab_scale,
    slab_df = slab_df
  )
  
# compile stan model
  bhfh_matform <- stan_model(
    file = here("Sparse_model_fits/Stan/BH_FH_Preliminary2_matrix_form.stan"),
    model_name = "bhfh_mf"
  )
  
# fit the model
  mfit <- sampling(
    bhfh_matform, data = data_prelim, 
    chains=3, iter = 3000, cores=3
    # control=list(adapt_delta = 0.9, max_treedepth=15)
  )

  
# diagnostics
  check_divergences(mfit)
  check_treedepth(mfit)
  check_energy(mfit)

# if there are no issues with model convergence, proceed to constraining species-specific
#   deviations from the means and refitting the model
  
# construct an array of posterior draws from alpha_hat matrix
  alpha_hat_post <- construct_alpha_hat(
    M,
    B_post = rstan::extract(mfit, pars="B")[[1]],
    covnames = c("Int", "Phos")
  )
  
# find which alpha_hat posteriors pulled away sufficiently from zero
#  function documentation in functions.R file
  Q <- non_generic(alpha_hat_post)
  
# refit the model with the new constraint matrix Q and std_normal() priors on remaining alpha_hats
  bh_matform <- stan_model(
    file = here("Sparse_model_fits/Stan/BH_Final2_matrix_form.stan"),
    model_name = "bh_mf"
  )
  
# data list for final fit
  data_final <- list(
    N = nrow(SpDataFocal),
    S = ncol(A),
    P_alpha = ncol(X_alpha),
    P_alpha_d = ncol(Z_alpha),
    P_eta = ncol(X_eta),
    P_lambda = ncol(X_lambda),
    Fecundity = Fecundity,
    A = A,
    C = Consp_A,
    X_alpha = X_alpha,
    Z_alpha = Z_alpha,
    M = M,
    Q = Q,
    X_lambda = X_lambda,
    X_eta = X_eta
  )
  
# refit the model  
  mfit_final <- sampling(
    bh_matform, data = data_final, 
    chains=3, iter = 3000, cores=3
    # control=list(adapt_delta = 0.9, max_treedepth=15)
  )
  
# save results to file
  FileName <- paste(here("Sparse_model_fits/"), FocalPrefix, "_", mod, "_FinalFit.rdata", sep = "")
  save(mfit_final, data_prelim, data_final, file = FileName)
  
  
  
  
    
  
  