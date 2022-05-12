# This script will run the empirical model fits for each focal species with
#   any number of environmental covariates. A description of the model formulation
#   can be found in the file "alt_model_form_ideas.pdf". 
#   A separate script will then make the empirical figures for the manuscript.

  rm(list = ls())
  library(rstan)
#  library(here)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)

  FocalLetter <- "T" # "W", "A", "T", "H"
  FocalPrefix <- "TRCY" # "WAAC", "ARCA", "HYGL", "TRCY"
  FocalSpecies <- "Trachymene.cyanopetala" # "Waitzia.acuminata", "Arctotheca.calendula", "Trachymene.cyanopetala", "Hypochaeris.glabra"
  
# Load in the data and subset out the current focal species.
  SpData <- read.csv("water_full_env.csv")
  SpData <- subset(SpData, select = -c(X.NA., Seedcount.extrapolated.integer))
  SpData <- na.omit(SpData) 
  FocalLetter
  SpDataFocal <- subset(SpData, Focal.sp.x == FocalLetter)

# Now calculate the total number of species to use for the model, discounting
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
    ~ Colwell.P_std, 
    data = env_data
  )
  
# model matrices for strength of intraspecific competition and LDGR
  # can differ from X_alpha if it makes sense
   X_lambda <- model.matrix(
     ~ Reserve.x*Colwell.P_std, 
     data = env_data
   )
   X_eta <- X_alpha
  
# create model matrix for species-level deviations from the mean
  Z_alpha <- model.matrix(
    ~ Reserve.x*Colwell.P_std - 1, 
    data = env_data
  )

# define hyperpriors for Finnish Horseshoe
  tau0 <- 1
  slab_scale <- sqrt(2)
  slab_df <- 4
  
# compile data into list for model fitting
  data_list <- list(
    N = nrow(SpDataFocal),
    S = ncol(A),
    P_alpha = ncol(X_alpha),
    P_alpha_d = ncol(Z_alpha),
    P_eta = ncol(X_eta),
    P_lambda = ncol(X_lambda),
    Fecundity = as.integer(SpDataFocal$Number.flowers.total),
    A = A,
    C = Consp_A,
    X_alpha = X_alpha,
    Z_alpha = Z_alpha,
    X_lambda = X_lambda,
    X_eta = X_eta,
    tau0 = tau0,
    slab_scale = slab_scale,
    slab_df = slab_df
  )
  
# compile stan model
  bhfh_matform <- stan_model(
    file = here("Sparse_model_fits/BH_FH_Preliminary2_matrix_form.stan"),
    model_name = "bhfh_mf"
  )
  
# fit the model
  mfit <- sampling(
    bhfh_matform, data = data_list, 
    chains=3, iter = 3000, cores=3
  )
  
# save model fit
  saveRDS(mfit, file = here("Sparse_model_fits/TRCY_Phos_bhfh_matform_rep2.rds"))
  


  
  
  