# check that we get the same answers using the new model formulation

library(rstan)
library(here)
library(tidyverse)
library(rlang)

# load the model fits
  trcy_mat <- readRDS(here("Sparse_model_fits/TRCY_Phos_peren_bhfh_matform.rds"))
  trcy_og <- readRDS(here("Sparse_model_fits/TRCY_phos_peren_og.rds"))
                      
# check that parameter posteriors are similar
  checks <- tibble(
    mfit = c("original", "reform")
  )
  
# extract posteriors
  lambdas_og <- as.data.frame(rstan::extract(trcy_og, pars = "lambdas"))
  lambdas_rf <- (rstan::extract(trcy_mat, pars = "beta_lambda_std"))[[1]]
  alphas_og <- as.data.frame(rstan::extract(trcy_og, pars = "alpha_generic"))
  alphas_rf <- rstan::extract(trcy_mat, pars = "beta_alpha")[[1]]
  intra_og <- as.data.frame(rstan::extract(trcy_og, pars = "alpha_intra"))
  intra_rf <- as.data.frame(rstan::extract(trcy_mat, pars = "beta_eta"))

# add posteriors
  checks <- checks %>% mutate(
    lambda_i0 = list(
      lambdas_og$lambdas.1,
      lambdas_rf[,1]
    ),
    lambda_iphos = list(
      lambdas_og$lambdas.2,
      lambdas_rf[,2]
    ),
    alpha_generic0 = list(
      alphas_og$alpha_generic.1,
      alphas_rf[,1]
    ),
    alpha_generic_phos = list(
      alphas_og$alpha_generic.2,
      alphas_rf[,2]
    ),
    alpha_intra0 = list(
      intra_og$alpha_intra.1,
      intra_rf$beta_eta.1
    ),
    alpha_intra_phos = list(
      intra_og$alpha_intra.2,
      intra_rf$beta_eta.2
    )
  )

  checks_long <- unnest(checks, cols = c(
    "lambda_i0", "lambda_iphos",
    "alpha_generic0", "alpha_generic_phos",
    "alpha_intra0", "alpha_intra_phos"
  ))

  check_posts <- function(df, c_name, vary_site = F) {
    if(isFALSE(vary_site)){
      ggplot(data = df, aes(x = {{c_name}}, linetype = mfit)) +
        geom_density() +
        theme_bw() +
        scale_linetype_manual(values = c(1,2))
    } else {
      ggplot(data = df, aes(x = {{c_name}}, color = site, linetype = mfit)) +
        geom_density() +
        theme_bw() +
        scale_color_manual(values = c("black", "grey")) +
        scale_linetype_manual(values = c(1,2))
    }
  }


gridExtra::grid.arrange(
  check_posts(checks_long, lambda_i0, vary_site = F),  
  check_posts(checks_long, lambda_iphos, vary_site = F),
  check_posts(checks_long, alpha_generic0),
  check_posts(checks_long, alpha_generic_phos),
  check_posts(checks_long, alpha_intra0),
  check_posts(checks_long, alpha_intra_phos),
  ncol = 2
)


# reconstruct lambdas
  X_new <- rbind(
    c(1,0,0,0),
    c(0,0,1,0),
    c(1,1,0,0),
    c(0,0,1,1)
  )
  lambda_mat <- matrix(nrow = dim(lambdas_rf)[1], ncol = 4)
  for(i in 1:nrow(lambda_mat)){
    lambda_mat[i,] <- t(X_new%*%lambdas_rf[i,])
  }
  
# what about numbers of non-generic interactions?
  B_post <- rstan::extract(trcy_mat, pars= "B")[[1]]
  
# reconstruct intercept and slope deviations
  alpha_hat_ij_rf <- array(dim = c(4500, 2, 37))
  alpha_hat_eij_rf <- array(dim = c(4500,2, 37))
  alpha_hat_ij_rf[,1,] = B_post[,1,]
  alpha_hat_ij_rf[,2,] = B_post[,1,] + B_post[,2,]
  alpha_hat_eij_rf[,1,] = B_post[,3,]
  alpha_hat_eij_rf[,2,] = B_post[,3,] + B_post[,4,]
  
# which pulled away sufficiently from the overall mean?
  non_g_ij_rf <- matrix(data = 0, nrow = 2, ncol = dim(B_post)[3])
  non_g_eij_rf <- matrix(data = 0, nrow = 2, ncol = dim(B_post)[3])
  
  IntLevel <- 0.5
  for(i in 1:2){
    for(s in 1:ncol(non_g_ij_rf)){
      Ints_ij <- HDInterval::hdi(alpha_hat_ij_rf[,i,s], credMass = IntLevel)
      Ints_eij <- HDInterval::hdi(alpha_hat_eij_rf[,i,s], credMass = IntLevel)
      if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
        non_g_ij_rf[i,s] <- 1
      }
      if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
        non_g_eij_rf[i,s] <- 1
      }
    }
  }
  
# repeat with original model fit
  alpha_hats_post <- rstan::extract(trcy_og, pars = c("alpha_hat_ij", "alpha_hat_eij"))
  non_g_ij <- matrix(data = 0, nrow = 2, ncol = dim(B_post)[3])
  non_g_eij <- matrix(data = 0, nrow = 2, ncol = dim(B_post)[3])
  
  for(i in 1:2){
    for(s in 1:ncol(non_g_ij)){
      Ints_ij <- HDInterval::hdi(alpha_hats_post$alpha_hat_ij[,i,s], credMass = IntLevel)
      Ints_eij <- HDInterval::hdi(alpha_hats_post$alpha_hat_eij[,i,s], credMass = IntLevel)
      if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
        non_g_ij[i,s] <- 1
      }
      if(Ints_eij[1] > 0 | Ints_eij[2] < 0){
        non_g_eij[i,s] <- 1
      }
    }
  }
  
  
# differing intercepts bendering 
  which(non_g_ij[1,] == 1)
  which(non_g_ij_rf[1,] == 1)
  
# differing intercepts perenjori
  which(non_g_ij[2,] == 1)
  which(non_g_ij_rf[2,] == 1)
  
# differing slopes bendering 
  which(non_g_eij[1,] == 1)
  which(non_g_eij_rf[1,] == 1)
  
# differing slopes perenjori
  which(non_g_eij[2,] == 1)
  which(non_g_eij_rf[2,] == 1)
  
  
  
  
# check against model fit with just perenjori data
  trcy_pmat <- readRDS(here("Sparse_model_fits/TRCY_Phos_peren_bhfh_matform.rds"))
  
  checks_p <- tibble(
    mfit = c("original", "reform"),
  )

# extract posteriors
  lambdas_rfp <- (rstan::extract(trcy_pmat, pars = "beta_lambda_std"))[[1]]
  
# add posteriors
  checks_p <- checks_p %>% mutate(
    lambda_i0 = list(
      lambdas_og$lambdas.2.1,
      lambdas_rfp[,1]
    ),
    lambda_iphos = list(
      lambdas_og$lambdas.2.2,
      lambdas_rfp[,2]
    )
  )
  
# unnest
  checks_plong <- unnest(checks_p, cols = c(
    "lambda_i0", "lambda_iphos"
  ))
  
  check_posts(checks_plong, lambda_i0)
  
  ggplot(data = checks_plong, aes(x = lambda_i0, linetype = mfit)) +
    geom_density() +
    theme_bw() +
    scale_linetype_manual(values = c(1,2))
  