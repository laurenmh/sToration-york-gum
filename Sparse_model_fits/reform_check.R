# check that we get the same answers using the new model formulation

library(rstan)
library(here)
library(tidyverse)
library(rlang)

# load the model fits
  trcy_mat <- readRDS(here("Sparse_model_fits/TRCY_Phos_bhfh_matform.rds"))
  trcy_og <- readRDS(here("Sparse_model_fits/TRCY_Phos_bhfh_og.rds"))
  
# check that parameter posteriors are similar
  checks <- tibble(
    mfit = rep(c("original", "reform"), each = 2),
    site = rep(c("Bendering", "Perenjori"), 2)
  )
  
# extract posteriors
  lambdas_og <- as.data.frame(rstan::extract(trcy_og, pars = "lambdas"))
  lambdas_rf <- (rstan::extract(trcy_mat, pars = "beta_lambda_std"))[[1]]
  alphas_og <- as.data.frame(rstan::extract(trcy_og, pars = "alphas"))
  alphas_rf <- rstan::extract(trcy_mat, pars = "beta_alpha")[[1]]

# add posteriors
  checks <- checks %>% mutate(
    lambda_i0 = list(
      lambdas_og$lambdas.1.1,
      lambdas_og$lambdas.2.1,
      lambdas_rf%*%c(1,0,0,0),
      lambdas_rf%*%c(1,1,0,0)
    ),
    lambda_iphos = list(
      lambdas_og$lambdas.1.2,
      lambdas_og$lambdas.2.2,
      lambdas_rf%*%c(0,0,1,0),
      lambdas_rf%*%c(0,0,1,1)
    ),
    alpha_0gen = list(
      alphas_og$alphas.1,
      alphas_og$alphas.1,
      alphas_rf[,1],
      alphas_rf[,1]
    ),
    alpha_1gen = list(
      alphas_og$alphas.2,
      alphas_og$alphas.2,
      alphas_rf[,2],
      alphas_rf[,2]      
    )
  )

  checks_long <- unnest(checks, cols = c(
    "lambda_i0", "lambda_iphos",
    "alpha_0gen", "alpha_1gen"
  ))

  check_posts <- function(df, c_name, vary_site = F) {
    if(isFALSE(vary_site)){
      ggplot(data = df, aes(x = {{c_name}}, linetype = mfit)) +
        geom_density(alpha=0.8) +
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

par(mfrow = c(2,2))
  check_posts(checks_long, lambda_i0, vary_site = T)  
  check_posts(checks_long, lambda_iphos, vary_site = T)  
  check_posts(checks_long, alpha_0gen)
  check_posts(checks_long, alpha_1gen)
dev.off()
  
  
  