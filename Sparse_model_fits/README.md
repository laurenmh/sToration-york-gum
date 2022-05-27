This folder contains code for and results from fitting Beverton-Holt models using Bayesian sparse methods to constrain some parameters to a global mean.

## Notes

As of May 27, 2022, models do not appear to capture all the variability and structure for TRCY based on posterior predictive checks. The model overestimates fecundity at some sites in which there are no competitors (even after accounting for variation among plots and blocks), indicating unmeasured variables also influence fecundity.

## Stan folder

Stan programs for different variants of the Beverton-Holt model. For full descriptions of the model encoded in the files `BH_FH_Preliminary.stan` and `BH_Final.stan`, see [Weiss-Lehman et al. (2022)](https://onlinelibrary.wiley.com/doi/10.1111/ele.13977). For descriptions of the models encoded in `BH_FH_Preliminary2_matrix_form.stan` and `BH_FH_Preliminary2_re_plot_block.stan` (and their `BH_Final` counterparts) see [alt_model_form_ideas.pdf](https://github.com/laurenmh/sToration-york-gum/blob/master/Sparse_model_fits/alt_model_form_ideas.pdf).

## .RData files

These files contain samples from the posterior distributions of parameters as well as the data structures used to fit the models. 

Naming conventions: "focal species code"_"model"_FinalFit.rdata

**Focal species codes**:
  * *Arctotheca calendula* (ARCA)
  * *Hypochaeris glabra* (HYGL)
  * *Trachymene cyanopetala* (TRCY)
  * *Waitzia acuminata* (WAAC)
  
**Model structure codes**:
  * Phos -- Log-linear model on low-density growth rate (LDGR) with phosphorous as a covariate. Log-linear model on competition coefficients (alpha) with phosphorous as a covariate. Some competing species are allowed to deviate from global average linear model [(Weiss-Lehman et al. 2022)](https://onlinelibrary.wiley.com/doi/10.1111/ele.13977). 
  * Shade -- Log-linear model on LDGR with canopy cover (shade) as a covariate. Log-linear model on competition coefficients (alpha) with canopy cover as a covariate. Some competing species are allowed to deviate from global average linear model [(Weiss-Lehman et al. 2022)](https://onlinelibrary.wiley.com/doi/10.1111/ele.13977).
  * `{PSI}-lambda_{PSI}-alpha` -- Log-linear model on LDGR (lambda) with phosphorous (P), shade (S), and/or the interaction between them (I) as covariates. Log-linear model on competition coefficients (alpha) with phosphorous (P), shade (S), and/or the interaction between them (I) as covariates. Some competing species are allowed to deviate from global average linear model.
  
## Functions (`functions.R`)

* `construct_M(n_env_covs, int = F)`: Construct a linear map from the linear model parameters to the deviations from the global mean effects (alpha_hat_eij) for each site (Bendering, Perenjori) based on 0, 1, or 2 environmental covariates, with and without the interaction between two environmental covariates.

* `construct_alpha_hat(M, B_post, ncovs, covnames, sitenames=c("bendering", "perenjori"))`: Construct array of draws from the posterior distribution of the alpha_hat vector for each competitor species based on the mapping `M` from B -> alpha_hat. 

* `non_generic(alpha_hat_draws, level = 0.5, sp_names=NA, all_devs = F)`: Determine which species' effects differ substantially from the overall mean effect based on which alpha_hat parameter posterior distributions pull away from zero.

* `ppcs(y_rep, obs, plot = T, pred_level = 0.95)`: Posterior predictive checks of the model fit.

* `plot_heatmap(m_formula, param_estims, sp_name, ncells = 100, dem_param = "lambda", ...)`: Create a heatmap of LDGR or another demographic parameter with phosphorous on one axis and canopy on the other.
  
  
  
  