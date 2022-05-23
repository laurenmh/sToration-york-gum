This folder contains code for and results from fitting Beverton-Holt models using Bayesian sparse methods to constrain some parameters to a global mean.

## Notes

As of May 23, 2022, models do not appear to capture all the variability and structure for TRCY and WAAC based on posterior predictive checks. Including effects to model the nested structure of the experiment may help, but no stan models are yet written to do so.

## Stan folder

Stan programs for different variants of the Beverton-Holt model. For full descriptions of the model encoded in the files `BH_FH_Preliminary.stan` and `BH_Final.stan`, see [Weiss-Lehman et al. (2022)](https://onlinelibrary.wiley.com/doi/10.1111/ele.13977). For descriptions of the models encoded in `BH_FH_Preliminary2_matrix_form.stan`, see [alt_model_form_ideas.pdf](https://github.com/laurenmh/sToration-york-gum/blob/combine-env-covs/Sparse_model_fits/alt_model_form_ideas.pdf).

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
  
## Functions

* `construct_M(n_env_covs, int = F)`: 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  