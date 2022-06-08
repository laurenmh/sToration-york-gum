This folder contains code for and results from fitting Beverton-Holt models using Bayesian sparse methods to constrain some parameters to a global mean.

## Notes

As of May 27, 2022, models do not appear to capture all the variability and structure for TRCY based on posterior predictive checks (although the fit seems reasonably good). The model overestimates fecundity at some sites in which there are no competitors (even after accounting for variation among plots and blocks), indicating unmeasured variables also influence fecundity.

## Stan folder

Stan programs for different variants of the Beverton-Holt model. For full descriptions of the model encoded in the files `BH_FH_Preliminary.stan` and `BH_Final.stan`, see [Weiss-Lehman et al. (2022)](https://onlinelibrary.wiley.com/doi/10.1111/ele.13977). For descriptions of the models encoded in `BH_FH_Preliminary2_matrix_form.stan` and `BH_FH_Preliminary2_re_plot_block.stan` (and their `BH_Final` counterparts) see [alt_model_form_ideas.pdf](https://github.com/laurenmh/sToration-york-gum/blob/master/Sparse_model_fits/alt_model_form_ideas.pdf).

## .RData files

These files contain samples from the posterior distributions of parameters as well as the data structures used to fit the models.

### Naming conventions

Each .RData file follows the convention: "focal species code"\_"model"\_FinalFit.rdata

**Focal species codes**:

-   *Arctotheca calendula* (ARCA)
-   *Hypochaeris glabra* (HYGL)
-   *Trachymene cyanopetala* (TRCY)
-   *Waitzia acuminata* (WAAC)

**Model structure codes**:

-   Phos -- Log-linear model on low-density growth rate (LDGR) with phosphorous as a covariate. Log-linear model on competition coefficients (alpha) with phosphorous as a covariate. Some competing species are allowed to deviate from global average linear model [(Weiss-Lehman et al. 2022)](https://onlinelibrary.wiley.com/doi/10.1111/ele.13977).
-   Shade -- Log-linear model on LDGR with canopy cover (shade) as a covariate. Log-linear model on competition coefficients (alpha) with canopy cover as a covariate. Some competing species are allowed to deviate from global average linear model [(Weiss-Lehman et al. 2022)](https://onlinelibrary.wiley.com/doi/10.1111/ele.13977).
-   `{PSI}-lambda_{PSI}-alpha` -- Log-linear model on LDGR (lambda) with phosphorous (P), shade (S), and/or the interaction between them (I) as covariates. Log-linear model on competition coefficients (alpha) with phosphorous (P), shade (S), and/or the interaction between them (I) as covariates. Some competing species are allowed to deviate from global average linear model.

### Objects included in each .RData file

-   `mfit_final`: Stan object from the final model fit (after using regularizing priors to decide which species-level deviations from global average effects should be set to zero). Parameters of interest include:
    -   ${\bf \beta}_\alpha$ (`beta_alpha`): Regression coefficients for modeling relationships between environmental covariates and the competitive effect of a generic neighbor.
    -   $\bf{\beta}_\lambda$ (`beta_lambda_std`): Regression coefficients for modeling relationships between environmental covariates and LDGR (received standard normal priors while all other regression coefficients had non-standard scales, hence the `_std`.)
    -   $\bf{\beta}_\eta$ (`beta_eta`): Regression coefficients for modeling relationships between environmental covariates and intra-specific competitive effects.
    -   ${\bf B}$ (`B`): $E \times S$ matrix of species-specific deviations from the average effects of $E$ environmental covariates for each of $S$ species.
-   `data_final`: List of data objects supplied to the Stan model
    -   `int<lower = 1> N`: Number of experimental units
    -   `int<lower = 1> K_plot`: Number of plots
    -   `int<lower = 1> K_block`: Number of blocks
    -   `int<lower = 1> S`: Number of species (not including the focal)
    -   `int<lower = 1> P_alpha`: Number of covariates determining competition coefficients (including the intercept)
    -   `int<lower = 1> P_alpha_d`: Number of variables that vary with species id (varying intercepts, slopes,...)
    -   `int<lower = 1> P_eta`: Number of covariates determining strength of intra-specific competition
    -   `int<lower = 1> P_lambda`: Number of covariates determining LDGR
    -   `int Fecundity[N]`: Fecundity of the focal species in each subplot
    -   `matrix[N, S] A`: Matrix of abundances for each species (not including non-focal conspecifics)
    -   `vector[N] C`: Vector of conspecific abundances in each subplot
    -   `matrix[N, P_alpha] X_alpha`: Model matrix for fixed effects
    -   `matrix[N, P_alpha_d] Z_alpha`: Model matrix for species-level deviations from the mean
    -   `matrix[P_alpha_d, P_alpha_d] M`: Mapping from b parameters to alpha_hat
    -   `matrix[P_alpha_d, S] Q`: Matrix of 0's and 1's to constrain some alpha_hats to zero
    -   `matrix[N, P_eta] X_eta`: Model matrix for modeling variability in the stength of intra-specific competition
    -   `matrix[N, P_lambda] X_lambda`: Model matrix for modeling variability in the LDGR
    -   `int plot_id[N]`: Index for the plot
    -   `int block_id[N]`: Index for the block
-   `data_prelim`: Same as `data_final` but without `Q` and includes:
    -   `real tau0`: Determines the scale of the global shrinkage parameter ($\tau$)
    -   `real slab_scale`: Scale for $\hat\alpha$ values that deviate from the mean
    -   `real slab_df`: Effective degrees of freedom for significant $\hat\alpha$ values
-   `ppc_list`: List posterior predictive model checks
    -   `plot`: Plot with 95% posterior-predictive intervals plotted beneath observed values.
    -   `prop_captured`: Proportion of observations captured in 95% posterior predictive intervals
    -   `df`: Dataframe with observed values and posterior predictive quantiles.

## Functions (`functions.R`)

-   `construct_M(n_env_covs, int = F)`: Construct a linear map from the linear model parameters to the deviations from the global mean effects (alpha_hat_eij) for each site (Bendering, Perenjori) based on 0, 1, or 2 environmental covariates, with and without the interaction between two environmental covariates.

-   `construct_alpha_hat(M, B_post, ncovs, covnames, sitenames=c("bendering", "perenjori"))`: Construct array of draws from the posterior distribution of the alpha_hat vector for each competitor species based on the mapping `M` from B -\> alpha_hat.

-   `non_generic(alpha_hat_draws, level = 0.5, sp_names=NA, all_devs = F)`: Determine which species' effects differ substantially from the overall mean effect based on which alpha_hat parameter posterior distributions pull away from zero.

-   `ppcs(y_rep, obs, plot = T, pred_level = 0.95)`: Posterior predictive checks of the model fit.

-   `plot_heatmap(m_formula, param_estims, sp_name, ncells = 100, dem_param = "lambda", ...)`: Create a heatmap of LDGR or another demographic parameter with phosphorous on one axis and canopy on the other.
