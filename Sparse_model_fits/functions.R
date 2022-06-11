### Functions to facilitate code readability ###






#' Construct mapping from b_j to alpha_hat_j
#'
#' @param n_env_covs Number of environmental covariates included in the model {0, 1, 2}
#' @param int Is there an interaction between the two environmental covariates in the model?
#'
#' @return Matrix of (nsites*ncovs) x (nsites*ncovs)
#' @export
#'
#' @examples
#' 
construct_M <- function(n_env_covs, int = F){
  opts <- c(0,1,2)
  n_env_covs %in% opts
  if(isFALSE(n_env_covs %in% opts)){
    stop("Sorry, not a very flexible function. n_env_covs must be in {0, 1, 2}.")
  }
  message("Constructing M based on default R dummy coding for Z_alpha.\n
          If Z_alpha was not constructed using dummy variables for the site (only two sites allowed),
          this function will give the wrong result.")
  
  # build upper left block
  M_sub <- diag(2)
  M_sub[,1] <- rep(1, 2)
  
  # intercept only
  if(n_env_covs == 0){ return(M_sub) }
  
  # for 1 env covariate
  if(n_env_covs == 1){
    # expand to block diagonal matrix
    Id <- diag(n_env_covs + 1)
    return(Id %x% M_sub)
  }
  
  # for 2 env covariates, no interaction
  if(n_env_covs == 2 & isFALSE(int)){
    return(
      rbind(
        c(1,0,0,0,0,0),
        c(1,1,0,0,0,0),
        c(0,0,1,0,0,0),
        c(0,0,1,0,1,0),
        c(0,0,0,1,0,0),
        c(0,0,0,1,0,1)
      )
    )
  }
  
  # for 2 env covariates plus an interaction between them
  if(n_env_covs == 2 & isTRUE(int)){
    return(
      rbind(
        c(1,0,0,0,0,0,0,0),
        c(1,1,0,0,0,0,0,0),
        c(0,0,1,0,0,0,0,0),
        c(0,0,1,0,1,0,0,0),
        c(0,0,0,1,0,0,0,0),
        c(0,0,0,1,0,1,0,0),
        c(0,0,0,0,0,0,1,0),
        c(0,0,0,0,0,0,1,1)
      )
    )
  }
  
}







#' Reconstruct an array of species-specific deviations from generic values
#' 
#' Reconstruct an R x (P*nsites) x S array of species-specific deviations from
#' generic model parameters, where R is number of posterior draws, P is number of parameters
#' of interest for each species in each site (e.g., 2 for intercept and one slope), and S is
#' number of species.
#'
#' @param M Matrix with one row for each parameter of interest (alpha_hat) mapping
#' b_j to alpha_hat. The matrix should be designed to construct b_0 for site 1 in row one, b_0 for site 2 in row 2, ...,
#' b_1 for site one in row nsites + 1, b_1 for site 2 in row nsites + 1, ...
#' @param B_post Array of posterior draws from deviations matrix (B) with dimensions R x P_alpha_d x S
#' @param covnames Names of covariates of interest (including intercept) in the same order as the
#' original model.matrix call for Z_alpha
#' @param sitenames Site/reserve names in same order as levels of factor from model.matrix call for Z_alpha
#'
#' @return Array with dimensions R x (P*nsites) x S
#' @export
#'
#' @examples
#' 
construct_alpha_hat <- function(
  M, B_post, ncovs, covnames, sitenames=c("bendering", "perenjori")
){
  npars <- ncovs
  nsites <- length(sitenames)
  
  # double-check constraints
  if(
    nrow(M) != npars * nsites |
    nrow(M) != npars * nsites
  ) {
    stop("Model matrix must be square with dimensions npars * nsites")
  }
    
  # empty array with same dimensions as B_post for filling in desired parameters
  alpha_hat_post <- array(
    dim = dim(B_post)
  )
  
  # construct a vector of parameter names for clarity
  alpha_hat_j <- paste("alpha_hat", covnames, sep = "_")
  parnames <- vector(mode="character", length(alpha_hat_j) * nsites)
  for(i in 1:length(alpha_hat_j)){
    for(j in 1:nsites){
      parnames[nsites*(i - 1) + j] <- paste(alpha_hat_j[i], sitenames[j], sep = "_")
    }
  }
  dimnames(alpha_hat_post)[[2]] <- parnames

  # compute parameters values of interest for each draw s from the posterior
  for(s in 1:dim(B_post)[1]){
    alpha_hat_post[s,,] <- M %*% B_post[s,,]
  }
  
  return(alpha_hat_post)
}









#' Construct constraints on alpha_hat matrix
#'
#' @param alpha_hat_draws R x P_alpha x S array of posterior draws from alpha_hat vector 
#' constructed using the function \code(construct_alpha_hat())
#' @param level Posterior credible level to construct interval
#' @param sp_names Optional vector of species names for tracking and checking purposes
#' @param all_devs If true, all deviations from generic effects will be possible for a species
#' with at least one alpha_hat that shows deviation from zero. If false, all but the alpha_hat 
#' whose posterior pulls away from zero will be constrained.
#'
#' @return Matrix of alpha_hats to include in the final model fit
#' @export
#'
#' @examples
#' 
non_generic <- function(alpha_hat_draws, level = 0.5, sp_names=NA, all_devs = F){
  
  # create Inclusion matrix
  # number of parameters
  P <- dim(alpha_hat_draws)[[2]]
  # number of heterospecifics
  S <- dim(alpha_hat_draws)[[3]]
  Inclusion <- matrix(
    data = 0,
    nrow = P,
    ncol = S
  )
  
  # row and column names
  rownames(Inclusion) <- dimnames(alpha_hat_draws)[[2]]
  if(!is.na(sp_names)){
    colnames(Inclusion) <- sp_names
  }
  
  # construct based on deviations from generic parameters
  for(i in 1:P){
    for(s in 1:S){
      Ints_ij <- HDInterval::hdi(alpha_hat_draws[,i,s], credMass = level)
      if(Ints_ij[1] > 0 | Ints_ij[2] < 0){
        Inclusion[i,s] <- 1
      }
    }
  }
  # Include all parameters for a species with at least one deviating from the mean?
  if(isFALSE(all_devs)){
    return(Inclusion)
  } else {
    key_sp <- which(colSums(Inclusion) > 0)
    ones <- rep(1, nrow(Inclusion))
    Inclusion[, key_sp] <- ones
    return(Inclusion)
  }
  
}
  


#' Posterior predictive checks
#'
#' @param y_rep Matrix of S predicted values sampled from the
#'  posterior predictive distribution for each of N observed values
#' @param obs Vector of N observed values
#' @param plot Should a PPC plot be returned
#' @param pred_level Level of posterior predictive interval
#'
#' @return Either a list or the single value of the proportion of observed values captured by
#' the posterior predictive intervals (should be around pred_level given a good fit)
#' @export
#'
#' @examples
#' 
ppcs <- function(y_rep, obs, plot = T, pred_level = 0.95){
  # build df
  df_ppc <- data.frame(
    y = obs,
    pred_low = apply(y_rep, 2, quantile, probs = (1 - pred_level)/2),
    pred_high = apply(y_rep, 2, quantile, probs = 1 - (1 - pred_level)/2)
  )
  
  # order the df by ascending fecundity
  df_ppc_ord <- df_ppc[order(df_ppc$y), ]
  df_ppc_ord$x <- 1:nrow(df_ppc_ord)
  df_ppc_ord$original_obs_id <- order(df_ppc$y)
  
  # create list with plot and summary stat
  if(isTRUE(plot)){
    ppc_plot <- ggplot2::ggplot(data = df_ppc_ord, aes(x = x))+
      geom_errorbar(
        aes(ymin = pred_low, ymax = pred_high),
        color = "grey",
        width = 0,
        size = 1.5
      )+
      geom_point(aes(y = y), color = "red", size = 0.5)+
      theme_bw()+
      xlab("")+
      ylab("Observed (red) and predicted fecundity")
    
    # calculate proportion of observations captured
    prop_captured <- mean(df_ppc$y >= df_ppc$pred_low & df_ppc$y <= df_ppc$pred_high)
    
    # find which observations fall outside the predictive intervals
    outside_obs <- which(df_ppc$y <= df_ppc$pred_low | df_ppc$y >= df_ppc$pred_high)
    
    # message about results
    print(paste(
      "Proportion of observations captured in ", 
      pred_level*100, 
      "% posterior predictive intervals: ", 
      round(prop_captured, 2),
      sep = ""
    ))
    
    return(list(
      plot = ppc_plot,
      prop_captured = prop_captured,
      df = df_ppc_ord
    ))
  }
  
  if(isFALSE(plot)){
    return(
      mean(df_ppc$y >= df_ppc$pred_low & df_ppc$y <= df_ppc$pred_high)
    )
  }
}






#' Posterior mode
#'
#' @param x Vector of draws from the marginal posterior of a parameter 
#'
#' @return Mode of the posterior distribution
#' @export
#'
#' @examples
post_mode <- function(x){
  max_x <- which.max(
    density(x)$y
  )
  return(x[max_x])
}








#' Plot a heatmap for a demographic parameter of interest
#'
#' @param m_formula Model formula used for the linear model component modeling the demographic
#' parameter of interest (should be in the same form as in the ModelFits2 script).
#' @param param_estims Bayes or other point estimator for the linear model coefficients
#' @param sp_name Name of the focal species
#' @param ncells Number of cells on a side for the heatmap
#' @param dem_param Demographic parameter of interest (e.g., lambda, alpha, etc.)
#' @param ... Additional options passed to the plotting function. Currently only takes a "viridis" color ramp option (see example).
#'
#' @return Two-panel heatmap plot, one for Bendering, the other for Perenjori.
#' @export
#'
#' @examples
#' 
plot_heatmap <- function(m_formula, param_estims, sp_name, ncells = 100, dem_param, ...){
  
  library(patchwork, quietly = T)
  arg_list <- list(...)
  if(is.null(arg_list$option)){arg_list$option <- "inferno"}
  # construct new df
  df_new <- data.frame(
    Reserve.x = as.factor(rep(c("Bendering", "Perenjori"), each = ncells^2)),
    Colwell.P_std = rep(rep(seq(-1, 1, length.out = ncells), ncells), 2),
    Canopy_std = rep(rep(seq(-1, 1, length.out = ncells), each = ncells), 2)
  )
  
  # construct model matrix
  X_new <- model.matrix(
    m_formula,
    data = df_new
  )
  
  # preds
  df_new$pred <- exp(X_new %*% param_estims)
  
  # plot theme
  hm_theme <- theme(
    legend.title = element_text(size = 14)
  )


  bend <- with(arg_list, {
    ggplot(data = subset(df_new, Reserve.x == "Bendering"))+
      geom_tile(
        aes(x = Colwell.P_std, y = Canopy_std, fill = pred)
      )+
      viridis::scale_fill_viridis(option = arg_list$option, discrete = F)+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      ggtitle("Bendering")+
      labs(fill = parse(text = paste("hat(", dem_param, ")")))+
      xlab("Standardized Phosphorous")+
      ylab("Standardized Canopy Cover")+
      hm_theme
  })
  
  pere <- with(arg_list, {
    ggplot(data = subset(df_new, Reserve.x == "Perenjori"))+
      geom_tile(
        aes(x = Colwell.P_std, y = Canopy_std, fill = pred)
      )+
      viridis::scale_fill_viridis(option = arg_list$option, discrete = F)+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      ggtitle("Perenjori")+
      labs(fill = parse(text = paste("hat(", dem_param, ")")))+
      xlab("Standardized Phosphorous")+
      ylab("")+
      hm_theme
  })
  
  # xlabel <- ggplot()+
  #   geom_text(aes(
  #     x = 1, y = 1,
  #     label = "Standardized Phosphorous"
  #   ))+
  #   theme_void()+
  #   coord_cartesian(clip = "off")
  
  (bend + pere)  +
    plot_annotation(sp_name)

}


get_lambda <- function(m_formula, param_estims, sp_name, ncells = 100, dem_param = "lambda", ...){
  
  library(patchwork, quietly = T)
  arg_list <- list(...)
  if(is.null(arg_list$option)){arg_list$option <- "inferno"}
  # construct new df
  df_new <- data.frame(
    Reserve.x = as.factor(rep(c("Bendering", "Perenjori"), each = ncells^2)),
    Colwell.P_std = rep(rep(seq(-1, 1, length.out = ncells), ncells), 2),
    Canopy_std = rep(rep(seq(-1, 1, length.out = ncells), each = ncells), 2)
  )
  
  # construct model matrix
  X_new <- model.matrix(
    m_formula,
    data = df_new
  )
  
  # preds
  df_new$pred <- exp(X_new %*% param_estims)
  
  return(df_new)
  
}


plotting_LDGR <- function(sp_name, df_new,  ...) {
  library(patchwork, quietly = T)
  arg_list <- list(...)
  if(is.null(arg_list$option)){arg_list$option <- "inferno"}
  # plot theme
  hm_theme <- theme(
    legend.title = element_text(size = 14)
  )
  
  
  bend <- with(arg_list, {
    ggplot(data = subset(df_new, Reserve.x == "Bendering"))+
      geom_tile(
        aes(x = Colwell.P_std, y = Canopy_std, fill = LDGR)
      )+
      viridis::scale_fill_viridis(option = arg_list$option, discrete = F)+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      ggtitle("Bendering")+
      labs(fill = parse(text = paste("LDGR")))+
      xlab("Standardized Phosphorous")+
      ylab("Standardized Canopy Cover")+
      hm_theme
  })
  
  pere <- with(arg_list, {
    ggplot(data = subset(df_new, Reserve.x == "Perenjori"))+
      geom_tile(
        aes(x = Colwell.P_std, y = Canopy_std, fill = LDGR)
      )+
      viridis::scale_fill_viridis(option = arg_list$option, discrete = F)+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      ggtitle("Perenjori")+
      labs(fill = parse(text = paste("LDGR")))+
      xlab("Standardized Phosphorous")+
      ylab("")+
      hm_theme
  })
  
  # xlabel <- ggplot()+
  #   geom_text(aes(
  #     x = 1, y = 1,
  #     label = "Standardized Phosphorous"
  #   ))+
  #   theme_void()+
  #   coord_cartesian(clip = "off")
  
  (bend + pere)  +
    plot_annotation(sp_name)
  
}


plotting_Diff <- function(sp_name, df_new, ...) {
  library(patchwork, quietly = T)
  arg_list <- list(...)
  if(is.null(arg_list$option)){arg_list$option <- "inferno"}
  # plot theme
  hm_theme <- theme(
    legend.title = element_text(size = 14)
  )
  
  
  bend <- with(arg_list, {
    ggplot(data = subset(df_new, Reserve.x == "Bendering"))+
      geom_tile(
        aes(x = Colwell.P_std, y = Canopy_std, fill = Diff)
      )+
      viridis::scale_fill_viridis(option = arg_list$option, discrete = F)+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      ggtitle("Bendering")+
      labs(fill = parse(text = paste("LDGR-log(lambda)")))+
      xlab("Standardized Phosphorous")+
      ylab("Standardized Canopy Cover")+
      hm_theme
  })
  
  pere <- with(arg_list, {
    ggplot(data = subset(df_new, Reserve.x == "Perenjori"))+
      geom_tile(
        aes(x = Colwell.P_std, y = Canopy_std, fill = Diff)
      )+
      viridis::scale_fill_viridis(option = arg_list$option, discrete = F)+
      scale_x_continuous(expand = c(0, 0))+
      scale_y_continuous(expand = c(0, 0))+
      ggtitle("Perenjori")+
      labs(fill = parse(text = paste("LDGR-log(lambda)")))+
      xlab("Standardized Phosphorous")+
      ylab("")+
      hm_theme
  })
  
  # xlabel <- ggplot()+
  #   geom_text(aes(
  #     x = 1, y = 1,
  #     label = "Standardized Phosphorous"
  #   ))+
  #   theme_void()+
  #   coord_cartesian(clip = "off")
  
  (bend + pere)  +
    plot_annotation(sp_name)
  
}








