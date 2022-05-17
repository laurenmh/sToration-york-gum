### Functions to facilitate code readability ###





#' Reconstruct an array of species-specific deviations from generic values
#' 
#' Reconstruct an R x (P*nsites) x S array of species-specific deviations from
#' generic model parameters, where R is number of posterior draws, P is number of parameters
#' of interest for each species in each site (e.g., 2 for intercept and one slope), and S is
#' number of species.
#'
#' @param mm_new Matrix with one row for each parameter of interest (alpha_hat) for each
#' site with 1's for model parameters necessary to reconstruct alpha_hat and zeros elsewhere.
#' The matrix should be designed to construct b_0 for site 1 in row one, b_0 for site 2 in row 2, ...,
#' b_1 for site one in row nsites + 1, b_1 for site 2 in row nsites + 1, ...
#' @param B_post Posterior draws from deviations matrix R x P_alpha_d x S
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
  mm_new, B_post, covnames = c("int", "phos", "shade"), sitenames=c("bendering", "perenjori")
){
  npars <- length(covnames)
  nsites <- length(sitenames)
  
  # double-check constraints
  if(
    nrow(mm_new) != npars * nsites |
    nrow(mm_new) != npars * nsites
  ) {
    stop("Model matrix must be square with dimensions npars * nsites")
  }
    
  # empty array with same dimensions as B_post for filling in desired parameters
  alpha_hat_post <- array(
    dim = dim(B_post)
  )
  
  # construct a vector of parameter names for clarity
  b_j <- paste("b", covnames, sep = "_")
  parnames <- vector(mode="character", length(b_j) * nsites)
  for(i in 1:length(b_j)){
    for(j in 1:nsites){
      parnames[nsites*(i - 1) + j] <- paste(b_j[i], sitenames[j], sep = "_")
    }
  }
  dimnames(alpha_hat_post)[[2]] <- parnames

  # compute parameters values of interest for each draw s from the posterior
  for(s in 1:dim(B_post)[1]){
    alpha_hat_post[s,,] <- mm_new %*% B_post[s,,]
  }
  
  return(alpha_hat_post)
}









#' Construct constraints on B matrix
#'
#' @param alpha_hat_draws R x P_alpha x S array of posterior draws from alpha_hat vector 
#' constructed using the function \code(construct_alpha_hat())
#' @param level Posterior credible level to construct interval
#' @param mm_new Matrix used to construct alpha_hat_draws in call to \code(construct_alpha_hat())
#' @param sp_names Optional vector of species names for tracking and checking purposes
#'
#' @return List with one matrix of alpha_hats to include and another with the corresponding constraints on B
#' @export
#'
#' @examples
#' 
non_generic <- function(alpha_hat_draws, level = 0.5, mm_new, sp_names=NA){
  
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
  
  # construct constrained parameter matrix
  B_constraints <- matrix(
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
  
  # construct constraint matrix for B
  for(s in 1:S){
    devs <- which(Inclusion[,s] == 1)
    if(length(devs) > 0){
      B_constraints[,s] <- colSums(mm_new[devs, , drop=FALSE])
    }
  }
  
  # reset occasions where param was counted twice because it 
  #  deviates from generic values in both sites
  B_constraints[B_constraints > 0] <- 1
  
  return(list(
    deviations = Inclusion,
    B_constraints = B_constraints
  ))
  
}
  










