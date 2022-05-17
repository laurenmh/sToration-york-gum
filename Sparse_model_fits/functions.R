### Functions to facilitate code readability ###






#' Construct mapping from b_j to alpha_hat_j
#'
#' @param nsites Number of sites (reserves) at which data were collected
#' @param ncovs Number of covariates (including the intercept) b_j
#'
#' @return Matrix of (nsites*ncovs) x (nsites*ncovs)
#' @export
#'
#' @examples
#' 
construct_M_dc <- function(nsites = 2, ncovs){
  message("Constructing M based on dummy coding for Z_alpha.\n
          If Z_alpha was not constructed using dummy variables for the site,
          this function will give the wrong result.")
  
  # build one block
  M_sub <- diag(nsites)
  M_sub[,1] <- rep(1, nsites)
  
  # expand to block diagonal matrix
  Id <- diag(ncovs)
  return(Id %x% M_sub)
  
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
  M, B_post, covnames = c("int", "phos", "shade"), sitenames=c("bendering", "perenjori")
){
  npars <- length(covnames)
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
#'
#' @return List with one matrix of alpha_hats to include and another with the corresponding constraints on B
#' @export
#'
#' @examples
#' 
non_generic <- function(alpha_hat_draws, level = 0.5, sp_names=NA){
  
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
  
  return(Inclusion)
  
}
  






