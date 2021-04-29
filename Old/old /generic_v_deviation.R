library(rstan)
library(coda)
rm(list=ls())

# load in the rdata posteriors - just do one focal species in one script 

load("Topher model fits/Waitzia_Phos_Finalfit.rdata")
waac.phos <- rstan::extract(FinalFit)
dim(waac.phos$alphas)
remove(FinalFit)


# load in s and g data 
load("SurvivalAndGermination/Germination.rdata")
germ <- rstan::extract(PrelimFit)
germ<-as.data.frame(germ)
remove(PrelimFit)

load("SurvivalAndGermination/Survival.rdata")
surv <- rstan::extract(PrelimFit)
surv<-as.data.frame(surv)
remove(PrelimFit)

germination <- data.frame(arca=germ$p.1, waac=germ$p.2)

survival <- data.frame(arca=surv$p.1, waac=surv$p.2)




# coexistence simulations - invasion growth rates -----------------------------------------

invader.abund <- 1

# things to set prior for hard coding 

germ=germination$waac
surv=survival$waac
lambda = waac.phos$lambdas
alphas = waac.phos$alphas

# BENDERING --------------------------------------------------------------------------------
# run the simulation - bendering generic only
posteriors=9000 # how many values from the posterior we want to sample?
plots=51 # run through each option for resident community population size 
b.waac.into.generic <- matrix(NA, nrow=plots, ncol=9000) # make an empty matrix (array for env??)
#e = # to run through increments for each environmental condition
for (r in 1:plots) {
  for (x in 1:posteriors) {
    #x <- sample(seq(1, 9000),1)
    log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]*0+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,]*0)*env[r] 
    #log_a_eij <- alphas[x,1]+ (alphas[x,2])*env[r] 
    # generic intercept effect + deviation intercept effect * whether it needs including + (generic slope effect + deviation effect*whether it needs including based on species presence) * environmental covariate
    
    #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
    a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
    
    
    y <- sample(seq(1,4500),1)
    # calculate resident abundance
    #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
    #y <- sample(seq(1,129),1)
    Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
    #germ_j <- germination$waac[y]
    Nj[45] <- 0 # setting WAAC community abundance to 0
    
    
    # calculate lambda 
    overall_lambda <- lambda[x,,1]+lambda[x,,2]*env[r] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
    lambdas <- exp(overall_lambda[reserve[r]])
    
    # invade WAAC
    waac_tp1 <- (1-germ[y])*surv[y]*1 + 
      1*germ[y]*lambdas/(1+sum(a_eij*Nj))
    # calculate LDGR of WAAC
    b.waac.into.generic[r,x] <- log(waac_tp1/1)
    
  }
  print(r)
}
hist(b.waac.into.generic, xlim=c(-2,2))

# Bendering - generic plus deviations
posteriors=9000 # how many values from the posterior we want to sample?
plots=51 # run through each option for resident community population size 
b.waac.into.generic.dev <- matrix(NA, nrow=plots, ncol=9000) # make an empty matrix (array for env??)
#e = # to run through increments for each environmental condition
for (r in 1:plots) {
  for (x in 1:posteriors) {
    #x <- sample(seq(1, 9000),1)
    #log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env[r] 
    log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env[r] 
    # generic intercept effect + deviation intercept effect * whether it needs including + (generic slope effect + deviation effect*whether it needs including based on species presence) * environmental covariate
    
    #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
    a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
    
    
    y <- sample(seq(1,4500),1)
    # calculate resident abundance
    #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
    #y <- sample(seq(1,129),1)
    Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
    #germ_j <- germination$waac[y]
    Nj[45] <- 0 # setting WAAC community abundance to 0
    
    
    # calculate lambda 
    overall_lambda <- lambda[x,,1]+lambda[x,,2]*env[r] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
    lambdas <- exp(overall_lambda[reserve[r]])
    
    # invade WAAC
    waac_tp1 <- (1-germ[y])*surv[y]*1 + 
      1*germ[y]*lambdas/(1+sum(a_eij*Nj))
    # calculate LDGR of WAAC
    b.waac.into.generic.dev[r,x] <- log(waac_tp1/1)
    
  }
  print(r)
}
hist(b.waac.into.generic.dev, xlim=c(-2,2))

# PERENJORI ----------------------------------------------------
# run the simulation -generic only
posteriors=9000 # how many values from the posterior we want to sample?
plots=129 # run through each option for resident community population size 
p.waac.into.generic <- matrix(NA, nrow=plots, ncol=9000) # make an empty matrix (array for env??)
#e = # to run through increments for each environmental condition
for (r in 52:plots) {
  for (x in 1:posteriors) {
    #x <- sample(seq(1, 9000),1)
    log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]*0+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,]*0)*env[r] 
    #log_a_eij <- alphas[x,1]+ (alphas[x,2])*env[r] 
    # generic intercept effect + deviation intercept effect * whether it needs including + (generic slope effect + deviation effect*whether it needs including based on species presence) * environmental covariate
    
    #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
    a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
    
    
    y <- sample(seq(1,4500),1)
    # calculate resident abundance
    #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
    #y <- sample(seq(1,129),1)
    Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
    #germ_j <- germination$waac[y]
    Nj[45] <- 0 # setting WAAC community abundance to 0
    
    
    # calculate lambda 
    overall_lambda <- lambda[x,,1]+lambda[x,,2]*env[r] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
    lambdas <- exp(overall_lambda[reserve[r]])
    
    # invade WAAC
    waac_tp1 <- (1-germ[y])*surv[y]*1 + 
      1*germ[y]*lambdas/(1+sum(a_eij*Nj))
    # calculate LDGR of WAAC
    p.waac.into.generic[r,x] <- log(waac_tp1/1)
    
  }
  print(r)
}
hist(p.waac.into.generic, xlim=c(-2,2))

# generic plus deviations
posteriors=9000 # how many values from the posterior we want to sample?
plots=129 # run through each option for resident community population size 
p.waac.into.generic.dev <- matrix(NA, nrow=plots, ncol=9000) # make an empty matrix (array for env??)
#e = # to run through increments for each environmental condition
for (r in 52:plots) {
  for (x in 1:posteriors) {
    #x <- sample(seq(1, 9000),1)
    #log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env[r] 
    log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env[r] 
    # generic intercept effect + deviation intercept effect * whether it needs including + (generic slope effect + deviation effect*whether it needs including based on species presence) * environmental covariate
    
    #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
    a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
    
    
    y <- sample(seq(1,4500),1)
    # calculate resident abundance
    #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
    #y <- sample(seq(1,129),1)
    Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
    #germ_j <- germination$waac[y]
    Nj[45] <- 0 # setting WAAC community abundance to 0
    
    
    # calculate lambda 
    overall_lambda <- lambda[x,,1]+lambda[x,,2]*env[r] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
    lambdas <- exp(overall_lambda[reserve[r]])
    
    # invade WAAC
    waac_tp1 <- (1-germ[y])*surv[y]*1 + 
      1*germ[y]*lambdas/(1+sum(a_eij*Nj))
    # calculate LDGR of WAAC
    p.waac.into.generic.dev[r,x] <- log(waac_tp1/1)
    
  }
  print(r)
}
hist(p.waac.into.generic.dev, xlim=c(-2,2))



