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


Nj <- list(SpMatrix, 10*SpMatrix)

# bendering 
posteriors=9000 
plots=51 # run through each option for resident community population size 
waac.into.observed <- array(data=NA, dim =c(plots, posteriors, 2))#e = # to run through increments for each environmental condition

for (c in 1:length(Nj)) {
for (r in 1:plots) {
  for (x in 1:posteriors) {
    #x <- sample(seq(1, 9000),1)
    log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env[1:51][r] 
    
    #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
    a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
    
    
    y <- sample(seq(1,4500),1)
    # calculate resident abundance
    #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
    #y <- sample(seq(1,129),1)
    #germ_j <- germination$waac[y]
    Nj[45] <- 0 # setting WAAC community abundance to 0
    
    
    # calculate lambda 
    overall_lambda <- lambda[x,,1]+lambda[x,,2]*env[1:51][r] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
    lambdas <- exp(overall_lambda[reserve[r]])
    
    # invade WAAC
    waac_tp1 <- (1-germ[y])*surv[y]*1 + 
      1*germ[y]*lambdas/(1+sum(a_eij*Nj[[c]]))
    # calculate LDGR of WAAC
    waac.into.observed[r,x,c] <- log(waac_tp1/1)
    
  }
  print(r)
}
}
hist(waac.into.observed[,,1], xlim=c(-2,2)) # im doing something wrong because these aren't different 
hist(waac.into.observed[,,2], xlim=c(-2,2))
