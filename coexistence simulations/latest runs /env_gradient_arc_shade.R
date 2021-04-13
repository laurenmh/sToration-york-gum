library(rstan)
library(coda)
rm(list=ls())

# load in the rdata posteriors - just do one focal species in one script 

load("Topher model fits/ARCA_Shade_Finalfit.rdata")
arca.phos <- rstan::extract(FinalFit)
dim(arca.phos$alphas)
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

germ=germination$arca
surv=survival$arca
lambda = arca.phos$lambdas
alphas = arca.phos$alphas


env_gradient <- cbind(env[1:28]-2, env[1:28]-1, env[1:28], 1+env[1:28], 2+env[1:28]) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation
#env_gradient <- cbind(rep(mean(env[1:51]-2),51), rep(mean(env[1:51]-1),51), rep(mean(env[1:51]),51), rep(mean(1+env[1:51]),51), rep(mean(2+env[1:51]),51))

# bendering 
posteriors=9000 
plots=28 # run through each option for resident community population size 
arca.into.observed <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      #x <- sample(seq(1, 9000),1)
      #log_a_eij <- array()
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
      
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
      #y <- sample(seq(1,129),1)
      Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      na.omit(Nj)
      #germ_j <- germination$waac[y]
      Nj[40] <- 0 # setting WAAC community abundance to 0
      
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[reserve[r]])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj))
      # calculate LDGR of WAAC
      arca.into.observed[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# with only generic competition 
arca.into.observed1 <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      #x <- sample(seq(1, 9000),1)
      #log_a_eij <- array()
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
      
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
      #y <- sample(seq(1,129),1)
      Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      na.omit(Nj)
      #germ_j <- germination$waac[y]
      Nj <- as.integer(rep(0,40)) # setting WAAC community abundance to 0
      
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[reserve[r]])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj))
      # calculate LDGR of WAAC
      arca.into.observed1[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}


# plots, posteriors, env gradient 

#waac.into.observed[1,,]
#plot(density(waac.into.observed[,,1]))
#plot(density(waac.into.observed[,,2]))

# alternative to heatmap - playing with ggridges (ridgeline package)
library(reshape2)
library(ggplot2)
#install.packages("ggridges")
library(ggridges)

data.melt <- melt(arca.into.observed, varnames = c("plot", "index", "e"), value.name = "ldgr")
data.melt$e <- as.factor(data.melt$e)
data.melt$plot <- as.factor(data.melt$plot)
data.melt1 <- melt(arca.into.observed1, varnames = c("plot", "index", "e"), value.name = "ldgr")
data.melt1$e <- as.factor(data.melt1$e)
data.melt1$plot <- as.factor(data.melt1$plot)
#str(data.melt)

ggplot(NULL) +
  stat_density_ridges(data=data.melt, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "comp")) +
  stat_density_ridges(data=data.melt1, geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      aes(x = ldgr, y = e, fill = "no comp")) +
  scale_fill_manual(values = c("#D55E0050", "#0072B250"), labels = c("comp", "no comp")) +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("-2 sd", "-1 sd", "sd", "+1 sd", "+2 sd"))+
  theme_ridges(center = T)

ggsave(filename = file.path("Figures","shade_grad_arca_bend.pdf"))
# Perenjori ----------------------------------------------------------------------------

env_gradient <- cbind(env[29:95]-2, env[29:95]-1, env[29:95], 1+env[29:95], 2+env[29:95]) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation
#env_gradient <- cbind(rep(mean(env[1:51]-2),51), rep(mean(env[1:51]-1),51), rep(mean(env[1:51]),51), rep(mean(1+env[1:51]),51), rep(mean(2+env[1:51]),51))

# pj 
posteriors=9000 
plots=67 # run through each option for resident community population size 
arca.into.observed <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      #x <- sample(seq(1, 9000),1)
      #log_a_eij <- array()
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
      #y <- sample(seq(1,129),1)
      Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      na.omit(Nj)
      #germ_j <- germination$waac[y]
      Nj[40] <- 0 # setting WAAC community abundance to 0
      
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj))
      # calculate LDGR of WAAC
      arca.into.observed[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# with only generic competition 
arca.into.observed1 <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      #x <- sample(seq(1, 9000),1)
      #log_a_eij <- array()
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
      #y <- sample(seq(1,129),1)
      Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      na.omit(Nj)
      #germ_j <- germination$waac[y]
      Nj <- as.integer(rep(0,40)) # setting WAAC community abundance to 0
      
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj))
      # calculate LDGR of WAAC
      arca.into.observed1[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}


# plots, posteriors, env gradient 

#waac.into.observed[1,,]
#plot(density(waac.into.observed[,,1]))
#plot(density(waac.into.observed[,,2]))

# alternative to heatmap - playing with ggridges (ridgeline package)
library(reshape2)
library(ggplot2)
#install.packages("ggridges")
library(ggridges)

data.melt <- melt(arca.into.observed, varnames = c("plot", "index", "e"), value.name = "ldgr")
data.melt$e <- as.factor(data.melt$e)
data.melt$plot <- as.factor(data.melt$plot)
data.melt1 <- melt(arca.into.observed1, varnames = c("plot", "index", "e"), value.name = "ldgr")
data.melt1$e <- as.factor(data.melt1$e)
data.melt1$plot <- as.factor(data.melt1$plot)
#str(data.melt)

ggplot(NULL) +
  stat_density_ridges(data=data.melt, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
                      aes(x = ldgr, y = e, fill = "comp")) +
  stat_density_ridges(data=data.melt1, geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      aes(x = ldgr, y = e, fill = "no comp")) +
  scale_fill_manual(values = c("#D55E0050", "#0072B250"), labels = c("comp", "no comp")) +
  scale_color_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Shade") +
  scale_y_discrete(labels =c("-2 sd", "-1 sd", "sd", "+1 sd", "+2 sd"))+
  theme_ridges(center = T)
ggsave(filename = file.path("Figures","shade_grad_arca_pj.pdf"))

