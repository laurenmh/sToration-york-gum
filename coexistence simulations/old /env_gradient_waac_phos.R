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


env_gradient <- cbind(env[1:51]-2, env[1:51]-1, env[1:51], 1+env[1:51], 2+env[1:51]) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation
#env_gradient <- cbind(rep(mean(env[1:51]-2),51), rep(mean(env[1:51]-1),51), rep(mean(env[1:51]),51), rep(mean(1+env[1:51]),51), rep(mean(2+env[1:51]),51))

# bendering 
posteriors=9000 
plots=51 # run through each option for resident community population size 
waac.into.observed <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
    #x <- sample(seq(1, 9000),1)
    #log_a_eij <- array()
    log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
    
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
    overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
    lambdas <- exp(overall_lambda[reserve[r]])
    
    # invade WAAC
    waac_tp1 <- (1-germ[y])*surv[y]*1 + 
      1*germ[y]*lambdas/(1+sum(a_eij*Nj))
    # calculate LDGR of WAAC
    waac.into.observed[r,x,e] <- log(waac_tp1/1)
    
  }
  }
}

# with no comp
waac.into.observed1 <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      #x <- sample(seq(1, 9000),1)
      #log_a_eij <- array()
      log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[reserve[r],]) #to get your alpha term for each environmental condition and each competitor
      
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
      #y <- sample(seq(1,129),1)
      Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #germ_j <- germination$waac[y]
      Nj <- as.integer(rep(0,45)) # setting WAAC community abundance to 0
      
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[reserve[r]])
      
      # invade WAAC
      waac_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj))
      # calculate LDGR of WAAC
      waac.into.observed1[r,x,e] <- log(waac_tp1/1)
      
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

data.melt <- melt(waac.into.observed, varnames = c("plot", "index", "e"), value.name = "ldgr")
data.melt$e <- as.factor(data.melt$e)
data.melt$plot <- as.factor(data.melt$plot)
data.melt1 <- melt(waac.into.observed1, varnames = c("plot", "index", "e"), value.name = "ldgr")
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
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("-2 sd", "-1 sd", "sd", "+1 sd", "+2 sd"))+
  theme_ridges(center = T)
ggsave(filename = file.path("Figures","phos_grad_waac_bend.pdf"))

# perenjori ---------------------------------------------------------------------------


invader.abund <- 1

# things to set prior for hard coding 

germ=germination$waac
surv=survival$waac
lambda = waac.phos$lambdas
alphas = waac.phos$alphas


env_gradient <- cbind(env[52:129]-2, env[52:129]-1, env[52:129], 1+env[52:129], 2+env[52:129]) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation
#env_gradient <- cbind(rep(mean(env[1:51]-2),51), rep(mean(env[1:51]-1),51), rep(mean(env[1:51]),51), rep(mean(1+env[1:51]),51), rep(mean(2+env[1:51]),51))

# bendering 
posteriors=9000 
plots=78 # run through each option for resident community population size 
waac.into.observed <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      #x <- sample(seq(1, 9000),1)
      #log_a_eij <- array()
      log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
      #y <- sample(seq(1,129),1)
      Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #germ_j <- germination$waac[y]
      Nj[45] <- 0 # setting WAAC community abundance to 0
      
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade WAAC
      waac_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj))
      # calculate LDGR of WAAC
      waac.into.observed[r,x,e] <- log(waac_tp1/1)
      
    }
  }
}

# with no comp
waac.into.observed1 <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      #x <- sample(seq(1, 9000),1)
      #log_a_eij <- array()
      log_a_eij <- alphas[x,1]+Inclusion_ij*waac.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*waac.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- arca.equilibrum[r] # do we still not actually need to simulate each species to equilibrium pop density first?
      #y <- sample(seq(1,129),1)
      Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #germ_j <- germination$waac[y]
      Nj <- as.integer(rep(0,45)) # setting WAAC community abundance to 0
      
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade WAAC
      waac_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj))
      # calculate LDGR of WAAC
      waac.into.observed1[r,x,e] <- log(waac_tp1/1)
      
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

data.melt <- melt(waac.into.observed, varnames = c("plot", "index", "e"), value.name = "ldgr")
data.melt$e <- as.factor(data.melt$e)
data.melt$plot <- as.factor(data.melt$plot)
data.melt1 <- melt(waac.into.observed1, varnames = c("plot", "index", "e"), value.name = "ldgr")
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
  ylab("Phosphorus") +
  scale_y_discrete(labels =c("-2 sd", "-1 sd", "sd", "+1 sd", "+2 sd"))+
  theme_ridges(center = T)
ggsave(filename = file.path("Figures","phos_grad_waac_pj.pdf"))


# OLD CODE 
# ggplot(NULL) +
#   stat_density_ridges(data=data.melt, geom = "density_ridges_gradient", calc_ecdf = TRUE, 
#                       aes(x = ldgr, y = e, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   stat_density_ridges(data=data.melt1, geom = "density_ridges_gradient", calc_ecdf = TRUE,
#                       aes(x = ldgr, y = e, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
#   #scale_fill_viridis_c(name = "Tail probability", direction = -1, option="magma") +
#   #scale_fill_viridis_c(name = "Tail probability", direction = -1, option="plasma") +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   ylab("Phosphorus") +
#   scale_y_discrete(labels =c("-2 sd", "-1 sd", "sd", "+1 sd", "+2 sd"))
# 
