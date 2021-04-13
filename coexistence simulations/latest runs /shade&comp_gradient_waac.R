library(rstan)
library(coda)
rm(list=ls())

# load in the rdata posteriors - just do one focal species in one script 

load("Topher model fits/Waitzia_Shade_Finalfit.rdata")
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

germ=germination$waac
surv=survival$waac
lambda = arca.phos$lambdas
alphas = arca.phos$alphas

# Bendering --------------------------------------------------------------------------------------------
colnames(SpMatrix) <- SpNames
# pick the natives
# subset(SpMatrix, select=c("Angianthus.tomentosus", "Brachyscome.iberidifolia", "Crassula.colorata.or.exserta", "Daucus.glochidiatus", "Gilberta.tenuifolia",
#                                   "Gonocarpus.nodulosus", "Goodenia.berardiana", "Haloragis.odontocarpa", "Hyalosperma.glutinosum.subsp.glutinosum", "Hydrocotyle.pilifera", 
#                                   "Lawrencella.rosea", "Nicotiana.rotundifolia", "Plantago.debilis", "Podolepis.canescens", "Podotheca.gnaphalioides", "Rhodanthe.laevis",
#                                   "Rhodanthe.manglesii", "Rhodanthe.polycephala", "Schoenia.cassiniana", "Schoenus.nanus", "Senecio.glossanthus", "Trachymene.cyanopetala",
#                                   "Velleia.cycnopotamica", "Velleia.rosea", "Waitzia.acuminata"))
# or the exotics 
# subset(SpMatrix, select=c("Arctotheca.calendula", "Calotis.hispidula", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
#                                   "Oxalis.sp", "Parentucellia.latifolia", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus", "Urospermum.picroides",
#                                   "Ursinia.anthemoides", "Vulpia.bromoides"))

# competitor gradient
#sort(colSums(SpMatrix))
Nj <- SpMatrix 
Nj2 <- SpMatrix 
Nj3 <- SpMatrix
Nj4 <- SpMatrix
Nj5 <- SpMatrix 
Nj6 <- SpMatrix


#colnames(SpMatrix)

# calculating average stem abundances or natives and exotics 
# nat.stems <- Nj[,c("Brachyscome.sp" , "Brachyscome.iberidifolia", "Crassula.colorata.or.exserta","Calotis.multicaulis", "Calotis.hispidula", 
#                    "Calandrinia.eremaea","Daucus.glochidiatus", "Gilberta.tenuifolia", "Phyllangium.sulcatum",
#                    "Gonocarpus.nodulosus", "Gilberta.tenuifolia", "Goodenia.berardiana", "Haloragis.odontocarpa", "Hyalosperma.glutinosum.subsp.glutinosum", "Hydrocotyle.pilifera", 
#                    "Lawrencella.rosea", "Nicotiana.rotundifolia", "Plantago.debilis", "Podolepis.canescens", "Podotheca.gnaphalioides", "Rhodanthe.laevis", "Rhodanthe.citrina",
#                    "Rhodanthe.manglesii", "Rhodanthe.polycephala", "Schoenia.cassiniana", "Schoenus.nanus", "Senecio.glossanthus", "Trachymene.cyanopetala",
#                    "Velleia.cycnopotamica", "Velleia.rosea", "Waitzia.acuminata")]
# mean(colSums(nat.stems))
# 
# ex.stems <- Nj3[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
#                    "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
#                    "Vulpia.bromoides")]
# mean(colSums(ex.stems))

# seeding natives (add any number) (all natives or can pick specific ones - eg deviation species from stan model)
Nj2[,c("Brachyscome.sp" , "Brachyscome.iberidifolia", "Crassula.colorata.or.exserta","Calotis.multicaulis", "Calotis.hispidula", 
       "Calandrinia.eremaea","Daucus.glochidiatus", "Gilberta.tenuifolia", "Phyllangium.sulcatum",
       "Gonocarpus.nodulosus", "Gilberta.tenuifolia", "Goodenia.berardiana", "Haloragis.odontocarpa", "Hyalosperma.glutinosum.subsp.glutinosum", "Hydrocotyle.pilifera", 
       "Lawrencella.rosea", "Nicotiana.rotundifolia", "Plantago.debilis", "Podolepis.canescens", "Podotheca.gnaphalioides", "Rhodanthe.laevis", "Rhodanthe.citrina",
       "Rhodanthe.manglesii", "Rhodanthe.polycephala", "Schoenia.cassiniana", "Schoenus.nanus", "Senecio.glossanthus", "Trachymene.cyanopetala",
       "Velleia.cycnopotamica", "Velleia.rosea", "Waitzia.acuminata")] <- Nj2[,c("Brachyscome.sp" , "Brachyscome.iberidifolia", "Crassula.colorata.or.exserta","Calotis.multicaulis", "Calotis.hispidula", 
                                                                                 "Calandrinia.eremaea","Daucus.glochidiatus", "Gilberta.tenuifolia", "Phyllangium.sulcatum",
                                                                                 "Gonocarpus.nodulosus", "Gilberta.tenuifolia", "Goodenia.berardiana", "Haloragis.odontocarpa", "Hyalosperma.glutinosum.subsp.glutinosum", "Hydrocotyle.pilifera", 
                                                                                 "Lawrencella.rosea", "Nicotiana.rotundifolia", "Plantago.debilis", "Podolepis.canescens", "Podotheca.gnaphalioides", "Rhodanthe.laevis", "Rhodanthe.citrina",
                                                                                 "Rhodanthe.manglesii", "Rhodanthe.polycephala", "Schoenia.cassiniana", "Schoenus.nanus", "Senecio.glossanthus", "Trachymene.cyanopetala",
                                                                                 "Velleia.cycnopotamica", "Velleia.rosea", "Waitzia.acuminata")]+160 # tripling native stem density
Nj3[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
       "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
       "Vulpia.bromoides")] <- Nj3[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
                                      "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
                                      "Vulpia.bromoides")]*0
Nj4[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
       "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
       "Vulpia.bromoides")] <- Nj4[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
                                      "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
                                      "Vulpia.bromoides")]*.5

Nj5[,c("Waitzia.acuminata", "Hypochaeris.glabra")] <- Nj5[,c("Waitzia.acuminata", "Hypochaeris.glabra")]*.5
Nj6[,c("Waitzia.acuminata", "Hypochaeris.glabra")] <- Nj6[,c("Waitzia.acuminata", "Hypochaeris.glabra")]*.75



Nj[,45]<-0
NjB <- Nj[1:28,]
Nj2[,45]<-0
Nj2B <- Nj2[1:28,]
Nj3[,45]<-0
Nj3B <- Nj3[1:28,]
Nj4[,45]<-0
Nj4B <- Nj4[1:28,]
Nj5[,45]<-0
Nj5B <- Nj5[1:28,]
Nj6[,45]<-0
Nj6B <- Nj6[1:28,]


#environmental gradient
bend.env <- env[1:28]
#env_gradient <- cbind(env[1:28]-2, env[1:28]-1, env[1:28], 1+env[1:28], 2+env[1:28]) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation
#env_gradient <- cbind(rep(mean(env[1:51]-2),51), rep(mean(env[1:51]-1),51), rep(mean(env[1:51]),51), rep(mean(1+env[1:51]),51), rep(mean(2+env[1:51]),51))
env_gradient <- cbind(rep(sort(bend.env)[1:14],2), rep(sort(bend.env)[15:28],2)) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation

# set up coexistence simulation
posteriors=9000 # we want to sample through all the posterior values  
plots=28 # run through each option for resident community population size 
arca.into.observedB <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)


# Nj ####
for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[1,]) #to get your alpha term for each environmental condition and each competitor [1,] for reserve 1
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #na.omit(Nj)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[1])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*NjB[r,]))
      #Nj needs to be array of plots by reserve by species and then subscripted to match a_eij [plot r and reserve 1]
      # calculate LDGR of WAAC
      arca.into.observedB[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# Nj2 ####
arca.into.observed1B <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[1,]) #to get your alpha term for each environmental condition and each competitor
      y <- sample(seq(1,4500),1)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[1])
      
      # invade ARCA
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj2B[r,]))
      
      # calculate LDGR of WAAC
      arca.into.observed1B[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# Nj3 ####
arca.into.observed2B <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[1,]) #to get your alpha term for each environmental condition and each competitor
      y <- sample(seq(1,4500),1)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[1])
      
      # invade ARCA
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj3B[r,]))
      
      # calculate LDGR of WAAC
      arca.into.observed2B[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# Nj4 ####
arca.into.observed3B <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[1,]) #to get your alpha term for each environmental condition and each competitor
      y <- sample(seq(1,4500),1)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[1])
      
      # invade ARCA
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj4B[r,]))
      
      # calculate LDGR of WAAC
      arca.into.observed3B[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}


#waac.into.observed[1,,]
#plot(density(waac.into.observed[,,1]))
#plot(density(waac.into.observed[,,2]))

# Nj5 ####
arca.into.observed5B <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[1,]) #to get your alpha term for each environmental condition and each competitor [1,] for reserve 1
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #na.omit(Nj)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[1])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj5B[r,]))
      #Nj needs to be array of plots by reserve by species and then subscripted to match a_eij [plot r and reserve 1]
      # calculate LDGR of WAAC
      arca.into.observed5B[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# Nj6 ####
arca.into.observed6B <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[1,]) #to get your alpha term for each environmental condition and each competitor [1,] for reserve 1
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #na.omit(Nj)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[1])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj6B[r,]))
      #Nj needs to be array of plots by reserve by species and then subscripted to match a_eij [plot r and reserve 1]
      # calculate LDGR of WAAC
      arca.into.observed6B[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# plot ####
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggridges)

data.melt <- melt(arca.into.observedB, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt$e <- as.factor(data.melt$e)
data.melt$plot <- as.factor(data.melt$plot)
data.melt$comp <- "seednat*2"

data.melt1 <- melt(arca.into.observed1B, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt1$e <- as.factor(data.melt1$e)
data.melt1$plot <- as.factor(data.melt1$plot)
data.melt1$comp <- "seednat*3"

data.melt2 <- melt(arca.into.observed2B, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt2$e <- as.factor(data.melt2$e)
data.melt2$plot <- as.factor(data.melt2$plot)
data.melt2$comp <- "weedex.75"

data.melt3 <- melt(arca.into.observed3B, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt3$e <- as.factor(data.melt3$e)
data.melt3$plot <- as.factor(data.melt3$plot)
data.melt3$comp <- "weedex.5"

data.melt4 <- melt(arca.into.observed5B, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt4$e <- as.factor(data.melt4$e)
data.melt4$plot <- as.factor(data.melt4$plot)
data.melt4$comp <- "Nj5"

data.melt5 <- melt(arca.into.observed6B, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt5$e <- as.factor(data.melt5$e)
data.melt5$plot <- as.factor(data.melt5$plot)
data.melt5$comp <- "Nj6"


dataB <- full_join(data.melt, data.melt1)
dataB <- full_join(dataB, data.melt2)
dataB <- full_join(dataB, data.melt3)
dataB <- full_join(dataB, data.melt4)
dataB <- full_join(dataB, data.melt5)

saveRDS(dataB, file = "Sim data/phos&comp_bend_bd.rds")


# PJ --------------------------------------------------------------
#environmental gradient
pj.env <- env[29:95]
#env_gradient <- cbind(env[1:28]-2, env[1:28]-1, env[1:28], 1+env[1:28], 2+env[1:28]) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation
#env_gradient <- cbind(rep(mean(env[1:51]-2),51), rep(mean(env[1:51]-1),51), rep(mean(env[1:51]),51), rep(mean(1+env[1:51]),51), rep(mean(2+env[1:51]),51))
env_gradient <- cbind(rep(sort(pj.env)[1:34],2), rep(sort(pj.env)[34:67],2)) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation


NjP <- Nj[29:95,]
Nj2P <- Nj2[29:95,]
Nj3P <- Nj3[29:95,]
Nj4P <- Nj4[29:95,]



# set up coexistence simulation
posteriors=9000 # we want to sample through all the posterior values  
plots=67 # run through each option for resident community population size 
arca.into.observed <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

# Nj ####
for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #na.omit(Nj)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade WAAC
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*NjP[r,]))
      
      # calculate LDGR of WAAC
      arca.into.observed[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# Nj2 ####
arca.into.observed1 <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      y <- sample(seq(1,4500),1)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade ARCA
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj2P[r,]))
      
      # calculate LDGR of WAAC
      arca.into.observed1[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# Nj3 ####
arca.into.observed2 <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      y <- sample(seq(1,4500),1)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade ARCA
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj3P[r,]))
      
      # calculate LDGR of WAAC
      arca.into.observed2[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# Nj4 ####
arca.into.observed3 <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient))) # make an empty matrix (array for env??)

for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      
      log_a_eij <- alphas[x,1]+Inclusion_ij*arca.phos$alpha_hat_ij[x,,]+ (alphas[x,2]+Inclusion_eij*arca.phos$alpha_hat_eij[x,,])*env_gradient[r,e] 
      
      #Gives a 2x45 matrix of log_alpha_e,i,j for speciesxlocation
      a_eij <- exp(log_a_eij[2,]) #to get your alpha term for each environmental condition and each competitor
      y <- sample(seq(1,4500),1)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade ARCA
      arca_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Nj4P[r,]))
      
      # calculate LDGR of WAAC
      arca.into.observed3[r,x,e] <- log(arca_tp1/1)
      
    }
  }
}

# plot ####
library(reshape2)
library(ggplot2)
library(ggridges)
library(dplyr)

data.melt <- melt(arca.into.observed, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt$e <- as.factor(data.melt$e)
data.melt$plot <- as.factor(data.melt$plot)
data.melt$comp <- "seednat*2"

data.melt1 <- melt(arca.into.observed1, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt1$e <- as.factor(data.melt1$e)
data.melt1$plot <- as.factor(data.melt1$plot)
data.melt1$comp <- "seednat*3"

data.melt2 <- melt(arca.into.observed2, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt2$e <- as.factor(data.melt2$e)
data.melt2$plot <- as.factor(data.melt2$plot)
data.melt2$comp <- "weedex.5"

data.melt3 <- melt(arca.into.observed3, varnames = c("plot", "post", "e"), value.name = "ldgr")
data.melt3$e <- as.factor(data.melt3$e)
data.melt3$plot <- as.factor(data.melt3$plot)
data.melt3$comp <- "weedex.75"


data <- full_join(data.melt, data.melt1)
data <- full_join(data, data.melt2)
data <- full_join(data, data.melt3)

saveRDS(data, file = "Sim data/phos&comp_gradient_waac_pj.rds")
