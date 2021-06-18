library(rstan)
library(coda)
rm(list=ls())

# load in the rdata posteriors - just do one focal species in one script 

load("Topher model fits/WAAC_Shade_Finalfit.rdata")
waac.phos <- rstan::extract(FinalFit)
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
lambda = waac.phos$lambdas
alphas = waac.phos$alphas

# Bendering --------------------------------------------------------------------------------------------
colnames(SpMatrix) <- SpNames
Nj <- SpMatrix 
Nj2 <- SpMatrix 
Nj3 <- SpMatrix
Nj4 <- SpMatrix
Nj5 <- SpMatrix 
Nj6 <- SpMatrix
Nj7 <- SpMatrix
Nj8 <- SpMatrix

#colnames(SpMatrix)

# calculating average stem abundances or natives and exotics 
#  nat.stems <- SpMatrix[,c("Brachyscome.sp" , "Brachyscome.iberidifolia", "Crassula.colorata.or.exserta","Calotis.multicaulis", "Calotis.hispidula", 
#                     "Calandrinia.eremaea","Daucus.glochidiatus", "Gilberta.tenuifolia", "Phyllangium.sulcatum",
#                     "Gonocarpus.nodulosus", "Gilberta.tenuifolia", "Goodenia.berardiana", "Haloragis.odontocarpa", "Hyalosperma.glutinosum.subsp.glutinosum", "Hydrocotyle.pilifera", 
#                     "Lawrencella.rosea", "Nicotiana.rotundifolia", "Plantago.debilis", "Podolepis.canescens", "Podotheca.gnaphalioides", "Rhodanthe.laevis", "Rhodanthe.citrina",
#                     "Rhodanthe.manglesii", "Rhodanthe.polycephala", "Schoenia.cassiniana", "Schoenus.nanus", "Senecio.glossanthus", "Trachymene.cyanopetala",
#                     "Velleia.cycnopotamica", "Velleia.rosea", "Waitzia.acuminata")]
#  mean(colSums(nat.stems))
# 
# ex.stems <- SpMatrix[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
 #                    "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
 #                    "Vulpia.bromoides")]
 # mean(colSums(ex.stems))

#double natives and remove ex
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
                                                                                 "Velleia.cycnopotamica", "Velleia.rosea", "Waitzia.acuminata")]+80 # double native stem density

Nj2[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
       "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
       "Vulpia.bromoides")] <- Nj2[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
                                      "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
                                      "Vulpia.bromoides")]*0

#remove ex
Nj3[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
       "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
       "Vulpia.bromoides")] <- Nj3[,c("Arctotheca.calendula", "Brassica.tournefortii", "Ehrharta.longiflora", "Hypochaeris.glabra", "Lysimachia.arvensis", "Monoculus.monstrosus",
                                      "Oxalis.sp", "Pentaschistis.airoides", "Petrorhagia.velutina", "Sonchus.oleraceus",
                                      "Vulpia.bromoides")]*0

#remove HYPO and thin common natives
Nj4[,c("Gonocarpus.nodulosus", "Trachymene.cyanopetala","Lawrencella.rosea")] <- Nj4[,c("Gonocarpus.nodulosus", "Trachymene.cyanopetala","Lawrencella.rosea")]*0.5 
Nj4[,c("Hypochaeris.glabra")] <- Nj4[,c("Hypochaeris.glabra")]*0

#remove HYPO and thin common natives
Nj5[,c("Gonocarpus.nodulosus", "Trachymene.cyanopetala","Lawrencella.rosea")] <- Nj5[,c("Gonocarpus.nodulosus", "Trachymene.cyanopetala","Lawrencella.rosea")]*0.5 
Nj5[,c("Hypochaeris.glabra")] <- Nj5[,c("Hypochaeris.glabra")]*0


Nj[,45]<-0
NjB <- Nj[1:28,]
Nj2[,45]<-0
Nj2B <- Nj2[1:28,]
Nj3[,45]<-0
Nj3B <- Nj3[1:28,]
Nj4[,45]<-0
Nj4B <- Nj4[1:28,]
#Nj5[,45]<-0
Nj5B <- Nj5[1:28,]


Neighhoods <- list(NjB, Nj2B, Nj3B, Nj4B*0)

#environmental gradient
bend.env <- env[1:51] 
env_gradient <- cbind(bend.env-2, bend.env-1, bend.env, bend.env+1, bend.env+2) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation

# set up coexistence simulation
posteriors=4500 # we want to sample through all the posterior values  
plots=28 # run through each option for resident community population size 
waac.into.observedB <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient), length(Neighhoods))) # make an empty matrix (array for env??)
res=1 # reserve bendering 

# Nj #### 
for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) { 
      for (n in 1:length(Neighhoods)) {
        log_a_eij <- (1-Intra) * waac.phos$alpha_generic[x,1] + Intra * waac.phos$alpha_intra[x,1] + 
          Inclusion_ij[res,] * waac.phos$alpha_hat_ij[x,res,] + ((1-Intra) * waac.phos$alpha_generic[x,2] + Inclusion_eij[res,] * waac.phos$alpha_hat_eij[x,res,] + Intra * waac.phos$alpha_intra[x,2])*env_gradient[r,e]
        a_eij <- exp(log_a_eij) 
        
      y <- sample(seq(1,4500),1)
      # calculate resident abundance
      #Nj <- SpMatrix[r,] # do this as a sample through values too? add y <- sample(seq(1,129),1) then should be able to change number of runs
      #na.omit(Nj)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[1])
      
      # invade WAAC
      waac_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Neighhoods[[n]][r,]))
      #Nj needs to be array of plots by reserve by species and then subscripted to match a_eij [plot r and reserve 1]
      # calculate LDGR of WAAC
      waac.into.observedB[r,x,e,n] <- log(waac_tp1/1)
      }
    }
  }
}

# plot ####
library(reshape2)
library(plyr)

dataB<- melt(waac.into.observedB, varnames = c("plot", "post", "e", "comp"), value.name = "ldgr")
dataB$e <- as.factor(dataB$e)
dataB$plot <- as.factor(dataB$plot)
dataB$comp <- as.factor(dataB$comp)
#dataB$comp <- revalue(dataB$comp, c("1"="observed", "2"="double.nat", "3"="remove.ex", "4"="double.ex", "5"="remove.dev.hypo", "6"="double.dev.hypo", "7"="seed.nat&weed.hypo", "8"="no.comp"))

saveRDS(dataB, file = "Sim data/shade&comp_waac_bd_new.rds")


# PJ --------------------------------------------------------------
#environmental gradient
pj.env <- env[52:129] # make is Perenjori specific 
env_gradient <- cbind(pj.env-2, pj.env-1, pj.env, pj.env+1, pj.env+2) # add or subtract (or a subset of all low values) adding one is shifting by one standard deviation


NjP <- Nj[29:95,]
Nj2P <- Nj2[29:95,]
Nj3P <- Nj3[29:95,]
Nj4P <- Nj4[29:95,]
Nj5P <- Nj5[29:95,]

Neighhoods <- list(NjP, Nj2P, Nj3P, Nj4P*0, Nj5P)


# set up coexistence simulation
posteriors=4500 # we want to sample through all the posterior values  
plots=67 # run through each option for resident community population size 
waac.into.observed <- array(data=NA, dim =c(plots, posteriors, ncol(env_gradient), length(Neighhoods))) # make an empty matrix (array for env??)
res=2
# Nj ####
for (r in 1:plots) {
  for (x in 1:posteriors) {
    for (e in 1:ncol(env_gradient)) {
      for (n in 1:length(Neighhoods)) {
        log_a_eij <- (1-Intra) * waac.phos$alpha_generic[x,1] + Intra * waac.phos$alpha_intra[x,1] + 
          Inclusion_ij[res,] * waac.phos$alpha_hat_ij[x,res,] + ((1-Intra) * waac.phos$alpha_generic[x,2] + Inclusion_eij[res,] * waac.phos$alpha_hat_eij[x,res,] + Intra * waac.phos$alpha_intra[x,2])*env_gradient[r,e]
        a_eij <- exp(log_a_eij) 
        
      y <- sample(seq(1,4500),1)
      
      # calculate lambda 
      overall_lambda <- lambda[x,,1]+lambda[x,,2]*env_gradient[r,e] # which order is the 2 and 2 for lambda? intercept/slope, reserve?
      lambdas <- exp(overall_lambda[2])
      
      # invade WAAC
      waac_tp1 <- (1-germ[y])*surv[y]*1 + 
        1*germ[y]*lambdas/(1+sum(a_eij*Neighhoods[[n]][r,]))
      
      # calculate LDGR of WAAC
      waac.into.observed[r,x,e,n] <- log(waac_tp1/1)
      }
    }
  }
}

# plot ####
data<- melt(waac.into.observed, varnames = c("plot", "post", "e", "comp"), value.name = "ldgr")
data$e <- as.factor(data$e)
data$plot <- as.factor(data$plot)
data$comp <- as.factor(data$comp)
#data$comp <- revalue(data$comp, c("1"="observed", "2"="double.nat", "3"="remove.ex", "4"="double.ex", "5"="no.comp"))

saveRDS(data, file = "Sim data/shade&comp_gradient_waac_pj.rds")

