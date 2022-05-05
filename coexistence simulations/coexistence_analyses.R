#update by Lauren Shoemaker Jan 16th

library(rstan)
library(coda)
library(HDInterval)
library(tidyverse)
library(hrbrthemes)
library(kableExtra)
library(ggridges)

rm(list=ls())

location <- 2 # set reserve to be Bendering 
             # a 2 will set this to be Perenjori

load("Topher model fits/WAAC_Shade_FinalFit.rdata")
Post <- rstan::extract(FinalFit)
Shade <- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity,
             reserve = reserve, SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij,
             Inclusion_eij = Inclusion_eij, Intra = Intra)
Shade_Params <- Shade$Post
remove(FinalFit)

load("Topher model fits/WAAC_Phos_FinalFit.rdata")
Post <- rstan::extract(FinalFit)
Phos <- list(Post = Post, SpNames = SpNames, N = N, S = S, Fecundity = Fecundity,
              reserve = reserve, SpMatrix = SpMatrix, env = env, Inclusion_ij = Inclusion_ij,
              Inclusion_eij = Inclusion_eij, Intra = Intra)
Phos_Params <- Phos$Post
remove(FinalFit)

env <- seq(-2, 2, by=0.1)
lambdas_shade <- Shade_Params$lambdas
lambda_shade_gradient <- array(NA, c(3, length(env)))

lambdas_phos <- Phos_Params$lambdas
lambda_phos_gradient <- array(NA, c(3, length(env)))
alphas_phos_gradient <- array(NA, c(3, length(env), Phos$S))

for (i in 1:length(env)) {
  lambda_shade_post <- exp(lambdas_shade[,location,1] + 
                             lambdas_shade[,location,2] * env[i])
  lambda_shade_gradient[1,i] <- mean(lambda_shade_post)
  lambda_shade_gradient[2:3,i] <- HDInterval::hdi(lambda_shade_post)
  
  lambda_phos_post <- exp(lambdas_phos[,location,1] + 
                             lambdas_phos[,location,2] * env[i])
  lambda_phos_gradient[1,i] <- mean(lambda_phos_post)
  lambda_phos_gradient[2:3,i] <- HDInterval::hdi(lambda_phos_post)
  
  alphas_phos_post <- matrix(NA, nrow=4500, ncol=Phos$S)
  for (xx in 1:4500) {
  alphas_phos_post[xx,] <- exp((1-Phos$Intra) * Phos_Params$alpha_generic[xx,1] + # generic, no env grad
                  Phos$Intra * Phos_Params$alpha_intra[xx,1] +     # intras, no env grad
                  Phos$Inclusion_ij[location,] * Phos_Params$alpha_hat_ij[xx,location,] + # non generic, no env grad
                  ((1-Phos$Intra) * Phos_Params$alpha_generic[xx,2] + # env grad, generic
                     Phos$Inclusion_eij[location,] * Phos_Params$alpha_hat_eij[xx,location,] + # slope, non generic 
                     Phos$Intra * Phos_Params$alpha_intra[xx,2]) * env[i])
  }
  
  for (s in 1:Phos$S) {
    alphas_phos_gradient[1,i,s] <-  mean(alphas_phos_post[,s])
    alphas_phos_gradient[2:3,i,s] <- HDInterval::hdi(alphas_phos_post[,s])
  }
  
}

env_polygon <- c(env, rev(env))
lambda_shade_polygon <- c(lambda_shade_gradient[2,], rev(lambda_shade_gradient[3,]))
lambda_phos_polygon <- c(lambda_phos_gradient[2,], rev(lambda_phos_gradient[3,]))

quartz()
plot(env, lambda_shade_gradient[1,], type="l", lwd=3, col="darkblue",
     xlim=c(-2,2), ylim=c(0,22), xlab="Environmental Gradient (Scaled)", ylab="Lambda (Density Independent Seed Production)")
polygon(env_polygon, lambda_shade_polygon, col=rgb(0,0,1,.2), border=rgb(0,0,1,.2))

lines(env, lambda_phos_gradient[1,], lty=1, lwd=3, col="darkgreen")
polygon(env_polygon, lambda_phos_polygon, col=rgb(0,1,0,.2), border=rgb(0,1,0,.2))
legend("topright", c("Shade", "Phosphorous"), col=c("darkblue", "darkgreen"), lwd=3)
legend("topright", c("Shade", "Phosphorous"), col=c(rgb(0,0,1,.2), rgb(0,1,0,.2)), lwd=15)

plot(env, alphas_phos_gradient[1,,1], type="l", lwd=2, col="darkblue",
     xlim=c(-2,2), ylim=c(0,.1), xlab="Environmental Gradient", ylab="Alphas")
for (s in 1:Phos$S) {
  lines(env, alphas_phos_gradient[1,,s], lty=1, lwd=2, col="darkgreen")
}


# load in s and g data 
load("SurvivalAndGermination/Germination.rdata")
germ <- rstan::extract(PrelimFit)
germ<-as.data.frame(germ)
remove(PrelimFit)

load("SurvivalAndGermination/Survival.rdata")
surv <- rstan::extract(PrelimFit)
surv<-as.data.frame(surv)
remove(PrelimFit)

g_waac<- germ$p.2
#germination <- data.frame(arca=germ$p.1, waac=germ$p.2)
s_waac <-surv$p.2 
#survival <- data.frame(arca=surv$p.1, waac=surv$p.2)

# coexistence simulations - invasion growth rates -----------------------------------------
# simulation of baseline community
plots <- which(Phos$reserve == location)
posteriors <- sample(4500, 1000, replac=FALSE)
obs_waac_ldgr <- data.frame(matrix(ncol=7, nrow=1000*length(plots)))
names(obs_waac_ldgr) <- c("env", "replicate", "obs_ldgr", "nocomp_ldgr", "lowp_ldgr", "nocomp_lowp_ldgr", "lowshade_ldgr")
lowp <- rep(-2, 129)
lowshade <- rep(-2, 129)
#obs_waac_ldgr <- matrix(NA, length(plots), length(posteriors))
neighbors <- Phos$SpMatrix[,]
neighbors[,45] <- 0
no_neighbors <- matrix(0, 129, 45)
counter <- 1
for (aa in plots) {
  post_counter <- 1
  for (bb in posteriors) {
      lambdas <- exp(Phos_Params$lambdas[bb,location,1] + Phos_Params$lambdas[bb,location,2]*Phos$env[aa])
      lowp_lambdas <- exp(Phos_Params$lambdas[bb,location,1] + Phos_Params$lambdas[bb,location,2]*lowp[aa])

      lowshade_lambdas <- exp(Shade_Params$lambdas[bb,location,1] + Shade_Params$lambdas[bb,location,2]*lowshade[aa])
      
      #sp. specific inter
      #intra
      alphas <- exp((1-Phos$Intra) * Phos_Params$alpha_generic[bb,1] + # generic, no env grad
                      Phos$Intra * Phos_Params$alpha_intra[bb,1] +     # intras, no env grad
                      Phos$Inclusion_ij[location,] * Phos_Params$alpha_hat_ij[bb,location,] + # non generic, no env grad
                      ((1-Phos$Intra) * Phos_Params$alpha_generic[bb,2] + # env grad, generic
                         Phos$Inclusion_eij[location,] * Phos_Params$alpha_hat_eij[bb,location,] + # slope, non generic 
                         Phos$Intra * Phos_Params$alpha_intra[bb,2]) * Phos$env[aa]) # slope intra
      
      lowp_alphas <- exp((1-Phos$Intra) * Phos_Params$alpha_generic[bb,1] + # generic, no env grad
                      Phos$Intra * Phos_Params$alpha_intra[bb,1] +     # intras, no env grad
                      Phos$Inclusion_ij[location,] * Phos_Params$alpha_hat_ij[bb,location,] + # non generic, no env grad
                      ((1-Phos$Intra) * Phos_Params$alpha_generic[bb,2] + # env grad, generic
                         Phos$Inclusion_eij[location,] * Phos_Params$alpha_hat_eij[bb,location,] + # slope, non generic 
                         Phos$Intra * Phos_Params$alpha_intra[bb,2]) * lowp[aa]) # slope intra
      
      lowshade_alphas <- exp((1-Shade$Intra) * Shade_Params$alpha_generic[bb,1] + # generic, no env grad
                           Shade$Intra * Shade_Params$alpha_intra[bb,1] +     # intras, no env grad
                           Shade$Inclusion_ij[location,] * Shade_Params$alpha_hat_ij[bb,location,] + # non generic, no env grad
                           ((1-Shade$Intra) * Shade_Params$alpha_generic[bb,2] + # env grad, generic
                              Shade$Inclusion_eij[location,] * Shade_Params$alpha_hat_eij[bb,location,] + # slope, non generic 
                              Shade$Intra * Shade_Params$alpha_intra[bb,2]) * lowshade[aa]) # slope intra
      
        
        # invade WAAC
        waac_growth <- (1-g_waac[bb])*s_waac[bb]*1 + 
          1*g_waac[bb]*lambdas/(1+sum(alphas*neighbors[aa,]))
        
        no_compwaac_growth <- (1-g_waac[bb])*s_waac[bb]*1 + 
          1*g_waac[bb]*lambdas/(1+sum(alphas*no_neighbors[aa,]))
        
       lowp_waac_growth <- (1-g_waac[bb])*s_waac[bb]*1 + 
          1*g_waac[bb]*lowp_lambdas/(1+sum(lowp_alphas*neighbors[aa,]))
       
       lowshade_waac_growth <- (1-g_waac[bb])*s_waac[bb]*1 + 
         1*g_waac[bb]*lowshade_lambdas/(1+sum(lowshade_alphas*neighbors[aa,]))
        
        nocomp_lowp_waac_growth <- (1-g_waac[bb])*s_waac[bb]*1 + 
          1*g_waac[bb]*lowp_lambdas/(1+sum(lowp_alphas*no_neighbors[aa,]))
        #Nj needs to be array of plots by reserve by species and then subscripted to match a_eij [plot r and reserve 1]
        # calculate LDGR of WAAC
        obs_waac_ldgr[counter,1] <- Phos$env[aa]
        obs_waac_ldgr[counter,2] <- post_counter
        obs_waac_ldgr[counter,3] <- log(waac_growth/1)
        obs_waac_ldgr[counter,4] <- log(no_compwaac_growth/1)
        obs_waac_ldgr[counter,5] <- log(lowp_waac_growth/1)
        obs_waac_ldgr[counter,6] <- log(nocomp_lowp_waac_growth/1)
        obs_waac_ldgr[counter,7] <- log(lowshade_waac_growth/1)

        post_counter <- post_counter + 1
        counter <- counter + 1
  }
  #counter <- counter + 1
}



quartz()
ggplot(obs_waac_ldgr, aes(x=ldgr, y=as.factor(1)))+
  geom_density_ridges(aes(x=obs_ldgr), alpha=0.5, fill="blue",scale=.5)+
  geom_density_ridges(aes(x=nocomp_ldgr), alpha=0.5, fill="grey",scale=.5)+
  geom_density_ridges(aes(x=lowp_ldgr), alpha=0.5, fill="green",scale=.5)+
  geom_density_ridges(aes(x=nocomp_lowp_ldgr), alpha=0.5, fill="purple",scale=.5) +
  theme_classic() +
  scale_x_continuous(name="Waitzia Low Density Growth Rate", breaks=c(-0.5, 0, .5, 1, 1.5, 2)) +
  geom_vline(xintercept = 0, color="black", size=2) +
  scale_y_discrete(name="Probability Density") 

cbp2 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

quartz(width=6, height=6)
ggplot(obs_waac_ldgr)+
  geom_density(aes(x=obs_ldgr, fill="red"))+
  geom_density(aes(x=nocomp_ldgr, fill="pink"))+
  geom_density(aes(x=lowp_ldgr, fill="green"))+
  geom_density(aes(x=lowshade_ldgr, fill="yellow"))+
  geom_density(aes(x=nocomp_lowp_ldgr, fill="blue")) +
  theme_classic() +
  scale_x_continuous(name="Waitzia Low Density Growth Rate", breaks=c(-0.5, 0, .5, 1, 1.5, 2)) +
  geom_vline(xintercept = 0, color="black", size=1.5) +
  scale_y_discrete(name="Probability Density") +
   scale_fill_manual(guide='legend', name = 'Scenario',
                       values=alpha(c("red"=cbp2[2],"pink"=cbp2[8],
                                "green"=cbp2[4], "yellow"=cbp2[1], "blue"=cbp2[6]),.5),
                        labels=c("Observed", "No Comp", "Low P",
                                 "Low Shade", "No Comp, Low P")) +
  theme(text = element_text(size = 14), legend.position = c(.8,.8))  


#theme(legend.position = c(1.0,.2))




ggplot(obs_waac_ldgr, aes(x=ldgr, y=as.factor(env)))+
  geom_density_ridges_gradient(aes(fill=..x..))+
  scale_fill_gradient(low="orange", high="navy") +
  scale_fill_viridis_c(name = "LDGR", direction = -1)


ggplot(obs_waac_ldgr, aes(x=ldgr, color=env, fill=env)) +
  geom_density_ridges_gradient(alpha=0.6) #+
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  geom_text( data=annot, aes(x=x, y=y, label=text, color=text), hjust=0, size=4.5) +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Assigned Probability (%)")

library(reshape2)
data.melt <- melt(obs_waac_ldgr, varnames = c("plot", "index"), value.name = "ldgr")
str(data.melt)
data.melt$index <- as.factor(data.melt$index)
data.melt$plot <- as.factor(data.melt$plot)
str(data.melt)

library(ggplot2)
#install.packages("ggridges")
library(ggridges)

ggplot(data.melt, aes(x = ldgr, y = index, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  #scale_fill_viridis_c(name = "Tail probability", direction = -1, option="plasma") +
  geom_vline(xintercept = 0, linetype = "dashed")








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

